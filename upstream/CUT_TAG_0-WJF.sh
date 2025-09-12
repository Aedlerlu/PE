#!/bin/sh
# 实验数据分析上游处理脚本
# 功能：处理双端测序数据，包括质量控制（fastp+fastqc）、序列比对、数据过滤和重复去除
# 输入：通过命令行参数$1指定的样本重复ID

    #命名： ${source}_${upload_date}_${tech}_${species}_${modification}_${sample}_${experiment}
    ###   5级目录定义方法： ~/cell/projects/测序技术/物种/修饰/细胞类型/单个项目/*

    #使用方法：
    #    1.改参考基因组路径   
    #    2.qsub run_job.sh  # 正确：qsub 直接跟 PBS 脚本路径
    #      qstat -u 你的用户名查看作业状态，确认Resources Used中的ncpus是否为你申请的核心数（如 32）。
    #      qdel 作业ID


# 实验基本信息配置
source="NHZY"
batch_num="X101SC24121131-Z01-J050"

save_date="."
seq_tech="CUT_TAG"               # 实验技术名称
species="human"                  # 物种类型
modification="undiff_blood"           # 目标修饰/蛋白名称
sample="blood"                      # 样本组织类型
rep="HC02599289K27A"                  # 从命令行接收的重复样本ID（子文件夹名）
experiment=$rep                  # 实验名称，与重复样本ID相同
fqgz_name="HC02599289K27A_SKDL250020185-1A_22WWNWLT4_L2" #（不带_1 _2,不带后缀）
tackle="0_assess"                  # 分析方法

# 计算资源配置
queen=core24                     # 指定计算队列
ppn=64                          # 进程数（每节点处理器核心数）
PPN=64                          # 冗余定义的进程数
mem=100gb                        # 所需内存大小

# 当前设备静态路径（基础路径）
cell_space="/home/data/weijianfan/cell" 

cell_biosoft="${cell_space}/biosoft"    
cell_ref_gene="${cell_space}/ref_gene"
cell_projects="${cell_space}/projects"
cell_rawdata="${cell_space}/rawdata"
#本项目名
#proj_folder_name="${source}_${batch_num}_${save_date}_${seq_tech}_${species}_${modification}_${sample}_${rep}"
proj_folder_name="${batch_num}"

cell_proj="${cell_projects}/$seq_tech/${species}/${modification}/${sample}/${proj_folder_name}"  # 核心工作目录变量
cell_proj_tackle="${cell_proj}/${rep}_${tackle}"  # 基于核心目录的通用子目录


# 原始测序数据路径定义
IN1="${cell_rawdata}/$proj_folder_name/01.RawData/${rep}/${fqgz_name}_1.fq.gz"
IN2="${cell_rawdata}/$proj_folder_name/01.RawData/${rep}/${fqgz_name}_2.fq.gz"

# 输入文件存在性检查
if [ ! -e $IN1 ] || [ ! -e $IN2 ] 
then
  echo -e "错误：未找到输入文件:\n $IN1 \n 或 \n $IN2:"
  exit  # 若文件不存在则终止脚本
fi

# 检查双端文件是否相同（不允许相同）
if [ $IN1 = $IN2 ]
then
  echo -e "错误：双端测序文件不能相同:\n IN1=$IN1 \n IN2=$IN2"
  exit  # 若文件相同则终止脚本
fi

# 路径配置
mydisk="$cell_space"      # 数据存储根目录
softbin="$cell_biosoft"      # 分析软件安装目录


# 创建输出目录并备份当前脚本
mkdir -p $cell_proj_tackle/
backup=`pwd`  # 记录当前工作目录路径
# 备份脚本到结果目录
cp "$backup/$(basename $0)" \
    "${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.sh"

# 参考基因组配置
ref="$mydisk/ref_gene/$species/hg38/bowtie2/hg38_index"  # 人类hg38基因组索引
chrsize="$mydisk/ref_gene/$species/hg38/hg38.chrom.sizes"  # 染色体大小信息文件

# 日志文件路径定义
trim_log="${cell_proj_tackle}/trim.log.txt"  # 质量修剪日志
log="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.mapping.log.txt"  # 主分析日志

# fastp质量修剪后的输出文件路径（替换Trimmomatic的输出）
# 配对输出
F_OUT_P="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.F.paired.trimed.fq"  # 正向配对reads
R_OUT_P="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.R.paired.trimed.fq"  # 反向配对reads
# 未配对输出
F_OUT_UP="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.F.unpaired.trimed.fq"  # 正向未配对reads
R_OUT_UP="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.R.unpaired.trimed.fq"  # 反向未配对reads
# fastp报告文件
fastp_html="${cell_proj_tackle}/fastp_report.html"
fastp_json="${cell_proj_tackle}/fastp_report.json"

### 新增：fastqc相关路径定义 ###
# 原始数据fastqc报告目录
fastqc_raw_dir="${cell_proj_tackle}/fastqc_raw"
# 质控后数据fastqc报告目录
fastqc_clean_dir="${cell_proj_tackle}/fastqc_clean"
# 创建fastqc输出目录
mkdir -p $fastqc_raw_dir $fastqc_clean_dir

# 定义后续分析使用的修剪后文件
TrimedIN1=$F_OUT_P  # 正向配对修剪后文件
TrimedIN2=$R_OUT_P  # 反向配对修剪后文件
OUT1="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.mapping.sam"  # 比对结果SAM文件

# 去重复后的最终BAM文件路径
treatment="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.deduplicated.bam"

# 创建作业提交脚本目录
mkdir "${cell_proj_tackle}/dsub_script"
cat >"${cell_proj_tackle}/dsub_script/${experiment}_nohup.sh" <<EOF
#!/bin/sh


# 创建日志文件
touch $log
touch $trim_log

### 新增：对原始数据进行fastqc评估 ###
echo "开始对原始数据进行fastqc质量评估..." >> $log
fastqc \
    --threads $ppn \
    -o $fastqc_raw_dir \
    $IN1 $IN2
echo "原始数据fastqc评估完成，结果在：$fastqc_raw_dir" >> $log
echo "" >> $log

########################################
## 数据质控（使用fastp）
########################################
fastp \
    --thread $ppn \
    -i "$IN1"  -I "$IN2" \
    -o "$F_OUT_P"  -O "$R_OUT_P" \
    --unpaired1 "$F_OUT_UP" \
    --unpaired2 "$R_OUT_UP" \
    --html "$fastp_html" \
    --json "$fastp_json" \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_front_mean_quality 10 \
    --cut_tail --cut_tail_mean_quality 10 \
    --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 \
    --length_required 40 \
    2> $trim_log  # 日志输出

##备用参数
    --adapter_sequence "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
    --adapter_sequence_r2 "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
    --max_len1 100 \
    --max_len2 100 \




### 新增：对质控后的数据进行fastqc评估 ###
echo "开始对fastp处理后的数据进行fastqc质量评估..." >> $log
fastqc \
    --threads $ppn \
    -o $fastqc_clean_dir \
    $F_OUT_P $R_OUT_P $F_OUT_UP $R_OUT_UP  # 评估所有质控后文件（配对+未配对）
echo "质控后数据fastqc评估完成，结果在：$fastqc_clean_dir" >> $log
echo "" >> $log


EOF

# 后台运行生成的作业脚本
nohup sh "${cell_proj_tackle}/dsub_script/${experiment}_nohup.sh" > \
        "${cell_proj_tackle}/${experiment}_output.log" 2>&1 &

# 备份作业脚本到结果目录
cp "${cell_proj_tackle}/dsub_script/${experiment}_nohup.sh" \
    "$cell_proj_tackle/"
echo "nohup已成功提交 ${experiment}_nohup.sh"  # 提示脚本提交成功