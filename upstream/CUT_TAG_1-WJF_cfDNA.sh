#!/bin/sh
# 实验数据分析上游处理脚本
# 功能：处理双端测序数据，包括质量控制、序列比对、数据过滤和重复去除
# 输入：通过命令行参数$1指定的样本重复ID

    #命名： ${source}_${upload_date}_${tech}_${species}_${modification}_${sample}_${experiment}
    ###   5级目录定义方法： ~/cell/projects/测序技术/物种/修饰/细胞类型/单个项目/*

    #使用方法：
    #    1.改参考基因组路径   
    #    2.qsub run_job.sh  # 正确：qsub 直接跟 PBS 脚本路径
    #      qstat -u 你的用户名查看作业状态，确认Resources Used中的ncpus是否为你申请的核心数（如 32）。
    #      qdel 作业ID


# 实验基本信息配置
# source="NHZY"
#batch_num="X101SC24121131-Z01-J083"
batch_num="$1"
save_date="."
seq_tech="CUT_TAG"               # 实验技术名称
species="human"                  # 物种类型
modification="undiff_blood"           # 目标修饰/蛋白名称
sample="blood"                      # 样本组织类型
rep="$2"  
#rep="PE6684654K4M"                # 从命令行接收的重复样本ID（子文件夹名）
experiment=$rep                  # 实验名称，与重复样本ID相同
fqgz_name="$rep" #（不带_1 _2,不带后缀）
tackle="1_cfDNA"                  # 分析方法

# 计算资源配置
queen=core24                     # 指定计算队列
ppn=16                          # 进程数（每节点处理器核心数）
PPN=16                          # 冗余定义的进程数
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
#IN1="/home/data/weijianfan/cell/projects/CUT_TAG/human/undiff_blood/blood/X101SC24121131-Z01-J050/01.RawData/HC02599289K4M/${fqgz_name}_1.fq.gz"
#IN2="/home/data/weijianfan/cell/projects/CUT_TAG/human/undiff_blood/blood/X101SC24121131-Z01-J050/01.RawData/HC02599289K4M/${fqgz_name}_2.fq.gz"
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

# 定义后续分析使用的修剪后文件
TrimedIN1=$F_OUT_P  # 正向配对修剪后文件
TrimedIN2=$R_OUT_P  # 反向配对修剪后文件
OUT1="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.mapping.sam"  # 比对结果SAM文件

# 去重复后的最终BAM文件路径
treatment="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.deduplicated.bam"

# 创建作业提交脚本目录(用于集群提交)
mkdir "${cell_proj_tackle}/dsub_script"
cat >"${cell_proj_tackle}/dsub_script/${experiment}_nohup.sh" <<EOF
#!/bin/sh

# 没有PBS集群
# # PBS作业调度参数
# # PBS -q core24                  # 指定队列名称
# # PBS -l walltime=10:00:00,nodes=$nodes:ppn=$ppn,mem=$mem  # 资源限制
# # PBS -e "${cell_proj_tackle}/${experiment}.err"  # 错误输出文件
# # PBS -o "${cell_proj_tackle}/${experiment}.out"  # 标准输出文件

# 创建日志文件
touch $log
touch $trim_log

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

########################################
## 序列比对（使用bowtie2.5.4）,conda
########################################
# 计算原始数据MD5值，用于数据完整性验证
echo "" >> $log
md5sum $IN1 $IN2 >>$log
echo "" >> $log

#conda安装
bowtie2 \
    -p $ppn \
    -x $ref \
    -1 $TrimedIN1 \
    -2 $TrimedIN2 \
    --very-sensitive \
    -X 2000 \
    -S $OUT1 2>> $log


# 比对成功后清理中间文件（修剪后的reads文件）
if [ -e $OUT1 ]
then
    rm $R_OUT_UP
    rm $R_OUT_P
    rm $F_OUT_P
    rm $F_OUT_UP
    rm $trim_log
fi
###################################################################################

### SAM文件处理：过滤低质量比对并转换为排序的BAM文件
cd $cell_proj_tackle  

# 过滤未比对和低质量比对的reads，并转换为BAM格式
samtools view  \
        -@ $ppn \
        -F 4 \
        -q 10 \
        -bS "${species}.${modification}.${sample}.${rep}.${experiment}.mapping.sam" \
        -o "${species}.${modification}.${sample}.${rep}.${experiment}.mapped.bam"

# 对BAM文件进行排序（按坐标排序）
samtools sort \
        -@ $ppn  \
        ${species}.${modification}.${sample}.${rep}.${experiment}.mapped.bam  \
        -o ${species}.${modification}.${sample}.${rep}.${experiment}.sorted.bam

# ---------------------------
# 15.1  新增：添加Read Group信息（解决Picard报错）
# ---------------------------
java -jar "${softbin}/picard/build/libs/picard.jar" AddOrReplaceReadGroups \
    I="${species}.${modification}.${sample}.${rep}.${experiment}.sorted.bam" \
    O="${species}.${modification}.${sample}.${rep}.${experiment}.sorted.rg.bam" \
    RGID="${experiment}_${rep}" \
    RGLB="${sample}_lib" \
    RGPL=ILLUMINA \
    RGPU="${experiment}_${rep}_unit1" \
    RGSM="${sample}" \
    CREATE_INDEX=true

### 使用picard工具去除PCR重复
mkdir -p "${cell_proj_tackle}/tmp_picard"
java -Xmx400G -jar "${softbin}/picard/build/libs/picard.jar" MarkDuplicates \
    --REMOVE_DUPLICATES true \
    --INPUT "${species}.${modification}.${sample}.${rep}.${experiment}.sorted.rg.bam" \
    --OUTPUT "${species}.${modification}.${sample}.${rep}.${experiment}.deduplicated.bam" \
    --METRICS_FILE "${species}.${modification}.${sample}.${rep}.${experiment}.deduplicate.metrics.txt" \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT \
    --MAX_RECORDS_IN_RAM 80000000 \
    --TMP_DIR "${cell_proj_tackle}/tmp_picard" \
    

# 去重成功后清理中间文件
if [ -e ${species}.${modification}.${sample}.${rep}.${experiment}.deduplicated.bam ]
then
    rm ${species}.${modification}.${sample}.${rep}.${experiment}.mapping.sam
    rm ${species}.${modification}.${sample}.${rep}.${experiment}.sorted.bam
    rm ${species}.${modification}.${sample}.${rep}.${experiment}.mapped.bam
fi

EOF

# 后台运行生成的作业脚本
nohup sh "${cell_proj_tackle}/dsub_script/${experiment}_nohup.sh" > \
        "${cell_proj_tackle}/${experiment}_output.log" 2>&1 &


echo "nohup已成功提交 ${experiment}_nohup.sh"  # 提示脚本提交成功