#!/bin/sh

# 实验基本信息配置
# source="NHZY"
#batch_num="X101SC24121131-Z01-J063"
batch_num="$1"
save_date="."
seq_tech="CUT_TAG"               # 实验技术名称
species="human"                  # 物种类型
modification="undiff_blood"           # 目标修饰/蛋白名称
sample="blood"                      # 样本组织类型
rep="$2"  
#rep="PE6684654K4M"                 # 从命令行接收的重复样本ID（子文件夹名）
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
# 本项目名
# proj_folder_name="${source}_${batch_num}_${save_date}_${seq_tech}_${species}_${modification}_${sample}_${rep}"
proj_folder_name="${batch_num}"

cell_proj="${cell_projects}/$seq_tech/${species}/${modification}/${sample}/${proj_folder_name}"  # 核心工作目录变量
cell_proj_tackle="${cell_proj}/${rep}_${tackle}"  # 基于核心目录的通用子目录
#-------------------------------------

mydisk="$cell_space"
softbin="$cell_biosoft"

# 输入BAM文件路径（处理后的比对文件）
treatment="${cell_proj_tackle}/${species}.${modification}.${sample}.${rep}.${experiment}.deduplicated.bam"

# 检查依赖工具是否存在
command -v macs3 >/dev/null 2>&1 || { echo >&2 "错误: macs3未安装或不在PATH中"; exit 1; }
command -v bamCoverage >/dev/null 2>&1 || { echo >&2 "错误: bamCoverage未安装或不在PATH中"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo >&2 "错误: samtools未安装或不在PATH中"; exit 1; }

# 检查输入文件是否存在
if [ ! -e "$treatment" ]
then
echo -e "error: File not found:$treatment"
exit
fi

# 创建必要的目录结构
# mkdir -p $mydisk/CUT-TAG/${species}/${modification}/${sample}/${rep}/commom/
mkdir -p ${cell_proj_tackle}/peakcalling_narrow

# 备份当前脚本到分析目录
backup=`pwd`
cp "$backup/$(basename $0)" \
        ${cell_proj_tackle}/peakcalling_narrow  ##To backup script for this run!!

# 参考基因组相关文件
#ref=$mydisk/ref/hg19_chip/hg19  ##for bowtie only
chrsize=${cell_ref_gene}/${species}/hg38/hg38.chrom.sizes  # 染色体大小文件，用于转换格式

# 创建并写入PBS作业脚本（用于集群提交）
# mkdir "${cell_proj_tackle}/dsub_script"
cat >"${cell_proj_tackle}/dsub_script/${experiment}_peakcalling_narrow_nohup.sh" <<EOF
#!/bin/sh

# 进入工作目录
cd $cell_proj_tackle  

# 生成标准化的覆盖度BigWig文件（用于可视化）
samtools index ${species}.${modification}.${sample}.${rep}.${experiment}.deduplicated.bam

# 生成标准化的信号轨道（使用spike-in比例）
# 注释掉spike-in标准化相关内容 - 标准ChIP-seq使用常规标准化
# bamCoverage --bam ${species}.${modification}.${sample}.${rep}.${experiment}.deduplicated.bam \
#   --normalizeUsing RPKM --binSize 50 -p $ppn --scaleFactor \$scale_factor \
#   --outFileName ${species}.${modification}.${sample}.${rep}.${experiment}.deeptools.bs50.scalefator\${scale_factor}.bw

生成未标准化的信号轨道（仅RPKM标准化）
bamCoverage --bam ${species}.${modification}.${sample}.${rep}.${experiment}.deduplicated.bam \
  --normalizeUsing RPKM --binSize 50 -p $ppn \
  --outFileName ${species}.${modification}.${sample}.${rep}.${experiment}.deeptools.bs50.bw

#################
# 峰识别分析流程
################

# 创建peak calling输出目录
# 使用统一工作目录cell_proj_tackle
mkdir -p ${cell_proj_tackle}/peakcalling_narrow
cd ${cell_proj_tackle}/peakcalling_narrow

# 使用MACS3进行peak calling (宽峰模式)

# 参数说明:
#     # -t: 处理的BAM文件
#     # --keep-dup all: 保留所有重复 reads
#     # -f BAM: 输入格式为BAM
#     # -g mm: 基因组大小(小鼠)mm     # -p 0.0001: p-value阈值
#     # --SPMR: 输出每百万reads的信号值
#     # -B: 输出bedGraph文件
#     # --broad: 宽峰模式
#     # --nomodel: 不构建片段大小模型
#     # --nolambda: 不计算局部lambda背景
#     # --extsize 200: 延伸片段大小为200bp

macs3 callpeak -t $treatment \
  --keep-dup all -f BAM -g hs -n ${species}.${modification}.${sample}.${rep}.${experiment} \
  -p 0.0001 --SPMR -B  --nomodel --nolambda --extsize 200

# 计算富集倍数(FE)
macs3 bdgcmp \
  -t ${species}.${modification}.${sample}.${rep}.${experiment}_treat_pileup.bdg \
  -c ${species}.${modification}.${sample}.${rep}.${experiment}_control_lambda.bdg \
  -o ${species}.${modification}.${sample}.${rep}.${experiment}_FE.bdg -m FE

# 将bedGraph文件转换为BigWig格式(更适合基因组浏览器可视化)
for i in *.bdg
do
  prefix=\${i%.bdg}
  Tembdg=\${prefix}tem.bdg
  Outbw=\${prefix}.bw
  
  # 系统自带sort，对 bedGraph 文件按染色体名称（第一列）和起始位置（第二列，数字排序）进行排序
  # BigWig 要求输入的 bedGraph 必须按基因组位置有序。
  sort -k1,1 -k2,2n \$i -o \$Tembdg

  #   # 核心命令 bedGraphToBigWig 接收 3 个关键参数：
  # 排序后的临时 bedGraph 文件（$Tembdg）
  # 染色体大小文件（$chrsize，如 hg38.chrom.sizes）
  # 输出的 BigWig 文件名（$Outbw）
  bedGraphToBigWig \$Tembdg $chrsize \$Outbw
  rm \$Tembdg
done
echo "[$(date +'%Y-%m-%d %H:%M:%S')] 所有bedGraph文件已成功转换为BigWig格式，输出路径：${cell_proj_tackle}/peakcalling_narrow"
echo "脚本结束"
EOF

# 在后台运行分析脚本并记录日志
nohup sh ${cell_proj_tackle}/dsub_script/${experiment}_peakcalling_narrow_nohup.sh > ${cell_proj_tackle}/peakcalling_narrow/${experiment}_peakcalling_narrow_output.log 2>&1 &

# 备份脚本并提示成功
cp ${cell_proj_tackle}/dsub_script/${experiment}_peakcalling_narrow_nohup.sh ${cell_proj_tackle}/peakcalling_narrow/
echo "nohup Successful for ${experiment}_peakcalling_narrow_nohup.sh"
