#!/bin/bash
spe="${1}"
name="${2}" 
exp="${3}"
type=${4}
bing=${5}
y=${6}
input=${7}
path="${HOME}/CHIPseq"
alignlog="${HOME}/CHIPseq/results/bowtie2/${name}align.log"
sam="${path}/results/samtools"
if [ ! -d "$sam"  ];then
 mkdir ${path}/results/samtools
 fi
mkdir -p ${path}/trimm/pe/${exp}
mkdir -p ${path}/trimm/se
softbin="${HOME}/software"
dirp="${path}/trimm/pe"
dirs="${path}/trimm/se" #这里应该加一个判断命令，如果文件夹已存在，就不执行,不过不加好像也不会影响
trimlog="${HOME}/CHIPseq/results/trim.log"
touch $trimlog
touch $alignlog
touch ${path}/${name}.log
mkdir -p ${HOME}/CHIPseq/results/sambamba/${exp}/filtered
 #处理CHIPseqrawdata文件
IN1="${path}/rawdata/${bing}/${exp}/${type}/${name}/${name}_1.fq.gz"
IN2="${path}/rawdata/${bing}/${exp}/${type}/${name}/${name}_2.fq.gz"
	r_outp="$dirp/${exp}/${name}_2.fq.gz"
	r_outup="$dirp/${exp}/${name}_R.unpaired.trimed.fq"
	f_outp="$dirp/${exp}/${name}_1.fq.gz"
	f_outup="$dirp/${exp}/${name}_F.unpaired.trimed.fq" 
	fastp \
  -i $IN1 \
  -I $IN2 \
  -o $f_outp \
  -O $r_outp \
  --detect_adapter_for_pe \
  --cut_front 6 --cut_tail 6 \
  --cut_window_size 4 \
  --cut_mean_quality 30 \
  --length_required 36 --trim_poly_g --trim_poly_x \
  --thread 16 \
  --html ${path}/report/fastp/${name}.fastp.html
 fastqc $f_outp -o ${HOME}/CHIPseq/results/fastqc&&fastqc $r_outp -o ${HOME}/CHIPseq/results/fastqc #修剪后文件再次进行质检
# 控制组处理
f_outIp="$dirp/${exp}/${input}_1.fq.gz"
r_outIp="$dirp/${exp}/${input}_2.fq.gz"
INP1="${path}/rawdata/${bing}/${exp}/${type}/${input}/${input}_1.fq.gz"
INP2="${path}/rawdata/${bing}/${exp}/${type}/${input}/${input}_2.fq.gz"
fastp \
  -i $INP1 \
  -I $INP2 \
  -o $f_outIp \
  -O $r_outIp \
  --detect_adapter_for_pe \
  --cut_front 6 --cut_tail 6 \
  --cut_window_size 4 \
  --cut_mean_quality 30 \
  --length_required 36 --trim_poly_g --trim_poly_x \
  --thread 16 \
  --html ${path}/report/fastp/${input}.fastp.html
 #修剪完成，bowtie2或者star对齐排序后picard和sambamba去重
 mkdir -p ${path}/results/bowtie2/${exp}
 mkdir -p ${path}/results/samtools/${exp}
 mkdir -p ${HOME}/CHIPseq/results/sambamba/${exp}
 bowtie2 -p 16 -q --no-contain --very-sensitive-local -X 2000 -x ${HOME}/CHIPseq/referencedata/index/${spe}_index -1 $f_outp -2 $r_outp -S ${path}/results/bowtie2/${exp}/${name}_${y}_aligned.sam 2>> ${path}/results/bowtie2/${exp}/${name}.log
samtools view -h -b -o ${path}/results/samtools/${exp}/${name}.${y}_aligned.bam ${path}/results/bowtie2/${exp}/${name}_${y}_aligned.sam 2>> ${path}/results/samtools/${exp}/${name}.${y}_stb.log
 bowtie2 -p 8 -q --very-sensitive-local --no-contain -X 2000 -x ${HOME}/CHIPseq/referencedata/index/${spe}_index -1 $f_outIp -2 $r_outIp -S ${path}/results/bowtie2/${exp}/${input}_${y}_aligned.sam 2>> ${path}/results/bowtie2/${exp}/${input}.log
 samtools view -h -b -o ${path}/results/samtools/${exp}/${input}.${y}_aligned.bam ${path}/results/bowtie2/${exp}/${input}_${y}_aligned.sam 2>> ${path}/results/sambamba/${exp}/${input}samtobam.log  #转成bam文件，用于macs callpeak
sambamba sort -t 16 -o ${HOME}/CHIPseq/results/sambamba/${exp}/${input}.${y}_sorted.bam ${path}/results/samtools/${exp}/${input}.${y}_aligned.bam 2>> ${path}/results/sambamba/${exp}/${input}sort.log 
##去重##
sambamba sort -t 16 -o ${HOME}/CHIPseq/results/sambamba/${exp}/${name}.${y}_sorted.bam ${path}/results/samtools/${exp}/${name}.${y}_aligned.bam 2>> $alignlog
mkdir -p ${HOME}/CHIPseq/results/picard/${exp}
java -jar ${softbin}/picard.jar AddOrReplaceReadGroups \
  I=${HOME}/CHIPseq/results/sambamba/${exp}/${name}.${y}_sorted.bam \
  O=${HOME}/CHIPseq/results/picard/${exp}/${name}.${y}_rg.bam \
  RGID=${name} \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=${name}
java -jar ${softbin}/picard.jar  MarkDuplicates REMOVE_DUPLICATES=true M=${HOME}/CHIPseq/results/picard/${name}expsample.metrics.txt I=${HOME}/CHIPseq/results/picard/${exp}/${name}.${y}_rg.bam O=${HOME}/CHIPseq/results/picard/${exp}/${name}.${y}_prededu.bam
java -jar ~/software/picard.jar AddOrReplaceReadGroups \
  I=${HOME}/CHIPseq/results/sambamba/${exp}/${input}.${y}_sorted.bam \
  O=${HOME}/CHIPseq/results/picard/${exp}/${input}.${y}_rg.bam \
  RGID=${input} \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=${input}
java -jar ${softbin}/picard.jar MarkDuplicates REMOVE_DUPLICATES=true M=${HOME}/CHIPseq/results/picard/${input}inputsample.metrics.txt I=${HOME}/CHIPseq/results/picard/${exp}/${input}.${y}_rg.bam O=${HOME}/CHIPseq/results/picard/${exp}/${input}.${y}_prededu.bam 2>> ${HOME}/CHIPseq/results/picard/${exp}/RUNNING.log
sambamba view -h -t 32 -f bam -F "[XS] == null and not unmapped  and not duplicate" ${HOME}/CHIPseq/results/picard/${exp}/${name}.${y}_prededu.bam > \
 ${HOME}/CHIPseq/results/sambamba/${exp}/filtered/${name}.${y}_filterd.bam 2>> $alignlog
sambamba view -h -t 32 -f bam -F "[XS] == null and not unmapped  and not duplicate" ${HOME}/CHIPseq/results/picard/${exp}/${input}.${y}_prededu.bam > \
 ${HOME}/CHIPseq/results/sambamba/${exp}/${input}.${y}_filterd.bam 2>> ${HOME}/CHIPseq/results/sambamba/${input}filter.log
#  sambamba index ${HOME}/CHIPseq/results/sambamba/${input}.${y}_filterd.bam
samtools index ${HOME}/CHIPseq/results/sambamba/${exp}/filtered/${name}.${y}_filterd.bam
mkdir -p ${path}/results/MACS3/${bing}/${exp}/narrow/bdg
#窄峰调用
mkdir -p ${HOME}/CHIPseq/results/MACS3/${bing}/${exp}
macs3 callpeak -t ${path}/results/sambamba/${exp}/filtered/${name}.${y}_filterd.bam -c ${HOME}/CHIPseq/results/sambamba/${exp}/${input}.${y}_filterd.bam -g hs -f BAM -p 0.0001 --SPMR -B --nomodel --nolambda --extsize 200 -n ${name} --outdir \
 ${HOME}/CHIPseq/results/MACS3/${bing}/${exp}/narrow #2>> ${HOME}/CHIPseq/results/MACS3/${name}.log #生成实验组文件
 #macs3 callpeak -t ${path}/results/sambamba/${name}input.${y}_filterd.bam -g hs -f BAM -p 0.0001 --SPMR -B --broad --nomodel --nolambda --extsize 200 -n ${name}input --outdir ${HOME}/CHIPseq/results/MACS3 2>> ${HOME}/CHIPseq/results/MACS3/${name}.log #生成控制组宽峰文件，没有必要执行，窄缝调用时有必要
#进行IDR分析，输出的结果是可信度较高的peak，适用于motif分析等等

macs3 bdgcmp -t ${HOME}/CHIPseq/results/MACS3/${bing}/${exp}/narrow/${name}_treat_pileup.bdg -c ${HOME}/CHIPseq/results/MACS3/${bing}/${exp}/narrow/${name}_control_lambda.bdg \
 -o ${HOME}/CHIPseq/results/MACS3/${bing}/${exp}/narrow/bdg/${name}.bdg -m FE
#bedtools intersect -v -a ${HOME}/CHIPseq/results/MACS3/${exp}/${name}_peaks.broadPeak -b ${HOME}/CHIPseq/referencedata/annotation/hg38-blacklist.v2.bed > ${HOME}/CHIPseq/results/MACS3/${exp}/${name}_peaks.broadPeak
#确定信号缩放的比例,因为后续要进行diffbind分析，这里不进行scalefactor基础上归一化
#将bam转为bw文件
mkdir -p ${path}/results/bedgraph/bwdoc/${bing}/${exp}/bam
echo "now we entered bwdoc"
#conda init
#conda activate base
#bamCoverage --bam ${HOME}/CHIPseq/results/sambamba/${name}.${y}_filterd.bam --normalizeUsing RPKM --binSize 50 -p 10 --outFileName ${name}.scaled.bw
#bamCoverage --bam ${HOME}/CHIPseq/results/sambamba/${exp}/filtered/${name}.${y}_filterd.bam --normalizeUsing RPKM --binSize 50 \
#  -p 16 --outFileName ${path}/results/bedgraph/bwdoc/${exp}/bam/${name}.bw
#将bdg文件转为bw文件
#conda activate idr39
mkdir -p  ${path}/results/MACS3/${bing}/${exp}/narrow/filtered
mkdir -p  ${path}/results/bedgraph/sort/${bing}/${exp}
awk '$4 > 0' ${HOME}/CHIPseq/results/MACS3/${bing}/${exp}/narrow/bdg/${name}.bdg > ${path}/results/MACS3/${bing}/${exp}/narrow/filtered/${name}.filtered.bdg
sort -k1,1 -k2,2n ${path}/results/MACS3/${bing}/${exp}/narrow/filtered/${name}.filtered.bdg > ${path}/results/bedgraph/sort/${bing}/${exp}/${name}_sorted.bdg #排序
bedGraphToBigWig ${path}/results/bedgraph/sort/${bing}/${exp}/${name}_sorted.bdg ${path}/referencedata/annotation/${spe}.chrom.sizes ${path}/results/bedgraph/bwdoc/${bing}/${exp}/narrow/${name}_sorted.bw
rm -f ${path}/trimm/pe/*.unpaired*
rm -f ${path}/results/samtools/*.sam
echo "好了！"