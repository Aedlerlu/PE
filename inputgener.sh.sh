mkdir /home/data/luzhuangzhuang/CHIPseq/results/input
type=$1
for sample in P8039838 P8020158 PE24022923; do
ls /home/data/share/cell_free/01_DB_YLab/tube_PE/PE_SD/*/*/01.RawData/${sample}*/${sample}*_1.fq.gz

i=$(ls /home/data/share/cell_free/01_DB_YLab/tube_PE/PE_SD/*/*/01.RawData/${sample}${type}/${sample}${type}_1.fq.gz)
I=$(ls /home/data/share/cell_free/01_DB_YLab/tube_PE/PE_SD/*/*/01.RawData/${sample}${type}/${sample}${type}_2.fq.gz)
o=/home/data/luzhuangzhuang/CHIPseq/results/input/${sample}${type}_1.fq
O=/home/data/luzhuangzhuang/CHIPseq/results/input/${sample}${type}_2.fq
fastp \
  -i $i \
  -I $I \
  -o $o \
  -O $O \
  --detect_adapter_for_pe \
  --cut_front 6 --cut_tail 6 \
  --cut_window_size 4 \
  --cut_mean_quality 30 \
  --length_required 36 --trim_poly_g --trim_poly_x \
  --thread 16 \
  --html /home/data/luzhuangzhuang/CHIPseq/results/input/${sample}.fastp.html
bowtie2 -p 64 -q --no-contain --very-sensitive-local -X 2000 -x ${HOME}/CHIPseq/referencedata/index/hg38_index -1 $o -2 $O -S /home/data/luzhuangzhuang/CHIPseq/results/input/${sample}${type}_aligned.sam 
samtools view -h -b -o /home/data/luzhuangzhuang/CHIPseq/results/input/${sample}${type}_aligned.bam /home/data/luzhuangzhuang/CHIPseq/results/input/${sample}${type}_aligned.sam 
sambamba sort -t 16 -o ${HOME}/CHIPseq/results/input/${sample}${type}_sorted.bam /home/data/luzhuangzhuang/CHIPseq/results/input/${sample}${type}_aligned.bam
picard AddOrReplaceReadGroups \
  I=${HOME}/CHIPseq/results/input/${sample}${type}_sorted.bam \
  O=${HOME}/CHIPseq/results/input/${sample}${type}_rg.bam \
  RGID=${sample}${type} \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=${sample}${type}
picard MarkDuplicates REMOVE_DUPLICATES=false M=${HOME}/CHIPseq/results/input/${sample}${type}.metrics.txt I=${HOME}/CHIPseq/results/input/${sample}${type}_rg.bam O=${HOME}/CHIPseq/results/input/${sample}${type}_prededu.bam
sambamba view -h -t 32 -f bam -F "[XS] == null and not unmapped  and not duplicate" ${HOME}/CHIPseq/results/input/${sample}${type}_prededu.bam > \
 ${HOME}/CHIPseq/results/input/${sample}${type}_filterd.bam
 sampletools flagstat /home/data/luzhuangzhuang/CHIPseq/results/input/${sample}${type}_aligned.bam > ${sample}${type}.txt
 samtools view -c ${HOME}/CHIPseq/results/input/${sample}${type}_filterd.bam > ${sample}${type}counts.txt
 done
