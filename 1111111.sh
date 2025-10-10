#!/bin/bash
set -e  # 遇到错误立即退出

outdir="/home/data/luzhuangzhuang/CHIPseq/results/input"
index="${HOME}/CHIPseq/referencedata/index/hg38_index"
mkdir -p "$outdir"

for sample in P8039838 P8020158 PE24022923; do
  for fq in /home/data/share/cell_free/01_DB_YLab/tube_PE/PE_SD/*/*/01.RawData/${sample}*/${sample}*_1.fq.gz; do
    example=$(basename "$(dirname "$fq")")

    echo "========== Processing ${example} =========="

    fq1=$(ls /home/data/share/cell_free/01_DB_YLab/tube_PE/PE_SD/*/*/01.RawData/${example}/${example}_1.fq.gz)
    fq2=$(ls /home/data/share/cell_free/01_DB_YLab/tube_PE/PE_SD/*/*/01.RawData/${example}/${example}_2.fq.gz)

    o=${outdir}/${example}_1.fq
    O=${outdir}/${example}_2.fq

    # --- fastp ---
    fastp \
      -i "$fq1" \
      -I "$fq2" \
      -o "$o" \
      -O "$O" \
      --detect_adapter_for_pe \
      --cut_front 6 --cut_tail 6 \
      --cut_window_size 4 \
      --cut_mean_quality 30 \
      --length_required 36 --trim_poly_g --trim_poly_x \
      --thread 16 \
      --html ${outdir}/${example}.fastp.html

    # --- bowtie2 ---
    bowtie2 -p 16 -q --no-contain --very-sensitive-local -X 2000 \
      -x "$index" -1 "$o" -2 "$O" \
      -S ${outdir}/${example}_aligned.sam

    # --- samtools + sambamba + picard ---
    samtools view -h -b -o ${outdir}/${example}_aligned.bam ${outdir}/${example}_aligned.sam
    sambamba sort -t 16 -o ${outdir}/${example}_sorted.bam ${outdir}/${example}_aligned.bam

    picard AddOrReplaceReadGroups \
      I=${outdir}/${example}_sorted.bam \
      O=${outdir}/${example}_rg.bam \
      RGID=${example} \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM=${example}

    picard MarkDuplicates \
      REMOVE_DUPLICATES=false \
      M=${outdir}/${example}.metrics.txt \
      I=${outdir}/${example}_rg.bam \
      O=${outdir}/${example}_prededu.bam

    sambamba view -h -t 16 -f bam \
      -F "[XS] == null and not unmapped and not duplicate" \
      ${outdir}/${example}_prededu.bam > ${outdir}/${example}_filtered.bam

    # --- QC statistics ---
    samtools flagstat ${outdir}/${example}_aligned.bam > ${outdir}/${example}.flagstat.txt
    samtools view -c ${outdir}/${example}_filtered.bam > ${outdir}/${example}.counts.txt

    echo "✅ Finished: ${example}"
    echo
  done
done