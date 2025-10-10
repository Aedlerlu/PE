#!/bin/bash

for sample in P8039838 P8020158 PE24022923; do
  for fq in /home/data/share/cell_free/01_DB_YLab/tube_PE/PE_SD/*/*/01.RawData/${sample}*/${sample}*_1.fq.gz; do
    example=$(basename "$(dirname "$fq")")
    outdir="/home/data/luzhuangzhuang/CHIPseq/results/input"
    macs2 callpeak --broad -t ${outdir}/${example}_filtered.bam -c ${outdir}/merge_extract.bam -n ${example} --outdir ${outdir} -f BAM -g hs
  done
done
