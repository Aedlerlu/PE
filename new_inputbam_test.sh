#！/bin/bash
for i in $(ls /home/data/luzhuangzhuang/CHIPseq/results/input/*filetered.bam);
do
  sample=$(basename $i|cut -d '_' -f1)
  mkdir /home/data/luzhuangzhuang/CHIPseq/results/input/peak/${sample}
  macs2 callpeak --broad -g hs -f BAM -t  /home/data/luzhuangzhuang/CHIPseq/results/input/${sample}_filtered.bam -c /home/data/luzhuangzhuang/CHIPseq/results/input/controlbam/merge_ext.bam -n $sample \
--outdir /home/data/luzhuangzhuang/CHIPseq/results/input/peak/${sample}
done
