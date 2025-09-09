#!/bin/bash
 sleep $((RANDOM % 60 + 60))
spe="${1}"
exp="${2}"
bing="${3}"
type="${4}"
path=${HOME}/CHIPseq
readname=${path}/name.txt
find ${path}/rawdata/${bing}/${exp} -type f \( \
  -name "*.fq" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fastq.gz" \
\) |\
xargs -I {} basename {} | \
cut -d'_' -f1|sort -u > "${readname}" #这里我们取的是双端测序_1之前的部分做输入变量，如果测序文件格式是这样的****_**_1.fq.gz,就得考虑更换提取方式
 #以字符长度判断控制组与实验组基于控制组文件命名更短的逻辑
for i in $(cat ${path}/name.txt);do
 if [[ "$i" == *INP ]];then
echo "$i" > input.txt
 else 
 echo "$i" > expname.txt
 fi
 for name in $(cat expname.txt);do
 input=$(cat input.txt)
  echo "Cleaning intermediate files for ${name}..."
  rm -f ${path}/results/fastp/${exp}/${name}*
  rm -f ${path}/results/samtools/${exp}/${name}*
  rm -f ${path}/results/sambamba/${exp}/${name}*
  rm -f ${path}/results/picard/${exp}/${name}*
  rm -f ${path}/results/picard/${exp}/${name}*
  rm -f ${path}/results/MACS3/${exp}/${name}*
  rm -f ${path}/results/MACS3/${exp}/${name}*
  rm -f ${path}/results/bedgraph/sort/${exp}/${name}*
  rm -f ${path}/results/bedgraph/bedoc/${exp}/${name}*
  rm -f ${path}/results/bedgraph/bedoc/${exp}/${name}*
  rm -f ${path}/results/bedgraph/bwdoc/${exp}/${name}*
  rm -f ${path}/results/bowtie2/${exp}/${name}*
  rm -f ${path}/trimm/pe/${exp}/${name}*
  rm -r ${path}/*.log
  rm -r ${path}/results/*.log
  echo "Clean complete."
  nohup bash ${HOME}/scripts/noinputchip.sh "${spe}"  "${name}" "${exp}" "${type}" "${bing}" y "${input}" 2>> ${path}/run.log
done
  

