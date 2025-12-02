#先shell把染色体分窗口
cd /home/data/luzhuangzhuang/CHIPseq/
awk '$1 ~ /^chr/ && $1 != "chrM"' "/home/data/luzhuangzhuang/CHIPseq/referencedata/annotation/hg38.chrom.sizes" > hg38.noM.sizes
bedtools makewindows -g hg38.noM.sizes -w 5000 > genome_5kb_bins.bed
########
awk '$5=="acen" {print $1"\t"$2"\t"$3}' cytoBand.txt \
> centromere.bed
###
#######bw列表获取
setwd("/home/data/luzhuangzhuang/CHIPseq/")
seqlengths <- read.table("hg38.noM.sizes",stringsAsFactors = FALSE)
seqlengths <- setNames(seqlengths$V2,seqlengths$V1)
bins <- tileGenome(seqlengths, tilewidth=5000, cut.last.tile.in.chrom=TRUE)
blacklist <- rtracklayer::import("/home/data/luzhuangzhuang/CHIPseq/ENCFF356LFX.bed")
centromere <- rtracklayer::import("centromere.bed")
bad_regions <- c(blacklist, centromere)
hits <- findOverlaps(bins, bad_regions)
bins_filtered <- bins[-unique(queryHits(hits))]
fk_bins_df <- as.data.frame(bins_filtered)[, c("seqnames","start","end")]
write.table(fk_bins_df, "5000_bins_filtered.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
######
ALL_K4M_bw_files <- list.files(pattern = "K4M_FE.*.bw$",
                                path = "/home/data/share/cell_free/05_results_YLab/01_PE_SD_results_qualified/",full.names = TRUE,recursive = TRUE)
writeLines(ALL_K4M_bw_files,"all_K4M_bw_list.txt")
save.image("5k_bed_treat.RData")
######
multiBigwigSummary BED-file \
--bwfiles $(cat all_K4M_bw_list.txt) \
--BED 5000_bins_filtered.bed \
--outRawCounts fk_bins_rawcounts.tsv \
--outFileName fk_bins_summary.npz \
--smartLabels \
--numberOfProcessors 64
##############
all_k4m_bw <- readLines("all_K4M_bw_list.txt")
sample_list <- list()
for (i in seq_along(all_k4m_bw)) { 
  group <- ifelse(grepl(pattern = "C",all_k4m_bw[i]),"HC","PE")
  samplename <- sub("_FE\\.bw","",basename(all_k4m_bw[i]))
  marker <-"H3K4M"
  df <- data.frame(samplename = samplename, group = group,marker = marker)
sample_list[[i]] <- df
}
sample_infor <- do.call(rbind,sample_list)
rownames(sample_infor)<- sample_infor$samplename
write.table(sample_infor,sep = "\t","sample_infor.tsv")
bins_file <- "bins_rawcounts.tsv"   # 输出自 multiBigwigSummary
#sample_info <- "sample_info.tsv"    # 两列: sample_name\tgroup\tmark (mark: K4M/K27A)
# sample_name 必须与 bins_rawcounts.tsv 的列名一致
out_prefix <- "bins5kb"
# ----- 读入 -----
bins <- fread(bins_file, stringsAsFactors = FALSE)
sum(is.na(bins))
sum(is.na(bins$chr))
sum(is.na(bins$start))
sum(is.na(bins$end))
bins[is.na(bins)] <- 0
sum(is.na(bins))
###########################
colnames(bins) <- gsub("'", "", colnames(bins))  # 先删除所有单引号
colnames(bins) <- gsub("^\\#","",colnames(bins))
colnames(bins) <- sub("_FE$", "", colnames(bins))  # 再删除模式
save.image("read_quene.RData")
# 第一 3 列是 chr start end
bins[, bin_id := paste0(chr, ":", start, "-", end)]
setnames(bins, c("chr","start","end"), c("chr","start","end"), skip_absent=TRUE) # robust
# 把列名规范化（multiBigwigSummary 可能在第一行写了 sample labels）
# 假设第4列开始是样本
count_mat <- as.matrix(bins[, -(1:3), with=FALSE])
rownames(count_mat) <- bins$bin_id

# ----- 读样本信息 -----
si <- sample_infor # sample_name group mark
# 检查样本顺序与 count_mat 列顺序一致
colnames(count_mat) <- colnames(count_mat)
if (!all(si$samplename %in% colnames(count_mat))) stop("sample names mismatch. Check sample_info and bins file.")
################
missing_samples <- setdiff(si$samplename, colnames(count_mat))
if (length(missing_samples) > 0) {
  warning("以下样本在count_mat中缺失: ", paste(missing_samples, collapse = ", "))
}
####################
# 调整列顺序以 match sample_info
count_mat <- count_mat[, si$samplename, drop=FALSE]

# ----- 初步过滤：去掉在多数样本中为 0 的 bins（减少计算量） -----
#min_frac_nonzero <- 0.1   # 在至少 10% 的样本中有非零
#keep <- rowSums(count_mat > 0) >= (min_frac_nonzero * ncol(count_mat))
#message("原始 bin 数: ", nrow(count_mat), " -> 过滤后: ", sum(keep))
#count_mat_f <- count_mat[keep, ]
sum(is.na(count_mat))
save.image("read_quene.RData")
count_mat_n <- matrix(as.numeric(count_mat), 
                            nrow = nrow(count_mat), 
                            ncol = ncol(count_mat))
rownames(count_mat_n) <- rownames(count_mat)
colnames(count_mat_n) <- colnames(count_mat)
count_mat <- count_mat_n
# ----- 建 DGEList 并归一化（TMM） -----
grp <- factor(si$group)   # e.g. "control" / "disease"
dge <- DGEList(counts = count_mat, group = grp)
dge <- calcNormFactors(dge, method="TMM")
save.image("K4M_dge.RData")
load("K4M_dge.RData")
# voom 变换
head(dge,n=2)
design <- model.matrix(~ 0 + grp)  # 使用无截距矩阵，方便对比
colnames(design) <- levels(grp)
v <- voom(dge, design = design, plot = TRUE)

# ----- 差异分析（limma） -----
# 例：如果只有 control vs disease 两组
if (nlevels(grp) == 2) {
  contrast.matrix <- makeContrasts(Diff = levels(grp)[2] - levels(grp)[1], levels=design)
  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  topbins <- topTable(fit2, coef = "Diff", number = Inf, sort.by = "P")
  fwrite(as.data.table(topbins), paste0(out_prefix, "_limma_topbins.tsv"), sep="\t")
} else {
  # 更多组的情况：按需求设 contrasts
  message("多个组，请自定义 contrast")
}


if (nlevels(grp) == 2) {
  # 构建对比字符串
  contrast_str <- paste(levels(grp)[2], levels(grp)[1], sep = "-")
  contrast.matrix <- makeContrasts(contrasts = contrast_str, levels = design)
  
  fit <- lmFit(v, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  topbins <- topTable(fit2, coef = contrast_str, number = Inf, sort.by = "P")
  fwrite(as.data.table(topbins), paste0(out_prefix, "_limma_topbins.tsv"), sep="\t")
}

