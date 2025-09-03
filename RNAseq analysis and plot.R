library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)

# 读取表达矩阵
counts_raw <- read.csv("biaodajuzhen.csv", stringsAsFactors = FALSE, check.names = FALSE)

# 提取行名（Ensembl ID），并设为唯一行名
ens_ids <- as.character(counts_raw[, 1])
counts <- counts_raw[, -1]  # 去除第一列（Ensembl ID）
rownames(counts) <- ens_ids

# 保存 gene_name 映射（用于后续合并）
gene_map <- data.frame(
  ens_id = ens_ids,
  gene_symbol = counts_raw$gene_name,
  stringsAsFactors = FALSE
)
# 去除 gene_name 列（避免出现在表达矩阵中）
counts$gene_name <- NULL

# 读取样本信息
coldata <- read.csv("sample_metadata.csv", row.names = 1)
stopifnot(all(colnames(counts) == rownames(coldata)))

# 可选：构建全局 dds 对象用于 VST 或 PCA 可视化等
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Condition)
dds <- DESeq(dds)

library(DESeq2)
library(apeglm)

# 创建一个函数，用于每一天的分析
run_deseq_for_day <- function(day, group1, group2, coef_name = NULL, use_apeglm = FALSE) {
  # 子集 coldata 和 counts
  col_sub <- coldata[coldata$Day == day & coldata$Condition %in% c(group1, group2), ]
  count_sub <- counts[, rownames(col_sub)]
  
  # 重新构建对象（以 Day 为单位）
  col_sub$Condition <- factor(col_sub$Condition, levels = c(group2, group1))  # 设置参考组为 Control
  dds_sub <- DESeqDataSetFromMatrix(countData = count_sub, colData = col_sub, design = ~ Condition)
  dds_sub <- DESeq(dds_sub)
  
  # 差异分析
  if (!is.null(coef_name) & use_apeglm) {
    res <- lfcShrink(dds_sub, coef = coef_name, type = "apeglm")
  } else {
    res <- lfcShrink(dds_sub, contrast = c("Condition", group1, group2), type = "normal")
  }
  
  return(res)
}

# Day4: ES vs Control（可以用 apeglm，因为是模型 coef）
res_day4 <- run_deseq_for_day(day = "Day4", group1 = "ES", group2 = "Control", 
                              coef_name = "Condition_ES_vs_Control", use_apeglm = TRUE)

# Day7: ES vs Control（使用 contrast，不能用 apeglm）
res_day7 <- run_deseq_for_day(day = "Day7", group1 = "ES", group2 = "Control")

# Day8: PostES vs Control
res_day8 <- run_deseq_for_day(day = "Day8", group1 = "PostES", group2 = "Control")

save_deg <- function(res, out_file) {
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$padj) & res_df$padj < 0.1 & abs(res_df$log2FoldChange) > 1, ]
  write.csv(res_df, out_file, row.names = FALSE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "clusterProfiler"))

# 加载注释和富集包
library(org.Hs.eg.db)
library(clusterProfiler)

# 函数：保存 DEG 结果，带基因名映射
save_deg_with_gene_name <- function(res, out_file) {
  res_df <- as.data.frame(res)
  res_df$ensembl_id <- sub("\\..*", "", rownames(res_df))
  
  # 🔍 过滤显著 DEG（padj < 0.05 且 log2FC > 1）
  res_df <- res_df[!is.na(res_df$padj) & res_df$padj < 0.1 & abs(res_df$log2FoldChange) > 1, ]
  
  # ✨ 使用 clusterProfiler::bitr 进行基因名注释
  gene_map <- bitr(res_df$ensembl_id,
                   fromType = "ENSEMBL",
                   toType = "SYMBOL",
                   OrgDb = org.Hs.eg.db)
  
  # 合并基因名
  res_merged <- merge(res_df, gene_map, by.x = "ensembl_id", by.y = "ENSEMBL", all.x = TRUE)
  
  # 📌 调整列顺序：gene_name, ensembl_id, log2FC, padj, ...
  res_final <- res_merged[, c("SYMBOL", "ensembl_id", setdiff(colnames(res_merged), c("SYMBOL", "ensembl_id")))]
  colnames(res_final)[1] <- "gene_name"  # 重命名 SYMBOL 为 gene_name
  
  # 保存为 CSV
  write.csv(res_final, file = out_file, row.names = FALSE)
}

save_deg_with_gene_name(res_day4, "DEG_Day4_ES_vs_Control_annotated.csv")
save_deg_with_gene_name(res_day7, "DEG_Day7_ES_vs_Control_annotated.csv")
save_deg_with_gene_name(res_day8, "DEG_Day8_PostES_vs_Control_annotated.csv")

# 加载必要库制作热图
library(DESeq2)
library(org.Hs.eg.db)
library(pheatmap)
library(clusterProfiler)

# 1. 创建 VST 表达矩阵
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

# 2. 计算每个基因的方差并选择 top 100
gene_var <- apply(vst_mat, 1, var)
top100_ids <- names(sort(gene_var, decreasing = TRUE))[1:100]

# 3. 去除 Ensembl ID 的版本号（.10）
ens_ids_clean <- sub("\\..*", "", top100_ids)

# 4. 映射为基因名 SYMBOL
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ens_ids_clean,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# 5. 提取 top 100 表达矩阵并替换行为 gene name
heat_mat <- vst_mat[top100_ids, ]
rownames(heat_mat) <- gene_symbols[ens_ids_clean]

# 6. 去除无法注释的行（NA）
heat_mat <- heat_mat[!is.na(rownames(heat_mat)), ]

# 7. 标准化表达（行 z-score）
heat_mat <- t(scale(t(heat_mat)))

# 8. 添加样本注释（按列）
annotation_col <- coldata[, c("Condition", "Day")]

# 9. 绘制热图
pheatmap(heat_mat,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 4,
         fontsize_col = 10,
         main = "Top 100 Most Variable Genes (Gene Symbols)")

# 提取 Day4 显著 DEG
deg4 <- as.data.frame(res_day4)
deg4$ensembl <- sub("\\..*", "", rownames(deg4))
deg4$entrez <- mapIds(org.Hs.eg.db, keys=deg4$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# 上调基因（log2FC > 1）
up_genes <- deg4$entrez[deg4$log2FoldChange > 1 & deg4$padj < 0.1]
down_genes <- deg4$entrez[deg4$log2FoldChange < -1 & deg4$padj < 0.1]

# GO 富集（上调）
ego_up <- enrichGO(gene = up_genes,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pvalueCutoff = 0.05)

dotplot(ego_up, showCategory = 20) +
  ggtitle("GO BP - Upregulated (Day4)") +
  theme(axis.text.y = element_text(size = 10),      # y轴路径名
        axis.text.x = element_text(size = 12),      # x轴数值
        plot.title = element_text(size = 14, face = "bold"))

# KEGG 富集
kk_up <- enrichKEGG(gene = up_genes,
                    organism = "hsa",
                    pvalueCutoff = 0.05)

dotplot(kk_up, showCategory = 20) +
  ggtitle("KEGG - Upregulated (Day4)") +
  theme(axis.text.y = element_text(size = 10),      # y轴路径名
        axis.text.x = element_text(size = 12),      # x轴数值
        plot.title = element_text(size = 14, face = "bold"))

# 提取 Day7 显著 DEG
deg7 <- as.data.frame(res_day7)
deg7$ensembl <- sub("\\..*", "", rownames(deg7))
deg7$entrez <- mapIds(org.Hs.eg.db, keys=deg7$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# 上调基因
up_genes7 <- deg7$entrez[deg7$log2FoldChange > 1 & deg7$padj < 0.1]

# GO 富集（上调）
ego_up7 <- enrichGO(gene = up_genes7,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.05)

dotplot(ego_up7, showCategory = 20) +
  ggtitle("GO BP - Upregulated (Day7)") +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))

# KEGG 富集（上调）
kk_up7 <- enrichKEGG(gene = up_genes7,
                     organism = "hsa",
                     pvalueCutoff = 0.05)

dotplot(kk_up7, showCategory = 20) +
  ggtitle("KEGG - Upregulated (Day7)") +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))

# 提取 Day8 显著 DEG
deg8 <- as.data.frame(res_day8)
deg8$ensembl <- sub("\\..*", "", rownames(deg8))
deg8$entrez <- mapIds(org.Hs.eg.db, keys=deg8$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# 上调基因
up_genes8 <- deg8$entrez[deg8$log2FoldChange > 1 & deg8$padj < 0.1]

# GO 富集（上调）
ego_up8 <- enrichGO(gene = up_genes8,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.05)

dotplot(ego_up8, showCategory = 20) +
  ggtitle("GO BP - Upregulated (Day8)") +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))

# KEGG 富集（上调）
kk_up8 <- enrichKEGG(gene = up_genes8,
                     organism = "hsa",
                     pvalueCutoff = 0.05)

dotplot(kk_up8, showCategory = 20) +
  ggtitle("KEGG - Upregulated (Day8)") +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))

if (!requireNamespace("fgsea", quietly = TRUE)) install.packages("fgsea")
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")

library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)

# 获取人类的 GO Biological Process (BP) gene sets
msig_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

# 提取成骨分化相关的 gene sets（关键词包含 bone/osteoblast/ossification）
osteo_sets <- msig_go[grep("OSTEOBLAST|OSSIFICATION|BONE", msig_go$gs_name, ignore.case = TRUE), ]

# 转换为 list 格式，fgsea 要求
osteo_pathways <- split(osteo_sets$gene_symbol, osteo_sets$gs_name)

# 加入 ENTREZ ID
res_day4$ens <- sub("\\..*", "", rownames(res_day4))
res_day4$symbol <- mapIds(org.Hs.eg.db, keys = res_day4$ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# 去除 NA 值
res_clean <- res_day4[!is.na(res_day4$log2FoldChange) & !is.na(res_day4$symbol), ]

# 创建命名排序向量（geneList）
geneList <- res_clean$log2FoldChange
names(geneList) <- res_clean$symbol
geneList <- sort(geneList, decreasing = TRUE)

library(dplyr)

# 确保 res_clean 是 data.frame 类型
res_clean_df <- as.data.frame(res_clean)

# 保留每个 gene symbol 最大的 log2FC（去重）
geneList_df <- res_clean_df %>%
  group_by(symbol) %>%
  summarize(log2FC = max(log2FoldChange, na.rm = TRUE)) %>%
  arrange(desc(log2FC))

# 创建 named vector 用于 fgsea
geneList <- geneList_df$log2FC
names(geneList) <- geneList_df$symbol

# 跑 fgsea
fgsea_res <- fgsea(pathways = osteo_pathways,
                   stats = geneList,
                   nperm = 10000)  # 建议至少 10000 次

# 查看显著结果
fgsea_sig <- fgsea_res[fgsea_res$padj < 0.1, ]

# 画 top 成骨通路的富集图
top_pathway <- fgsea_sig$pathway[1]  # 选第一个
plotEnrichment(osteo_pathways[[top_pathway]], geneList) +
  ggtitle(paste("GSEA Enrichment(Day4):", top_pathway))
fgsea_sig[fgsea_sig$pathway == "GOBP_NEGATIVE_REGULATION_OF_BONE_REMODELING", ]

library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(dplyr)

# 获取 GO Biological Process (BP) gene sets
msig_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

# 提取成骨相关通路（含关键词）
osteo_sets <- msig_go[grep("OSTEOBLAST|OSSIFICATION|BONE", msig_go$gs_name, ignore.case = TRUE), ]

# 转为 list 格式供 fgsea 使用
osteo_pathways <- split(osteo_sets$gene_symbol, osteo_sets$gs_name)

# 对 Day7 结果加入基因 symbol 注释
res_day7$ens <- sub("\\..*", "", rownames(res_day7))
res_day7$symbol <- mapIds(org.Hs.eg.db,
                          keys = res_day7$ens,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# 清洗缺失值
res_clean_day7 <- as.data.frame(res_day7[!is.na(res_day7$log2FoldChange) & !is.na(res_day7$symbol), ])

# 保留每个 symbol 最大 log2FC，避免重复
geneList_df_day7 <- res_clean_day7 %>%
  group_by(symbol) %>%
  summarize(log2FC = max(log2FoldChange, na.rm = TRUE)) %>%
  arrange(desc(log2FC))

# 创建排序向量
geneList_day7 <- geneList_df_day7$log2FC
names(geneList_day7) <- geneList_df_day7$symbol

# 跑 fgsea
fgsea_res_day7 <- fgsea(pathways = osteo_pathways,
                        stats = geneList_day7,
                        nperm = 10000)

# 查看显著通路
fgsea_sig_day7 <- fgsea_res_day7[fgsea_res_day7$padj < 0.1, ]
print(fgsea_sig_day7)

# 画 enrichment 图（第一个显著的通路）
if (nrow(fgsea_sig_day7) > 0) {
  top_pathway_day7 <- fgsea_sig_day7$pathway[1]
  plotEnrichment(osteo_pathways[[top_pathway_day7]], geneList_day7) +
    ggtitle(paste("GSEA Enrichment (Day7):", top_pathway_day7))
}
fgsea_sig_day7[fgsea_sig_day7$pathway == "GOBP_BONE_GROWTH", ]

# 加载必要的包
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(dplyr)

# 获取 GO BP 成骨相关 gene sets（关键词含 bone/osteoblast/ossification）
msig_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
osteo_sets <- msig_go[grep("OSTEOBLAST|OSSIFICATION|BONE", msig_go$gs_name, ignore.case = TRUE), ]
osteo_pathways <- split(osteo_sets$gene_symbol, osteo_sets$gs_name)

# 准备 res_day8 数据
res_day8$ens <- sub("\\..*", "", rownames(res_day8))
res_day8$symbol <- mapIds(org.Hs.eg.db,
                          keys = res_day8$ens,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# 清洗缺失值
res_clean8 <- res_day8[!is.na(res_day8$log2FoldChange) & !is.na(res_day8$symbol), ]
res_clean8_df <- as.data.frame(res_clean8)

# 每个 gene symbol 只保留最大 log2FC
geneList_df8 <- res_clean8_df %>%
  group_by(symbol) %>%
  summarize(log2FC = max(log2FoldChange, na.rm = TRUE)) %>%
  arrange(desc(log2FC))

# 创建命名排序向量供 fgsea 使用
geneList8 <- geneList_df8$log2FC
names(geneList8) <- geneList_df8$symbol

# 跑 fgsea
fgsea_res8 <- fgsea(pathways = osteo_pathways,
                    stats = geneList8,
                    nperm = 10000)

# 提取显著通路
fgsea_sig8 <- fgsea_res8[fgsea_res8$padj < 0.1, ]

# 可视化 top1 成骨通路
if (nrow(fgsea_sig8) > 0) {
  top_pathway8 <- fgsea_sig8$pathway[1]
  plotEnrichment(osteo_pathways[[top_pathway8]], geneList8) +
    ggtitle(paste("GSEA Enrichment (Day8):", top_pathway8))
} else {
  message("Day8 中未发现 padj < 0.05 的成骨通路")
}

# PCA 图
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("Condition", "Day"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = Condition, shape = Day)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of all samples")
