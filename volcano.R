# ===== 1. 加载包 =====
library(DESeq2)
library(readxl)
library(dplyr)

# ===== 2. 读取数据 =====
# 修改为你的文件路径
data <- readxl::read_excel("~/电刺激成骨RNA-seq/day8.xlsx")
data <- data[!is.na(data$gene_name), ]
# 暂不去除重复，保留 gene_name 用于后续合并
expr_counts <- data[, -1]
rownames(expr_counts) <- make.unique(data$gene_name)  # 自动生成唯一 rownames（如 gene.1, gene.2）

# ===== 3. 构建样本信息（control 拼写修正）=====
colData <- data.frame(
  row.names = colnames(expr_counts),
  condition = factor(c("Control", "Control", "Control", "ES", "ES", "ES"),
                     levels = c("Control", "ES"))
)

# ===== 4. 构建 DESeq2 对象并分析 =====
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(expr_counts)),
                              colData = colData,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "ES", "Control"))

# ===== 5. 整理结果（合并原始 gene_name）=====
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)  # rowname 是 unique id（可能是 gene, gene.1 等）
res_df$gene_name <- data$gene_name  # 保留原始重复的 gene_name 顺序

# ===== 6. 保留最显著的重复 gene（p 值最小）=====
final_df <- res_df %>%
  dplyr::select(gene_name, log2FoldChange, pvalue) %>%
  dplyr::filter(!is.na(log2FoldChange), !is.na(pvalue), !is.na(gene_name)) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice_min(order_by = pvalue, with_ties = FALSE) %>%
  dplyr::ungroup()

# ===== 7. 导出结果 =====
write.csv(final_df, "DEG_day8_ES_vs_Control.csv", row.names = FALSE)
