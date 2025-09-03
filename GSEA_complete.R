rm(list = ls())

library(dplyr)
library(org.Hs.eg.db)#物种注释包(Homo sapiens)
library(clusterProfiler)#富集分析
library(readxl)
library(writexl)
library(tidyverse)
library(DESeq2)
# library(org.Mm.eg.db)#物种注释包(小鼠)
library(enrichplot)
# BiocManager::install("ReactomePA")
# install.package("ReactomePA")
# library(ReactomePA)
setwd("C:\\Users\\Administrator\\Desktop\\转录组下游")

# 加载所需的包
library("pacman") 
p_load(
  colorspace, stringi, data.table, ggvenn, DOSE,
  enrichplot, fgsea, msigdbr, cowplot, GseaVis
)
####mRNA4####

DESeq2 <- read_excel("./1.GSEA/mRNA4.xlsx")  ##从4.DEGs/DESeq2.csv这个文件里复制基因列和log2FoldChange列就可以
length(rownames(DESeq2)) 

# 选择差异基因数据中的"SYMBOL"和"logFC"两列
diff_genes <- DESeq2[, c("ID", "log2FoldChange")]
names(diff_genes)[names(diff_genes) == "ID"] <- "SYMBOL"
#添加entrez ID列：
##symbol转entrez ID：
entrez <- bitr(diff_genes$SYMBOL,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

df_all <- merge(entrez, diff_genes,by = "SYMBOL", all = FALSE)

# Glu=read_excel("C:/Users/Administrator/Desktop/diff.xlsx")
# Glu=merge(Glu,df_all,by="ENTREZID")
# write_xlsx(Glu,"C:/Users/Administrator/Desktop/Glu_TATB.xlsx")

# 排序合并后的数据框
df_all_sort <- df_all[order(df_all$log2FoldChange, decreasing = TRUE), ]

# 提取foldchange
gene_fc <- df_all_sort$log2FoldChange

# 将gene_fc对象的命名改为df_all_sort数据框中ENTREZID列的值
names(gene_fc) <- df_all_sort$ENTREZID

####基于KEGG基因集的GSEA富集
KEGG_ges <- gseKEGG(
  geneList = gene_fc,
  organism = "hsa",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#重新转化为symbol基因
KEGG_ges <- setReadable(KEGG_ges,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID")
# 将结果转换为数据框
KEGG_table <- as.data.frame(KEGG_ges)
write_xlsx(KEGG_table,"./1.GSEA/KEGG_table4.xlsx")

pdf("./1.GSEA/MAPK_4.pdf",width =6,height =5)
# 使用gseaNb函数绘制GSEA图
gseaNb(
  object = KEGG_ges,          # GSEA结果的对象，这里是KEGG_ges
  geneSetID = "hsa04010",    # 想要绘制的基因集ID，这里是KEGG pathway的ID
  subPlot = 3,               # 绘制的子图的选择，这里是第三个子图
  # kegg = T,                  #这样基因ID就是symbol
  # addGene = TRUE,            # 是否在曲线图上添加基因名，这里是添加
  lineSize = 0.8,            # 曲线的线宽度
  rmSegment = FALSE,         # 是否移除曲线图上的段，这里是不移除
  segCol = "red",            # 曲线图上的段的颜色
  pvalX = 0.55, pvalY = 0.8,  # P值标签的位置坐标
  pHjust = 0,                # P值标签的水平对齐方式
  # geneCol = "black",         # 基因名的颜色
  curveCol = c("#76BA99", "#EB4747", "#996699"),  # 曲线的颜色
  addPval = TRUE,            # 是否添加P值和NES值，确保它们在GSEA结果中已经计算
  pCol = "grey30",           # P值标签的颜色
  # markTopgene = TRUE,       # 是否突出显示前n个基因
  # topGeneN = 5               # 要突出显示的前n个基因的数量
  # ... 其他需要的参数
)

dev.off()

####mRNA7####

DESeq2 <- read_excel("./1.GSEA/mRNA7.xlsx")  ##从4.DEGs/DESeq2.csv这个文件里复制基因列和log2FoldChange列就可以
length(rownames(DESeq2)) 

# 选择差异基因数据中的"SYMBOL"和"logFC"两列
diff_genes <- DESeq2[, c("ID", "log2FoldChange")]
names(diff_genes)[names(diff_genes) == "ID"] <- "SYMBOL"
#添加entrez ID列：
##symbol转entrez ID：
entrez <- bitr(diff_genes$SYMBOL,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

df_all <- merge(entrez, diff_genes,by = "SYMBOL", all = FALSE)

# Glu=read_excel("C:/Users/Administrator/Desktop/diff.xlsx")
# Glu=merge(Glu,df_all,by="ENTREZID")
# write_xlsx(Glu,"C:/Users/Administrator/Desktop/Glu_TATB.xlsx")

# 排序合并后的数据框
df_all_sort <- df_all[order(df_all$log2FoldChange, decreasing = TRUE), ]

# 提取foldchange
gene_fc <- df_all_sort$log2FoldChange

# 将gene_fc对象的命名改为df_all_sort数据框中ENTREZID列的值
names(gene_fc) <- df_all_sort$ENTREZID

####基于KEGG基因集的GSEA富集
KEGG_ges <- gseKEGG(
  geneList = gene_fc,
  organism = "hsa",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#重新转化为symbol基因
KEGG_ges <- setReadable(KEGG_ges,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID")
# 将结果转换为数据框
KEGG_table <- as.data.frame(KEGG_ges)
write_xlsx(KEGG_table,"./1.GSEA/KEGG_table7.xlsx")

pdf("./1.GSEA/Calcium_7.pdf",width =6,height =5)
# 使用gseaNb函数绘制GSEA图
gseaNb(
  object = KEGG_ges,          # GSEA结果的对象，这里是KEGG_ges
  geneSetID = "hsa04020",    # 想要绘制的基因集ID，这里是KEGG pathway的ID
  subPlot = 3,               # 绘制的子图的选择，这里是第三个子图
  # kegg = T,                  #这样基因ID就是symbol
  # addGene = TRUE,            # 是否在曲线图上添加基因名，这里是添加
  lineSize = 0.8,            # 曲线的线宽度
  rmSegment = FALSE,         # 是否移除曲线图上的段，这里是不移除
  segCol = "red",            # 曲线图上的段的颜色
  pvalX = 0.55, pvalY = 0.8,  # P值标签的位置坐标
  pHjust = 0,                # P值标签的水平对齐方式
  # geneCol = "black",         # 基因名的颜色
  curveCol = c("#76BA99", "#EB4747", "#996699"),  # 曲线的颜色
  addPval = TRUE,            # 是否添加P值和NES值，确保它们在GSEA结果中已经计算
  pCol = "grey30",           # P值标签的颜色
  # markTopgene = TRUE,       # 是否突出显示前n个基因
  # topGeneN = 5               # 要突出显示的前n个基因的数量
  # ... 其他需要的参数
)

dev.off()



####mRNA8####

DESeq2 <- read_excel("./1.GSEA/mRNA8.xlsx")  ##从4.DEGs/DESeq2.csv这个文件里复制基因列和log2FoldChange列就可以
length(rownames(DESeq2)) 

# 选择差异基因数据中的"SYMBOL"和"logFC"两列
diff_genes <- DESeq2[, c("ID", "log2FoldChange")]
names(diff_genes)[names(diff_genes) == "ID"] <- "SYMBOL"
#添加entrez ID列：
##symbol转entrez ID：
entrez <- bitr(diff_genes$SYMBOL,
               fromType = "SYMBOL",#现有的ID类型
               toType = "ENTREZID",#需转换的ID类型
               OrgDb = "org.Hs.eg.db")
head(entrez)

df_all <- merge(entrez, diff_genes,by = "SYMBOL", all = FALSE)

# Glu=read_excel("C:/Users/Administrator/Desktop/diff.xlsx")
# Glu=merge(Glu,df_all,by="ENTREZID")
# write_xlsx(Glu,"C:/Users/Administrator/Desktop/Glu_TATB.xlsx")

# 排序合并后的数据框
df_all_sort <- df_all[order(df_all$log2FoldChange, decreasing = TRUE), ]

# 提取foldchange
gene_fc <- df_all_sort$log2FoldChange

# 将gene_fc对象的命名改为df_all_sort数据框中ENTREZID列的值
names(gene_fc) <- df_all_sort$ENTREZID

####基于KEGG基因集的GSEA富集
KEGG_ges <- gseKEGG(
  geneList = gene_fc,
  organism = "hsa",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

#重新转化为symbol基因
KEGG_ges <- setReadable(KEGG_ges,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID")
# 将结果转换为数据框
KEGG_table <- as.data.frame(KEGG_ges)
write_xlsx(KEGG_table,"./1.GSEA/KEGG_table8.xlsx")

