####时序分析####
#用差异代谢物做表达模式的探索

# BiocManager::install("Mfuzz")
# 加载
library(Mfuzz)
library(readxl)
library(ClusterGVis)
library(RColorBrewer)
library(tidyverse)

setwd("D:\\project\\7.转录组数据\\66.时序分析")

color <- colorRampPalette(rev(c("#E02401", "Yellow","#0F52BA", "#3E7C17")))(1000)

mfuzz=read_excel("./mfuzz.xlsx") %>% column_to_rownames(var="ID")

Control_D4=apply(mfuzz[,1:3],1,FUN = mean)
Control_D7=apply(mfuzz[,4:6],1,FUN = mean)
Control_D8=apply(mfuzz[,7:9],1,FUN = mean)
ES_D4=apply(mfuzz[,10:12],1,FUN = mean)
ES_D7=apply(mfuzz[,13:15],1,FUN = mean)
PostES_D8=apply(mfuzz[,16:18],1,FUN = mean)


Mfuzz= data.frame(Con4=Control_D4,
                  Con7=Control_D7,
                  Con8=Control_D8,
                  ES4=ES_D4,
                  ES7=ES_D7,
                  PostES8=PostES_D8)

Mfuzz=apply(Mfuzz,2,as.numeric)
rownames(Mfuzz)=rownames(mfuzz)


getClusters(exp = mfuzz)   ##9个左右


#构建 mfuzz 对象并处理缺失及异常值
mfuzz_class <- new('ExpressionSet',exprs = Mfuzz)
# mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
# mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)

mfuzz_class <- standardise(mfuzz_class)

# 基于 fuzzy c-means 的算法进行聚类，详情 ?mfuzz
set.seed(123)
cluster_num <- 9


mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))

# 作图，详情 ?mfuzz.plot2
color.2 <- colorRampPalette(rev(c("#ff0000", "Yellow", "OliveDrab1")))(1000)



# 导出为PDF矢量图（适合高质量印刷）
pdf("./mfuzz_clusters.pdf", width =10, height = 8)
mfuzz.plot2(mfuzz_class, cl = mfuzz_cluster, mfrow = c(3, 3), 
            time.labels = colnames(Mfuzz), xlab="", x11 = F)
dev.off()





#查看每个cl中的差异基因数量
cluster_size <- mfuzz_cluster$size
names(cluster_size) <- 1:cluster_num
cluster_size

#查看每个差异基因所属的cl
head(mfuzz_cluster$cluster)
# Mfuzz 通过计算一个叫 membership 的统计量判断差异基因所属的聚类群，以最大的 membership 值为准
# 查看每个基因的 membership 值
head(mfuzz_cluster$membership)

# 最后，提取所有差异代谢物所属的聚类群
gene_cluster <- mfuzz_cluster$cluster
gene_cluster <- cbind(mfuzz[names(gene_cluster), ], gene_cluster)

head(gene_cluster)

write.table(gene_cluster, './cluster.xls', sep = '\t', col.names = NA, quote = FALSE)
