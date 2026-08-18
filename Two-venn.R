rm(list = ls())
library(tidyverse)
library(ggvenn)
gene_DEGs <- readRDS('output/GEO_DEGs.Rds')
# 内质网应激表型基因
gene_ERS <- read.table('input/Endoplasmic_reticulum_stress.txt', header = T)$Endoplasmic_reticulum_stress

# 交集韦恩图 ---------------------------------------------------------------

p <- ggvenn(
  list(DEGs = gene_DEGs, `ERS` = gene_ERS),
  # 下面要与列表中的命名一致
  c('DEGs', 'ERS'),
  # 不展示比例
  show_percentage = F,
  fill_alpha = 0.5,
  stroke_color = NA,
  fill_color = c('#4DBBD5', '#E64B35')
)

ggsave(file = 'output/geo_venn.pdf', p, width = 6, height = 6)


# 交集情况 ---------------------------------------------------------------
gene_geoERSRDEGs <- intersect(gene_DEGs, gene_ERS)
gene_geoERSRDEGs


# 数据保存 ---------------------------------------------------------------
saveRDS(gene_geoERSRDEGs, file = 'output/geoERSRDEGs.Rds')
write.csv(gene_geoERSRDEGs, file = 'output/geoDEGs-ERS.CSV')
