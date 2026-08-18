rm(list = ls())
library(tidyverse)
library(limma)
library(ggplot2)
dat_expr <- readRDS('input/GEO_combined_dataset.Rds')
dat_group <- readRDS('input/GEO_combined_group.Rds')
# write.csv(dat_expr, file='input/GEO_dataset.csv')
# write.csv(dat_group, file='input/GEO_group.csv')

# 差异分析 ---------------------------------------------------------------

# 创建模型矩阵
design <- model.matrix( ~ 0 + factor(dat_group$group))
colnames(design) <- levels(factor(dat_group$group))
rownames(design) <- dat_group$sample

# 使用线性模型拟合表达数据
fit <- lmFit(dat_expr, design)

# 为组间差异创建对照
contrasts <- c('PCOS-Control')
# 一定要疾病组在前，对照组在后

# 创建线性对比矩阵
cont_matrix <- makeContrasts(contrasts = contrasts, levels = design)

# 对拟合的线性模型进行再次拟合
fit2 <- contrasts.fit(fit, cont_matrix)

# Bayes调整
fit2 <- eBayes(fit2)

# 提取差异分析结果
dat_res_diff <- topTable(fit2, coef = contrasts, n = Inf) %>% na.omit()

# 将为0的p值改成除0外最小的p值，以防止p值为0取不到对数
dat_res_diff$P.Value[dat_res_diff$P.Value == 0] <- min(dat_res_diff$P.Value[dat_res_diff$P.Value != 0])
# dat_res_diff$adj.P.Val[dat_res_diff$adj.P.Val == 0] <- min(dat_res_diff$adj.P.Val[dat_res_diff$adj.P.Val != 0])

# 筛选差异基因的阈值
logFC_cutoff <- 0.5
dat_res_diff$change <- ifelse(
  dplyr::between(dat_res_diff$logFC, -logFC_cutoff, logFC_cutoff) |
    dat_res_diff$P.Value >= 0.05,
  'Not',
  ifelse(dat_res_diff$logFC >= logFC_cutoff, 'Up', 'Down')
)
table(dat_res_diff$change)

# 差异基因
gene_DEGs <- dat_res_diff %>%
  dplyr::filter(change != 'Not') %>%
  dplyr::arrange(desc(logFC)) %>%
  rownames()


# 火山图 ---------------------------------------------------------------

# 火山图数据并对p值取负对数
dat_vocano <- dat_res_diff[, c('logFC', 'P.Value', 'change')]
dat_vocano$logP <- -log10(dat_vocano$P.Value)
dat_vocano$change <- factor(dat_vocano$change, levels = c('Down', 'Not', 'Up'))

# 火山图
p <- ggplot(dat_vocano, aes(logFC, logP, color = change)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  labs(x = 'LogFC', y = '-Log10(adj. p)', color = 'Change') +
  # 添加水平线和垂直线，线型均为虚线，位置分别在p值阈值（的负对数）和正负logFC阈值
  geom_hline(yintercept = -log10(0.05), lty = 2) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 2) +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = c('#4DBBD5', 'grey', '#E64B35'))

ggsave(file = 'output/vocano_GEO_limma.pdf', p, width = 6, height = 6)

# 数据保存 ---------------------------------------------------------------
saveRDS(dat_res_diff, file = 'output/GEO_diff_res.Rds')
saveRDS(gene_DEGs, file = 'output/GEO_DEGs.Rds')
