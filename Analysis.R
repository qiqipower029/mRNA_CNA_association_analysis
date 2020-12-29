#------------------------Association Test--------------------------------#
#--------------------------Jieqi Tu--------------------------------------#

library(tidyverse)
library(ggpubr)

#---------------------HNSC dataset---------------------------#
# Tidy data for HNSC
CNA_HNSC = read.delim("./Downloaded Data/CNA HNSC.txt")
mRNA_HNSC = read.delim("./Downloaded Data/mRNA expression z-scores HNSC.txt")
HNSC = cbind.data.frame(CNA_HNSC, mRNA_HNSC) 
HNSC = HNSC[-c(5,6)]

#-------FAT1 analysis in HNSC---------#
FAT1_HNSC = HNSC[c(1, 2, 3, 5)] 
colnames(FAT1_HNSC) = c("Study ID", "Sample ID", "Copy Number Alterations", "mRNA Expression")  

# Making linear regression plot
P1 = 
FAT1_HNSC %>% 
  ggplot(aes(x = `Copy Number Alterations`, y = `mRNA Expression`)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm") + theme_bw() + 
  labs(caption = "FAT1 in HNSCC")

# Pearson Correlation
cor.test(FAT1_HNSC$`Copy Number Alterations`, FAT1_HNSC$`mRNA Expression`, method="pearson")

# Spearman Correlation
cor.test(FAT1_HNSC$`Copy Number Alterations`, FAT1_HNSC$`mRNA Expression`, method="spearman")

# Kendall Correlation
cor.test(FAT1_HNSC$`Copy Number Alterations`, FAT1_HNSC$`mRNA Expression`, method="kendall")

# Anova analysis
FAT1_HNSC$`Copy Number Alterations` = as.character(FAT1_HNSC$`Copy Number Alterations`)
analysis1 = aov(FAT1_HNSC$`mRNA Expression`~FAT1_HNSC$`Copy Number Alterations`)
summary.aov(analysis1)


#-------EGFR analysis in HNSC---------#
EGFR_HNSC = HNSC[c(1, 2, 4, 6)] 
colnames(EGFR_HNSC) = c("Study ID", "Sample ID", "Copy Number Alterations", "mRNA Expression")  

# Making linear regression plot
P2 = 
EGFR_HNSC %>% 
  ggplot(aes(x = `Copy Number Alterations`, y = `mRNA Expression`)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm") + theme_bw() + 
  labs(caption = "EGFR in HNSCC")

# Pearson Correlation
cor.test(EGFR_HNSC$`Copy Number Alterations`, EGFR_HNSC$`mRNA Expression`, method="pearson")

# Spearman Correlation
cor.test(EGFR_HNSC$`Copy Number Alterations`, EGFR_HNSC$`mRNA Expression`, method="spearman")

# Kendall Correlation
cor.test(EGFR_HNSC$`Copy Number Alterations`, EGFR_HNSC$`mRNA Expression`, method="kendall")

# Anova analysis
EGFR_HNSC$`Copy Number Alterations` = as.character(EGFR_HNSC$`Copy Number Alterations`)
analysis2 = aov(EGFR_HNSC$`mRNA Expression`~EGFR_HNSC$`Copy Number Alterations`)
summary.aov(analysis2)


#---------------------LSCC dataset---------------------------#
# Tidy data for LSCC
CNA_LSCC = read.delim("./Downloaded Data/CNA LSCC.txt")
mRNA_LSCC = read.delim("./Downloaded Data/mRNA expression z-scores LSCC.txt")
LSCC = cbind.data.frame(CNA_LSCC, mRNA_LSCC) 
LSCC = LSCC[-c(5,6)]

#-------FAT1 analysis in LSCC---------#
FAT1_LSCC = LSCC[c(1, 2, 3, 5)] 
colnames(FAT1_LSCC) = c("Study ID", "Sample ID", "Copy Number Alterations", "mRNA Expression")  

# Making linear regression plot
P3 = 
FAT1_LSCC %>% 
  ggplot(aes(x = `Copy Number Alterations`, y = `mRNA Expression`)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm") + theme_bw() + 
  labs(caption = "FAT1 in Lung SCC")

# Pearson Correlation
cor.test(FAT1_LSCC$`Copy Number Alterations`, FAT1_LSCC$`mRNA Expression`, method="pearson")

# Spearman Correlation
cor.test(FAT1_LSCC$`Copy Number Alterations`, FAT1_LSCC$`mRNA Expression`, method="spearman")

# Kendall Correlation
cor.test(FAT1_LSCC$`Copy Number Alterations`, FAT1_LSCC$`mRNA Expression`, method="kendall")

# Anova analysis
FAT1_LSCC$`Copy Number Alterations` = as.character(FAT1_LSCC$`Copy Number Alterations`)
analysis3 = aov(FAT1_LSCC$`mRNA Expression`~FAT1_LSCC$`Copy Number Alterations`)
summary.aov(analysis3)


#-------EGFR analysis in LSCC---------#
EGFR_LSCC = LSCC[c(1, 2, 4, 6)] 
colnames(EGFR_LSCC) = c("Study ID", "Sample ID", "Copy Number Alterations", "mRNA Expression")  

# Making linear regression plot
P4 = 
EGFR_LSCC %>% 
  ggplot(aes(x = `Copy Number Alterations`, y = `mRNA Expression`)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm") + theme_bw() + 
  labs(caption = "EGFR in Lung SCC")

# Pearson Correlation
cor.test(EGFR_LSCC$`Copy Number Alterations`, EGFR_LSCC$`mRNA Expression`, method="pearson")

# Spearman Correlation
cor.test(EGFR_LSCC$`Copy Number Alterations`, EGFR_LSCC$`mRNA Expression`, method="spearman")

# Kendall Correlation
cor.test(EGFR_LSCC$`Copy Number Alterations`, EGFR_LSCC$`mRNA Expression`, method="kendall")

# Anova analysis
EGFR_LSCC$`Copy Number Alterations` = as.character(EGFR_LSCC$`Copy Number Alterations`)
analysis4 = aov(EGFR_LSCC$`mRNA Expression`~EGFR_LSCC$`Copy Number Alterations`)
summary.aov(analysis4)


#---------------------CSCC dataset---------------------------#
# Tidy data for CSCC
CNA_CSCC = read.delim("./Downloaded Data/CNA CSCC.txt")
mRNA_CSCC = read.delim("./Downloaded Data/mRNA expression z-scores CSCC.txt")
CSCC = cbind.data.frame(CNA_CSCC, mRNA_CSCC) 
CSCC = CSCC[-c(5,6)]

#-------FAT1 analysis in CSCC---------#
FAT1_CSCC = CSCC[c(1, 2, 3, 5)] 
colnames(FAT1_CSCC) = c("Study ID", "Sample ID", "Copy Number Alterations", "mRNA Expression")  

# Making linear regression plot
P5 = 
FAT1_CSCC %>% 
  ggplot(aes(x = `Copy Number Alterations`, y = `mRNA Expression`)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm") + theme_bw() + 
  labs(caption = "FAT1 in Cervical SCC")

# Pearson Correlation
cor.test(FAT1_CSCC$`Copy Number Alterations`, FAT1_CSCC$`mRNA Expression`, method="pearson")

# Spearman Correlation
cor.test(FAT1_CSCC$`Copy Number Alterations`, FAT1_CSCC$`mRNA Expression`, method="spearman")

# Kendall Correlation
cor.test(FAT1_CSCC$`Copy Number Alterations`, FAT1_CSCC$`mRNA Expression`, method="kendall")

# Anova analysis
FAT1_CSCC$`Copy Number Alterations` = as.character(FAT1_CSCC$`Copy Number Alterations`)
analysis5 = aov(FAT1_CSCC$`mRNA Expression`~FAT1_CSCC$`Copy Number Alterations`)
summary.aov(analysis5)


#-------EGFR analysis in CSCC---------#
EGFR_CSCC = CSCC[c(1, 2, 4, 6)] 
colnames(EGFR_CSCC) = c("Study ID", "Sample ID", "Copy Number Alterations", "mRNA Expression")  

# Making linear regression plot
P6 = 
EGFR_CSCC %>% 
  ggplot(aes(x = `Copy Number Alterations`, y = `mRNA Expression`)) + geom_point(alpha = 0.3) + geom_smooth(method = "lm") + theme_bw() + 
  labs(caption = "EGFR in Cervical SCC")

# Pearson Correlation
cor.test(EGFR_CSCC$`Copy Number Alterations`, EGFR_CSCC$`mRNA Expression`, method="pearson")

# Spearman Correlation
cor.test(EGFR_CSCC$`Copy Number Alterations`, EGFR_CSCC$`mRNA Expression`, method="spearman")

# Kendall Correlation
cor.test(EGFR_CSCC$`Copy Number Alterations`, EGFR_CSCC$`mRNA Expression`, method="kendall")

# Anova analysis
EGFR_CSCC$`Copy Number Alterations` = as.character(EGFR_CSCC$`Copy Number Alterations`)
analysis6 = aov(EGFR_CSCC$`mRNA Expression`~EGFR_CSCC$`Copy Number Alterations`)
summary.aov(analysis6)

ggarrange(P1, P2, P3,
          P4, P5, P6,
          ncol = 2, nrow = 3)
