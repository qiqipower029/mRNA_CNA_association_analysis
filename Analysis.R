#------------------------Association Test--------------------------------#
#--------------------------Jieqi Tu--------------------------------------#



#---------------------HNSC dataset---------------------------#
# Tidy data
library(tidyverse)
CNA_HNSC = read.delim("./Downloaded Data/CNA HNSC.txt")
mRNA_HNSC = read.delim("./Downloaded Data/mRNA expression z-scores HNSC.txt")
HNSC = cbind.data.frame(CNA_HNSC, mRNA_HNSC) 
HNSC = HNSC[-c(5,6)]
