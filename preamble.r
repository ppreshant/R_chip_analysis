setwd("E:/R files_array comparision")
# f <- 'F_100181_pmt-950_5um.gpr'
# con <- 'C_100182_pmt-950_5um.gpr'
# g <- 'gst/F_pilot1_gst.gpr'
# c <- "pilot4/C_pilot4_9049232.gpr"
# f <- "pilot4/F_pilot4_9049233.gpr"
# c <- 'pilot5/C_pilot5.gpr'
# j <- 'pilot5/J_pilot5.gpr'
# F <- hits(f)
c <- 'C_F1.gpr'
f <- 'F2.gpr'

fl <- list.files('results/F', pattern = '*.gpr$', recursive= FALSE, include.dirs = TRUE)
#fl <- list.files('results/pilot3', pattern = '*.gpr$', recursive= F, include.dirs = TRUE)

# f2 <- 'F_100185_pilot2.gpr'
# c2 <- 'C_100183_pilot2.gpr'
library(ggplot2)
# library(cowplot)

# source('E:/R files_array comparision/simple_read.r')
source('E:/R files_array comparision/simpler.r')
source('E:/R files_array comparision/check_hits_g.r')
source('E:/R files_array comparision/make_scatterplots.r')
# Fl <- lapply(fl,make_checklist)
# fc <- make_sctr(6,3,Fl)

# C <- make_checklist(c)
# J <- make_checklist(j)
# jc <- make_sctr(1,2,list(J,C))