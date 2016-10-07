setwd("E:/R files_array comparision") # change it to the dropbox location of the R files_array comparision folder

# give the path of the input files (it should be the same unless changed)
folo <- 'results/Kinase/F_late_2dayexp'  # for Fgfr2 arrays
# folo <- 'results/Kinase/J_6dayexp' # for Jak2 arrays

# reads the list of all 12 files in the folder (Control, WT, 4 mutants all in duplicate = 6 * 2 = 12)
fl <- list.files(folo, pattern = '*.gpr$', recursive= FALSE, include.dirs = TRUE)

# defining path for output of figures etc.
oufl <- 'results/latest plots/min of 4 spots/'

# sourcing ggplot2 and cowplot - to install the packages use install.packages(ggplot2) and install.packages(cowplot)
library(ggplot2)
library(cowplot)

# sourcing the scripts being used

source('E:/R files_array comparision/readarray.r') # for the radioactive assay - 
source('E:/R files_array comparision/check_hits_g.r')
source('E:/R files_array comparision/make_scatterplots.r')
# source('E:/R files_array comparision/simple_read.r') # for pilot 4,5 arrays - not tested to work

# Fl <- lapply(fl,make_checklist)
# fc <- make_sctr(6,3,Fl)

# C <- make_checklist(c)
# J <- make_checklist(j)
# jc <- make_sctr(1,2,list(J,C))

Flf <- lapply(fl,make_checklist) # Reads all files in the folder, makes a checklist column for each file- for known hit(from Heng's paper), known Protein protein interaction, etc.
Fl2f <- all_min4() # Takes all 12 files in the folder and Merges data from the 2 duplicate arrays retaining the minimum of the 2 values for each spot

# draws all the required scatterplots
# a <- graphitall(Fl2,'F')


# --------------------IGNORE BELOW THIS LINE-------------------------------------------

# f <- 'F_100181_pmt-950_5um.gpr'
# con <- 'C_100182_pmt-950_5um.gpr'
# g <- 'gst/F_pilot1_gst.gpr'
# c <- "pilot4/C_pilot4_9049232.gpr"
# f <- "pilot4/F_pilot4_9049233.gpr"
# c <- 'pilot5/C_pilot5.gpr'
# j <- 'pilot5/J_pilot5.gpr'
# F <- hits(f)
# c <- 'C_F1.gpr'
# f <- 'F2.gpr'
# 

# folo <- 'C:/Dropbox (Zhu Lab)/prashant_zhu lab/results/Huprot scans/Kinase/J' 
# fl <- list.files('results/F/w photoshop', pattern = '*.gpr$', recursive= FALSE, include.dirs = TRUE)
#fl <- list.files('results/pilot3', pattern = '*.gpr$', recursive= F, include.dirs = TRUE)

# f2 <- 'F_100185_pilot2.gpr'
# c2 <- 'C_100183_pilot2.gpr'

# source('E:/R files_array comparision/simpler.r')

# Fl2 <- all_merg()
# # Al <- rbind(cbind(x[[1]], Kinase = 'Control'),cbind(x[[2]], Kinase = 'Fgfr2'),cbind(x[[3]], Kinase = 'F1'),cbind(x[[4]], Kinase = 'F2'),cbind(x[[5]], Kinase = 'F3'),cbind(x[[6]], Kinase = 'F4'))
# # p <- ggplot(data = Alf, aes(B22.Median, fill = Kinase, alpha = .5)) + geom_density() + theme_classic() + ggtitle('Full background density')
# Flr2 <- all_rmcntrl(Fl2)
# # plot on multiple rows - for sig spots - int comparision
# # plt <- ggplot(data = Als, aes(Name, `Kinase Intensity`, fill = Kinase)) + geom_bar(stat = 'identity', position = 'dodge') + ggtitle('Significant spots : Intensity comparision') + scale_fill_brewer(palette = 'Dark2') + theme(axis.text = element_text(angle = 90)) + facet_wrap(~dum, scales = 'free')
# # > plt + theme(strip.background = element_blank(), strip.text.x = element_blank())
