# # run simple read and append phosphonetworks hits
# source check_hits_g for check_list function

# f <- cbind(all_hits(fl[9]), check = 0)
# # ft <- give_list(tst_F,f)
# # ft$check = 1
# # f <- rbind(f,ft)
# f <- check_list(tst_F,f)
# # c <- all_hits(fl[5])
# c <- cbind(all_hits(fl[5]), check = 0)
# # ct <-  give_list(tst_F,c)
# # c <- rbind(c,ct)
# # qplot(c[[5]],f[[5]], xlim = c(0,5), ylim = c(0,5), colour = as.factor(f$check), size = f$check) + scale_color_manual(breaks = c("0", "1"), values=c("blue", "red"))
# qplot(c[[5]],f1[[5]], xlim = c(0,5), ylim = c(0,5), colour = as.factor(f1$check)) + scale_color_manual(breaks = c("0", "1"), values=c("grey", "black")) + ggtitle ('Pilot3 Comparision') + xlab('control') + ylab('Kinase')

make_checklist <- function(a)
{
  f <- cbind(all_hits(a), status = 'normal')
  if (grepl('F_*',a)) 
  {
#     f <- check_list(F_PPI_string,f,'PPI')
#     f <- check_list(F_PPI_biogrid,f,'PPI-expt')
    f <- check_list(raw_F,f,'raw_hits')
    f <- check_list(tst_F,f,'known hit')
    
  }
  
  if (grepl('J_*',a)) 
    {
    # f <- check_list(tst_J,f)
#     f <- check_list(J_PPI_string,f,'PPI')
#     f <- check_list(J_PPI_biogrid,f,'PPI-expt')
    f <- check_list(raw_F,f,'raw_hits')
    f <- check_list(tst_J,f,'known hit')
    }
  # write.table(f, file = paste('results/no threshold/',a,'.txt', sep = ''),sep = "\t", row.names = FALSE)
  
  f
} 

make_sctr <- function(k,c, Fl)
{ # takes 2 numbers for kinase, control and Fl of the make_checklist lapplied
  # kin <- readline(prompt = 'which Kinase is this?')
  fa <- Fl[[k]][order(Fl[[k]]$status),]
  ca <- Fl[[c]][order(Fl[[k]]$status),]
  dat <- data.frame(Name = fa$Name, Jak2.Intensity = fa$`F/B n`, Control.Intensity = ca$`F/B n`, status = fa$status, Ratio = fa$`F/B n`/ca$`F/B n` , stringsAsFactors = FALSE)
  s <- ggplot(data = dat, aes(Control.Intensity,Jak2.Intensity, color = status, alpha = .5)) + geom_point() + scale_color_manual(breaks = c("normal","PPI","PPI-expt", "known hit"), values=c("yellow","green",'blue', "black")) + geom_abline(intercept = 0, slope = 1) + geom_hline(yintercept = median(dat$Jak2.Intensity)) #+ xlim(0,10) + ylim(0,10)
  print(s)
  dat
}

# f <- Fl[[4]][order(Fl[[4]]$status),]
# c <- Fl[[1]][order(Fl[[4]]$status),]

# some other plot
## ggplot(data = cc, aes(index, `F/B n`, color = scanner)) + geom_point() + coord_cartesian(ylim = c(0,2))
find_thr <- function(s,x)
{
  # make dummy distribution of I <= 1
  d <- s[5][s[[5]] <= 1, ]
  d <- rbind (d, 2 - d ) # complete a mirror image of distribution
  thr <- mean(d) + x * sd(d)  # threshold mean + 6 sd of dummy distribution
  # # filtering the rows/ spots
  # t <- s[[8]] > thr  #mean(s[[7]]) + 1 * sd(s[[7]])   # Hit criterea - Greater than 3 standard deviations
  # s <- s[t,]
}