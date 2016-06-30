comb_plot <- function(f1,f2)
{ # plots two variables from different data frames in same graph by combining them
  a <- 'phosphatase'
  f1[paste(a)] = 'no phosphatase'
  f2[paste(a)] = '+ phosphatase'
#   names(f1)[3] <- names(f2)[3] <- 'F.Median'
#   names(f1)[4] <- names(f2)[4] <- 'B.Median'
#   
  f <- rbind(f1,f2)
#   p <- ggplot( data = f, aes(B.Median, fill = phosphatase, alpha = .5)) + geom_density()
#   print(p)
#   f
}