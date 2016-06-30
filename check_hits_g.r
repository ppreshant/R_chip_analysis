 # Many random functions here
 # gives all hits in the list from phosphonetworks with their data 
check <- function(x,a)
{  
  t <- x == as.character(a[['Name']])  # check the names
  #a[t,]
  
  #check <- 1
}

tst_J <- read.table(paste('input files/','J_substrates.txt', sep = ''), header = FALSE)
tst_F <- read.table(paste('input files/','F_substrates.txt', sep = ''), header = FALSE)
F_PPI_string <- read.table(paste('input files/','F_PPI_string.txt', sep = ''), header = FALSE)
F_PPI_biogrid <- read.table(paste('input files/','F_PPI_biogrid.txt', sep = ''), header = FALSE)
J_PPI_string <- read.table(paste('input files/','J_PPI_string.txt', sep = ''), header = FALSE)
J_PPI_biogrid <- read.table(paste('input files/','J_PPI_biogrid.txt', sep = ''), header = FALSE)
raw_J <- read.table(paste('input files/','raw_J.txt', sep = ''), header = FALSE)
raw_F <- read.table(paste('input files/','raw_F.txt', sep = ''), header = FALSE)

#s <- 0
#ss <- lapply(c("APITD1","BCL3"), check, a = a)#, simplify = 'array')

# 
#  give_list <- function(tst,fl)
# {
# p <- 0
# #t <- FALSE
# for (x in as.character(tst[[1]]))
# {
#   tt <-  check(x,fl) # check same name as in list
#   #t <- as.logical(t + tt)
#   if(sum(tt) > 0)
#   {
#   e <- fl[tt,]
#   #e <- e[c(1:8)]
#   p <-  rbind(p,e)
#    }
# }
# p <- p[2:dim(p)[1],]
# #b <- fl[t,]
# #write.table(p, file = 'results/summary_significance/Fgfr2.txt',sep = "\t", row.names = FALSE)
#  }
#  

 # checks a list and flags the spot with same names as variable 'ms' under status column 
 #  tst <- array file handle
 #  fl  <- data frame of the array
 #  ms  <- name of the flag
 check_list <- function(tst,fl, ms)
 {
   levels(fl$status) <- c('normal','PPI','PPI-expt', 'raw_hits', 'known hit')
   p <- 0
   #t <- FALSE
   for (x in as.character(tst[[1]]))
   {
     tt <-  check(x,fl) # check same name as in list
     if(sum(tt) > 0)
     {
       fl$status[tt] <- ms
     }
   }
   fl 
 }