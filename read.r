# reads an array and makes file of intensities - can apply threshold
# set unik = 1 to avg two spots side by side

hit <- function(file_name)
{  # finding classes of each column to make it easier for reading the full table
  unik = 1 # should i merge unique entries?
  thresholding = 1 # should I remove points below threshold?
  a <- read.table('input files/test.txt',skip= 31, header = TRUE,nrows = 4)
  x <- y <- sapply(a,class)
  # Read in only reqired columns
  x[c(1:3,6:8,10:13,15:40)] <- y[c(5:40)] <- 'NULL'  # we need rows 4 name, 5 ID, 9 F median, 14 B median
  x[4] <- y[4]  <- 'character'
  a <- read.table(paste('results/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 46080,colClasses = x)
  b <- read.table(paste('results/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 46080,colClasses = y)
  
  # remove controls  : buffer, control, empty
  p <- grepl('empty',a[[2]]) | (grepl('buffer',a[[2]])) | (grepl('CONTROL',a[[2]]))
  a <- a[!p,c(1,3:dim(a)[2])]
  b <- b[!p,]
  # index unique spots
  #ix <- c(1:23040,1:23040)
  # index <- (1:46080 - 1)%/%2 + 1
  # a <- cbind(a,index)
  
  a <- cbind(a, 'F/B n' = a[[2]]/a[[3]]) # normalize foreground to background
  a <- cbind(a, 'F-B' = a[[2]]-a[[3]]) # foreground - Background 
  
  s <- a
  # look for block medians
  block_med <- 0
  for(i in 1:max(b$Block))
  {
    tmp <- b[[1]] == i  # all spots within block-i
    st <- s[tmp,]
    block_med[i] <- median(st[[4]]) 
    s[4][tmp,] <- st[[4]] - block_med[i] + 1   # normalize : within chip norm. median Ij
  }
  
  if (thresholding == 1)
  { index <- (1:46080 - 1)%/%2 + 1
    s <- cbind(s,index)
    # make dummy distribution of I <= 1
    d <- s[4][s[[4]] <= 1, ]
    d <- rbind (d, 2 - d ) # complete a mirror image of distribution
    thr <- mean(d) + 3 * sd(d)  # threshold mean + 6 sd of dummy distribution
    # filtering the rows/ spots
    t <- s[[4]] > thr  #mean(s[[7]]) + 1 * sd(s[[7]])   # Hit criterea - Greater than 3 standard deviations
    # s <- s[t,]
    s[!t,4] <- 0
  }
  
  if (unik == 1)
  {# merge duplicates
    s1 <- s[c(TRUE,FALSE), 2:dim(s)[2]]  # extract alternate elements to be avg
    s2 <- s[c(FALSE,TRUE), 2:dim(s)[2]]
    Name <- s[c(TRUE,FALSE),1]            # names
    index <- (1:dim(s1)[1])
    # ss <- s[c(TRUE,FALSE),1:3]          # other const elements
    s <- cbind( Name, index, (s1 + s2)/2 ) 
    s[[1]] <- as.character(s[[1]])
    s[c(3,4,6)] <-  round(s[c(3,4,6)])  # rounding off
    s[5] <- round(s[5],2)
  }
  s
  
  #write.table(s, file = paste('results/summary_significance/sum_',file_name,'.txt', sep = ''),sep = "\t", row.names = FALSE)
  #write(paste('\n Number of Hits :',dim(s)[1]), paste('results/summary_significance/sum_',file_name,'.txt', sep = ''), append = TRUE)
  #hits <- s  # output
}
