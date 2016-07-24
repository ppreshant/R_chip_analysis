# applicable for radioactive film scans (Kinase J,F)
# reads an array and makes file of intensities - can't apply threshold
# set unik = 1 to avg two spots - diagonal
# pTYR = wavelength 22 => manually set number 

all_hitsR <- function(file_name)
{  # finding classes of each column to make it easier for reading the full table
  
  unik = 1 # should i merge unique entries?
  badspo = 0 # should I remove bad spots - high background critireon (1000) and ratio < 3?
  thresholding <- -1 #  threshold above how many sd. (-1 => no thresholding)
  
  a <- read.table('input files/test.txt',skip= 31, header = TRUE,nrows = 4) # read a similar file to get column classes
  x <- y <- sapply(a,class)
  # Read in only reqired columns
  x[c(1:3,6:8,10:13,15:length(x))] <- y[c(5:37,39:length(x))] <- 'NULL'  # we need rows 4 name, 5 ID,21 F median, 26 B median) 
  x[4] <- y[4] <- 'character'
  a <- read.table(paste('results/F/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 46080,colClasses = x)
  b <- read.table(paste('results/F/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 46080,colClasses = y) # dummy file to identify blocks
  
  # remove controls  : buffer, control, empty
  p <- grepl('empty',a[[2]]) | (grepl('buffer',a[[2]])) | (grepl('Control',a[[2]]))
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
  
  if(thresholding != -1)
  {thr <- (s$`F/B n` >= find_thr(s,thresholding))   # thresholding - 3 sd above mean for dummy dist
  s[!thr,4] <- 0                                    # if spot < threshold, make it 0
  }
  
  if (unik == 1)
  {# merge duplicates
    s1 <- s[c(rep(T,32),rep(FALSE,32)), 2:dim(s)[[2]]]
    s2 <- s[c(rep(FALSE,32),rep(T,32)), 2:dim(s)[[2]]]
    # s1 <- s[c(TRUE,FALSE), 2:dim(s)[2]]  # extract alternate elements to be avg
    # s2 <- s[c(FALSE,TRUE), 2:dim(s)[2]]
    
    Name <- s[c(rep(T,32),rep(FALSE,32)), 1]            # names
    index <- (1:dim(s1)[1])
    # ss <- s[c(TRUE,FALSE),1:3]          # other const elements
   
     # filtering steps - remove Inf, flags, thresholding etc.
    s1[(s1[3] == Inf),3] <- 1    # remove Inf
    s2[(s2[3] == Inf),3] <- 1
    
    tt <- (s1[3] == 0) | (s2[3] == 0)  # If one spot is < threshold (it will be 0), make the other one also o
    s1[tt,3] <- s2[tt,3] <- 0
    
    
    s <- cbind( Name, index, (s1 + s2)/2 ) 
    s[[1]] <- as.character(s[[1]])
    s[c(3,4,6)] <-  round(s[c(3,4,6)])  # rounding off
    s[5] <- round(s[5],2)
  }
  if(badspo == 1) s <- s[(s[[4]] < 1000 |  s[[5]] > 3), ] # arbitrary criterea backg > 1000 and ratio <3 excluded
  s <- s[s[5] > 0,]              # only retain spots above threshold 
  
  
  
  #write.table(s, file = paste('results/summary_significance/sum_',file_name,'.txt', sep = ''),sep = "\t", row.names = FALSE)
  #write(paste('\n Number of Hits :',dim(s)[1]), paste('results/summary_significance/sum_',file_name,'.txt', sep = ''), append = TRUE)
  #hits <- s  # output
}

