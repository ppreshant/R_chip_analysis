# applicable for radioactive film scans (Kinase J,F)
# reads an array and makes file of intensities - can apply threshold
# set unik = 1 to avg two spots - diagonal - 32 spots apart
# pTYR = wavelength 22 => manually set number 
# Caution: F/(B+1) is being used instead of F/B - to cover case where B = 0
hitsR <- function(file_name)
{ # parameters 
  
  unik = 4 # USEFUL # should i merge duplicates of every spot? 0 = no , 1 = Average, 4 = select the minimum of 2 values
  thresholding <- 3 # USEFUL #  threshold above how many standard deviations. (-1 => no thresholding)
  badspo = 0 # IGNORE # should I remove bad spots - high background critireon (1000) and ratio < 3?
  deflag <- 0 # IGNORE # get rid of bad spots manually flagged? 0 = NO, (-50 for Fgfr2 arrays, -100 for Jak2 arrays) = YES : not flagged in the data being used
  # norm_median <- 0 # should I normalize the median F/B to 1 in each block?
  
  # finding classes of each column to make it easier for reading only the required columns in the table 
  a <- read.table('input files/test.txt',skip= 31, header = TRUE,nrows = 4) # read a similar file to get column classes
  x <- y <- sapply(a,class)
  # Assign the reqired columns      - the vector y - it is to read blocks, flags and other auxillary data.
  x[c(1:3,6:8,10:13,15:length(x))] <- y[c(5:37,39:length(x))] <- 'NULL'  # classes marked NULL will not be read into the program : we need rows 4 name, 5 ID,9 F median, 14 B median for x, y reads the block number any flags etc) 
  x[4] <- y[4] <- 'character' 
  
  # reading the actual file 
  a <- read.table(paste(folo,'/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 48384,colClasses = x)       # 'results/F/w photoshop/'
  b <- read.table(paste(folo,'/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 48384,colClasses = y) # dummy file to identify blocks, flags etc.
  
  # remove controls  : buffer, control, empty
  p <- grepl('empty',a[[2]]) | (grepl('buffer',a[[2]])) | (grepl('Control',a[[2]]))
  a <- a[!p,c(1,3:dim(a)[2])]
  b <- b[!p,]

  a <- cbind(a, 'F/B n' = a[[2]]/(a[[3]]+1)) # to normalize foreground to background : into a new column
  a <- cbind(a, 'F-B' = a[[2]]-a[[3]]) # foreground - Background : into a new column
  
  s <- a
  # look for block medians - used when normalizing all block medians to be = 1
  # block_med <- 0
  # for(i in 1:max(b$Block))
  # {
  #   tmp <- b[[1]] == i  # all spots within block-i
  #   st <- s[tmp,]
  #   block_med[i] <- median(st[[4]]) 
  #   s[4][tmp,] <- st[[4]] - block_med[i] + 1   # normalize : within chip norm. median Ij
  # }
  
  # filtering steps - remove Inf, flags, thresholding etc. 
  if(thresholding != -1)
  {thr <- (s$`F/B n` >= find_thr(s,thresholding))   # thresholding - 3 sd above mean for dummy dist
  s[!thr,4] <- 0                                    # if spot < threshold, make it 0
  }
  
  # to remove manually flagged items - spots around saturated areas etc.
  if(deflag == 1)
  { 
    if (grepl('F_*',file_name)) flg <- -50  else if (grepl('J_*',file_name)) flg <- -100
    sf <- b$Flags == flg
    s[sf,4] <- s[sf,2] <- s[sf,3] <- 0
  }
  
  # merge duplicates - average of two data points
  if (unik == 1)
  {
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
  }
  
  # Dont merge duplicates - returns all points
  else if (unik == 0)
  { 
    s1 <- s[c(rep(T,32),rep(FALSE,32)), 2:dim(s)[[2]]]
    s2 <- s[c(rep(FALSE,32),rep(T,32)), 2:dim(s)[[2]]]
    # s1 <- s[c(TRUE,FALSE), 2:dim(s)[2]]  # extract alternate elements to be avg
    # s2 <- s[c(FALSE,TRUE), 2:dim(s)[2]]
    
    Name <- s$Name            # names
    index <- (1:dim(s)[1])
    # ss <- s[c(TRUE,FALSE),1:3]          # other const elements
    
    # filtering steps - remove Inf, flags, thresholding etc.
    # s1[(s1[3] == Inf),3] <- 1    # remove Inf
    # s2[(s2[3] == Inf),3] <- 1
    
    tt <- (s[c(rep(T,32),rep(FALSE,32)),3] == 0) | (s[c(rep(FALSE,32),rep(T,32)),3] == 0)  # If one spot is < threshold (it will be 0), make the other one also o
    s[tt,3] <- s[tt,3] <- 0
    
    
    s <- cbind(Name, index, rbind(s1, s2)) 
  }
  
  # Dont merge duplicates - return minimum values of the 2 duplicate spots
  else if (unik == 4)
  { 
    s1 <- s[c(rep(T,32),rep(FALSE,32)), 2:dim(s)[[2]]]
    s2 <- s[c(rep(FALSE,32),rep(T,32)), 2:dim(s)[[2]]]
    # s1 <- s[c(TRUE,FALSE), 2:dim(s)[2]]  # extract alternate elements to be avg
    # s2 <- s[c(FALSE,TRUE), 2:dim(s)[2]]
    
    Name <- s[c(rep(T,32),rep(FALSE,32)), 1]            # names
    index <- (1:dim(s1)[1])

    # make filtering on both spots - used when removing spots below the threshold
    tt <- (s[c(rep(T,32),rep(FALSE,32)),3] == 0) | (s[c(rep(FALSE,32),rep(T,32)),3] == 0)  # If one spot is < threshold (it will be 0), make the other one also o
    s[tt,3] <- s[tt,3] <- 0
    
    s <- cbind(Name, index, pmin(s1, s2)) 
  }
  
  # formatting, rounding off. 
  s[[1]] <- as.character(s[[1]])
  s[c(3,4,6)] <-  round(s[c(3,4,6)])  # rounding off
  s[5] <- round(s[5],2)
  
  
  if(badspo == 1) s <- s[(s[[4]] < 1000 |  s[[5]] > 3), ] #IGNORE: Testing:  arbitrary criterea backg > 1000 and ratio <3 excluded
  
  # outputs the data frame  - Data points above the threshold
  if(thresholding != -1) s <- s[s[5] > 0,]   else s           # only retain spots above threshold 
  
  #write.table(s, file = paste('results/summary_significance/sum_',file_name,'.txt', sep = ''),sep = "\t", row.names = FALSE)
  #write(paste('\n Number of Hits :',dim(s)[1]), paste('results/summary_significance/sum_',file_name,'.txt', sep = ''), append = TRUE)
  #hits <- s  # output
}

