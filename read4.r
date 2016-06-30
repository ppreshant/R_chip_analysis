# applicable for dual scan two wavelength(s) file (Pilot 4,5)
# reads an array and makes file of intensities - can apply threshold
# set unik = 1 to avg two spots side by side
# pTYR = wavelength 532 => Cy3 

hit4 <- function(file_name)
{  # finding classes of each column to make it easier for reading the full table
  
  unik = 0 # should i merge unique entries?
  badspo = 0 # should I remove bad spots - high background critireon (1000) and ratio < 3?
  
  a <- read.table('results/pilot4/C_pilot4_9049232.gpr',skip= 31, header = TRUE,nrows = 4)
  x <- y <- sapply(a,class)
  # Read in only reqired columns
  x[c(1:3,6:20,22:25,27:length(x))] <- y[5:length(x)] <- 'NULL'  # we need rows 4 name, 5 ID,21 F median, 26 B median) of 532 - ptyr
  x[4] <- y[4] <- 'character'
  a <- read.table(paste('results/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 46080,colClasses = x)
  b <- read.table(paste('results/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 46080,colClasses = y)
  
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
  
  if (unik == 1)
  {# merge duplicates
    s1 <- s[c(rep(T,32),rep(F,32)), 2:dim(s)[2]]
    s2 <- s[c(rep(F,32),rep(T,32)), 2:dim(s)[2]]
    # s1 <- s[c(TRUE,FALSE), 2:dim(s)[2]]  # extract alternate elements to be avg
    # s2 <- s[c(FALSE,TRUE), 2:dim(s)[2]]
    
    Name <- s[c(rep(T,32),rep(F,32)), 1]            # names
    index <- (1:dim(s1)[1])
    # ss <- s[c(TRUE,FALSE),1:3]          # other const elements
    s <- cbind( Name, index, (s1 + s2)/2 ) 
    s[[1]] <- as.character(s[[1]])
    s[c(3,4,6)] <-  round(s[c(3,4,6)])  # rounding off
    s[5] <- round(s[5],2)
  }
  if(badspo == 1) s <- s[(s[[4]] < 1000 |  s[[5]] > 3), ] # arbit criterea backg > 1000 and ratio <3 excluded
  s
  
  
  
  #write.table(s, file = paste('results/summary_significance/sum_',file_name,'.txt', sep = ''),sep = "\t", row.names = FALSE)
  #write(paste('\n Number of Hits :',dim(s)[1]), paste('results/summary_significance/sum_',file_name,'.txt', sep = ''), append = TRUE)
  #hits <- s  # output
}
