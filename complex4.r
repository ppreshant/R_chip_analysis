# removing flagged spots....

# reads an array and makes file of intensities - can apply threshold
# set unik = 1 to avg two spots side by side
# ptyr is wavelength 532 => Cy3

call_hits4 <- function(file_name)
{  # finding classes of each column to make it easier for reading the full table
  
  unik = 1 # should i merge unique entries?
  # should I remove bad spots?
  
  a <- read.table('results/pilot4/C_pilot4_9049232.gpr',skip= 31, header = TRUE,nrows = 4)
  x <- sapply(a,class)
  # Read in only reqired columns
  x[c(1:3,6:20,22:25,27:53,55:length(x))] <- 'NULL'  # we need rows 4 name, 5 ID,21 F median, 26 B median) of 532 - ptyr
  x[4] <- 'character'
  a <- read.table(paste('results/',file_name, sep = ''),skip= 31, header = TRUE,nrows = 46080,colClasses = x)
  
  # remove controls  : buffer, control, empty
  p <- grepl('empty',a[[2]]) | (grepl('buffer',a[[2]])) | (grepl('Control',a[[2]]))
  a <- a[!p,c(1,3:dim(a)[2])]
  
  
  # index unique spots
  #ix <- c(1:23040,1:23040)
  # index <- (1:46080 - 1)%/%2 + 1
  # a <- cbind(a,index)
  
  a <- cbind(a, 'F/B n' = a[[2]]/a[[3]]) # normalize foreground to background
  a <- cbind(a, 'F-B' = a[[2]]-a[[3]]) # foreground - Background 
  
  s <- a
  # look for block medians
  block_med <- 0
  for(i in 1:44)
  {
    tmp <- s[[1]] == i  # all blocks with number i
    st <- s[tmp,]
    block_med[i] <- median(st[[5]]) 
    s[5][tmp,] <- st[[5]] - block_med[i] + 1   # normalize : within chip norm. median Ij
  }
  # s
  # to select for desired controls and remove duplicates etc.
  # p <- !grepl('empty',s[[4]]) # grepl('anti',s[[4]]) | (grepl('IgG',s[[4]])) #(duplicated(s[4]) | duplicated(s[4],fromLast = TRUE)) & !(grepl('rabbit',s[[4]])) & !(grepl('IgG555',s[[4]]))  # have duplicated ones, remove rabbit-anti-biotin etc.
  # s <- s[p,]
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
    # remove bad spots
    p <- s[[5]] == 0
    s <- s[p, c(1:4,6:dim(s)[2])]
    
    s[[1]] <- as.character(s[[1]])
    s[c(3,4,6)] <-  round(s[c(3,4,6)])  # rounding off
    s[5] <- round(s[5],2)
  }
  s
  
  
  
  #write.table(s, file = paste('results/summary_significance/sum_',file_name,'.txt', sep = ''),sep = "\t", row.names = FALSE)
  #write(paste('\n Number of Hits :',dim(s)[1]), paste('results/summary_significance/sum_',file_name,'.txt', sep = ''), append = TRUE)
  #hits <- s  # output
}
