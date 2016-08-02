pull_nd <- function(el = NULL)
{
  root <- 'https://collection.cdi-lab.com/public/home/search?utf8=%E2%9C%93&key=symbol&value='
  
  u <- paste(root,'&block=', el$Block, '&row=', el$Row, '&column=', el$Column, '&version=', 3, '&huprot=Search+HuProt+Array', sep = '')
  url <- URLencode(u)
  su <- readLines(url)
  data = as.character(readHTMLTable(su)[[1]]$Symbol);
}

all_pull <- function(ndfl = ndf)
{
  # sapply(ndfl, pull_nd)
  ou <- 0
  for (i in 1:dim(ndfl)[1])
  {
    ou[i] <- pull_nd(ndfl[i,])
  }
  ou
}

rem_ND <- function(file_name)
{
  a <- read.table(paste('results/F/with old settings/with NDs/',file_name, sep = ''),nrows = 31, sep = '?')
  b <- read.table(paste('results/F/with old settings/with NDs/',file_name, sep = ''), skip= 31, header = TRUE, nrows = 50000, sep = '\t')
  b$Name <- xgal$Name
  write.table(a, paste('results/F/',file_name, sep = ''), row.names = F, col.names = F)
  write.table(b, paste('results/F/',file_name, sep = ''), append = T, sep = '\t', row.names = F)
}