#' Call netMHC for a peptide
#' @param peptide Peptide sequence
#' @param hla HLA type
#' @export
netmhc <- function(peptide, netmhc.dir, hla){

  flag <- FALSE
  if(!dir.exists('tmp')){
    flag <- TRUE
    dir.create('tmp')
  }
  if(!dir.exists(netmhc.dir)) stop(paste0(netmhc.dir, ' does not exist'))

  fasta <- data.frame(peptide)
  colnames(fasta) <- '> pept'
  write.table(fasta, file = 'tmp/x.fasta', row.names=F, col.names=T,quote=F,sep='\t')
  header <- paste0('#!/bin/bash\nexport NETMHCpan=',netmhc.dir,'\n')
  write(header, file='tmp/cmd.sh')
  hla <- gsub('[*]','',hla)
  cmd <- paste0(netmhc.dir,
                 '/bin/netMHCpan -s -l 9,10,11 -f tmp/x.fasta -inptype 0 -xls 1 -BA 1 -a ',
                hla, ' >> tmp/output')
  write(cmd, file='tmp/cmd.sh', append=TRUE)
  system('bash tmp/cmd.sh')
  b1 <- readLines('tmp/output')
  system('rm -rf ./tmp/output ./tmp/x..fasta')
  if(flag) system('rmdir tmp')

  idx <- grep('---------',b1)
  re <- b1[seq(idx[2]+1,idx[3]-1)]
  res <- strsplit(re,split=' ')[[1]]
  res <- res[res!='']
  res <- res[res!='<=']

  col <- b1[idx[1]+1]
  colu <- strsplit(col,split=' ')[[1]]
  colum <- colu[colu!='']

  if(length(res)==16) res <- c(res,'')
  dat <- data.frame(t(res))

  colnames(dat) <- colum
  dat <- dat[,colnames(dat)!='Pos']

  return(dat)
}
