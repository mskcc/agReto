#' Wrapper Script for netMHCpan
#' @param input Data frame with columns \code{peptide}, \code{hla}
#' @param netmhc.dir netMHCpan source directory
#' @param progress.bar Display progress bar
#' @export
run.netmhc <- function(input, netmhc.dir, progress.bar = TRUE){

  cols <- c('peptide','hla')
  if(any(!cols %in% colnames(input)))
    stop('Some required input columns missing')
  if(!dir.exists(netmhc.dir))
    stop(paste0(netmhc.dir, ' does not exist'))

  m <- NROW(input)
  if(progress.bar) pb <- txtProgressBar(style=3)
  b <- do.call(rbind, lapply(seq(m), function(i){
    if(is.na(input[i,'peptide']) | is.na(input[i,'hla']))
      return(NULL)
    y <- netmhc(peptide = input[i, 'peptide'], hla = input[i, 'hla'],
                netmhc.dir = netmhc.dir)
    if(progress.bar) setTxtProgressBar(pb, i/m)
    return(y)
  }))
  if(progress.bar) close(pb)

  return(b)
}
