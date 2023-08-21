#' Finds the wild-type sequence of a peptide
#' @param neoag Data.frame with columns \code{peptide}, \code{gene}, and \code{mutation}
#' @param z List of known sequences for all genes to be considered
#' @export
wildtype <- function(neoag, z = NA){

  if(is.na(z)) z <- system.file(package = 'agReto','extdata','z.rds')
  if(is.character(z)){
    if(!file.exists(z)) stop(paste0(z, ' does not exist'))
    z <- readRDS(z)
  }

  cols <- c('peptide','gene','mutation')
  if(any(!cols %in% colnames(neoag)))
    stop(paste0('Input column(s) in neoag missing'))

  m <- NROW(neoag)
  all.genes <- names(z)

  wt <- lapply(seq(m), FUN = function(i){
    igene <- neoag[i, 'gene']
    if(!igene %in% all.genes) return(NA)
    prot <- z[[igene]]
    x <- wt.pept(prot = prot, hgsp = neoag[i, 'mutation'],
            neoag = neoag[i, 'peptide'])
    return(x)
  })
  wt <- unlist(wt)
  v <- cbind(neoag, data.frame(wt=wt))

  x <- v[match(neoag$peptide, v$peptide), ]
  return(x)
}

#' Finds wild-type peptide from protein
wt.pept <- function(prot, hgsp, neoag){

  m <- nchar(neoag)
  mut <- hgsp
  smut <- strsplit(mut,split='')[[1]]
  zmut <- smut[!is.na(as.integer(smut))]
  imut <- as.integer(paste0(zmut,collapse=''))
  wta <- substr(mut,start=3,stop=which(smut==zmut[1])-1)

  pos <- imut
  flag <- FALSE
  for(i in seq(NROW(prot))){
    p1 <- prot[i,1]
    if(!is.character(p1)) next()
    a1 <- strsplit(p1,split='')[[1]]
    if(length(a1)<pos) next()
    a1pos <- paste(a1[seq(pos,pos+nchar(wta)-1)],collapse='')
    if(a1pos==wta & pos + m - 1 <= nchar(p1)){
      flag <- TRUE
      break()
    }
  }
  if(!flag){
    warning(paste0('Mutated AA in ',mut,'does not match peptide sequences'))
    return(NA)
  }

  aneo <- strsplit(neoag,split='')[[1]]
  dist <- rep(0,m)
  for(i in seq(m)){
    start <- max(c(1,pos-i+nchar(wta)))
    stop <- min(c(pos-i+m+nchar(wta)-1,length(a1)))
    pept <- a1[seq(start,stop)]
    aa <- paste(pept,collapse='')
    dista <- sum(pept!=aneo)
    dist[i] <- dista
    names(dist)[i] <- aa
  }
  wt <- names(which.min(dist))
  return(wt)
}
