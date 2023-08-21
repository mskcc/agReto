library(agReto)
x <- read.table('data/input_data.txt',header=TRUE,sep='\t') # example input data
z0 <- wildtype(neoag=x)  # find WT peptide sequences
z <- z0[!is.na(z0$wt),]  # remove failed entries

dir <- '/juno/work/ccs/wooh/pkgs/netMHCpan-4.1/Linux_x86_64'  # main directory of netMHCpan

# single-peptide netMHC runs
b1 <- netmhc(peptide=z[1,'wt'], hla=hla[1,'hla'],netmhc.dir = dir)
b2 <- netmhc(peptide=z[2,'wt'], hla=hla[2,'hla'],netmhc.dir = dir)

# using a wrapper for multi-row data.frame input
input <- data.frame(peptide=z$wt, hla=z$hla)  # necessary to have two input columns as named
b <- run.netmhc(input, netmhc.dir=dir)
write.table(b, file='output.txt',row.names=F,col.names=T,quote=F,sep='\t')
