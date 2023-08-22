# agReto
This package contains utilities for finding wild-type sequences of neoantigens and running netMHCpan on them. The affinity data obtained can be used to estimate agretopicity (ratio of mutant vs. wild-type peptide binding affinities).

## Regenerating wild-type sequences of neoantigens

The function [wildtype](man/wildtype.Rd) attemps to regenerate wild-type sequences of a set of neoantigens solely based on neoantigen sequences along with the gene name and identity of mutations (expressed in `HGVSp_short` notation; hg19 genome assumed). The input file is taken into the argument `neoag`, where the expected column names are `peptide`, `gene`, and `mutation`:

| peptide     | gene  | mutation | hla         |
| ----------- | ----- | -------- | ----------- |
| DYNWLIYHL   | DHX36 | p.H882D  | HLA-A*24:02 |
| YGLFVTRAL   | EIF3H | p.S127L  | HLA-C*12:03 |
| ...         | ...   | ...      | ...         |

An example file can be seen in [input_data.txt](example/data/input_data.txt). The column `hla` above is not necssary for this call but will just be carried over into the output. The other argument to `wildtype` is a list of all known protein sequences with each element named by gene name (Hugo symbol). It defaults to an object shipped within the package [z.rds](inst/extdata/z.rds). Any queries for genes not in this list will return `NA`. 

A call like the following:
```
library(agReto)
x <- read.table('input_data.txt', header = TRUE, sep = '\t')
z <- wildtype(neoag = x)
head(z)
```

would produce

| peptide     | gene  | mutation | hla         |     wt    |     
| ----------- | ----- | -------- | ----------- | --------- |
| DYNWLIYHL   | DHX36 | p.H882D  | HLA-A*24:02 | HYNWLIYHL | 
| YGLFVTRAL   | EIF3H | p.S127L  | HLA-C*12:03 | YGSFVTRAL |
| ...         | ...   | ...      | ...         | ...       |

Note that any mutations that are non-SNVs would also generate `NA`s. It is best to filter them out first.

## Binding affinity prediction for peptides

For convenience in estimating agretopicity (relative affinity of neoantigens with respect to their wild-type counterparts), a function [netmhc](man/netmhc.Rd) for calling [netMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) on a single peptide (wild-type in this context but could be arbitrary peptides) and a wrapper [run.netmhc](man/run.netmhc.Rd) for multiple peptides are provided. 

The single-peptide function can be called, e.g., by

```
dir <- '/juno/work/ccs/wooh/pkgs/netMHCpan-4.1/Linux_x86_64'
b1 <- netmhc(peptide = 'HYNWLIYHL', hla = 'HLA-A*24:02', netmhc.dir = dir)
b1
```
The path `netmhc.dir` is the main directory path of the [netMHCpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/). The path shown above works on `juno` (needs to be linux).

The output `b1` is a formated data frame of the standard netMHCpan output.

To call the wrapper, put together a data frame with two columns `peptide` and `hla`:

| peptide   | hla         |
| --------- | ----------- |
| HYNWLIYHL | HLA-A*24:02 |
| YGSFVTRAL | HLA-C*12:03 |
| ...       | ...         |

and call it by

```
b <- run.netmhc(input, netmhc.dir = dir)
```

The output `b` is a multi-row data frame version of `b1` above.

See [example](example/run.R) for further information.
