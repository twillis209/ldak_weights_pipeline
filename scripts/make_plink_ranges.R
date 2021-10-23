library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Produce PLINK SNP ID files to extract designated SNPs from PLINK files')
parser$add_argument('-i', '--gwas_file', type = 'character', help = 'Path to merged GWAS summary statistics file')
parser$add_argument('-b', '--bim_dir', type = 'character', help = 'Path to directory containing bim files')
parser$add_argument('-r', '--bim_regex', type = 'character', help = 'Regex for bim files. Must contain \'%d\' for autosome numbers.', default = 'chr%d.bim')
parser$add_argument('-chr', type = 'character', help = 'Label of chromosome column in GWAS file')
parser$add_argument('-bp', type = 'character', help = 'Label of BP column in GWAS file')
parser$add_argument('-ref', type = 'character', help = 'Label of reference allele column in GWAS file')
parser$add_argument('-alt', type = 'character', help = 'Label of alternative allele column in GWAS file')
parser$add_argument('-prin', type = 'character', help = 'Label of principal p-value column')
parser$add_argument('-aux', type = 'character', help = 'Label of auxiliary p-value column')
parser$add_argument('-o', '--output_dir', type = 'character', help = 'Path to output directory', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args<-parser$parse_args()

setDTthreads(threads=args$no_of_threads)

gwas_dat <- fread(args$gwas_file, sep = '\t', header = T, select = c(args$chr, args$bp, args$ref, args$alt, args$prin, args$aux))

gwas_dat[ , args$chr := as.character(get(args$chr))]

for(i in 1:22) {
  bim_dat <- fread(file.path(args$bim_dir, sprintf(args$bim_regex, i)), sep = '\t', header = F, col.names = c('CHR38', 'ID', 'Cm', 'BP38', 'A1', 'A2'))

  bim_dat[, CHR38 := as.character(CHR38)]

  bim_join <- merge(bim_dat, gwas_dat, by.x = c(args$chr, args$bp), by.y = c('CHR38', 'BP38'), sort = F)

  # Make sure alleles match
  bim_join <- bim_join[(REF == A1 & ALT == A2) | (REF == A2 & ALT == A1)]

  bim_join <- bim_join[, .(ID)]

  fwrite(bim_join, file = file.path(args$output_dir, sprintf("chr%d.txt", i)), row.names = F, sep = ' ', col.names = F, quote = F)
}
