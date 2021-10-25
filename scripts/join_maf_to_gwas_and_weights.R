library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Join merged GWAS and weight file to MAF file')
parser$add_argument('-gf', '--gwas_file', type = 'character', help = 'Path to merged GWAS summary statistics file')
parser$add_argument('-mf', '--maf_file', type = 'character', help = 'Path to MAF file')
parser$add_argument('-chr_g', type = 'character', help = 'Label of chromosome column in GWAS file')
parser$add_argument('-chr_m', type = 'character', help = 'Label of chromosome column in MAF file')
parser$add_argument('-bp_g', type = 'character', help = 'Label of BP column in GWAS file')
parser$add_argument('-bp_m', type = 'character', help = 'Label of BP column in MAF file')
parser$add_argument('-ref_g', type = 'character', help = 'Label of reference allele column in GWAS file')
parser$add_argument('-ref_m', type = 'character', help = 'Label of reference allele column in MAF file')
parser$add_argument('-alt_g', type = 'character', help = 'Label of alternative allele column in GWAS file')
parser$add_argument('-alt_m', type = 'character', help = 'Label of alternative allele column in MAF file')
parser$add_argument('-maf', type = 'character', help = 'Label of MAF column in MAF file')
parser$add_argument('-o', '--output_path', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

gwas_dat <- fread(args$gwas_file, sep = '\t', header = T)

maf_dat <- fread(args$maf_file, sep = ' ', header = T)
maf_cols <- c(args$chr_m, args$bp_m, args$ref_m, args$alt_m, args$maf)

maf_dat <- maf_dat[, ..maf_cols]

merged_dat <- merge(gwas_dat, maf_dat, all.x = T, by.x = c(args$chr_g, args$bp_g), by.y = c(args$chr_m, args$bp_m), sort = F)

merged_dat <- merged_dat[( get(args$ref_g) == get(args$ref_m) & get(args$alt_g) == get(args$alt_m) ) | ( get(args$ref_g) == get(args$alt_m) & get(args$alt_g) == get(args$ref_m) )]

merged_dat[, c(args$ref_m, args$alt_m) := NULL]

fwrite(merged_dat, file = args$output_path, sep = '\t', col.names = T, row.names = F, quote = F)
