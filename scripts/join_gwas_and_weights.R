library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Join LDAK weightings file to merged GWAS summary statistics')
parser$add_argument('-gf', '--gwas_file', type = 'character', help = 'Path to merged GWAS summary statistics file')
parser$add_argument('-wf', '--weights_file', type = 'character', help = 'Path to LDAK weightings file')
parser$add_argument('-chr_g', type = 'character', help = 'Label of chromosome column in GWAS file')
parser$add_argument('-chr_w', type = 'character', help = 'Label of chromosome column in weights file')
parser$add_argument('-bp_g', type = 'character', help = 'Label of BP column in GWAS file')
parser$add_argument('-bp_w', type = 'character', help = 'Label of BP column in weights file')
parser$add_argument('-ref_g', type = 'character', help = 'Label of reference allele column in GWAS file')
parser$add_argument('-ref_w', type = 'character', help = 'Label of reference allele column in weights file')
parser$add_argument('-alt_g', type = 'character', help = 'Label of alternative allele column in GWAS file')
parser$add_argument('-alt_w', type = 'character', help = 'Label of alternative allele column in weigts file')
parser$add_argument('-of', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

gwas_dat <- fread(args$gwas_file, sep = '\t', header = T)
gwas_dat[, (args$chr_g) := as.character(get(args$chr_g))]

weights_dat <- fread(args$weights_file, sep = ' ', header = T)
weights_cols <- c('Weight', args$chr_w, args$bp_w, args$ref_w, args$alt_w)
weights_dat[, (args$chr_w) := as.character(get(args$chr_w))]

# Weights ref/alt pair is actually arbitrary so we check both
merged_dat <- merge(gwas_dat, weights_dat[, ..weights_cols], by.x = c(args$chr_g, args$bp_g), by.y = c(args$chr_w, args$bp_w), all.x = T, sort = F)

merged_dat <- merged_dat[( get(args$ref_g) == get(args$ref_w) & get(args$alt_g) == get(args$alt_w) ) | ( get(args$ref_g) == get(args$alt_w) & get(args$alt_g) == get(args$ref_w) )]

merged_dat[, c(args$ref_w, args$alt_w) := NULL]

fwrite(merged_dat, file = args$output_file, sep = '\t', col.names = T, row.names = F, quote = F)
