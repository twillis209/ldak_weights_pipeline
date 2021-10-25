library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Add basepair column from bim file to combined MAF data')
parser$add_argument('-mf', '--maf_file', type = 'character', help = 'Path to MAF file')
parser$add_argument('-bf', '--bim_file', type = 'character', help = 'Path to combined bim file')
parser$add_argument('-chr_m', type = 'character', help = 'Label of chromosome column in MAF file')
parser$add_argument('-chr_b', type = 'character', help = 'Label of chromosome column in bim file')
parser$add_argument('-bp_b', type = 'character', help = 'Label of BP column in bim file')
parser$add_argument('-ref_m', type = 'character', help = 'Label of reference allele column in MAF file')
parser$add_argument('-alt_m', type = 'character', help = 'Label of alternative allele column in MAF file')
parser$add_argument('-ref_b', type = 'character', help = 'Label of reference allele column in bim file')
parser$add_argument('-alt_b', type = 'character', help = 'Label of alternative allele column in bim file')
parser$add_argument('-id_m', type = 'character', help = 'Label of ID column in MAF file')
parser$add_argument('-id_b', type = 'character', help = 'Label of ID column in bim file')
parser$add_argument('-maf', type = 'character', help = 'Label of MAF column in MAF file')
parser$add_argument('-of', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(args$no_of_threads)

maf_dat <- fread(args$maf_file, sep = ' ', header = T)
bim_dat <- fread(args$bim_file, sep = '\t', header = F, col.names = c(args$chr_b, args$id_b, 'Cm', args$bp_b, args$ref_b, args$alt_b))

bim_cols <- c(args$id_b, args$bp_b)

merged_dat <- merge(maf_dat, bim_dat[, ..bim_cols], by.x = args$id_m, by.y = args$id_b, sort = F)

merged_cols <- c(args$id_m, args$chr_m, args$bp_b, args$ref_m, args$alt_m, args$maf)
merged_dat <- merged_dat[, ..merged_cols]

fwrite(merged_dat, file = args$output_file, sep = ' ', col.names = T, row.names = F, quote = F)
