library(data.table)
library(argparse)

parser <- ArgumentParser(description = 'Add alleles to output weights files')
parser$add_argument('-bf', '--bim_file', type = 'character', help = 'Path to bim file', required = T)
parser$add_argument('-wf', '--weights_file', type = 'character', help = 'Path to weights file', required = T)
parser$add_argument('-of', '--output_file', type = 'character', help = 'Path to output file', required = T)
parser$add_argument('-nt', '--no_of_threads', type = 'integer', help = 'Number of threads to use', default = 1)

args <- parser$parse_args()

setDTthreads(threads=args$no_of_threads)

bim_dat <- fread(args$bim_file, sep = '\t', header = F, col.names = c('Chr', 'ID', 'Cm', 'BP19', 'A1', 'A2'))

weights_dat <- fread(args$weights_file, sep = ' ', header = T)

# Drop rows with NA values
weights_dat <- na.omit(weights_dat, cols = c('Predictor', 'Weight'))

join_dat <- merge(weights_dat, bim_dat[, .(ID, Chr, BP19, A1, A2)], all.x = T, by.x = c('Predictor', 'Chr'), by.y = c('ID', 'Chr'), sort = F)

# Check that we have not lost any SNPs
if(any(is.na(join_dat$BP19))) {
  stop("Failed to match at least one SNP in weights file for %s")
}

fwrite(join_dat, file = args$output_file, sep = ' ', col.names = T, row.names = F, quote = F)
