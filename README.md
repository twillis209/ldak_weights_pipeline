# A `snakemake` pipeline for generating LDAK weightings for use in `fcfdr`

## Software and R package dependencies

The pipeline requires that the following binaries be on the path:

* `ldak` (pipeline was developed using `Version 5.1`)
* `plink` (pipeline was developing using `PLINK v1.90b6.21 64-bit (19 Oct 2020)`

If you don't want to put these on the path, you can instead specify the full path to each binary in the shell invocations in the `snakefile`.

In addition, the R scripts depend on the following R packages:

* `data.table`
* `argparse`

## Things to beware

* It may be necessary to customise the arguments passed to the pipeline's R scripts. 
* 
