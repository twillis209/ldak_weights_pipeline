wildcard_constraints:
  chr = "chr[0-9X]{1,2}"

def get_mem_mb(wildcards, threads):
    return threads * 3420

CHROMS = ["chr%d" % i for i in range(1,23)]

rule download_1000g_genotype_data:
     output:
      "resources/1000g/{chr}.vcf.gz"
     run:
      if wildcards.chr == 'chrX':
         shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz")
      else:
        shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz")

rule download_1000g_sample_metadata:
     output:
      "resources/1000g/20130606_g1k_3202_samples_ped_population.txt"
     shell:
      "wget -O resources/1000g/20130606_g1k_3202_samples_ped_population.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

rule vcf_to_bed:
     input:
      "resources/1000g/{chr}.vcf.gz"
     output:
      "resources/1000g/{chr}.bed",
      "resources/1000g/{chr}.bim",
      "resources/1000g/{chr}.fam"
     threads: 8
     resources:
         mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --vcf resources/1000g/{wildcards.chr}.vcf.gz --make-bed --out resources/1000g/{wildcards.chr}"

rule make_euro_fam:
     input:
      "resources/1000g/20130606_g1k_3202_samples_ped_population.txt"
     output:
      "resources/1000g/euro.fam"
     shell:
      "Rscript scripts/get_euro_fam_file.R"

rule get_euro_samples:
     input:
      "resources/1000g/{chr}.bed",
      "resources/1000g/{chr}.bim",
      "resources/1000g/{chr}.fam",
      "resources/1000g/euro.fam"
     output:
      "resources/1000g/euro/{chr}_euro.bed",
      "resources/1000g/euro/{chr}_euro.bim",
      "resources/1000g/euro/{chr}_euro.fam"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
         "plink --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/{wildcards.chr} --keep resources/1000g/euro.fam --make-bed --silent --out resources/1000g/euro/{wildcards.chr}_euro"

rule qc:
     input:
      "resources/1000g/euro/{chr}_euro.bed",
      "resources/1000g/euro/{chr}_euro.bim",
      "resources/1000g/euro/{chr}_euro.fam"
     output:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam"
     threads: 8
     resources:
       mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/euro/{wildcards.chr}_euro --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --make-bed --silent --out resources/1000g/euro/qc/{wildcards.chr}_qc"

rule maf:
    input:
        "resources/1000g/euro/qc/{chr}_qc.bed",
        "resources/1000g/euro/qc/{chr}_qc.bim",
        "resources/1000g/euro/qc/{chr}_qc.fam"
    output:
        "resources/1000g/euro/qc/maf/{chr}.frq"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    shell:
        "plink --memory {resources.mem_mb} --threads {threads} --bfile resources/1000g/euro/qc/{wildcards.chr}_qc --freq --out resources/1000g/euro/qc/maf/{wildcards.chr}"

rule join_gwas:
     input:
      A = "resources/gwas/{imd_a}.tsv.gz",
      B = "resources/gwas/{imd_b}.tsv.gz"
     output:
      AB = "resources/gwas/{imd_a}_{imd_b}/{imd_a}_{imd_b}.tsv.gz"
     threads: 8
     shell:
      "Rscript scripts/join_gwas_stats.R -a {input.A} -b {input.B} -chr_a CHR38 -chr_b CHR38 -bp_a BP38 -bp_b BP38 -ref_a REF -ref_b REF -alt_a ALT -alt_b ALT -p_a P -p_b P -o {output.AB} -nt {threads}"

rule make_plink_ranges:
     input:
      ("resources/1000g/euro/qc/chr%d_qc.bim" % x for x in range(1,23)),
      gwas_file = "resources/gwas/{imd_a}_{imd_b}/{imd_a}_{imd_b}.tsv.gz"
     output:
      ("resources/gwas/{imd_a}_{imd_b}/matching_ids/chr%d.txt" % x for x in range(1,23))
     params:
      input_dir = "resources/1000g/euro/qc",
      output_dir = "resources/gwas/{imd_a}_{imd_b}/matching_ids"
     threads: 2
     shell:
      "Rscript scripts/make_plink_ranges.R -i {input.gwas_file} -b {params.input_dir} -r chr%d_qc.bim -chr CHR38 -bp BP38 -ref REF -alt ALT -prin P.A -aux P.B -o {params.output_dir} -nt {threads}"

rule subset_reference:
     input:
      "resources/1000g/euro/qc/{chr}_qc.bed",
      "resources/1000g/euro/qc/{chr}_qc.bim",
      "resources/1000g/euro/qc/{chr}_qc.fam",
      range_file = "resources/gwas/{imd_a}_{imd_b}/matching_ids/{chr}.txt"
     output:
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.bed",
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.bim",
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.fam"
     params:
      bfile = "resources/1000g/euro/qc/{chr}_qc",
      out = "resources/gwas/{imd_a}_{imd_b}/plink/{chr}"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     shell:
      "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --extract {input.range_file} --make-bed --out {params.out}"

rule cut_weights:
     input:
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.bed",
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.bim",
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.fam"
     output:
      predictors = "resources/ldak/{imd_a}_{imd_b}/{chr}/weights.predictors",
      log = "resources/ldak/{imd_a}_{imd_b}/{chr}/cut_weights.log"
     params:
      input_dir = "resources/gwas/{imd_a}_{imd_b}/plink/{chr}",
      output_dir = "resources/ldak/{imd_a}_{imd_b}/{chr}"
     shell:
      "ldak --cut-weights {params.output_dir} --bfile {params.input_dir} > {output.log}"

rule calc_weights_all:
     input:
      "resources/ldak/{imd_a}_{imd_b}/{chr}/weights.predictors",
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.bed",
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.bim",
      "resources/gwas/{imd_a}_{imd_b}/plink/{chr}.fam"
     output:
      weights = "resources/ldak/{imd_a}_{imd_b}/{chr}/weights.all",
      log = "resources/ldak/{imd_a}_{imd_b}/{chr}/calc_weights_all.log"
     params:
      input_dir = "resources/gwas/{imd_a}_{imd_b}/plink/{chr}",
      output_dir = "resources/ldak/{imd_a}_{imd_b}/{chr}"
     shell:
      "ldak --calc-weights-all {params.output_dir} --bfile {params.input_dir} > {output.log}"

rule combine_weights:
     input:
      ["resources/ldak/{imd_a}_{imd_b}/%s/weights.all" % x for x in CHROMS]
     output:
      "results/{imd_a}_{imd_b}/ldak/combined_weights.all"
     shell:
      """
      echo -e "Predictor Weight Neighbours Tagging Info Check" >>{output}
      for x in {input}; do tail -n +2 $x >>{output}; done
      """

rule combine_bim_files:
     input:
      ["resources/gwas/{imd_a}_{imd_b}/plink/%s.bim" % x for x in CHROMS]
     output:
      "resources/gwas/{imd_a}_{imd_b}/plink/chr_all.bim"
     shell:
      "for x in {input}; do cat $x >>{output}; done"

rule add_bp_alleles_to_combined_weights:
     input:
      bim_file = "resources/gwas/{imd_a}_{imd_b}/plink/chr_all.bim",
      weights_file = "results/{imd_a}_{imd_b}/ldak/combined_weights.all"
     output:
      "results/{imd_a}_{imd_b}/ldak/combined_weights_meta.all"
     threads: 4
     shell:
      "Rscript scripts/add_bp_alleles_to_combined_weights.R -bf {input.bim_file} -wf {input.weights_file} -of {output} -nt {threads}"
