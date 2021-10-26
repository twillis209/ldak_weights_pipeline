configfile: 'config.yaml'

wildcard_constraints:
  chr = "chr[0-9X]{1,2}"

def get_mem_mb(wildcards, threads):
    return threads * 3420

CHROMS = ["chr%d" % i for i in range(1,23)]

if not(config["reference"] in ["hg38", "hg19"]):
    raise Exception("Invalid \"reference\" value: %s" % config["reference"])

if config["reference"] == "hg38":
    kg_bp_label = "BP38"
else:
    kg_bp_label = "BP19"

rule download_1000g_genotype_data:
     output:
      "resources/1000g/{chr}.vcf.gz"
     run:
        if config["reference"] == "hg38":
            if wildcards.chr == 'chrX':
                shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz")
            else:
                shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz")
        else:
            if wildcards.chr == 'chrX':
                shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz")
            else:
                shell("wget -O resources/1000g/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

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

rule combine_qc_bim_files:
    input:
        ["resources/1000g/euro/qc/%s_qc.bim" % x for x in CHROMS]
    output:
        "resources/1000g/euro/qc/chr_all.bim"
    shell:
        """
        for x in {input}; do cat $x >>{output}; done
        """

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

rule combine_maf_files:
    input:
        ["resources/1000g/euro/qc/maf/%s.frq" % x for x in CHROMS]
    output:
        "resources/1000g/euro/qc/maf/chr_all.frq"
    shell:
        """
        echo -e "CHR SNP A1 A2 MAF NCHROBS" >>{output}
        for x in {input}; do tail -n +2 $x >>{output}; done
        """

rule add_bp_to_maf_file:
    input:
        maf_file = "resources/1000g/euro/qc/maf/chr_all.frq",
        bim_file = "resources/1000g/euro/qc/chr_all.bim"
    output:
        "resources/1000g/euro/qc/maf/chr_all_bp.frq"
    threads: 4
    resources:
        mem_mb=get_mem_mb
    shell:
        f"Rscript scripts/add_bp_to_maf_file.R -mf {{input.maf_file}} -bf {{input.bim_file}} -chr_m CHR -chr_b Chr -bp_b {kg_bp_label} -ref_m A1 -alt_m A2 -ref_b A1 -alt_b A2 -id_m SNP -id_b ID -maf MAF -of {{output}} -nt {{threads}}"

rule join_gwas:
     input:
      A = "resources/gwas/{imd_a}.tsv.gz",
      B = "resources/gwas/{imd_b}.tsv.gz"
     output:
      AB = "resources/gwas/{imd_a}_{imd_b}/{imd_a}_{imd_b}.tsv.gz"
     threads: 8
     shell:
         f"Rscript scripts/join_gwas_stats.R -a {{input.A}} -b {{input.B}} -chr_a {config[gwas_chr_label]} -chr_b {config[gwas_chr_label]} -bp_a {config[gwas_bp_label]} -bp_b {config[gwas_bp_label]} -ref_a {config[gwas_ref_label]} -ref_b {config[gwas_ref_label]} -alt_a {config[gwas_alt_label]} -alt_b {config[gwas_alt_label]} -p_a {config[gwas_pvalue_label]} -p_b {config[gwas_pvalue_label]} -o {{output.AB}} -nt {{threads}}"

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
      f"Rscript scripts/make_plink_ranges.R -i {{input.gwas_file}} -b {{params.input_dir}} -r chr%d_qc.bim -chr {config[gwas_chr_label]} -bp {config[gwas_bp_label]} -ref {config[gwas_ref_label]} -alt {config[gwas_alt_label]} -prin {config[gwas_prin_pvalue_label]} -aux {config[gwas_aux_pvalue_label]} -o {{params.output_dir}} -nt {{threads}}"

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

rule combine_subsetted_bim_files:
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

rule join_gwas_and_weights:
    input:
     weights_file = "results/{imd_a}_{imd_b}/ldak/combined_weights_meta.all",
     merged_gwas_file = "resources/gwas/{imd_a}_{imd_b}/{imd_a}_{imd_b}.tsv.gz"
    output:
     "results/{imd_a}_{imd_b}/gwas/{imd_a}_{imd_b}_with_weights.tsv.gz"
    threads: 4
    shell:
        f"Rscript scripts/join_gwas_and_weights.R -gf {{input.merged_gwas_file}} -wf {{input.weights_file}} -chr_g {config[gwas_chr_label]} -chr_w Chr -bp_g {config[gwas_chr_label]} -bp_w {kg_bp_label} -ref_g {config[gwas_ref_label]} -ref_w A1 -alt_g {config[gwas_alt_label]} -alt_w A2 -of {{output}} -nt {{threads}}"

rule join_maf_to_gwas_and_weights:
    input:
     merged_gwas_file = "results/{imd_a}_{imd_b}/gwas/{imd_a}_{imd_b}_with_weights.tsv.gz",
     maf_file = "resources/1000g/euro/qc/maf/chr_all_bp.frq"
    output:
     "results/{imd_a}_{imd_b}/gwas/{imd_a}_{imd_b}_with_weights_and_maf.tsv.gz"
    threads: 4
    shell:
        f"Rscript scripts/join_maf_to_gwas_and_weights.R -gf {{input.merged_gwas_file}} -mf {{input.maf_file}} -chr_g {config[gwas_chr_label]} -chr_m CHR -bp_g {config[gwas_bp_label]} -bp_m {kg_bp_label} -ref_g {config[gwas_ref_label]} -ref_m A1 -alt_g {config[gwas_alt_label]} -alt_m A2 -maf MAF -o {{output}} -nt {{threads}}"
