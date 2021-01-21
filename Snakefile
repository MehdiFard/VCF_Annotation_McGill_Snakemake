import os
configfile:"config.yaml"
VCF_DIR=config["vcf_directory"]
VCF_EXT=config["vcf_extension"]
FILES={f[:-len(VCF_EXT)] for f in os.listdir(VCF_DIR)}


rule all:
    input:
        "intermediate/vcfnames.list",
        "intermediate/qcthr.csv",
        "intermediate/merged.vcf.gz",
        "intermediate/merged_qcok.vcf.gz",
        "intermediate/filter.king.cutoff.in.id",
        "intermediate/filter.king.cutoff.out.id",
        "intermediate/merged_qc_unrel_kin.kin0",
        "results/list_of_unrelated_individuals.txt",
        "intermediate/merged_clean.vcf.gz",
        "intermediate/merged_clean_annotated.vcf.bgz",
        "intermediate/merged_clean_annotated_filtered.vcf.gz",
        "intermediate/merged_clean_annotated_filtered.vcf.gz.tbi",
        expand("results/{file}.out.vcf.gz", file=FILES)


rule parsing_metadata:
    input:
        info="data/info_files/metadata.csv"
    output:
        vcflist="intermediate/vcfnames.list",
        thrs="intermediate/qcthr.csv"
    conda:
        "envs/environment.yaml",
    log:
        "log/parsing.log"
    shell:
        """
        cut -d"," -f1 {input.info} | sed '1d' | ts {VCF_DIR}/ | tr -d ' ' > {output.vcflist} 2> {log} # VCF file names.
        cut -d"," -f1 --complement {input.info} | awk 'NR==1 || NR==2' > {output.thrs} 2>> {log}  # Thresholds
        """


rule mergevcfs:
    input:
        "intermediate/vcfnames.list",
    output:
        mvcf="intermediate/merged.vcf.gz",
        mvcf_inx="intermediate/merged.vcf.gz.tbi"
    conda:
        "envs/environment.yaml",
    log:
        "log/merge.log"
    shell:
        "PicardCommandLine MergeVcfs I={input} O={output.mvcf} 2> {log}"


rule qc:
    input:
        thr="intermediate/qcthr.csv",
        mvcfs="intermediate/merged.vcf.gz"
    output:
        mvcfsqc="intermediate/merged_qcok.vcf.gz"
    conda:
        "envs/environment.yaml",
    log:
        "log/qc.log"
    shell:
        """
        bash scr/thrs_to_variables.sh {input.thr} {input.mvcfs} {output.mvcfsqc} 2> {log}
        """

rule identify_related_inds:
    input:
        "intermediate/merged_qcok.vcf.gz"
    output:
        out_in="intermediate/filter.king.cutoff.in.id",
        out_out="intermediate/filter.king.cutoff.out.id",
        out_kin="intermediate/merged_qc_unrel_kin.kin0",
        out_list="results/list_of_unrelated_individuals.txt"
    conda:
        "envs/environment.yaml"
    log:
        "log/unrelated.log"
    shell:
        """
        plink2 --vcf {input} --king-cutoff 0.0442 --out intermediate/filter &> {log}
        grep -v "IID" {output.out_in} > {output.out_list}
        echo "# Check remained individuals for kinship ..."
        plink2 --vcf {input} --remove {output.out_out} --make-king-table --out intermediate/merged_qc_unrel_kin  &> {log} # Making pairwase kinship table with unrelated individuals
        Rscript scr/kinship_check.R
        """

rule remove_related_inds:
    input:
        list="results/list_of_unrelated_individuals.txt",
        mvcfqc="intermediate/merged_qcok.vcf.gz"
    output:
        clean="intermediate/merged_clean.vcf.gz",
    conda:
        "envs/environment.yaml"
    log:
        "log/filter.log"
    shell:
        """
        bcftools view -S {input.list} {input.mvcfqc} | bgzip -c > {output.clean} 2> {log}  # Removing related individuals
        """

rule annotation_Rscript:
    input:
        "intermediate/merged_clean.vcf.gz"
    output:
        "intermediate/merged_clean_annotated.vcf.bgz"
    conda:
        "envs/environment.yaml"
    log:
        "log/annotation.log"
    shell:
        """
        Rscript scr/info_field_cal_opt.R 2> {log}
        """


rule filter_AF:
    input:
        "intermediate/merged_clean_annotated.vcf.bgz"
    output:
        vcf="intermediate/merged_clean_annotated_filtered.vcf.gz",
        vcf_ind="intermediate/merged_clean_annotated_filtered.vcf.gz.tbi"
    conda:
        "envs/environment.yaml"
    log:
        "log/filter_AF.log"
    shell:
        """
        bcftools view -e "AF=0 | AF=1"  {input} | bgzip -c > {output.vcf} # Filtering sites with AF=1 and AF=0 (keeping sitis with 0<AF<1)
        tabix -f -p vcf {output.vcf}
        """


rule split_reheader:
    input:
        "intermediate/merged_clean_annotated_filtered.vcf.gz"
    output:
        "results/{file}.out.vcf.gz"
    conda:
        "envs/environment.yaml"
    log:
        "log/chr{file}_split_reheader.log"
    shell:
        """
        bcftools view {input} --regions chr{wildcards.file} | sed -e "/ID=chr{wildcards.file}/p" -e "/contig/d" -e "/##bcftools/d" | bgzip -c > results/{wildcards.file}.out.vcf.gz 2> {log}
        """

### END
