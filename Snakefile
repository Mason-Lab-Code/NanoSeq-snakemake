configfile: "config.yaml"

rule all:
    input:
        qc="03_qc/multiqc_report.html",
        analysis=expand("10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}/", sample=config["SAMPLES"], duplex_type=config["DUPLEX_TYPES"], undiluted_type=config["UNDILUTED_TYPES"]),
        contamination=expand("11_contamination_check/{sample}_{type}.out", sample=config["SAMPLES"], type=config["TYPES"]),
        efficiency=expand("12_efficiency_estimate/{sample}_{type}.RBs", sample=config["SAMPLES"], type=config["TYPES"])

rule fastqc:
    input:
        "00_raw/{sample}_{type}_{readnum}.fastq.gz",
    output:
        directory("03_qc/{sample}_{type}_{readnum}_fastqc")
    params:
        threads=8
    resources:
        runtime=60,
        mem_mb=2000,
        cpus_per_task=9
    shell:
        r"""
        module load bio/FastQC

        mkdir {output}
        fastqc -t {params.threads} -o {output} {input}
        """

rule multiqc:
    input:
        expand("03_qc/{sample}_{type}_{readnum}_fastqc", sample=config["SAMPLES"], type=config["TYPES"], readnum=["read1", "read2"])
    output:
        "03_qc/multiqc_report.html"
    resources:
        runtime=30,
        cpus_per_task=9
    shell:
        r"""
        module purge
        module load bio/MultiQC

        multiqc 03_qc/*_fastqc/ -o 03_qc/
        """

rule extract_tags:
    input:
        fq1="00_raw/{sample}_{type}_read1.fastq.gz",
        fq2="00_raw/{sample}_{type}_read2.fastq.gz"
    output:
        fq1="04_fastq_extract_tags/{sample}_{type}_read1.extract_tags.fastq.gz",
        fq2="04_fastq_extract_tags/{sample}_{type}_read2.extract_tags.fastq.gz"
    params:
        read_length=150
    resources:
        runtime=300,
        mem_mb=1000,
        cpus_per_task=1
    shell:
        r"""
        module purge
        module load bio/NanoSeq

        extract_tags.py -a {input.fq1} -b {input.fq2} -c {output.fq1} -d {output.fq2} -m 3 -s 4 -l {params.read_length}
        """

rule bwa_index:
    input:
        "01_ref/GRCh38.primary_assembly.genome.fa"
    output:
        amb="01_ref/GRCh38.primary_assembly.genome.fa.amb",
        ann="01_ref/GRCh38.primary_assembly.genome.fa.ann",
        bwt="01_ref/GRCh38.primary_assembly.genome.fa.bwt",
        fai="01_ref/GRCh38.primary_assembly.genome.fa.fai",
        pac="01_ref/GRCh38.primary_assembly.genome.fa.pac",
        sa="01_ref/GRCh38.primary_assembly.genome.fa.sa"
    resources:
        runtime=360,
        mem_mb=32000,
        cpus_per_task=8
    shell:
        r"""
        module load bio/BWA

        bwa index {input}
        """

rule bwa_mem:
    input:
        fq1="04_fastq_extract_tags/{sample}_{type}_read1.extract_tags.fastq.gz",
        fq2="04_fastq_extract_tags/{sample}_{type}_read2.extract_tags.fastq.gz",
        fa="01_ref/GRCh38.primary_assembly.genome.fa",
        amb="01_ref/GRCh38.primary_assembly.genome.fa.amb",
        ann="01_ref/GRCh38.primary_assembly.genome.fa.ann",
        bwt="01_ref/GRCh38.primary_assembly.genome.fa.bwt",
        fai="01_ref/GRCh38.primary_assembly.genome.fa.fai",
        pac="01_ref/GRCh38.primary_assembly.genome.fa.pac",
        sa="01_ref/GRCh38.primary_assembly.genome.fa.sa"
    output:
        temp("05_sam/{sample}_{type}.sam")
    params:
        threads=8
    resources:
        runtime=2880, # 48h
        mem_mb=48000,
        cpus_per_task=8
    shell:
        r"""
        module load bio/BWA

        bwa mem -C -t {params.threads} \
          {input.fa} \
          {input.fq1} \
          {input.fq2} > {output}
        """

rule sort_and_add_rc_mc_tags:
    input:
        "05_sam/{sample}_{type}.sam"
    output:
        "06_bam_rc_mc_tags/{sample}_{type}.rc_mc_tags.bam"
    params:
        threads=4
    resources:
        runtime=540,
        mem_mb=24000,
        cpus_per_task=4
    shell:
        r"""
        module purge
        module load bio/NanoSeq

        bamsormadup inputformat=sam rcsupport=1 threads={params.threads} < {input} > {output}
        """

rule mark_dups:
    input:
        "06_bam_rc_mc_tags/{sample}_{type}.rc_mc_tags.bam"
    output:
        bam="07_bam_mark_dups/{sample}_{type}.mark_dups.bam",
        bai="07_bam_mark_dups/{sample}_{type}.mark_dups.bam.bai"
    params:
        threads=4
    resources:
        runtime=540,
        mem_mb=3000,
        cpus_per_task=4
    shell:
        r"""
        module purge
        module load bio/NanoSeq

        bammarkduplicatesopt inputthreads={params.threads} optminpixeldif=2500 I={input} O={output.bam}

        module purge
        module load bio/SAMtools/1.16.1-GCC-11.3.0

        samtools index -o {output.bai} {output.bam}
        """

rule filter_and_append_rb_tags:
    input:
        "07_bam_mark_dups/{sample}_{type}.mark_dups.bam"
    output:
        bam="08_bam_rb_tags/{sample}_{type}.rb_tags.bam",
        bai="08_bam_rb_tags/{sample}_{type}.rb_tags.bam.bai"
    resources:
        runtime=480,
        mem_mb=1000,
        cpus_per_task=1
    shell:
        r"""
        module purge
        module load bio/NanoSeq

        bamaddreadbundles -I {input} -O {output.bam}

        module purge
        module load bio/SAMtools/1.16.1-GCC-11.3.0

        samtools index -o {output.bai} {output.bam}
        """

rule keep_random_read:
    input:
        "08_bam_rb_tags/{sample}_{type}.rb_tags.bam"
    output:
        bam="09_bam_random_read/{sample}_{type}.random_read.bam",
        bai="09_bam_random_read/{sample}_{type}.random_read.bam.bai"
    resources:
        runtime=420,
        mem_mb=4000,
        cpus_per_task=1
    shell:
        r"""
        module purge
        module load bio/NanoSeq

        randomreadinbundle -I {input} -O {output.bam}

        module purge
        module load bio/SAMtools/1.16.1-GCC-11.3.0

        samtools index -o {output.bai} {output.bam}
        """

rule run_analysis:
    input:
        bam_undiluted="09_bam_random_read/{sample}_{undiluted_type}.random_read.bam",
        bam_duplex="08_bam_rb_tags/{sample}_{duplex_type}.rb_tags.bam",
        genome="01_ref/GRCh38.primary_assembly.genome.fa",
        mask_SNP="01_ref/SNP.sorted.GRCh38.bed.gz",
        mask_NOISE="01_ref/NOISE.sorted.GRCh38.bed.gz"
    output:
        dir=directory("10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}"),
        muts="10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}/tmpNanoSeq/post/results.muts.vcf.gz",
        indel="10_analysis/{sample}_{duplex_type}-vs-{undiluted_type}/tmpNanoSeq/post/results.indel.vcf.gz"
    resources:
        runtime=840,
        mem_mb=96000,
        cpus_per_task=21
    shell:
        r"""
        module purge
        module load bio/NanoSeq

        cd {output.dir}

        # cov
        runNanoSeq.py -t 10 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          cov \
          -Q 0
        # part
        runNanoSeq.py -t 1 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          part \
          -n 20
        # dsa
        runNanoSeq.py -t 20 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          dsa \
          -C ../../{input.mask_SNP} \
          -D ../../{input.mask_NOISE} \
          -d 2 \
          -q 30
        # var
        runNanoSeq.py -t 20 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          var \
          -a 50 \
          -b 0 \
          -c 0.02 \
          -d 2 \
          -f 0.9 \
          -i 1 \
          -m 8 \
          -n 3 \
          -p 0 \
          -q 60 \
          -r 144 \
          -v 0.01 \
          -x 8 \
          -z 15
        # indel
        runNanoSeq.py -t 20 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          indel \
          -s {wildcards.sample}_{wildcards.duplex_type}-vs-{wildcards.undiluted_type} \
          --rb 2 \
          --t3 136 \
          --t5 8 \
          -a 50 \
          -c 0.02 \
          -z 15 \
          -v 0.01
        # post
        runNanoSeq.py -t 2 \
          -A ../../{input.bam_undiluted} \
          -B ../../{input.bam_duplex} \
          -R ../../{input.genome} \
          post
        """

rule check_contamination:
    input:
        bam="09_bam_random_read/{sample}_{type}.random_read.bam",
        hla_fa="01_ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    output:
        out="11_contamination_check/{sample}_{type}.out",
        selfsm="11_contamination_check/{sample}_{type}.selfSM"
    resources:
        runtime=60, 
        mem_mb=2000,
        cpus_per_task=4
    conda:
        "VerifyBamID2-conda-env.yaml"
    shell:
        r"""
        mkdir -p 11_contamination_check
        verifybamid2 1000g.phase3 100k b38 --Reference {input.hla_fa} --BamFile {input.bam} --Output 11_contamination_check/{wildcards.sample}_{wildcards.type}
        """

rule estimate_efficiency:
    input:
        bam_rbtags="08_bam_rb_tags/{sample}_{type}.rb_tags.bam",
        bam_randomread="09_bam_random_read/{sample}_{type}.random_read.bam",
        genome="01_ref/GRCh38.primary_assembly.genome.fa"
    output:
        rbs="12_efficiency_estimate/{sample}_{type}.RBs",
        gc_inserts="12_efficiency_estimate/{sample}_{type}.RBs.GC_inserts.tsv",
        pdf="12_efficiency_estimate/{sample}_{type}.RBs.pdf"
    params:
        threads=2
    resources:
        runtime=60, 
        mem_mb=10000,
        cpus_per_task=2
    shell:
        r"""
        module purge
        module load bio/NanoSeq

        mkdir -p 12_efficiency_estimate
        cd 12_efficiency_estimate
        efficiency_nanoseq.pl -t {params.threads} -d ../{input.bam_randomread} -x ../{input.bam_rbtags} -r ../{input.genome} -o {wildcards.sample}_{wildcards.type}
        """
