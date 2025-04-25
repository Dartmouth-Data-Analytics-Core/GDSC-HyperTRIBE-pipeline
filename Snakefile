#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GDSC-HyperTribe Pipeline
# Adapted from: HyperTribe (https://hypertribe.readthedocs.io)
#
# Authors: Shannon Soucy, Mike Martinez
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##########---------- TO DO ----------##########


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
import pandas as pd 

#----- Set config file
configfile: "config.yaml"

#----- Read in the sample data
samples_df = pd.read_table(config["sample_csv"], delimiter=",").set_index("Sample_ID", drop=False)
sample_list = list(samples_df["Sample_ID"])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Final Rule
rule all:
    input:
        #----- Rule trim outputs
        expand("trimming/{sample}_R1.paired.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}_R2.paired.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}_R1.unpaired.trim.fastq.gz", sample=sample_list),
        expand("trimming/{sample}_R2.unpaired.trim.fastq.gz", sample=sample_list),
       
        #----- Rule align outputs
        expand("alignment/{sample}.sam", sample=sample_list),
        
        #----- Rule sam2bam outputs
        expand("alignment/{sample}.bam", sample=sample_list),
        expand("alignment/{sample}.sort.bam", sample=sample_list),

        #----- Rule mark_duplicates outputs
        expand("noDups/{sample}.nodup.bam", sample=sample_list),
        expand("noDups/{sample}.dup.txt", sample=sample_list),
        expand("noDups/{sample}.nodup.sort.bam", sample=sample_list),
        expand("noDups/{sample}.nodup.sort.bam.bai", sample=sample_list),

        #----- Rule sort sam outputs
        expand("noDups/{sample}.nodup.sort.sam", sample=sample_list),

        #----- Rule sam2matrix outputs
        expand("matrix/{sample}.matrix", sample=sample_list),

    output: "multiqc_report.html"
    conda: "rnaseq1"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    shell: """
        
        #----- Run MultiQC
        multiqc trimming nodup
    """

#----- Rule to execute trimming
rule trim:
    output:
        trimR1 = "trimming/{sample}_R1.paired.trim.fastq.gz",
        trimR2 = "trimming/{sample}_R2.paired.trim.fastq.gz",
        unpairedR1 = "trimming/{sample}_R1.unpaired.trim.fastq.gz",
        unpairedR2 = "trimming/{sample}_R2.unpaired.trim.fastq.gz"
    conda: "align"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        fastq1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"],
        avgquality = config["avgquality"]
    shell: """
        
        trimmomatic PE -phred33 \
            {params.fastq1} {params.fastq2} \
            {output.trimR1} {output.unpairedR1} \
            {output.trimR2} {output.unpairedR2} \
            HEADCROP:6 \
            LEADING:25 \
            TRAILING:25 \
            AVGQUAL:{params.avgquality} \
            MINLEN:19
		
        
        
    """

rule align:
    input:
        trimR1 = "trimming/{sample}_R1.paired.trim.fastq.gz",
        trimR2 = "trimming/{sample}_R2.paired.trim.fastq.gz"
    output:
        aligned = "alignment/{sample}.sam"
    conda: "align"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        starIndex = config["starIndex"],
    shell: """
        #----- align with STAR
        STAR  \
            --runThreadN 8 \
            --readFilesCommand zcat \
            --outFilterMismatchNoverLmax 0.07 \
        	--outFileNamePrefix alignment/{params.sample}"_" \
            --outFilterMatchNmin 16 \
        	--outFilterMultimapNmax 1 \
            --genomeDir {params.starIndex} \
        	--readFilesIn {input.trimR1} {input.trimR2}

        #----- Rename the file
        mv alignment/{params.sample}_Aligned.out.sam {output.aligned}
    
    """

#----- Rule to convert sam to bam
rule sam2bam:
    input:
        sam = "alignment/{sample}.sam"
    output:
        bam = "alignment/{sample}.bam",
        sortedBam = "alignment/{sample}.sort.bam"
    conda: "samtools"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
    shell: """
        
        #----- View and sort sam tile
        samtools view -@ 4 -Sh -q 10 {input.sam} > alignment/{params.sample}.highquality.sam
        
        #----- Convert to bam and sort
        samtools view -@ 4 -bhS alignment/{params.sample}.highquality.sam > {output.bam}
        samtools sort -@ 6 {output.bam} -o {output.sortedBam}
    """

#----- Rule to mark duplicates
rule mark_duplicates:
    input:
        sortedBam = "alignment/{sample}.sort.bam"
    output:
        noDup = "noDups/{sample}.nodup.bam",
        metrics = "noDups/{sample}.dup.txt",
        noDupSorted = "noDups/{sample}.nodup.sort.bam",
        index = "noDups/{sample}.nodup.sort.bam.bai"
    conda: "rnaseq1"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
    shell:"""
    
        #----- Run picard markduplictes
        picard MarkDuplicates \
        INPUT={input.sortedBam} \
        OUTPUT={output.noDup} \
        METRICS_FILE={output.metrics} \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=true \
        TMP_DIR=tmp \
        ASSUME_SORTED=true

        #----- Sort the new bam file
        samtools sort -@ 6 {output.noDup} -o {output.noDupSorted}

        #----- Index the sorted bam file
        samtools index {output.noDupSorted}
    """

#----- Rule to sort sam file
rule sort_sam:
    input:
        noDupSorted = "noDups/{sample}.nodup.sort.bam"
    output:
        samSorted = "noDups/{sample}.nodup.sort.sam"
    conda: "samtools"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
    shell: """

        #----- Create sam from bam file
        samtools view -@ 4 -h {input.noDupSorted} > {output.samSorted}
    
    """

#----- Rule to generate mySql Table (might need to add script argument to sample file)
rule sam2matrix:
    input:
        samSorted = "noDups/{sample}.nodup.sort.sam"
    output:
        matrixFile = "matrix/{sample}.matrix"
    conda: "perl"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        expName = config["experimentName"],
        replicate = lambda wildcards: samples_df.loc[wildcards.sample, "replicate"],
        loadTable = config["loadTableScript"],
        samToMatrix = config["samToMatrixScript"],
        loadMatrix = config["loadMatrix"]
    shell: """
        #----- Run the sam_to_matrix.pl script
        perl {params.samToMatrix} \
            {input.samSorted} \
            {params.sample} \
            {params.expName} \
            {params.replicate}
       
       #----- Move perl output to new file name
       mv noDups/{params.sample}.nodup.sort.sam.matrix.wig {output.matrixFile}

       #----- Create SQL database
       #perl {params.loadMatrix} -t {params.sample} -d {output.matrixFile}
    """






