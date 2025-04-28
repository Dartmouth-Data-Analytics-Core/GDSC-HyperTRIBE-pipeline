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
        expand("trimming/{sample}_R2.paired.trim.fastq.gz", sample=sample_list) if config["layout"] == "paired" else [],
        expand("trimming/{sample}_R1.unpaired.trim.fastq.gz", sample=sample_list) if config["layout"] == "paired" else [],
        expand("trimming/{sample}_R2.unpaired.trim.fastq.gz", sample=sample_list) if config["layout"] == "paired" else [],
       
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

        #----- Rule build_SQL outputs
        expand("mySQL_logs/{sample}.connected.txt", sample=sample_list),

        #----- Rule diagnose_SQL outputs
        expand("mySQL_logs/{sample}.diagnostic.txt", sample=sample_list),

        #----- Rule Fine_RNA_edit_sites outputs
        expand("RNA_Edits/{sample}.a2g.txt", sample=sample_list),

        #----- Rule edits2bedgraph outputs
        expand("bedgraphs/{sample}.bedgraph", sample=sample_list),

    output: "multiqc_report.html"
    conda: "rnaseq1"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    shell: """
        
        #----- Run MultiQC
        multiqc trimming noDups
    """

#----- Rule to execute trimming
rule trim:
    output:
        trimR1 = "trimming/{sample}_R1.paired.trim.fastq.gz",
        trimR2 = "trimming/{sample}_R2.paired.trim.fastq.gz" if config["layout"]=="paired" else [],
        unpairedR1 = "trimming/{sample}_R1.unpaired.trim.fastq.gz" if config["layout"]=="paired" else [],
        unpairedR2 = "trimming/{sample}_R2.unpaired.trim.fastq.gz" if config["layout"]=="paired" else []
    conda: "align"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        layout = config["layout"],
        fastq1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"] if config["layout"]=="paired" else [],
        avgquality = config["avgquality"]
    shell: """
        
        if [ "{params.layout}" == "paired" ]
        then
            trimmomatic PE -phred33 \
                {params.fastq1} {params.fastq2} \
                {output.trimR1} {output.unpairedR1} \
                {output.trimR2} {output.unpairedR2} \
                HEADCROP:6 \
                LEADING:25 \
                TRAILING:25 \
                AVGQUAL:{params.avgquality} \
                MINLEN:19
        else
            trimmomatic SE -phred33 \
                {params.fastq1} \
                {output.trimR1} \
                HEADCROP:6 \
                LEADING:25 \
                TRAILING:25 \
                AVGQUAL:{params.avgquality} \
                MINLEN:19
        fi     
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

#----- Rule to generate matrix needed to build SQL table
rule sam2matrix:
    input:
        samSorted = "noDups/{sample}.nodup.sort.sam"
    output:
        matrixFile = "matrix/{sample}.matrix",
    conda: "perl"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        replicate = lambda wildcards: samples_df.loc[wildcards.sample, "replicate"],
        samToMatrix = config["samToMatrixScript"],
    shell: """
        #----- Run the sam_to_matrix.pl script
        perl {params.samToMatrix} \
            {input.samSorted} \
            {params.sample} \
            {params.replicate}
       
        #----- Move perl output to new file name
        mv noDups/{params.sample}.nodup.sort.sam.matrix.wig {output.matrixFile}
    """

#----- Rule to build SQL table 
rule build_SQL:
    input:
        matrixFile = "matrix/{sample}.matrix"
    output:
        sqlLog = "mySQL_logs/{sample}.connected.txt"
    conda: "perl"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        loadMatrix = config["loadMatrix"],
    shell: """
        
        #----- Load matrix
        touch "{output.sqlLog}"
        echo "Starting connection for table {params.sample} in dmseq database..." >> {output.sqlLog}
        perl {params.loadMatrix} -t {params.sample} -d {input.matrixFile}  && \
        echo "Connection successful!" >> {output.sqlLog}
    """

#----- Rule to check table entry is same size as matrix
rule diagnose_SQL:
    input:
        sqlLog = "mySQL_logs/{sample}.connected.txt"
    output:
        diagnostic = "mySQL_logs/{sample}.diagnostic.txt"
    conda: "perl"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        diagnoseDatabase = config["diagnoseDatabase"]
    shell: """
        
        #----- Diagnose that table is loaded correctly
        touch "{output.diagnostic}"
        echo "Diagnosing {params.sample} database..." >> {output.diagnostic}
        perl {params.diagnoseDatabase} -t {params.sample} "matrix/{params.sample}.matrix" && \
        echo "Data table exists!" >> {output.diagnostic}


    """

#----- Rule to identify RNAedit sites
# Eventually, will need to add support for replicates...
rule find_RNA_edit_sites:
    input:
        diagnostic = "mySQL_logs/{sample}.diagnostic.txt",
        sqlLog = "mySQL_logs/{sample}.connected.txt",
        matrixFile = "matrix/{sample}.matrix"
    output:
        a2g = "RNA_Edits/{sample}.a2g.txt"
    conda: "perl"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        findRNA = config["findRNA"],
        annotations = config["annotationFile"],
        refSample = config["refSample"],
        timepoint = config["timepoint"],
        wtTimepoint = config["wtTimepoint"],
    shell: """

        #-----HELP MENU
        # -a = annotation file
        # -t = RNA table name
        # -e = RNA experiment name (sample name)
        # -c = RNA timepoint (from sample sheet)
        # -o = output file name
        # -g = reference RNA sample
        # -j = reference RNA table name
        # -k = reference RNA timepoint (from sample sheet)

        #----- Check that the current sample is not the reference sample
        if [ {params.sample} == {params.refSample} ]
        then 
            echo "Skipping reference sample {params.sample}"
            touch RNA_Edits/{params.refSample}.a2g.txt
        else
            #----- Run perl script to find RNA edit sites
            perl {params.findRNA} \
                -a {params.annotations} \
                -t {params.sample} \
                -e {params.sample} \
                -c {params.timepoint} \
                -o {output.a2g} \
                -g {params.refSample} \
                -j {params.refSample} \
                -k {params.wtTimepoint}
        fi 
    """

#----- Rule to convert RNA Edits to bedgraph
rule edits2bedgraph:
    input:
        a2g = "RNA_Edits/{sample}.a2g.txt"
    output:
        bedgraph = "bedgraphs/{sample}.bedgraph"
    conda: "rnaseq1"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        bedgraphScript = config["bedgraphScript"]
    shell: """

        #----- Run the python script
        python {params.bedgraphScript} {input.a2g}

        #----- Move python output to new file name
        mv RNA_Edits/{params.sample}.a2g.txt.bedgraph {output.bedgraph}
    
    """




    






