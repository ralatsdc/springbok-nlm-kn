#!/usr/bin/env nextflow

// salmon.nf

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --reads {*_1, *_2}.fq.gz
    
    Inputs Options:
    --reads               This is where the fq.gz files are entered
    """.stripIndent()
}

params.help = ""

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

reads = Channel.fromFilePairs(params.reads, size: 2)

// Download the transcriptome
process obtain_human_cdna {
   tag "$name"
   publishDir "results", mode: 'copy'
   container 'pgc-images.sbgenomics.com/deslattesmaysa2/curl

   output:
   file "hg38.cdna.all.fa.gz" into "cdna_results"

   script:
   """
   curl https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o hg38.cdna.all.fa.gz
   """
}

// Create an index from the latest download
process salmon_index {

    tag "$name"
    publishDir "results", mode: 'copy'
    container 'pgc-images.sbgenomics.com/deslattesmaysa2/salmon:v1.9'

    input:
    file ('hg38_cdna_gz') from cdna_results.collect()

    output:
    file ('hg38_cdna_index') into salmon_index

    script:
    """
    salmon index -t hg38_cdna_gz -i hg38_cdna_index
    """
}
process salmon-quant {

    publishDir "results", mode: 'copy'
    container 'pgc-images.sbgenomics.com/deslattesmaysa2/salmon:v1.9'

    input:
    set file (reads), val (hg38_cdna_index) from salmon_index

    output:
    file "*_quant" into quants

    script:
    """
    samp=`basename ${fn}`
    salmon quant -i athal_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o ${samp}_quant
    """
}
