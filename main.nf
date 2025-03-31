nextflow.enable.dsl=2

params.genome = 'GRCh38'
params.csv_file = "samples.csv"
params.reference_genome = params.genomes[params.genome].fasta
params.index_bowtie2    = params.genomes[params.genome].bowtie2
params.reference_gtf    = params.genomes[params.genome].gtf
params.adapters = "adapters.fa"

Channel
    .fromPath(params.csv_file)
    .splitCsv(header: true, sep: '\t')
    .map { row -> 
        return tuple(row.gsm_id, row.sample_name)
    }
    .set { sample_run_ch }

process EXTRACT_SRR {
    tag "Extract SRR for ${gsm_id}"

    input:
    tuple val(gsm_id), val(sample_name)

    output:
    tuple val(gsm_id), val(sample_name), 
    path("${gsm_id}_srr_list.txt"),
    path("${gsm_id}_read_type.txt")

    script:
    """
    singularity exec ${params.img_edirect} esearch -db gds -query ${gsm_id} | elink -target sra | efetch -format runinfo > ${gsm_id}_runinfo.txt

    # Extract SRR IDs
    grep ${gsm_id} ${gsm_id}_runinfo.txt | cut -d ',' -f 1 | grep SRR > ${gsm_id}_srr_list.txt

    # Determine read type
    if grep ${gsm_id} ${gsm_id}_runinfo.txt | grep -q 'PAIRED'; then
        read_type='pe'
        echo "pe" > ${gsm_id}_read_type.txt
    else
        read_type='se'
        echo "se" > ${gsm_id}_read_type.txt
    fi

    """
    
}


process DOWNLOAD_SRA {                                                          
                                                                                
    tag "Download SRA for ${sample_name}"                                       
                                                                                
    input:                                                                      
    tuple val(gsm_id),                                                          
    val(sample_name),                                                           
    path(srr_file),                                                             
    path(readtype_file)                                                         
                                                                                
    output:                                                                     
    tuple val(sample_name), path("${sample_name}_{1,2}.fastq.gz"), val('auto')

    script: 
    """                                                                         
    mkdir -p ${sample_name}_fastq                                                  
                                                                                
    read_type=\$(cat "${readtype_file}")                                        
                                                                                
    while read srr_id; do                                                       
        singularity exec ${params.img_sratool} /opt/sratoolkit.3.1.0-ubuntu64/bin/prefetch "\$srr_id"

        if [ "\$read_type" = "pe" ]; then                                     
            singularity exec ${params.img_sratool} /opt/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump --split-files "\$srr_id" -O ${sample_name}_fastq     
        else                                                                    
            singularity exec ${params.img_sratool} /opt/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump "\$srr_id" -O ${sample_name}_fastq                   
            mv ${sample_name}_fastq/\${srr_id}.fastq ${sample_name}_fastq/\${srr_id}_1.fastq 
            touch ${sample_name}_fastq/\${srr_id}_2.fastq 
        fi                                                                      
    done < "${srr_file}"                                                          

    pattern="${sample_name}_fastq/*_1.fastq"
    cat \$(ls \$pattern | sort) > ${sample_name}_1.fastq
    gzip ${sample_name}_1.fastq


    pattern="${sample_name}_fastq/*_2.fastq"
    cat \$(ls \$pattern | sort) > ${sample_name}_2.fastq
    gzip ${sample_name}_2.fastq

    """ 
}  

process GENERATE_SAMPLESHEET {
    input:
    tuple val(sample_name), path(fastq_files), val(strandedness)

    output:
    path("samplesheet.csv")

    script:
    """
    echo "sample,fastq_1,fastq_2,strandedness" > samplesheet.csv
    echo "${sample_name},${fastq_files[0]},${fastq_files[1]},${strandedness}" >> samplesheet.csv
    """
}

workflow {

    parsed_samples = sample_run_ch

    srr_samples = parsed_samples | EXTRACT_SRR

    fastq_samples = srr_samples | DOWNLOAD_SRA

    sample_sheet = fastq_samples | GENERATE_SAMPLESHEET

}
