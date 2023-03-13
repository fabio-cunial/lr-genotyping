version 1.0


# 
#
workflow Sv2Igv {
    input {
        File vcf_file
        File bam_addresses
        String output_bucket_dir
        Int n_cpus
        File reference_fa
    }
    parameter_meta {
        bam_addresses: "File containing a list of bucket addresses."
    }
    
    call Sv2IgvImpl {
        input:
            vcf_file = vcf_file,
            bam_addresses = bam_addresses,
            output_bucket_dir = output_bucket_dir,
            n_cpus = n_cpus,
            reference_fa = reference_fa
    }
    output {
        File report = Sv2IgvImpl.report
    }
}


# 
#
task Sv2IgvImpl {
    input {
        File vcf_file
        File bam_addresses
        String output_bucket_dir
        Int n_cpus
        File reference_fa
    }
    parameter_meta {
    }
    
    Int ram_size_gb = n_cpus*4  # Arbitrary
    Int disk_size_gb = 128  # Arbitrary

    command <<<
        set -euxo pipefail
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        CHR21_LENGTH="46709983"
        CHR22_LENGTH="50818468"
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))  # 2: Trying to use hyperthreading as well.
        
        function downloadThread() {
            local CHUNK_FILE=$1
            local REGION=$2
            
            while read REMOTE_BAM; do
                SAMPLE_ID=$(basename -s .bam ${REMOTE_BAM})
                TEST=$(samtools view --threads 1 -h --bam ${REMOTE_BAM} ${REGION} > ${SAMPLE_ID}.bam && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
                    ${TIME_COMMAND} samtools view --threads 1 -h --bam ${REMOTE_BAM} ${REGION} > ${SAMPLE_ID}.bam
                fi
                samtools index ${SAMPLE_ID}.bam
            done < ${CHUNK_FILE}
        }
        
        REGION="chr22:50674415-50733298"
        echo -e "chr22\t50674415\t50733298" > region.bed
        
        BAM_NAME=$(basename -s .bam ~{bam_addresses})
        N_ROWS=$(wc -l < ~{bam_addresses})
        N_ROWS_PER_CHUNK=$(( ${N_ROWS}/${N_THREADS} ))
        split -l ${N_ROWS_PER_CHUNK} ~{bam_addresses} chunk_
        for CHUNK in chunk_*; do
            downloadThread ${CHUNK} ${REGION} &
        done
        wait
        ${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o all.bam *.bam
        
        # IGV-REPORT
        ${TIME_COMMAND} create_report region.bed ~{reference_fa} --flanking 1000 --exclude-flags 0 --sort BASE --tracks all.bam --output report.html --sequence 1 --begin 2 --end 3 --standalone
    >>>

    output {
        File report = "report.html"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: n_cpus
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
