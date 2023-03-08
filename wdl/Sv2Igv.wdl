version 1.0


# 
#
workflow Sv2Igv {
    input {
        File vcf_file
        File bam_addresses
        String output_bucket_dir
        Int n_cpus
    }
    parameter_meta {
        bam_addresses: "File containing a list of bucket addresses."
    }
    
    call Sv2IgvImpl {
        input:
            vcf_file = vcf_file,
            bam_addresses = bam_addresses,
            output_bucket_dir = output_bucket_dir,
            n_cpus = n_cpus
    }
    output {
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
    }
    parameter_meta {
    }
    
    Int ram_size_gb = n_cpus*4  # Arbitrary
    Int disk_size_gb = 128  # Arbitrary

    command <<<
        set -euxo pipefail
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        HORIZONTAL_SLACK="5000"
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
                TEST=$(samtools view --threads 1 ${REMOTE_BAM} ${REGION} > ${SAMPLE_ID}.sam && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
                    ${TIME_COMMAND} samtools view --threads 1 ${REMOTE_BAM} ${REGION} > ${SAMPLE_ID}.sam
                    samtools view --bam ${SAMPLE_ID}.sam > ${SAMPLE_ID}.bam
                    samtools index ${SAMPLE_ID}.bam
                fi
            done < ${CHUNK_FILE}
        }
        
        tail -n 1 ~{vcf_file} | tr '\t' ',' > sv.txt
        CHR=$(cut -d , -f 1 sv.txt)
        if [ ${CHR} = "chr21" ]; then
            CHR_LENGTH=${CHR21_LENGTH}
        elif [ ${CHR} = "chr22" ]; then
            CHR_LENGTH=${CHR22_LENGTH}
        fi
        REGION=$(java -cp / Vcf2Region ~{vcf_file} ${HORIZONTAL_SLACK} ${CHR_LENGTH})
        BAM_NAME=$(basename -s .bam ~{bam_addresses})
        N_ROWS=$(wc -l < ~{bam_addresses})
        N_ROWS_PER_CHUNK=$(( ${N_ROWS}/${N_THREADS} ))
        split -l ${N_ROWS_PER_CHUNK} ~{bam_addresses} chunk_
        for CHUNK in chunk_*; do
            downloadThread ${CHUNK} ${REGION} &
        done
        wait
        ${TIME_COMMAND} java -cp / -Xmx$((~{ram_size_gb}-4))g Pigv ~{vcf_file} . ${HORIZONTAL_SLACK} ${CHR_LENGTH} image.png
        
        # Samplot
        CHR=${REGION%":*"}; 
        REGION=${REGION#"*:"}; START=${REGION%"-*"}; END=${REGION#"*-"}
        ${TIME_COMMAND} samplot plot --bams *.bam --max_depth 10000 --output_file samplot.png --chrom ${CHR} --start $(( ${START}-${HORIZONTAL_SLACK} )) --end $(( ${END}+${HORIZONTAL_SLACK} )) 
        
        while : ; do
            TEST=$(gsutil -m cp './*.png' ~{output_bucket_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading PNG files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>

    output {
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: n_cpus
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
