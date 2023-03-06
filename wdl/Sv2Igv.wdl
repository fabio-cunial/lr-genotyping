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
        HORIZONTAL_SLACK="2000"
        CHR21_LENGTH="46709983"
        CHR22_LENGTH="50818468"
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))  # 2: Trying to use hyperthreading as well.
        
        function downloadThread() {
            local CHUNK_FILE=$1
            local CHR=$2
            local START=$3
            local END=$4
            
            while read REMOTE_BAM; do
                SAMPLE_ID=$(basename -s .bam ${REMOTE_BAM})
                TEST=$(samtools view --threads 1 ${REMOTE_BAM} ${CHR}:${START}-${END} > ${SAMPLE_ID}.sam && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
                    ${TIME_COMMAND} samtools view --threads 1 ${REMOTE_BAM} ${CHR}:${START}-${END} > ${SAMPLE_ID}.sam
                fi
            done < ${CHUNK_FILE}
        }
        
        tr '\t' ',' < ~{vcf_file} > sv.txt
        CHR=$(cut -d , -f 1 sv.txt)
        START=$(cut -d , -f 2 sv.txt)
        START=$(( ${START}-${HORIZONTAL_SLACK} ))
        if [ ${START} -lt 1 ]; then
            START="1"
        fi
        END=$(cut -d , -f 3 sv.txt)
        END=$(( ${END}+${HORIZONTAL_SLACK} ))
        if [ ${END} -le ${START} ]; then
            END=$(( ${START}+1 ))
        fi
        if [ ${CHR} = "chr21" -a ${END} -gt ${CHR21_LENGTH} ]; then
            END=$(( ${CHR21_LENGTH} ))
        elif [ ${CHR} = "chr22" -a ${END} -gt ${CHR22_LENGTH} ]; then
            END=$(( ${CHR22_LENGTH} ))
        fi
        BAM_NAME=$(basename -s .bam ~{bam_addresses})
        N_ROWS=$(wc -l < ~{bam_addresses})
        N_ROWS_PER_CHUNK=$(( ${N_ROWS}/${N_THREADS} ))
        split -l ${N_ROWS_PER_CHUNK} ~{bam_addresses} chunk_
        for CHUNK in chunk_*; do
            downloadThread ${CHUNK} ${CHR} ${START} ${END} &
        done
        wait
        ${TIME_COMMAND} java -Xmx$((~{ram_size_gb}-4))g /Pigv ~{vcf_file} . ${HORIZONTAL_SLACK} ${CHR22_LENGTH} image.png
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
