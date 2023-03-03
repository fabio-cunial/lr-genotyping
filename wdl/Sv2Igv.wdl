version 1.0


# 
#
workflow Sv2Igv {
    input {
        File bam_addresses
        File regions_bed
        File reference_fa
        File reference_fai
        String bucket_dir
        Int n_cpus
        Int bam_size_gb
    }
    parameter_meta {
        bam_addresses: "File containing a list of bucket addresses."
        bam_size_gb: "Upper bound on the size of a single BAM."
    }
    
    call Sv2IgvImpl {
        input:
            bam_addresses = bam_addresses,
            regions_bed = regions_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            bucket_dir = bucket_dir,
            n_cpus = n_cpus,
            bam_size_gb = bam_size_gb
    }
    output {
    }
}


# 
#
task Sv2IgvImpl {
    input {
        File bam_addresses
        File regions_bed
        File reference_fa
        File reference_fai
        String bucket_dir
        Int n_cpus
        Int bam_size_gb
    }
    parameter_meta {
        bam_size_gb: "Upper bound on the size of a single BAM."
    }
    
    Int ram_size_gb = n_cpus*2  # Arbitrary
    Int disk_size_gb = 128  # Arbitrary

    command <<<
        set -euxo pipefail
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))  # 2: Hoping that hyperthreading can speed parallel downloads up...
        CHR21_LENGTH="46709983"
        CHR22_LENGTH="50818468"
        
        function downloadThread() {
            local CHUNK_FILE=$1
            local BAM_LIST="${CHUNK_FILE}.bamlist.txt"
            
            rm -f ${BAM_LIST}; touch ${BAM_LIST}
            ID=${CHUNK_FILE#chunk_}
            i="0"
            while read REMOTE_BAM; do
                TEST=$(samtools view --threads 1 --target-file ~{regions_bed} --reference ~{reference_fa} --fai-reference ~{reference_fai} --bam --output alignments_${ID}_${i}.bam ${REMOTE_BAM} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
                    ${TIME_COMMAND} samtools view --threads ${N_THREADS} --target-file ~{regions_bed} --reference ~{reference_fa} --fai-reference ~{reference_fai} --bam --output alignments_${ID}_${i}.bam ${REMOTE_BAM}
                fi
                samtools index alignments_${ID}_${i}.bam
                echo -e "alignments_${ID}_${i}.bam\t$(basename -s .bam ${REMOTE_BAM})" >> ${BAM_LIST}
                i=$(( ${i}+1 ))
            done < ${CHUNK_FILE}
        }
        
        # Downloading the local BAM of every sample across all the SVs
        BAM_NAME=$(basename -s .bam ~{bam_addresses})
        BED_NAME=$(basename -s .bed ~{regions_bed})
        TEST=$(gsutil -q stat ~{bucket_dir}/${BAM_NAME}_${BED_NAME}/report.tar.gz && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            echo "The report already exists in the bucket"
            exit 0
        fi
        N_ROWS=$(wc -l < ~{bam_addresses})
        N_ROWS_PER_CHUNK=$(( ${N_ROWS}/${N_THREADS} ))
        split -l ${N_ROWS_PER_CHUNK} ~{bam_addresses} chunk_
        for CHUNK in chunk_*; do
            downloadThread ${CHUNK} &
        done
        wait
        cat *.bamlist.txt > bamlist.txt; rm *.bamlist.txt
        
        # Printing a report with distinct sections for each individual
        HORIZONTAL_SLACK="10000"
        REPORT_DIR="report"
        rm -rf ${REPORT_DIR}; mkdir ${REPORT_DIR}
        BAMSNAP_BED="bamsnap.bed"
        tr '\t' ',' < ~{regions_bed} > regions_prime.txt
        i="0"
        while read SV; do
            CHR=$(echo ${SV} | cut -d , -f 1)
            START=$(echo ${SV} | cut -d , -f 2)
            START=$(( ${START}-${HORIZONTAL_SLACK} ))
            if [ ${START} -lt 1 ]; then
                START="1"
            fi
            END=$(echo ${SV} | cut -d , -f 3)
            END=$(( ${END}+${HORIZONTAL_SLACK} ))
            if [ ${END} -le ${START} ]; then
                END=$(( ${START}+1 ))
            fi
            if [ ${CHR} = "chr21" -a ${END} -gt ${CHR21_LENGTH} ]; then
                END=$(( ${CHR21_LENGTH} ))
            elif [ ${CHR} = "chr22" -a ${END} -gt ${CHR22_LENGTH} ]; then
                END=$(( ${CHR22_LENGTH} ))
            fi
            echo -e "${CHR}\t${START}\t${END}" >> ${BAMSNAP_BED}
            i=$(( ${i}+1 ))
        done < regions_prime.txt
        cat ${BAMSNAP_BED}; cat bamlist.txt
        ${TIME_COMMAND} bamsnap -process ${N_THREADS} -ref ~{reference_fa} -bamlist bamlist.txt -bed ${BAMSNAP_BED} -out ${REPORT_DIR} \
            -separated_bam \
            -bamplot read -read_thickness 2 -read_gap_height 0 -read_gap_width 1 \
            -show_soft_clipped
        tar -czvf report.tar.gz ${REPORT_DIR}
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp report.tar.gz ~{bucket_dir}/${BAM_NAME}_${BED_NAME}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading file <~{bucket_dir}/${BAM_NAME}_${BED_NAME}/report.tar.gz>. Trying again..."
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
