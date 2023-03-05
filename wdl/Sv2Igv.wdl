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
    
    Int ram_size_gb = n_cpus*4  # Arbitrary
    Int disk_size_gb = 128  # Arbitrary

    command <<<
        set -euxo pipefail
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        Xvfb :1 &
        export DISPLAY=":1"
        IMAGE_HEIGHT="20"; HORIZONTAL_SLACK="10000"
        CHR21_LENGTH="46709983"
        CHR22_LENGTH="50818468"
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))  # 2: Trying to use hyperthreading as well.
        
        function downloadThread() {
            local CHUNK_FILE=$1
            
            while read REMOTE_BAM; do
                SAMPLE_ID=$(basename -s .bam ${REMOTE_BAM})
                TEST=$(samtools view --threads 1 --target-file ~{regions_bed} --reference ~{reference_fa} --fai-reference ~{reference_fai} --bam --output ${SAMPLE_ID}.bam ${REMOTE_BAM} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
                    ${TIME_COMMAND} samtools view --threads ${N_THREADS} --target-file ~{regions_bed} --reference ~{reference_fa} --fai-reference ~{reference_fai} --bam --output ${SAMPLE_ID}.bam ${REMOTE_BAM}
                fi
                samtools index ${SAMPLE_ID}.bam
                IGV_SCRIPT="${SAMPLE_ID}.txt";
                echo "new" > ${IGV_SCRIPT}
                echo "snapshotDirectory ." >> ${IGV_SCRIPT}
                echo "load ${SAMPLE_ID}.bam" >> ${IGV_SCRIPT}
                echo "genome ~{reference_fa}" >> ${IGV_SCRIPT}
                echo "maxPanelHeight ${IMAGE_HEIGHT}" >> ${IGV_SCRIPT}
                tr '\t' ',' < ~{regions_bed} > regions_prime.txt
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
                    echo "goto ${CHR}:${START}-${END}" >> ${IGV_SCRIPT}
                    echo "sort base" >> ${IGV_SCRIPT}
                    echo "squish" >> ${IGV_SCRIPT}
                    echo "snapshot ${CHR}_${START}_${END}_${SAMPLE_ID}.png" >> ${IGV_SCRIPT}
                done < regions_prime.txt
                echo "exit" >> ${IGV_SCRIPT}
                ${TIME_COMMAND} python /IGV-snapshot-automator/make_IGV_snapshots.py -mem 2000 -onlysnap ${IGV_SCRIPT} ${SAMPLE_ID}.bam
            done < ${CHUNK_FILE}
        }
        
        BAM_NAME=$(basename -s .bam ~{bam_addresses})
        BED_NAME=$(basename -s .bed ~{regions_bed})
        N_ROWS=$(wc -l < ~{bam_addresses})
        N_ROWS_PER_CHUNK=$(( ${N_ROWS}/${N_THREADS} ))
        split -l ${N_ROWS_PER_CHUNK} ~{bam_addresses} chunk_
        for CHUNK in chunk_*; do
            downloadThread ${CHUNK} &
        done
        wait
        while : ; do
            TEST=$(gsutil -m cp './*.png' ~{bucket_dir}/${BAM_NAME}_${BED_NAME}/ && echo 0 || echo 1)
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
