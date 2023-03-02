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
        #File report = Sv2IgvImpl.report
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
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))  # 2: Trying to use hyperthreading as well.
        
        function downloadThread() {
            local CHUNK_FILE=$1
            
            ID=${CHUNK_FILE#chunk_}
            i="0"
            while read REMOTE_BAM; do 
                TEST=$(samtools view --threads 1 --target-file ~{regions_bed} --reference ~{reference_fa} --fai-reference ~{reference_fai} --bam --output alignments_${ID}_${i}.bam ${REMOTE_BAM} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
                    ${TIME_COMMAND} samtools view --threads ${N_THREADS} --target-file ~{regions_bed} --reference ~{reference_fa} --fai-reference ~{reference_fai} --bam --output alignments_${ID}_${i}.bam ${REMOTE_BAM}
                fi
                i=$(( ${i}+1 ))
            done < ${CHUNK_FILE}
        }
        
        # Building the merged BAM across all samples and SVs
        BAM_NAME=$(basename -s .bam ~{bam_addresses})
        BED_NAME=$(basename -s .bed ~{regions_bed})
        TEST=$(gsutil -q stat ~{bucket_dir}/${BAM_NAME}_${BED_NAME}.bam && echo 0 || echo 1)
        if [ ${TEST} -eq 0 ]; then
            while : ; do
                TEST2=$(gsutil -m cp ~{bucket_dir}/${BAM_NAME}_${BED_NAME}.bam ./all.bam && echo 0 || echo 1)
                if [ ${TEST2} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/${BAM_NAME}_${BED_NAME}.bam>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            N_ROWS=$(wc -l < ~{bam_addresses})
            N_ROWS_PER_CHUNK=$(( ${N_ROWS}/${N_THREADS} ))
            split -l ${N_ROWS_PER_CHUNK} ~{bam_addresses} chunk_
            for CHUNK in chunk_*; do
                downloadThread ${CHUNK} &
            done
            wait
            ${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o all.bam alignments_*.bam
            rm -f alignments_*.bam
            while : ; do
                TEST=$(gsutil -m cp ./all.bam ~{bucket_dir}/${BAM_NAME}_${BED_NAME}.bam && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <all.bam>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} all.bam
        
        # Printing one image per SV
        IMAGE_HEIGHT="10500"
        HORIZONTAL_SLACK="10000"
        FIGURES_DIR="${pwd}/figures"
        rm -rf ${FIGURES_DIR}; mkdir ${FIGURES_DIR}
        IGV_SCRIPT="script.txt"
        echo "new" > ${IGV_SCRIPT}
        echo "snapshotDirectory ${FIGURES_DIR}" >> ${IGV_SCRIPT}
        echo "load all.bam" >> ${IGV_SCRIPT}
        echo "genome ~{reference_fa}" >> ${IGV_SCRIPT}
        echo "maxPanelHeight ${IMAGE_HEIGHT}" >> ${IGV_SCRIPT}
        i="0"
        while read SV; do
            CHR=$(echo ${SV} | cut -f 1)
            START=$(echo ${SV} | cut -f 2)
            START=$(( ${START}-${HORIZONTAL_SLACK} ))
            END=$(echo ${SV} | cut -f 3)
            END=$(( ${END}+${HORIZONTAL_SLACK} ))
            echo "goto ${CHR}:${START}-${END}" >> ${IGV_SCRIPT}
            echo "sort base" >> ${IGV_SCRIPT}
            echo "squish" >> ${IGV_SCRIPT}
            echo "snapshot sv-${i}.png" >> ${IGV_SCRIPT}
            i=$(( ${i}+1 ))
        done < ~{regions_bed}
        echo "exit" >> ${IGV_SCRIPT}
        cat ${IGV_SCRIPT}
        python ./IGV-snapshot-automator/make_IGV_snapshots.py -mem $(( ~{ram_size_gb}-4 )) -onlysnap ${IGV_SCRIPT}
        tar -czvf report.tar.gz ${FIGURES_DIR}
        
        #${TIME_COMMAND} create_report ~{regions_bed} ~{reference_fa} --flanking 10000 --exclude-flags 0 --sort BASE --tracks all.bam --output report.html --sequence 1 --begin 2 --end 3 --standalone
    >>>

    output {
        File report = "report.tar.gz"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: n_cpus
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
