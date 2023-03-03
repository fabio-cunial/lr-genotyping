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
        File report_merged_1 = Sv2IgvImpl.report_merged_1
        File report_distinct = Sv2IgvImpl.report_distinct
        File report_merged_2 = Sv2IgvImpl.report_merged_2
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
                TEST2=$(gsutil -m cp '~{bucket_dir}/${BAM_NAME}_${BED_NAME}.bam' ./all.bam && echo 0 || echo 1)
                if [ ${TEST2} -eq 1 ]; then
                    echo "Error downloading BAM files at <~{bucket_dir}/${BAM_NAME}_${BED_NAME}/>. Trying again..."
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
            while : ; do
                TEST=$(gsutil -m cp './*.bam' ~{bucket_dir}/${BAM_NAME}_${BED_NAME}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading BAM files. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} all.bam
        
        # Printing the union of reads from all individuals (one image per SV).
        IMAGE_HEIGHT="10500"
        HORIZONTAL_SLACK="10000"
        FIGURES_DIR_MERGED="figures_merged_1"
        rm -rf ${FIGURES_DIR_MERGED}; mkdir ${FIGURES_DIR_MERGED}
        IGV_SCRIPT="script.txt"; BAMSNAP_BED="bamsnap.bed"
        echo "new" > ${IGV_SCRIPT}
        echo "snapshotDirectory ${FIGURES_DIR_MERGED}" >> ${IGV_SCRIPT}
        echo "load all.bam" >> ${IGV_SCRIPT}
        echo "genome ~{reference_fa}" >> ${IGV_SCRIPT}
        echo "maxPanelHeight ${IMAGE_HEIGHT}" >> ${IGV_SCRIPT}
        tr '\t' ',' < ~{regions_bed} > regions_prime.txt
        i="0"
        while read SV; do
            CHR=$(echo ${SV} | cut -d , -f 1)
            START=$(echo ${SV} | cut -d , -f 2)
            START=$(( ${START}-${HORIZONTAL_SLACK} ))
            END=$(echo ${SV} | cut -d , -f 3)
            END=$(( ${END}+${HORIZONTAL_SLACK} ))
            echo "goto ${CHR}:${START}-${END}" >> ${IGV_SCRIPT}
            echo "sort base" >> ${IGV_SCRIPT}
            echo "squish" >> ${IGV_SCRIPT}
            echo "snapshot sv-${i}.png" >> ${IGV_SCRIPT}
            echo "${CHR}\t${START}\t${END}" >> ${BAMSNAP_BED}
            i=$(( ${i}+1 ))
        done < regions_prime.txt
        echo "exit" >> ${IGV_SCRIPT}
        cat ${IGV_SCRIPT}
        python /IGV-snapshot-automator/make_IGV_snapshots.py -mem $(( ~{ram_size_gb}-4 )) -onlysnap ${IGV_SCRIPT} all.bam
        tar -czvf report_merged_1.tar.gz ${FIGURES_DIR_MERGED}
        
        # Printing reads from different individuals separately in an HTML report
#        FIGURES_DIR_DISTINCT="figures_distinct"
#        bamsnap -process ${N_THREADS} -ref ~{reference_fa} -bam alignments_*.bam -bed ${BAMSNAP_BED} -out ${FIGURES_DIR_DISTINCT} \
#            -separated_bam \
#            -bamplot read -read_thickness 2 -read_gap_height 0 -read_gap_width 1 \
#            -show_soft_clipped \
#        tar -czvf report_distinct.tar.gz ${FIGURES_DIR_DISTINCT}
        
        # Printing reads from different individuals (one image per SV).
#        FIGURES_DIR_MERGED="figures_merged_2"
#        bamsnap -process ${N_THREADS} -ref ~{reference_fa} -bam alignments_*.bam -bed ${BAMSNAP_BED} -out ${FIGURES_DIR_MERGED} \
#            -bamplot read -read_thickness 2 -read_gap_height 0 -read_gap_width 1 \
#            -show_soft_clipped \
#        tar -czvf report_merged_2.tar.gz ${FIGURES_DIR_MERGED}
        
        
        
        
        
        
        touch report_distinct.tar.gz report_merged_2.tar.gz
        
        #${TIME_COMMAND} create_report ~{regions_bed} ~{reference_fa} --flanking 10000 --exclude-flags 0 --sort BASE --tracks all.bam --output report.html --sequence 1 --begin 2 --end 3 --standalone
    >>>

    output {
        File report_merged_1 = "report_merged_1.tar.gz"
        File report_distinct = "report_distinct.tar.gz"
        File report_merged_2 = "report_merged_2.tar.gz"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: n_cpus
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
