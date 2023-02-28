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
        File report = Sv2IgvImpl.report
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
    
    Int ram_size_gb = 64  # Arbitrary
    Int disk_size_gb = 100*bam_size_gb + ceil(size(reference_fa,"GB"))

    command <<<
        set -euxo pipefail
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        TEST1=$(gsutil -q stat ~{bucket_dir}/all.bam && echo 0 || echo 1)
        if [ ${TEST1} -eq 0 ]; then
            while : ; do
                TEST2=$(gsutil -m cp ~{bucket_dir}/all.bam . && echo 0 || echo 1)
                if [ ${TEST2} -eq 1 ]; then
                    echo "Error downloading file <~{bucket_dir}/all.bam>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        else
            i="0"; FILES=""
            while read BAM_FILE; do                
                samtools view --threads ${N_THREADS} --target-file ~{regions_bed} --reference ~{reference_fa} --fai-reference ~{reference_fai} --bam --output alignments_${i}.bam ${BAM_FILE}
                FILES="${FILES} alignments_${i}.bam"
                i=$(( $i + 1 ))
            done < ~{bam_addresses}
            ${TIME_COMMAND} samtools merge -@ ${N_THREADS} -o all.bam ${FILES}
            rm -f ${FILES}
            while : ; do
                TEST2=$(gsutil -m cp all.bam ~{bucket_dir}/all.bam && echo 0 || echo 1)
                if [ ${TEST2} -eq 1 ]; then
                    echo "Error uploading file <all.bam>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
        fi
        samtools index -@ ${N_THREADS} all.bam
        ${TIME_COMMAND} create_report ~{regions_bed} ~{reference_fa} --flanking 10000 --exclude-flags 0 --sort BASE --tracks all.bam --output report.html --sequence 1 --begin 2 --end 3 --standalone
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
