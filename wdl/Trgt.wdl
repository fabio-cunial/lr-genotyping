version 1.0


# 
#
workflow Trgt {
    input {
        String region
        File repeats_bed
        String output_dir
        File bam_addresses
        File reference_fa
        File reference_fai
        Int n_nodes
        Int n_cpus
    }
    parameter_meta {
        output_dir: "Remote directory where results are stored."
        repeats_bed: "The BED file needed by TRGT."
        bam_addresses: "File containing a list of bucket addresses."
        n_nodes: "Use this number of nodes in parallel."
        n_cpus: "Lower bound on the number of CPUs per node."
    }
    
    call GetChunks {
        input:
            bam_addresses = bam_addresses,
            n_chunks = n_nodes
    }
    scatter(chunk_file in GetChunks.chunks) {
        call TrgtImpl {
            input:
                chunk = chunk_file,
                region = region,
                repeats_bed = repeats_bed,
                bucket_dir = output_dir,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                n_cpus = n_cpus
        }
    }
    Array[Int] values = [0, 1,]
    scatter(split_multiallelics in values) {
        call MergeIndividuals {
            input:
                bucket_dir = output_dir,
                split_multiallelics = split_multiallelics,
                n_cpus = n_cpus,
                artificial_input = TrgtImpl.artificial_output
        }
    }
}


# Creates an array of balanced lists of BAM addresses.
#
task GetChunks {
    input {
        File bam_addresses
        Int n_chunks
    }
    parameter_meta {
        bam_addresses: "File containing a list of bucket addresses."
    }
    
    command <<<
        set -euxo pipefail
        
        N_LINES=$(wc -l < ~{bam_addresses})
        N_LINES_PER_CHUNK=$(( (${N_LINES} + ~{n_chunks} - 1) / ~{n_chunks} ))
        split -d -l ${N_LINES_PER_CHUNK} ~{bam_addresses} chunk-
    >>>
    output {
        Array[File] chunks = glob("chunk-*")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


# Builds a VCF for each individual in a chunk.
#
task TrgtImpl {
    input {
        File chunk
        String region
        File repeats_bed
        String bucket_dir
        File reference_fa
        File reference_fai
        Int n_cpus
    }
    parameter_meta {
        chunk: "A list of remote BAM file addresses."
    }
    
    Int ram_size_gb = 64
    Int disk_size_gb = 64
    String docker_dir = "/"
    String work_dir = "/cromwell_root/lr-genotyping"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_DELAY_S="600"
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export RUST_BACKTRACE="full"
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        
        while read REMOTE_FILE; do
            INDIVIDUAL=$( basename ${REMOTE_FILE} .bam )
            TEST=$(samtools view --threads ${N_THREADS} --with-header --bam --output ${INDIVIDUAL}.bam ${REMOTE_FILE} ~{region} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
                samtools view --threads ${N_THREADS} --with-header --bam --output ${INDIVIDUAL}.bam ${REMOTE_FILE} ~{region}
            fi
            samtools index -@ ${N_THREADS} ${INDIVIDUAL}.bam
            ${TIME_COMMAND} ~{docker_dir}trgt --genome ~{reference_fa} --repeats ~{repeats_bed} --reads ${INDIVIDUAL}.bam --output-prefix tmp
            bcftools sort --output-type z --output ${INDIVIDUAL}.vcf.gz tmp.vcf.gz
            rm -f tmp.vcf.gz
            bcftools index --threads ${N_THREADS} ${INDIVIDUAL}.vcf.gz
            bcftools norm --multiallelics -any --multi-overlaps 0 --output-type z --output tmp.vcf.gz ${INDIVIDUAL}.vcf.gz
            bcftools sort --output-type z --output ${INDIVIDUAL}_split.vcf.gz tmp.vcf.gz
            rm -f tmp.vcf.gz
            bcftools index --threads ${N_THREADS} ${INDIVIDUAL}_split.vcf.gz
            samtools sort -@ ${N_THREADS} -o ${INDIVIDUAL}.spanning.bam tmp.spanning.bam
            rm -f tmp.spanning.bam
            samtools index ${INDIVIDUAL}.spanning.bam
            ${TIME_COMMAND} ~{docker_dir}trvz --genome ~{reference_fa} --repeats ~{repeats_bed} --vcf ${INDIVIDUAL}.vcf.gz --spanning-reads ${INDIVIDUAL}.spanning.bam --repeat-id id --image ${INDIVIDUAL}.svg || echo "trvz failed"
            rm -f ${INDIVIDUAL}.spanning.bam ${INDIVIDUAL}.bam
            while : ; do
                TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ${INDIVIDUAL}'*' ~{bucket_dir}/ && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error uploading file <${INDIVIDUAL}*>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            rm -f ${INDIVIDUAL}*
        done < ~{chunk}
    >>>

    output {
        Int artificial_output = 1
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: n_cpus
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Merges all the <.vcf.gz> files produced by trgt.
#
task MergeIndividuals {
    input {
        String bucket_dir
        Int split_multiallelics
        Int n_cpus
        Array[Int] artificial_input
    }
    parameter_meta {
        bucket_dir: "Containing all the <.vcf.gz> files produced by trgt."
        split_multiallelics: "1=split, 0=do not split."
        artificial_input: "Just to force sequentiality."
    }
    
    Int ram_size_gb = 64
    Int disk_size_gb = 128

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Downloading VCFs
        if [ ~{split_multiallelics} -eq 1 ]; then
            SUFFIX="_split"
        else 
            SUFFIX=""
        fi
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp ~{bucket_dir}/'*'${SUFFIX}'.vcf.gz*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading VCF files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        ls *${SUFFIX}.vcf.gz > list.txt
        
        # Merging VCFs
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --missing-to-ref --file-list list.txt --output-type z --output merged${SUFFIX}.vcf.gz --write-index
        while : ; do
            TEST=$(gsutil -m ${GSUTIL_UPLOAD_THRESHOLD} cp merged${SUFFIX}.vcf.gz'*' ~{bucket_dir}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading merged VCF. Trying again..."
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
