version 1.0


# 
#
workflow LocalAssembly {
    input {
        String region
        File bam_addresses
        Int k
        Int n_cpus
    }
    parameter_meta {
    }
    
    call LocalAssemblyImpl {
        input:
            bams_list = bam_addresses,
            region = region,
            k = k,
            n_cpus = n_cpus
    }
}


#
task LocalAssemblyImpl {
    input {
        File bams_list
        String region
        Int k
        Int n_cpus
    }
    parameter_meta {
        bams_list: "A list of remote BAM file addresses."
    }
    
    Int ram_size_gb = 64
    Int disk_size_gb = 256
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
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        BIFROST_COMMAND="/bifrost/build/src/Bifrost"
        
        # Downloading all BAMs
        while read REMOTE_FILE; do
            INDIVIDUAL=$( basename ${REMOTE_FILE} .bam )
            TEST=$(samtools view --threads ${N_THREADS} --with-header --bam --output ${INDIVIDUAL}.bam ${REMOTE_FILE} ~{region} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
                samtools view --threads ${N_THREADS} --with-header --bam --output ${INDIVIDUAL}.bam ${REMOTE_FILE} ~{region}
            fi
            samtools fastq ${INDIVIDUAL}.bam > ${INDIVIDUAL}.fastq
            rm -f ${INDIVIDUAL}.bam
        done < ~{bams_list}
        cat *.fastq > all.fastq
        
        # Assembling all BAMs
        ${TIME_COMMAND} ${BIFROST_COMMAND} build --threads ${N_THREADS} --kmer-length ~{k} --input-seq-file all.fastq --output-file bifrost
        ls -laht
        #${TIME_COMMAND} hifiasm -t ${N_THREADS} -o all all.fastq
        #tar -czf all.tar.gz *.gfa
    >>>

    output {
        File assembly = work_dir + "/bifrost.gfa.gz"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: n_cpus
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
