version 1.0


# 
#
workflow BAMtracks {
    input {
        String region
        File bam_addresses
        Int window_length
        Int window_step
        Int kmer_length
        Int n_nodes
    }
    parameter_meta {
        bam_addresses: "File containing a list of bucket addresses."
        n_nodes: "Use this number of nodes to regenotype in parallel."
    }
    
    call GetChunks {
        input:
            bam_addresses = bam_addresses,
            n_chunks = n_nodes
    }
    scatter(chunk_file in GetChunks.chunks) {
        call BAMtracksImpl { 
            input:
                chunk = chunk_file,
                region = region,
                window_length = window_length,
                window_step = window_step,
                kmer_length = kmer_length
        }
    }
    call PasteTracks {
        input:
            prefix_file = BAMtracksImpl.prefix[0],
            table_files = BAMtracksImpl.table
    }
    output {
        File out = PasteTracks.out
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


task BAMtracksImpl {
    input {
        File chunk
        String region
        Int window_length
        Int window_step
        Int kmer_length
    }
    parameter_meta {
        chunk: "A list of remote BAM file addresses."
    }

    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        
        touch prefix.txt table.txt
        while read REMOTE_FILE; do
            INDIVIDUAL=$( basename ${REMOTE_FILE} .bam )
            TEST=$(samtools view --threads ${N_THREADS} ${REMOTE_FILE} ~{region} > ${INDIVIDUAL}.sam && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
                ${TIME_COMMAND} samtools view --threads ${N_THREADS} ${REMOTE_FILE} ~{region} > ${INDIVIDUAL}.sam
            fi
            ${TIME_COMMAND} java -cp / BAMtracks ./${INDIVIDUAL}.sam ./${INDIVIDUAL}.tracks ~{window_length} ~{window_step} ~{kmer_length}
            rm -f ${INDIVIDUAL}.sam
            cut -d , -f 1,2 ./${INDIVIDUAL}.tracks > prefix.txt
            cut -d , -f 3,7 ./${INDIVIDUAL}.tracks | paste -d , table.txt -
        done < ~{chunk}
    >>>

    output {
        File prefix = "prefix.txt"
        File table = "table.txt"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        preemptible: 0
    }
}


task PasteTracks {
    input {
        File prefix_file
        Array[File] table_files
    }

    command <<<
        set -euxo pipefail
        
        cat ~{prefix_file} > out.txt
        for FILE in ~{sep=' ' table_files}; do
            paste -d , out.txt ${FILE} > tmp.txt
            mv tmp.txt out.txt
        done
    >>>

    output {
        File out = "out.txt"
    }
    runtime {
        docker: "ubuntu:latest"
    }
}
