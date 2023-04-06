version 1.0


# 
#
workflow JediGraph {
    input {
        String merged_vcf_gz
        String region
        File bam_addresses
        String reference_fa
        Int n_nodes
        Int n_cpus
    }
    parameter_meta {
        merged_vcf_gz: "Remote address. Might cover more than $region$."
        bam_addresses: "File containing a list of bucket addresses."
        n_nodes: "Use this number of nodes to regenotype in parallel."
        n_cpus: "Lower bound on the number of CPUs per regenotype node."
    }
    
    call GetVcfToGenotype {
        input:
            merged_vcf_gz = merged_vcf_gz,
            reference_fa = reference_fa,
            region = region
    }
    call GetChunks {
        input:
            bam_addresses = bam_addresses,
            n_chunks = n_nodes
    }
    scatter(chunk_file in GetChunks.chunks) {
        call RegenotypeChunk { 
            input:
                chunk = chunk_file,
                vcf_to_genotype = GetVcfToGenotype.vcf_to_genotype,
                region = region,
                graph = GetVcfToGenotype.gfa,
                n_cpus = n_cpus
        }
    }
    call PasteGenotypedChunks {
        input:
            vcf_to_genotype = GetVcfToGenotype.vcf_to_genotype,
            genotypes = RegenotypeChunk.genotypes,
            format = RegenotypeChunk.format[0]
    }
    output {
        File vcf_gz = PasteGenotypedChunks.vcf_gz
        File tbi = PasteGenotypedChunks.tbi
    }
}


# Outputs: (1) a .vcf file that equals $merged_vcf_gz$ but lacks any sample
# information and is limited to a given region of the genome; (2) .fa and .fai
# files that contain just $region$; (3) the .gfa graph built by svjedi.
#
# Remark: the program that transforms a VCF into a GFA distributed with
# SVJedi-graph is sequential.
#
task GetVcfToGenotype {
    input {
        String merged_vcf_gz  # Remote address
        String reference_fa  # Remote address
        File reference_fai
        String region
    }
    
    Int disk_size_gb = 2*ceil(size(merged_vcf_gz, "GB")) + 6
    Int ram_size_gb = disk_size_gb

    command <<<
        set -euxo pipefail
        
        SVJEDI_PATH="/opt/conda/bin"
        TIME_COMMAND="/usr/bin/time --verbose"
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        
        bcftools view --with-header --output-type v --output filtered.vcf --regions-overlap variant --regions ~{region} ~{merged_vcf_gz}
        grep '#' filtered.vcf > vcf_header.txt
        N_ROWS=$(wc -l < vcf_header.txt)
        head -n $(( ${N_ROWS} - 1 )) vcf_header.txt > variants_only.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> variants_only.vcf
        tail -n +$((${N_ROWS} + 1)) filtered.vcf | cut -f 1,2,3,4,5,6,7,8 >> variants_only.vcf
        samtools faidx ~{reference_fa} ~{region} > new_reference.fa
        samtools faidx new_reference.fa
        ${TIME_COMMAND} python ${SVJEDI_PATH}/construct-graph.py --vcf variants_only.vcf --ref new_reference.fa -o graph.gfa
    >>>
    
    output {
        File vcf_to_genotype = "variants_only.vcf"
        File new_reference_fa = "new_reference.fa"
        File new_reference_fai = "new_reference.fai"
        File gfa = "graph.gfa"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: 1  # <construct-graph.py> is sequential
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
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


# Genotypes $vcf_to_genotype$ using every remote address of a BAM in list
# $chunk$. The task returns a file with all and only the genotype columns that
# will have to be added to $vcf_to_genotype$ to create the final VCF (such
# columns include the heading with sample ID).
#
task RegenotypeChunk {
    input {
        File chunk
        File vcf_to_genotype
        String region
        File graph
        Int n_cpus
    }
    parameter_meta {
        chunk: "A list of remote BAM file addresses."
    }
    
    Int ram_size_gb = 10*ceil( size(vcf_to_genotype,"GB") + 3 + size(graph,"GB") )
    Int disk_size_gb = ram_size_gb

    command <<<
        set -euxo pipefail
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        
        mv ~{graph} graph.gfa
        touch format.txt genotypes.txt
        i="0"
        while read REMOTE_FILE; do
            # Computing genotypes
            INDIVIDUAL=$( basename ${REMOTE_FILE} .bam )
            TEST=$(samtools view --threads ${N_THREADS} --with-header --bam --output ${INDIVIDUAL}.bam ${REMOTE_FILE} ~{region} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
                samtools view --threads ${N_THREADS} --with-header --bam --output ${INDIVIDUAL}.bam ${REMOTE_FILE} ~{region}
            fi
            samtools fastq --threads ${N_THREADS} ${INDIVIDUAL}.bam > ${INDIVIDUAL}.fastq
            rm -f ${INDIVIDUAL}.bam
            ${TIME_COMMAND} svjedi-graph.py --threads ${N_THREADS} --vcf ~{vcf_to_genotype} --ref unused --reads ${INDIVIDUAL}.fastq --prefix graph
            rm -f ${INDIVIDUAL}.fastq
            mv graph_genotype.vcf genotypes.vcf
            rm -f *.gaf *.json
            
            # Appending genotypes
            N_LINES=$(grep '#' genotypes.vcf | wc -l)
            echo "FORMAT" > format.txt
            tail -n +$(( ${N_LINES} + 1 )) genotypes.vcf | cut -f 9 >> format.txt
            echo ${INDIVIDUAL} > new_genotypes.txt
            GENOTYPE_COLUMN="10"
            tail -n +$(( ${N_LINES} + 1 )) genotypes.vcf | cut -f ${GENOTYPE_COLUMN} >> new_genotypes.txt
            rm -f genotypes.vcf
            if [ $i -eq 0 ]; then
                mv new_genotypes.txt genotypes.txt
                i="1"
            else
                paste genotypes.txt new_genotypes.txt > tmp.txt
                mv tmp.txt genotypes.txt; rm -f tmp.txt
            fi
            head genotypes.txt
        done < ~{chunk}
    >>>

    output {
        File format = "format.txt"
        File genotypes = "genotypes.txt"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: n_cpus
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Adds the new genotypes of every chunk as columns of a new VCF file.
#
task PasteGenotypedChunks {
    input {
        File vcf_to_genotype
        Array[File] genotypes
        File format
    }
    parameter_meta {
        vcf_to_genotype: "The .vcf file that was used as input to $RegenotypeChunk$."
        format: "Any of the format files returned by $RegenotypeChunk$ (they are assumed to be all equal)."
    }

    Int disk_size_gb = 10*( ceil(size(vcf_to_genotype,"GB")) + ceil(size(genotypes,"GB")) )
    String first_genotyped_file = genotypes[0]

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Building the VCF header (excluding the column title line).
        grep '#' ~{vcf_to_genotype} > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > out_header.txt
        
        # Building the VCF body (including the column title line).
        tail -n 1 header.txt > out_body.txt
        tail -n +$(( ${N_ROWS} + 1 )) ~{vcf_to_genotype} >> out_body.txt
        
        # Adding the FORMAT column (assumed to be the same for all chunks).
        paste out_body.txt ~{format} > tmp.txt
        mv tmp.txt out_body.txt
        
        # Adding the genotype columns
        for FILE in ~{sep=' ' genotypes}; do
            paste out_body.txt ${FILE} > tmp.txt
            mv tmp.txt out_body.txt
        done
        
        # Building the output VCF
        cat out_header.txt out_body.txt > out.vcf
        bgzip --threads ${N_THREADS} out.vcf
        tabix out.vcf.gz
    >>>

    output {
        File vcf_gz = "out.vcf.gz"
        File tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: 1
        memory: "8 GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
