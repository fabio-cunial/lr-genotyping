version 1.0


# Refines the genotypes of a population.
#
workflow Regenotype {
    input {
        File merged_vcf
        File bam_addresses
        Int use_lrcaller
        Int use_cutesv
        File reference_fa
        File reference_fai
        Int n_nodes
        Int n_cpus
        Int bam_size_gb
        String backup_address
    }
    parameter_meta {
        merged_vcf: "The output of the merging step, whose genotypes must be refined. Can be either .vcf or .vcf.gz."
        bam_addresses: "File containing a list of bucket addresses."
        n_nodes: "Use this number of nodes to regenotype in parallel."
        n_cpus: "Lower bound on the number of CPUs per regenotype node."
        bam_size_gb: "Upper bound on the size of a single BAM."
    }
    
    call GetVcfToGenotype {
        input:
            merged_vcf = merged_vcf
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
                use_lrcaller = use_lrcaller,
                use_cutesv = use_cutesv,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                n_cpus = n_cpus,
                bam_size_gb = bam_size_gb,
                backup_address = backup_address
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


# Outputs a .vcf file that equals $merged_vcf$ but lacks any sample information.
#
task GetVcfToGenotype {
    input {
        File merged_vcf
    }
    parameter_meta {
        merged_vcf: "The output of the merging step, whose genotypes must be refined. Can be either .vcf or .vcf.gz."
    }
    
    Int disk_size_gb = 2*ceil(size(merged_vcf, "GB"))

    command <<<
        set -euxo pipefail
        
        bcftools view -h ~{merged_vcf} > vcf_header.txt
        N_ROWS=$(wc -l < vcf_header.txt)
        head -n $(( ${N_ROWS} - 1 )) vcf_header.txt > variants_only.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> variants_only.vcf
        bcftools view -H ~{merged_vcf} | cut -f 1,2,3,4,5,6,7,8 >> variants_only.vcf
    >>>
    
    output {
        File vcf_to_genotype = "variants_only.vcf"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
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
        Int use_lrcaller
        Int use_cutesv
        File reference_fa
        File reference_fai
        Int n_cpus
        Int bam_size_gb
        String backup_address
    }
    parameter_meta {
        bam_size_gb: "Upper bound on the size of a single BAM."
    }
    
    Int ram_size_gb = n_cpus*8 + 2*( ceil(size(vcf_to_genotype,"GB")) + bam_size_gb + ceil(size(reference_fa,"GB")) )
    Int disk_size_gb = bam_size_gb + 10*ceil(size(vcf_to_genotype,"GB")) + ceil(size(reference_fa,"GB"))

    command <<<
        set -euxo pipefail
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        touch format.txt genotypes.txt
        i="0"
        while read BAM_FILE; do
            while : ; do
                TEST=$(gsutil -m cp ${BAM_FILE} ${BAM_FILE}.bai . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${BAM_FILE}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            if [ ~{use_lrcaller} -eq 1 ]; then
                ${TIME_COMMAND} LRcaller --number_of_threads ${N_THREADS} -fa ~{reference_fa} --genotyper joint $(basename ${BAM_FILE}) ~{vcf_to_genotype} genotypes.vcf
            fi
            if [ ~{use_cutesv} -eq 1 ]; then
                mkdir ./cutesv_tmp
                ${TIME_COMMAND} cutesv --threads ${N_THREADS} -Ivcf ~{vcf_to_genotype} --max_ cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.8 -mi 500 -md 500 -s 3 --genotype -L -1 $(basename ${BAM_FILE}) ~{reference_fa} genotypes.vcf ./cutesv_tmp
            fi
            
            
            gsutil -m cp genotypes.vcf ~{backup_address}/genotypes.vcf
            
            
            
            rm -f $(basename ${BAM_FILE}) $(basename ${BAM_FILE}).bai
            echo "FORMAT" > format.txt
            bcftools view -H genotypes.vcf | cut -f 9 >> format.txt
            INDIVIDUAL=$(basename ${BAM_FILE} -s .bam)
            echo ${INDIVIDUAL} > new_genotypes.txt
            bcftools view -H genotypes.vcf | cut -f 10 >> new_genotypes.txt
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

    Int disk_size_gb = 2*( ceil(size(vcf_to_genotype,"GB")) + ceil(size(genotypes,"GB")) )

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Building the VCF header (excluding the column title line).
        bcftools view -h ~{vcf_to_genotype} > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > out_header.txt
        
        # Building the VCF body (including the column title line).
        tail -n 1 header.txt > out_body.txt
        bcftools view -H ~{vcf_to_genotype} >> out_body.txt
        
        # Adding the FORMAT column (assumed to be the same for all chunks).
        paste out_body.txt ~{format} > tmp.txt
        mv tmp.txt out_body.txt
        
        # Adding the genotype columns
        echo ~{sep='-' genotypes} | tr '-' '\n' > genotypes.txt
        while read FILE; do
            paste out_body.txt ${FILE} > tmp.txt
            mv tmp.txt out_body.txt
        done < genotypes.txt
        
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
