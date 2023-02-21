version 1.0


workflow Standardize {
    input {
        File vcf
        File tbi
        File ref_fai
        String caller
        String prefix
    }
    call Standardize_impl {
        input:
            vcf = vcf,
            tbi = tbi,
            ref_fai = ref_fai,
            caller = caller,
            prefix = prefix
    }
    output {
        File standardized_vcf = Standardize_impl.standardized_vcf
        File standardized_tbi = Standardize_impl.standardized_tbi
    }
}


task Standardize_impl {
    input {
        File vcf
        File tbi
        File ref_fai
        String caller
        String prefix
    }

    Int disk_size = 20*ceil(size([vcf, tbi, ref_fai], "GB")) + 1

    command <<<
        set -euxo pipefail
        
        svtk standardize \
            --include-reference-sites \
            --contigs ~{ref_fai} \
            --prefix ~{prefix} ~{vcf} - ~{caller} > standardized.vcf || echo "SVTK exited with an error code"
        tail -n 100 standardized.vcf
        bcftools sort standardized.vcf -o ~{prefix}.truvari.std.vcf.gz -O z
        tabix ~{prefix}.truvari.std.vcf.gz
    >>>

    output {
        File standardized_vcf = "~{prefix}.truvari.std.vcf.gz"
        File standardized_tbi = "~{prefix}.truvari.std.vcf.gz.tbi"
    }
    runtime {
        cpu: 4
        memory: "24 GiB"
        disks: "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        docker: "fcunial/svtk"
    }
}
