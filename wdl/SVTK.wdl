version 1.0


task Standardize {
    input {
        File vcf
        File tbi
        File ref_fai

        String caller
        String prefix
    }

    Int disk_size = 8*ceil(size([vcf, tbi, ref_fai], "GB")) + 1

    command <<<
        set -euxo pipefail

        svtk standardize \
            --include-reference-sites \
            --contigs ~{ref_fai} \
            --prefix ~{prefix} ~{vcf} - ~{caller} | \
            bcftools sort /dev/stdin -o ~{prefix}.truvari.std.vcf.gz -O z

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