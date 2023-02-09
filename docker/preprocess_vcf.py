# Converts DUP to INS and retrieves sequence for ref/alt
# If DEL ref/alt is symbolic (as in Sniffles2 output), retrieves sequence for that as well

import sys
from pysam import VariantFile, FastaFile

def main():
    if len(sys.argv)!=3:
        sys.stderr.write("Usage: program.py vcf ref_fasta\n")
        sys.exit(1)
    vcf_in = VariantFile(sys.argv[1])
    vcf_out = VariantFile('-', 'w', header=vcf_in.header)
    ref_fa = FastaFile(sys.argv[2])

    for rec in vcf_in.fetch():
        # check whether SVLEN is reasonable length for inserting sequence data (<=100kb)
        if 'SVLEN' in rec.info.keys(): # missing for BNDs
            if isinstance(rec.info['SVLEN'], tuple): # this is to work around a confusing bug where sometimes SVLEN is a tuple even though there's only one (e.g. SVLEN=500 becomes (500,))
                svlen = rec.info['SVLEN'][0]
            else:
                svlen = rec.info['SVLEN']
            if abs(svlen) > 100000: # don't fill in sequences longer than 100 kb
                getseq = False
            else:
                getseq = True
        else:
            getseq = False

        # DUP-->INS
        if rec.info['SVTYPE']=='DUP':
            # retrieve sequences from indexed reference fasta
            # doing -1 because pysam coordinates (i.e. input to fetch command) are 0-based, vcf is 1-based (tested this to confirm)
            if getseq:
                rec.alts = (ref_fa.fetch(rec.chrom, rec.pos-1, rec.stop-1),)
                rec.ref = ref_fa.fetch(rec.chrom, rec.pos-1, rec.pos)
            else:
                rec.alts = (rec.alts[0].replace('DUP','INS'),)

            rec.info['SVTYPE']='INS'
            rec.stop = rec.pos+1

        # add sequence for symbolic DELs
        if rec.info['SVTYPE']=='DEL':
            if rec.alts[0]=="<DEL>" and getseq:
                rec.ref = ref_fa.fetch(rec.chrom, rec.pos-1, rec.stop)
                rec.alts = (ref_fa.fetch(rec.chrom, rec.pos-1, rec.pos),)

        vcf_out.write(rec) # write every variant

if __name__=="__main__":main()
