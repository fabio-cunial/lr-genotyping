#!/usr/bin/env python3

"""*******************************************************************************
    Name: SVJedi-graph
    Description: SVjedi-graph aims to genotype structural variant with long reads data using a variation graph.
    Authors: Lolita Lecompte, Sandra Romain
    Contact: sandra.romain@inria.fr, Inria/Univ Rennes/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
    
    Copyright (C) 2019 Inria
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************"""

import sys
#import numpy as np
import argparse
import math
from decimal import *
import json
import scipy.special as special

def main(args):
    """ Parsing arguments """
    parser = argparse.ArgumentParser(description="Structural variations genotyping using long reads")
    
    parser.add_argument("-d", "--aln", metavar="<alndict>", nargs=1, required=True)
    
    parser.add_argument("-v", "--vcf", metavar="<vcffile>", help="vcf format", required=True)
    
    parser.add_argument("-o", "--output", metavar="<output>", nargs=1, help="output file")

    # parser.add_argument("-i", "--minid", metavar="<min identity>", nargs=1, type=float, help="minimum pct of identity")

    parser.add_argument("-e", "--err", nargs=1, type=float, help="allele error probability")

    parser.add_argument("-ladj", metavar="<allele_size>", nargs=1, type=int, default=[5000], help="Sequence allele adjacencies at each side of the SV")

    parser.add_argument(
        "-ms",
        "--minsupport",
        metavar="<minNbAln>",
        type=int,
        default=3,
        help="Minimum number of alignments to genotype a SV (default: 3>=)",
    )

    args = parser.parse_args()

    # check if outputfile
    if args.output is None:
        output = "genotype_results.txt"
    else:
        output = args.output[0]
    
    # if args.minid is not None:
    #     min_id = args.minid[0]
    # else:
    #     min_id = 0
    
    if args.err is not None:
        e = args.err[0]
    else:
        e = 0.00005

    vcffile = args.vcf
    aln_dict_file = args.aln[0]
    min_support = args.minsupport
    l_adj = args.ladj[0]

    with open(aln_dict_file, 'r') as file:
        dict_of_informative_aln = json.load(file)

    missing_id = []

    missing_id = decision_vcf(dict_of_informative_aln, vcffile, output, min_support, e, l_adj, missing_id)

    # print(len(missing_id))
    # print(missing_id[:10])

def get_info(info, label):
    #Info label first in info field
    if info.split(";")[0].startswith(label+"="):
        return info.split(label+"=")[1].split(";")[0]

    #Info label last in info field
    elif info.split(";")[-1].startswith(label+"="):
        return info.split(label.join([";", "="]))[1]

    else:
        return info.split(label.join([";", "="]))[1].split(";")[0]

def decision_vcf(dictReadAtJunction, inputVCF, outputDecision, minNbAln, e, l_adj, missing_id):
    """ Output in VCF format and take genotype decision """
    getcontext().prec = 28
    outDecision = open(outputDecision, "w")
    
    genotype_format = "GT:DP:AD:PL"
    genotyped_svs = 0
    ungenotyped_svs = [0, []]

    with open(inputVCF) as inputFile:
        for line in inputFile:
            if line.startswith("##FORMAT"):
            	continue
            
            if line.startswith("##"):
                outDecision.write(line)

            elif line.startswith("#C"):
                outDecision.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                outDecision.write('##FORMAT=<ID=DP,Number=1,Type=Float,Description="Total number of informative read alignments across all alleles (after normalization for unbalanced SVs)">\n')
                outDecision.write('##FORMAT=<ID=AD,Number=2,Type=Float,Description="Number of informative read alignments supporting each allele (after normalization by breakpoint number for unbalanced SVs)">\n')
                
                outDecision.write('##FORMAT=<ID=PL,Number=3,Type=Integer,Description="Phred-scaled likelihood for each genotype">\n')
                #outDecision.write(line.rstrip("\n") + "\t" + "\t".join(["FORMAT", "SAMPLE"]) + "\n")
                outDecision.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n")
                
            else:
                in_chrom, in_start, _, __, in_type, ___, ____, in_info, *_ = line.rstrip("\n").split("\t")

                ### get SVTYPE ###
                if 'SVTYPE' in in_info:
                    if in_info.split(';')[-1].startswith('SVTYPE='):
                        svtype = in_info.split('SVTYPE=')[1]
                    else:
                        svtype = in_info.split('SVTYPE=')[1].split(';')[0]
                else:
                    svtype = ''
                
                ### get END ###
                #end = line.split(";END=")[1].split(";")[0]
                #if svtype != "BND":
                #    end = get_info(in_info, "END")

                
                ### get LENGTH for DELETION ###                 
                if svtype == 'DEL':
                    if "SVLEN=FALSE" in in_info:
                        if in_info.startswith("END="):
                            end = in_info.split("END=")[1].split(";")[0]
                        else:
                            end = in_info.split(";END=")[1].split(";")[0]

                        in_length = int(end) - int(in_start)

                    elif "SVLEN=" in in_info:
                        if in_info.split(';')[-1].startswith('SVLEN='):
                            in_length = abs(int(in_info.split("SVLEN=")[1]))
                        else:
                            in_length = abs(int(in_info.split("SVLEN=")[1].split(";")[0]))
                
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    # in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV
                    in_sv = in_chrom + ":" + in_start + "-" + str(int(in_start)+in_length-1)
                    
                
                ### get LENGTH for INSERTION ###                    
                elif svtype == 'INS': 
                    if "SVLEN=" in in_info:
                        if in_info.split(';')[-1].startswith('SVLEN='):
                            in_length = int(in_info.split("SVLEN=")[1])
                        else:
                            in_length = int(in_info.split("SVLEN=")[1].split(";")[0])

                    else:
                        in_length = len(in_type) 
                    
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    # in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV
                    in_sv = in_chrom + ":" + in_start + "-" + in_start
                
                
                ### get LENGTH for INVERSION ###                
                elif svtype == 'INV':
                    if in_info.startswith("END="):
                        end = in_info.split("END=")[1].split(';')[0]
                    else:
                        end = in_info.split(";END=")[1].split(';')[0]               
                    #in_length = int(end) - int(in_start)
                    in_length = int(in_info.split("SVLEN=")[1])
                    
                    #if abs(in_length) < 50: continue #focus on svlength of at least 50 bp
                    # in_sv = in_chrom + "_" + in_start + "-" + str(in_length) #define sv id for DEL, INS, INV
                    # in_sv = in_chrom + "_" + in_start + "-" + str(end)
                    #in_sv = in_chrom + ":" + in_start + "-" + str(end)
                    in_sv = in_chrom + ":" + in_start + "-" + str(int(in_start)+in_length-1)


                ### get sv id for TRANSLOCATION ###             
                elif svtype == 'BND': 

                    in_length = 50
                    
                    # if in_info.startswith("END="):
                    #     end = in_info.split("END=")[1].split(';')[0]
                    # elif ";END=" in in_info:
                    #     end = in_info.split(";END=")[1].split(';')[0]   
                    # elif "[" in  in_type:
                    #     end = in_type.split(':')[1].split('[')[0]
                    # elif "]" in  in_type:
                    #     end = in_type.split(':')[1].split(']')[0]
                    
                    # if 'CHR2=' in in_info: 
                    #     chr2 = in_info.split('CHR2=')[1].split(';')[0]
                    # elif '[' in in_type:
                    #         chr2 = in_type.split(':')[0].split('[')[1]     #ALT[CHR2:POSTION[
                    # elif ']' in in_type:
                    #         chr2 = in_type.split(':')[0].split(']')[1]     #ALT]CHR2:POSTION]

                    if "[" in in_type:
                        alt = list(filter(bool, in_type.split("[")))

                        if ":" in alt[1]:
                            chr2 = alt[1].split(":")[0]
                            end = alt[1].split(":")[1]
                            in_sv = ":".join([in_chrom, in_start]) + "[" + alt[1] + "["
                        else:
                            chr2 = alt[0].split(":")[0]
                            end = alt[0].split(":")[1]
                            in_sv = "[" + alt[1] + "[" + ":".join([in_chrom, in_start])
                    
                    elif "]" in in_type:
                        alt = list(filter(bool, in_type.split("]")))

                        if ":" in alt[1]:
                            chr2 = alt[1].split(":")[0]
                            end = alt[1].split(":")[1]
                            in_sv = ":".join([in_chrom, in_start]) + "]" + alt[1] + "]"
                        else:
                            chr2 = alt[0].split(":")[0]
                            end = alt[0].split(":")[1]
                            in_sv = "]" + alt[1] + "]" + ":".join([in_chrom, in_start])
                    
                    else:
                        in_sv = "wrong_format"

                            
                    # in_sv = in_chrom + "_" + in_start + "-" + chr2 + "-" + end #define sv id for TRANS
                
                else:
                    in_sv = in_chrom + ":" + in_start + "-" + str(end)
                
                #######################################################################################
                #Asign genotype 

                # in_sv = in_chrom + ":" + in_start + "-" + str(end)
                if svtype in ('DEL', 'INS', 'INV', 'BND') and in_sv in list(dictReadAtJunction.keys()) and abs(in_length) >= 50:

                    #-------------------#
                    allele_alns = [[], []]
                    for a in [0, 1]:
                        for aln in dictReadAtJunction[in_sv][a]:
                            # if float(aln.split("\tid:f:")[1].split("\t")[0]) >= min_id:
                            allele_alns[a].append(aln)
                    #-------------------#

                    nbAln = [len(x) for x in allele_alns]
                    genotype, proba = likelihood(nbAln, svtype, in_length, minNbAln, l_adj, e)

                    genotyped_svs += 1
            
                else: #if svtype different from DEL, INS, INV, BND or if sv not in supported by alignment

                    if in_sv not in list(dictReadAtJunction.keys()):

                        ungenotyped_svs[1].append(in_sv)
                        
                        #Try find version pb###############################################################
                        # corr_in_sv = in_chrom + "_" + in_start + "-" + in_start
                        # if corr_in_sv in list(dictReadAtJunction.keys()):
                        #     #-------------------#
                        #     #Sup. filter alns (to be moved to filter script later)
                        #     allele_alns = [[], []]
                        #     for a in [0, 1]:
                        #         for aln in dictReadAtJunction[corr_in_sv][a]:
                        #             if float(aln.split("\tid:f:")[1].split("\t")[0]) >= min_id:
                        #                 allele_alns[a].append(aln)
                        #     #-------------------#

                        #     nbAln = [len(x) for x in allele_alns]
                        #     genotype, proba = likelihood(nbAln, svtype, in_length, minNbAln, l_adj)
                        ###################################################################################

                        # else:
                        #     missing_id.append(in_sv)

                    nbAln = [0,0]
                    genotype = "./."
                    proba = [".",".","."]

                    ungenotyped_svs[0] += 1

                    
                #######################################################################################
                #Output genotype in VCF
                
                numbers = ",".join(str(y) for y in nbAln)
                
                if len(line.split("\t")) <= 8:    
                	# input VCF does not contain any genotype information, just adding our genotype at the end of the input line
                    line_before_genotype = line.rstrip("\n")
                else:
                	# input VCF does contain at least one genotype information, we replace the old genotype information by ours
                    fields_without_genotype = line.split("\t")[0:8]
                    line_before_genotype = "\t".join(fields_without_genotype)
                    
                new_line = (
                        line_before_genotype
                        + "\t"
                        + genotype_format
                        + "\t"
                        + genotype
                        + ":"
                        + str(round(sum(nbAln), 3))
                        + ":"
                        + str(numbers)
                        + ":"
                        + str(','.join(proba))
                    )
                outDecision.write(new_line + "\n")

    outDecision.close()

    print("Genotyped svs: " + str(genotyped_svs))
    #print("Ungenotyped svs: " + str(ungenotyped_svs[0]))
    #print(ungenotyped_svs[1][:10])

    return missing_id

def likelihood(all_count, svtype, svlength, minNbAln, l_adj, e):
    """ Compute likelihood """
    
    unbalanced_sv = ("DEL", "INS")
    if svtype in unbalanced_sv:
        c1, c2 = allele_normalization(all_count, svtype)  # allelic sequence normalization for unbalanced SV
        # c1, c2 = all_count
    else:
        c1, c2 = all_count
    
    rc1 = int(round(c1,0))
    rc2 = int(round(c2,0))
    # e = 0.00005 #sequencing err

    lik0 = Decimal(c1*math.log10(1-e)) + Decimal(c2*math.log10(e)) 
    lik1 = Decimal((c1+c2)*math.log10(1/2)) 
    lik2 = Decimal(c2*math.log10(1-e)) + Decimal(c1*math.log10(e))
    
    L = [lik0, lik1, lik2]
    
    index_of_L_max = [i for i, x in enumerate(L) if x == max(L)]
    if len(index_of_L_max) == 1: 
        geno_not_encoded = str(index_of_L_max[0])
        geno = encode_genotype(geno_not_encoded)
        
    else:
        geno = "./."    #no genotype estimation since likelihood are not conclusive
    
    #Check for minimum number of alignment to assign a genotype
    if not sum(all_count) >= minNbAln: # minimum support
        geno = "./."
                
    combination = Decimal(math.log10(special.comb(rc1 + rc2, rc1)))
    lik0 += combination 
    lik1 += combination
    lik2 += combination

    #phred scaled score
    prob0 = -10*lik0
    prob1 = -10*lik1
    prob2 = -10*lik2                                 
                 
    prob = [str(int(prob0)), str(int(prob1)), str(int(prob2))]
            
    return geno, prob

# def allele_normalization(nb_aln_per_allele, svtype, svlength, l_adj):
#     ''' Allele length normalization '''
#     if svlength > (2*l_adj): svlength = 2*l_adj #for upper bound, if case of sv size > 2XLadj, only 2 sequences of 2Ladj are represented
    
#     if svtype == "DEL":
#         nb_aln_longest_allele_seq = nb_aln_per_allele[0]
#         if nb_aln_longest_allele_seq > 0:
#             nb_aln_per_allele[0] = round(nb_aln_longest_allele_seq * (2*l_adj) / ((2*l_adj) + svlength), 3)
    
#     elif svtype == "INS":
#         nb_aln_longest_allele_seq = nb_aln_per_allele[1]
#         if nb_aln_longest_allele_seq > 0:
#             nb_aln_per_allele[1] = round(nb_aln_longest_allele_seq * (2*l_adj) / ((2*l_adj) + svlength), 3)
    
#     return nb_aln_per_allele

def allele_normalization(nb_aln_per_allele, svtype):
    ''' Allele breakpoint number normalization '''
    if svtype == "DEL":
        allele_2bkpt = 0
    elif svtype == "INS":
        allele_2bkpt = 1

    nb_aln_2bkpt_allele = nb_aln_per_allele[allele_2bkpt]
    if nb_aln_2bkpt_allele > 0:
        nb_aln_per_allele[allele_2bkpt] = round(nb_aln_2bkpt_allele / 2, 1)
    
    return nb_aln_per_allele

def encode_genotype(g):
    ''' Encode genotype from 0, 1, 2 to 0/0, 0/1, 1/1 '''
    if g == '0': genotype = "0/0"
    elif g == '1': genotype = "0/1"
    elif g == '2': genotype = "1/1"
    
    return genotype

if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")
    else:
        main(sys.argv[1:])
