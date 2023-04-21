#!/usr/bin/env python
import numpy
import argparse
import scipy.stats
import gzip

def non_param_test(total_count,alt_count,n,all_variants):
    p=0
    selected_vars=numpy.random.randint(low=0,high=len(all_variants),size=n)
    snp_alt_ratio=alt_count/float(total_count)
    ref_count=total_count-alt_count
    diff=abs(alt_count-ref_count)/float(total_count)

    for var in selected_vars:
        alt_ratio=all_variants[var][0]/float(all_variants[var][1])
        sim_ref=all_variants[var][1]-all_variants[var][0]
        sim_diff=abs(sim_ref-all_variants[var][0])/float(all_variants[var][1])

        if diff < sim_diff:
            p+=1
    p=p/float(n)
    return(p)

parser = argparse.ArgumentParser("""Add bootstrap and binomial P values to a vcf file, the output vcf is printed to stdout""")
parser.add_argument('--vcf'        ,required = True, type=str, help="vcf file")
parser.add_argument('--ase'        ,required = True, type=str, help="GATK-ASE file (needs to be generated using the input vcf)")
parser.add_argument('-m'        ,type=int, default=20, help="minimum variant support")
parser.add_argument('-r'        ,type=float, default=0.15, help="minimum normal alt/ref ratio")
parser.add_argument('-R'        ,type=float, default=0.85, help="maximum normal alt/ref ratio")
parser.add_argument('-n'        ,type=int,default=1000, help="permutations")
args = parser.parse_args()

ase_list={}
all_variants=[]


first=True
for line in open(args.ase):
    if first:
        first=False
        continue
    content=line.strip().split()
    if not content[0] in ase_list:
        ase_list[content[0]] = {}
    if not content[1] in ase_list[content[0]]:
        ase_list[content[0]][content[1]]={}

    p_bin=scipy.stats.binom_test(int(content[6]), n=int(content[7]), p=0.5)
    ase_list[content[0]][content[1]][content[4]]={ "ref_count":int(content[5]),"alt_count":int(content[6]),"tot_count":int(content[7]),"p_bin":p_bin,"non_param":-1}

    if not int(content[7]):
       continue
    if int(content[6])/float(content[7]) < args.r or int(content[6])/float(content[7]) > args.R or args.m > int(content[7]):
       continue

    all_variants.append([int(content[6]),int(content[7])])

all_variants=numpy.array(all_variants)

for chromosome in ase_list:
    for pos in ase_list[chromosome]:
        for alt in ase_list[chromosome][pos]:
            total_count=ase_list[chromosome][pos][alt]["tot_count"]
            alt_count=ase_list[chromosome][pos][alt]["alt_count"]

            if total_count < args.m:
               continue

            ase_list[chromosome][pos][alt]["non_param"]=non_param_test(total_count,alt_count,args.n,all_variants)

if args.vcf.endswith('.gz'):
    reader = gzip.open
else:
    reader = open

for line in reader(args.vcf):
    line = line.decode('utf-8')
    if line[0] == "#":
        if not line[1] == "#":
            print ("##FORMAT=<ID=BT,Number=4,Type=Float,Description=\"BootstrapAnn p-values and GATK-ASEcounter stats (alt_count,total_count,binomial,nonparametric(-1 if not tested) )\">")
        print (line.strip())
        continue

    content=line.strip().split()
    if content[0] in ase_list:
        if content[1] in ase_list[content[0]]:
            if content[4] in ase_list[content[0]][content[1]]:
                content[8]+= ":BT"

                bin_p="{}".format(round(float(ase_list[content[0]][content[1]][content[4]]["p_bin"]),5))
                non_param_p="{}".format(round( float(ase_list[content[0]][content[1]][content[4]]["non_param"]) ,5))

                content[9]+=":{},{},{},{}".format(ase_list[content[0]][content[1]][content[4]]["alt_count"],ase_list[content[0]][content[1]][content[4]]["tot_count"],bin_p,non_param_p)
                print ("\t".join(content))

            else:
                print (line.strip())
        else:
            print (line.strip())
    else:
        print (line.strip())