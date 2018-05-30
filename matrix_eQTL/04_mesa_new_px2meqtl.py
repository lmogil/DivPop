#!/usr/bin/python
'''
converts PrediXcan dosages to MatrixEQTL dosages and snplocs
filters individuals to include those with expression data
dose allele is Allele2(B), see https://github.com/hakyimlab/PrediXcan/
creates 3 files
1) predixcan dosage file with individuals with expresision data
2)meqtl dosage file
3)meqtl location file
'''

#each chr takes 1-12 hrs depending on cluster use

import gzip
import sys

chr = sys.argv[1]
pop=sys.argv[2]


#chr='22'
#pop='AFA'

dir='/home/lauren/imputation_mesa_2/'+pop.lower()+'_imp/UMich_dosages/'
d_file='chr'+chr+'.maf0.01.r20.8noambig.dosage.txt.gz'



#expidfile = "/home/lauren/mesa_dosages/samples_gen_"+pop+".txt"
#expidlist = open(expidfile).read().split()


outdosage = gzip.open(dir+pop+"chr"+chr+"dosage_w_expression.txt.gz","wb")


# get a list of ids that have expression                                                                                                                   
exfile = "/home/lauren/mesa_dosages/samples_gen_"+pop+".txt"
exidlist = open(exfile).read().rstrip().split('\n')
exidlist1 = [sam.split(' ')[0] for sam in exidlist]
 #exclude FID  

# get list of predixcan dosage sample ids
samplefile = dir+'samples.txt'
samplelist = open(samplefile).read().rstrip().split('\n')
samplelist = [sam.split(' ')[1] for sam in samplelist]


# get sampleid header for output
outsamples = open(dir+'samples_'+pop+'.txt',"w")
ids = [sam for sam in exidlist1 if sam in samplelist]
#outsamples.write('\n'.join(ids) + '\n')


# get dosage file data
#SNP ID, RS #, base-pair position, allele coded A, and allele coded B
dosagefile=dir+d_file
outndosage = gzip.open(dir+pop+"chr"+chr+"dosage_w_expression.txt.gz","wb")
outdosage = gzip.open(dir+pop+"chr"+chr+"meqtl_w_expression.txt.gz","wb")
outsnplist = open(dir+pop+"_pop_dosage_chr" + chr + "_snpsloc.txt", "w")


head = 'rsid'+'\t' + '\t'.join(ids) + '\n'

outdosage.write(head)
outsnplist.write('rsid\tchr\tpos\n')

for line in gzip.open(dosagefile):
    arr = line.strip().split()
    (c, rsid, pos, a1, a2, maf) = arr[:6]
    dosagerow = arr[6:]
    #keep dosages with PC data, dosagerow and samplelist are same length/order  
    dosagerowid = [samplelist[i] for i in range(len(dosagerow)) if samplelist[i] in exidlist1]                                               
    dosagerow1 = [dosagerow[i] for i in range(len(dosagerow)) if samplelist[i] in ids]
    dosages = '\t'.join(map(str,dosagerow1))
    dosages2=' '.join(map(str,dosagerow1))
    outnewdosage=c+' '+rsid+' '+pos+' '+a1+' '+a2+' '+maf+' '+dosages2+'\n'
    outndosage.write(outnewdosage)
    output = rsid + '\t' + dosages + '\n'
    outdosage.write(output)
    snpout = rsid + '\t' + c + '\t' + pos + '\n'
    outsnplist.write(snpout)

outsamples.write('\n'.join(dosagerowid) + '\n') 


outdosage.close()
outndosage.close()
outsnplist.close()
outsamples.close()
