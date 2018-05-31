###Lauren S. Mogil###
###make MESA data in the correct format for predictDB###

import gzip
import sys

pop=sys.argv[1]
chr=sys.argv[2]

#pop='CAU'
#chr='22'

######change file paths to fit your data 
dosage = "/home/lauren/files_for_revisions_plosgen/px_dosages_expression_ppl/"+pop+"chr"+chr+"dosage_w_expression.txt.gz"
samples="/home/lauren/files_for_revisions_plosgen/px_dosages_expression_ppl/samples_"+pop+".txt"


anot_file=open("/home/lauren/files_for_revisions_plosgen/en_v7/prepare_data/genotypes/"+pop+"_"+chr+"_annot.txt","wb")
snp_file=open("/home/lauren/files_for_revisions_plosgen/en_v7/prepare_data/genotypes/"+pop+"_"+chr+"_snp.txt","wb")


########
anot_file.write('chr'+'\t'+'pos'+'\t'+'varID'+'\t'+'refAllele'+'\t'+'effectAllele'+'\t'+'rsid'+'\n')

samp=[]
for sam in open(samples):
	s=sam.strip().split()
	s1=s[0]
	samp.append(s1)
	
	
samp1 = str(samp) #convert to string
splitheader = samp1.replace("'", "") #format and remove extra characters
splitheader = splitheader.replace("]", "")
splitheader = splitheader.replace("[", "")
splitheader = splitheader.replace(" ", "")
splitheader = splitheader.replace(",", "\t") #make tab delim file 
head = "id\t" + splitheader #add label for ids
snp_file.write(head + "\n")
	
for line in gzip.open(dosage):
	arr=line.strip().split()
	(c, rs, pos, a1, a2, maf)=arr[0:6]
	dose=arr[6:]
	varid= chr+'_'+pos+'_'+a1+'_'+a2+'_'+'b37'
	snpid='snp'+'_'+chr+'_'+pos
	dosages = '\t'.join(map(str,dose))
	dosages=dosages.replace(",","")
	snp= varid+'\t'+dosages+'\n'
	annot= chr+'\t'+pos+'\t'+varid+'\t'+a1+'\t'+a2+'\t'+rs+'\n'
	anot_file.write(annot)
	snp_file.write(snp)
	
	
	
anot_file.close()
snp_file.close()



