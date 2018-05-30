######replace affyID with RSID in plink files
###Lauren S. Mogil

import gzip
from collections import defaultdict



ref_file_rs= "/home/lauren/MESA_dbGaP_55081/phg000071.v1.p1.NHLBI_SHARE_MESA.marker-info.MULTI/GenomeWideSNP_6.na24.annot.csv"
###"Probe Set ID","Affy SNP ID","dbSNP RS ID"
#make dictionary of affy id and rsid
snprsdict = {}
for line in open(ref_file_rs):
	if line.startswith('#') == False:
		arr = line.strip().split(',')
		affy = arr[0]
		affy = affy.replace('"', "")
		rs=arr[2]
		rs = rs.replace('"', "")
		snprsdict[affy]=rs


#open file to write out new bim file with rsids
out=open('/home/lauren/MESA_dbGaP_55081/all_mesa_merged/MESA_all_merged_rsid_new.bim','wb')

bim_file_rs= "/home/lauren/MESA_dbGaP_55081/all_mesa_merged/MESA_all_merged.bim"
###"Probe Set ID","Affy SNP ID","dbSNP RS ID"
for line in open(bim_file_rs):
	arr = line.strip().split('\t')
	chr = arr[0]
	aff = arr[1]
	loc=arr[2]
	bp=arr[3]
	a1=arr[4]
	a2=arr[5]
	#if the affy id is in the dictionary pull rsid
	if aff in snprsdict:
		rsid =snprsdict[aff]
		outf= chr+'\t'+rsid+'\t'+loc+'\t'+bp+'\t'+a1+'\t'+a2+'\n'
		out.write(outf)
		
		
out.close()


