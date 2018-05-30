import gzip
from collections import defaultdict

pop=sys.argv[1]

pop="CAU"
#afa samples file path
afa_ids="/home/lauren/mesa_dosages/samples_gen_"+pop+".txt"

#initiate list of afa sample ids
afa=[]
for i in open(afa_ids): #open sample file
	afai = i.strip().split() #split by white space
	afaid=afai[0]#assign first 5 indices
	afa.append(afaid) #put sample ids in list called afa

#MESA fam file that contains all of the populatiuons
afafam="/home/lauren/MESA_dbGaP_55081/all_mesa_merged/"+pop.lower()+"_samples.txt"

#open new file that will have FID/IID columns
ofile1=open("/home/lauren/MESA_dbGaP_55081/all_mesa_merged/"+pop.lower()+"_wexp_samples.txt","wb")

#read and open fam file
for a in open(afafam):
	afaf = a.strip().split() #split by white space
	fid=afaf[0] #first column is the FID
	iid=afaf[1] #second column is the IID
	if iid in afa: #if the IID is in the sample list from previous step
		out= fid+' '+iid+'\n' #make column 1 FID and column 2 IID and then skip to a new line
		ofile1.write(out) #write out this line to file created

ofile1.close() #close file
