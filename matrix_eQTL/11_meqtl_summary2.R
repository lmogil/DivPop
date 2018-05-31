"%&%" = function(a,b) paste(a,b,sep="")
library(data.table)
library(dplyr)
my.dir <- '/home/lauren/'
res.dir <- my.dir %&% 'files_for_revisions_plosgen/meqtl_results/MESA/'
date <- Sys.Date()

poplist <- c('AFA','HIS','CAU')
nklist <- c(0,10,20,30,50)
pclist <- c(0,3,5,10)

fdrarray <- array(0,dim = c(length(pclist),7,length(nklist),length(poplist)))

for(popnum in 1:length(poplist)){
  pop <- poplist[popnum]
  for(nknum in 1:length(nklist)){
    Nk <- nklist[nknum]
    for(pcnum in 1:length(pclist)){
    	pc <- pclist[pcnum]
    	fullres <- data.frame()
    	for(i in 1:22){
      		infile <- res.dir %&% pop %&%'_Nk_'%&% Nk %&%'_PFs_chr' %&% i %&% 'pcs_'%&% pc %&%'.meqtl.cis.2018-04-*.txt.gz'
      		res <- fread(sprintf("zcat %s", infile))
      		#res <- fread(sprintf("%s",infile))
      		fullres <- rbind(fullres, res)
    	}
    		fullres <- mutate(fullres, allFDR = p.adjust(pvalue,"BH"))
    		p <- dim(dplyr::filter(fullres,pvalue<5e-8))[1]
    		f1 <- dim(dplyr::filter(fullres,allFDR<0.1))[1]
    		f5 <- dim(dplyr::filter(fullres,allFDR<0.05))[1]
    		tot <- dim(fullres)[1]
    		fdrarray[pcnum,,nknum,popnum] <- c(pop,pc,Nk,p,f1,f5,tot)
    		png(filename = res.dir %&% pop %&% "_Nk_" %&% Nk %&% "pcs"%&% pc %&%"_PFs_ALL_chr1-22_" %&% date %&% ".png")
    		plot(hist(fullres$pvalue))
    		dev.off()
  }
}
}
fdrmat <- apply(fdrarray, 2, rbind)
colnames(fdrmat) <- c("pop","pc","Nk","P_5e-8","FDR_0.1","FDR_0.05","cisTested")
write.table(fdrmat, file=res.dir %&% "MESA_meqtl.cis_summary_" %&% date %&% ".txt",quote=F,row.names=F)
