########new GWAS QC for mesa dbgap data (not imputed)
$DIR = /home/lauren/MESA_dbGaP_55081/phg000071.v2.NHLBI_SHARE_MESA.genotype-calls-matrixfmt.c1/plink_qc_files/
cd $DIR
#1. merge two mesa data sets 

#2. change affy ID to rsid...new bim file


#3. use newly formated bim file rsids...need to liftover to hg19 from hg18 (later)
#make sure you are in the correct directory and/or have the correct path to your files while using plink flags

plink --bed MESA_all_merged.bed --bim MESA_all_merged_rsid_new.bim  --fam MESA_all_merged.fam --make-bed  --out /home/lauren/MESA_dbGaP_55081/all_mesa_merged/MESA_all_merged_1

#515881 MB RAM detected; reserving 257940 MB for main workspace.
#909622 variants loaded from .bim file.
#7377 people (3440 males, 3937 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 5450 founders and 1927 nonfounders present.
#Calculating allele frequencies... done.
#Warning: 447174 het. haploid genotypes present (see
#/home/lauren/MESA_dbGaP_55081/phg000071.v2.NHLBI_SHARE_MESA.genotype-calls-matrixfmt.c1/plink_qc_files/MESA_c1_start.hh
#); many commands treat these as missing.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
#treat these as missing.
#Total genotyping rate is 0.996968.
#909622 variants and 7377 people pass filters and QC.
#Note: No phenotypes present.

#4. Check sex and calculate call rates for flagging poorly called SNPs and individuals

plink --bfile MESA_all_merged_1 --check-sex --missing --out MESA_all_merged_1

#515881 MB RAM detected; reserving 257940 MB for main workspace.
#909622 variants loaded from .bim file.
#7377 people (3440 males, 3937 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 5450 founders and 1927 nonfounders present.
#Calculating allele frequencies... done.
#Warning: 447174 het. haploid genotypes present (see MESA_c1_start.sexmiss.hh );
#many commands treat these as missing.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
#treat these as missing.
#Total genotyping rate is 0.996968.
#--missing: Sample missing data report written to MESA_c1_start.sexmiss.imiss,
#and variant-based missing data report written to MESA_c1_start.sexmiss.lmiss.
#909622 variants and 7377 people pass filters and QC.
#Note: No phenotypes present.
#--check-sex: 36538 Xchr and 0 Ychr variant(s) scanned, 689 problems detected.
#Report written to MESA_c1_start.sexmiss.sexcheck .

#5. use autosome flag to remove all non-autosomal chromosomes

plink --bfile MESA_all_merged_1 --autosome --make-bed --out MESA_all_merged_autosome


#6. Recalculate individual call rates after removing SNPs with call rates <99%
plink --bfile MESA_all_merged_autosome --geno 0.01 --make-bed --out MESA_all_merged_autosome.geno0.01

#515881 MB RAM detected; reserving 257940 MB for main workspace.
#909622 variants loaded from .bim file.
#7377 people (3440 males, 3937 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 5450 founders and 1927 nonfounders present.
#Calculating allele frequencies... done.
#Warning: 447174 het. haploid genotypes present (see MESA_c1_start.geno0.01.hh
#); many commands treat these as missing.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
#treat these as missing.
#Total genotyping rate is 0.996968.
#53838 variants removed due to missing genotype data (--geno).
#855784 variants and 7377 people pass filters and QC.
#Note: No phenotypes present.

plink --bfile MESA_all_merged_autosome.geno0.01 --missing --out MESA_all_merged_autosome.geno0.01

#515881 MB RAM detected; reserving 257940 MB for main workspace.
#855784 variants loaded from .bim file.
#7377 people (3440 males, 3937 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 5450 founders and 1927 nonfounders present.
#Calculating allele frequencies... done.
#Warning: 404203 het. haploid genotypes present (see MESA_c1_start.geno0.01.hh
#); many commands treat these as missing.
#Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
#treat these as missing.
#Total genotyping rate is 0.998716.
#--missing: Sample missing data report written to MESA_c1_start.geno0.01.imiss,
#and variant-based missing data report written to MESA_c1_start.geno0.01.lmiss.



#5. LD prune (rm 1 SNP if r2>0.3 in 50 SNP window) for relationship check and heterozygosity calculation
plink --bfile MESA_all_merged_autosome.geno0.01 --indep-pairwise 50 5 0.3 --out MESA_all_merged_autosome.geno0.01.LD3

#6. Relationship check
plink --bfile MESA_all_merged_autosome.geno0.01 --extract MESA_all_merged_autosome.geno0.01.LD3.prune.in --genome --min 0.05 --out MESA_all_merged_autosome.geno0.01.LD3.0.05

#7. Separate pops
**make new *pop*_samples file from python script that pulls FID/IID pairs from fam file in /home/lauren/scripts/make_samples_keep.py

### location of sample files here /home/lauren/MESA_dbGaP_55081/all_mesa_merged
# example: plink --bfile MESA_all_merged_autosome.geno0.01 --keep /home/lauren/MESA_dbGaP_55081/all_mesa_merged/afa_samples.txt --make-bed --out MESA_AFA_merged_autosome.geno0.01


plink --bfile MESA_all_merged_autosome.geno0.01 --keep afa_samples.txt --make-bed --out MESA_AFA_merged_autosome.geno0.01

plink --bfile MESA_all_merged_autosome.geno0.01 --keep his_samples.txt --make-bed --out MESA_HIS_merged_autosome.geno0.01

plink --bfile MESA_all_merged_autosome.geno0.01 --keep cau_samples.txt --make-bed --out MESA_CAU_merged_autosome.geno0.01


#8. LD prune (rm 1 SNP if r2>0.3 in 50 SNP window) for relationship check and heterozygosity calculation
plink --bfile MESA_AFA_merged_autosome.geno0.01 --indep-pairwise 50 5 0.3 --out MESA_AFA_merged_autosome.geno0.01.LD3

plink --bfile MESA_HIS_merged_autosome.geno0.01 --indep-pairwise 50 5 0.3 --out MESA_HIS_merged_autosome.geno0.01.LD3

plink --bfile MESA_CAU_merged_autosome.geno0.01 --indep-pairwise 50 5 0.3 --out MESA_CAU_merged_autosome.geno0.01.LD3


#9.Relationship check for individual pops
plink --bfile MESA_AFA_merged_autosome.geno0.01 --extract MESA_AFA_merged_autosome.geno0.01.LD3.prune.in --genome --min 0.05 --out MESA_AFA_merged_autosome.geno0.01.LD3.0.125

plink --bfile MESA_HIS_merged_autosome.geno0.01 --extract MESA_HIS_merged_autosome.geno0.01.LD3.prune.in --genome --min 0.05 --out MESA_HIS_merged_autosome.geno0.01.LD3.0.05

plink --bfile MESA_CAU_merged_autosome.geno0.01 --extract MESA_CAU_merged_autosome.geno0.01.LD3.prune.in --genome --min 0.05 --out MESA_CAU_merged_autosome.geno0.01.LD3.0.05

#####repeat step 6.1 with individual pops

#10.Check heterozygosity (across all autosomal SNPs) -- look at that distribution across individuals to check for and rm outliers (F: mean +/-3 sd)
plink --bfile MESA_AFA_merged_autosome.geno0.01 --het --out MESA_AFA_merged_autosome.geno0.01

plink --bfile MESA_HIS_merged_autosome.geno0.01 --het --out MESA_HIS_merged_autosome.geno0.01

plink --bfile MESA_CAU_merged_autosome.geno0.01 --het --out MESA_CAU_merged_autosome.geno0.01

#11. liftover steps---change from genome build hg18 to hg19
#make ped/map files of pop plink files
plink --bfile MESA_AFA_merged_autosome.geno0.01 --recode --out MESA_AFA_4lift_all

plink --bfile MESA_HIS_merged_autosome.geno0.01 --recode --out MESA_HIS_4lift_all

plink --bfile MESA_CAU_merged_autosome.geno0.01 --recode --out MESA_CAU_4lift_all

#run liftover script with new ped/mapfiles --example with AFA
python /home/lauren/MESA_dbGaP_55081/all_mesa_merged/liftover/LiftMap.py -m MESA_AFA_4lift_all.map -p MESA_AFA_4lift_all.ped -o mesa_afa_new_all

plink --file MESA_AFA_4lift_all --exclude mesa_afa_new_all.bed.unlifted --recode --out AFA_tolift_all

python /home/lauren/MESA_dbGaP_55081/all_mesa_merged/liftover/LiftMap.py -m AFA_tolift_all.map -p AFA_tolift_all.ped -o mesa_afa_use_all

plink --file mesa_afa_use_all --make-bed --out MESA_AFA_lifted_final_all

#HIS
python /home/lauren/MESA_dbGaP_55081/all_mesa_merged/liftover/LiftMap.py -m MESA_HIS_4lift_all.map -p MESA_HIS_4lift_all.ped -o mesa_his_new_all

plink --file MESA_HIS_4lift_all --exclude mesa_his_new_all.bed.unlifted --recode --out HIS_tolift_all

python /home/lauren/MESA_dbGaP_55081/all_mesa_merged/liftover/LiftMap.py -m HIS_tolift_all.map -p HIS_tolift_all.ped -o mesa_his_use_all

plink --file mesa_his_use_all --make-bed --out MESA_HIS_lifted_final_all

#CAU
python /home/lauren/MESA_dbGaP_55081/all_mesa_merged/liftover/LiftMap.py -m MESA_CAU_4lift_all.map -p MESA_CAU_4lift_all.ped -o mesa_cau_new_all

plink --file MESA_CAU_4lift_all --exclude mesa_cau_new_all.bed.unlifted --recode --out CAU_tolift_all

python /home/lauren/MESA_dbGaP_55081/all_mesa_merged/liftover/LiftMap.py -m CAU_tolift_all.map -p CAU_tolift_all.ped -o mesa_cau_use_all

plink --file mesa_cau_use_all --make-bed --out MESA_CAU_lifted_final_all
###########
############
###########
###########
#13. filter by genotyping rate and maf (unmerged)
plink --bfile MESA_AFA_lifted_final_all --geno 0.05 --maf 0.05 --make-bed   --out MESA_AFA_lifted_filt_all

plink --bfile MESA_HIS_lifted_final_all --geno 0.05 --maf 0.05 --make-bed  --out MESA_HIS_lifted_filt_all

plink --bfile MESA_CAU_lifted_final_all --geno 0.05 --maf 0.05 --make-bed  --out MESA_CAU_lifted_filt_all


#14. LD prune MESA files with removed individuals and make ped/map files (unmerged)
plink --bfile MESA_AFA_lifted_filt_all --indep-pairwise 50 5 0.3 --recode  --out MESA_AFA_autosome.geno0.01.sd.pruned

plink --bfile  MESA_HIS_lifted_filt_all --indep-pairwise 50 5 0.3 --recode --out MESA_HIS_autosome.geno0.01.sd.pruned

plink --bfile MESA_CAU_lifted_filt_all --indep-pairwise 50 5 0.3 --recode --out MESA_CAU_autosome.geno0.01.sd.pruned



#20.calculate allele freq. 
plink --bfile  MESA_AFA_lifted_final_all --freq --out  MESA_AFA_lifted_final_freq_all

plink --bfile  MESA_HIS_lifted_final_all --freq --out  MESA_HIS_lifted_final_freq_all

plink --bfile  MESA_CAU_lifted_final_all --freq --out  MESA_CAU_lifted_final_freq_all

#21. run perl file that checks samples against 1000g
perl HRC-1000G-check-bim.pl -b MESA_AFA_lifted_final_all.bim -f MESA_AFA_lifted_final_freq_all.frq -r 1000GP_Phase3_combined.legend -g

perl HRC-1000G-check-bim.pl -b MESA_HIS_lifted_final_all.bim -f MESA_HIS_lifted_final_freq_all.frq -r 1000GP_Phase3_combined.legend -g

perl HRC-1000G-check-bim.pl -b MESA_CAU_lifted_final_all.bim -f MESA_CAU_lifted_final_freq_all.frq -r 1000GP_Phase3_combined.legend -g

#22. use generated files from #21 to format for imputation
plink --bfile MESA_AFA_lifted_final_all --exclude Exclude-MESA_AFA_lifted_final_all-1000G.txt --make-bed --out TEMP1
plink --bfile TEMP1 --update-map Chromosome-MESA_AFA_lifted_final_all-1000G.txt --update-chr --make-bed --out TEMP2
plink --bfile TEMP2 --update-map Position-MESA_AFA_lifted_final_all-1000G.txt --make-bed --out TEMP3
plink --bfile TEMP3 --flip Strand-Flip-MESA_AFA_lifted_final_all-1000G.txt --make-bed --out TEMP4
plink --bfile TEMP4 --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --make-bed --out MESA_AFA_lifted_finalALL-updated

plink --bfile MESA_HIS_lifted_final_all --exclude Exclude-MESA_HIS_lifted_final_all-1000G.txt --make-bed --out TEMP1
plink --bfile TEMP1 --update-map Chromosome-MESA_HIS_lifted_final_all-1000G.txt --update-chr --make-bed --out TEMP2
plink --bfile TEMP2 --update-map Position-MESA_HIS_lifted_final_all-1000G.txt --make-bed --out TEMP3
plink --bfile TEMP3 --flip Strand-Flip-MESA_HIS_lifted_final_all-1000G.txt --make-bed --out TEMP4
plink --bfile TEMP4 --reference-allele Force-Allele1-MESA_HIS_lifted_final_all-1000G.txt --make-bed --out MESA_HIS_lifted_finalALL-updated

plink --bfile MESA_CAU_lifted_final_all --exclude Exclude-MESA_CAU_lifted_final_all-1000G.txt --make-bed --out TEMP1
plink --bfile TEMP1 --update-map Chromosome-MESA_CAU_lifted_final_all-1000G.txt --update-chr --make-bed --out TEMP2
plink --bfile TEMP2 --update-map Position-MESA_CAU_lifted_final_all-1000G.txt --make-bed --out TEMP3
plink --bfile TEMP3 --flip Strand-Flip-MESA_CAU_lifted_final_all-1000G.txt --make-bed --out TEMP4
plink --bfile TEMP4 --reference-allele Force-Allele1-MESA_CAU_lifted_final_all-1000G.txt --make-bed --out MESA_CAU_lifted_finalALL-updated


##make a Run_plink_makevcf.sh file that contains these lines ans run
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 1 --out MESA_AFA_lifted_final_all_imputationchr1
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 2 --out MESA_AFA_lifted_final_all_imputationchr2
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 3 --out MESA_AFA_lifted_final_all_imputationchr3
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 4 --out MESA_AFA_lifted_final_all_imputationchr4
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 5 --out MESA_AFA_lifted_final_all_imputationchr5
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 6 --out MESA_AFA_lifted_final_all_imputationchr6
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 7 --out MESA_AFA_lifted_final_all_imputationchr7
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 8 --out MESA_AFA_lifted_final_all_imputationchr8
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 9 --out MESA_AFA_lifted_final_all_imputationchr9
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 10 --out MESA_AFA_lifted_final_all_imputationchr10
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 11 --out MESA_AFA_lifted_final_all_imputationchr11
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 12 --out MESA_AFA_lifted_final_all_imputationchr12
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 13 --out MESA_AFA_lifted_final_all_imputationchr13
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 14 --out MESA_AFA_lifted_final_all_imputationchr14
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 15 --out MESA_AFA_lifted_final_all_imputationchr15
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 16 --out MESA_AFA_lifted_final_all_imputationchr16
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 17 --out MESA_AFA_lifted_final_all_imputationchr17
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 18 --out MESA_AFA_lifted_final_all_imputationchr18
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 19 --out MESA_AFA_lifted_final_all_imputationchr19
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 20 --out MESA_AFA_lifted_final_all_imputationchr20
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 21 --out MESA_AFA_lifted_final_all_imputationchr21
plink --bfile MESA_AFA_lifted_finalALL-updated --reference-allele Force-Allele1-MESA_AFA_lifted_final_all-1000G.txt --recode vcf --chr 22 --out MESA_AFA_lifted_final_all_imputationchr22


#23. sort and out put for each chromosome
for i in {1..22}; do vcf-sort MESA_AFA_lifted_final_all_imputationchr${i}.vcf | bgzip -c > MESA_AFA_lifted_final_all_imputation_chr${i}.vcf.gz; done

#24 Keep people with expression data 
plink --bfile MESA_AFA_lifted_final_all --keep afa_wexp_samples.txt --make-bed --out MESA_AFA_w_expression
plink --bfile MESA_HIS_lifted_final_all --keep all_mesa_merged/his_wexp_samples.txt --make-bed --out MESA_HIS_w_expression
plink --bfile MESA_CAU_lifted_final_all --keep cau_wexp_samples.txt --make-bed --out MESA_CAU_w_expression
