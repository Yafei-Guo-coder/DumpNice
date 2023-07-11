###Berg and Coop method: new with full data: Actually only Europe so, check the old script, optimize and re-run the necessary parts
###part1: create the input files
###Part2: run the method

####input files
#1 gwasdatafile (associated SNPs with the traits and the effect size and frequency)
#2 freqsfile (frequency of the associated SNPs in each population)
#3 EnvVar files (several files with environmental traits)
#4 Fulldataset (full data sets of all SNPs)
#5 matchpop (match data to one pop)

###############
#1 gwasdatafile
###############
#I need the phenotype data for this
load("/media/bene/Data/Phenotypic_data_output/HLLL_exp_2016/New_Analysis201718/corrected_genotype_means_GWAS_HL_274acc.R")#corrected_data_hl
head(corrected_data_hl)
load("/media/bene/Data/Phenotypic_data_output/HLLL_exp_2016/New_Analysis201718/corrected_genotype_means_GWAS_LL_277acc.R")
head(corrected_data_ll)

#load genomic data:
load("/media/bene/Genomic_data/Genotype_data/SNP_data/233_2061001_27china/snpmat_full_beagle_corrected_lowNA_highMAF.R")
dim(snpmat_beagle_corrected_filtered_lowNA_highMAF)


#create a gwasdatafile for each GWAS trait:
#Trait 1: Final Size in HL
load(file="/media/bene/Data/Output_files/Gwas/GWAS_EurandChi/Final_size_HL_205acc_Europe.R")
head(gwas_FSHL)
gwasdatafile_FSHL <- gwas_FSHL
#filter by -log10(P) > 4
gwasdatafile_FSHL_lowP <- gwasdatafile_FSHL[which(-log10(gwasdatafile_FSHL$Pval)>4),]
dim(gwasdatafile_FSHL_lowP)#69 SNPs
#are there any SNPs with too low MAF?
nrow(gwasdatafile_FSHL_lowP[which(gwasdatafile_FSHL_lowP$MAF>0.05),])#2 SNPs with too low MAF
gwasdatafile_FSHL_lowP_MAFfiltered <-  gwasdatafile_FSHL_lowP[which(gwasdatafile_FSHL_lowP$MAF>0.05),]#67 SNPs

#####find the effect allele
pheno_FSHL <- corrected_data_hl[,c("name.accession","ID.genotype","population","FinalSize")]
pheno_FSHL$ID.genotype <- as.numeric(pheno_FSHL$ID.genotype)
pheno_FSHL <- rbind(pheno_FSHL[which(pheno_FSHL$population=="Northern_Europe"),],pheno_FSHL[which(pheno_FSHL$population=="Spain"),])
#205 accession without france
pheno_FSHL <- pheno_FSHL[order(pheno_FSHL$ID),]
rownames(pheno_FSHL) <- pheno_FSHL$ID
head(pheno_FSHL)
ids <- as.character(pheno_FSHL$ID.genotype)
#only 206 ids are shared
ids <- ids[which(ids %in% rownames(snpmat_beagle_corrected_filtered_lowNA_highMAF))]
allele_state_FSHL_lowP <- snpmat_beagle_corrected_filtered_lowNA_highMAF[ids,gwasdatafile_FSHL_lowP_MAFfiltered$SNP]
#bind allele state and phenotype
FSHL_lowP_effallele <- merge(x=pheno_FSHL,y=allele_state_FSHL_lowP, by = 0)

#find the effect allele and output the FRQ for every allele
phenomeans_0 <- NULL
phenomeans_1 <- NULL
#allele state
A1 <- NULL
A2 <- NULL
eff_alleles <- NULL
noeff_alleles <- NULL

for (i in 1:ncol(allele_state_FSHL_lowP)){
  allele <- allele_state_FSHL_lowP[,i]
  pheno_mean_0 <- mean(pheno_FSHL[which(allele==0),4])#mean of all alleles of state 0 at one locus
  pheno_mean_1 <- mean(pheno_FSHL[which(allele==1),4])
  #save the means into a matrix
  phenomeans_0 <- cbind(phenomeans_0,pheno_mean_0)
  phenomeans_1 <- cbind(phenomeans_1,pheno_mean_1)
  #find the effect alleles
  if (pheno_mean_1 < pheno_mean_0){ ##smaller value means a lower t50 -> faster (both ways can be defined)
    eff_allele <- "AC_1"
    noeff_allele <- "AC_0"
  } else {
    eff_allele <- "AC_0"
    noeff_allele <- "AC_1"
  }
  #save all the effect/non-effect alleles in 2 vectors
  eff_alleles <- cbind(eff_alleles,eff_allele)
  noeff_alleles <- cbind(noeff_alleles,noeff_allele)
  
  A1_ <- gwasdatafile_FSHL_lowP_MAFfiltered[i,eff_allele]
  A2_ <- gwasdatafile_FSHL_lowP_MAFfiltered[i,noeff_allele]
  #bind all A1/A2 together
  A1 <- c(A1,A1_)
  A2 <- c(A2,A2_)
}

#convert into data frame with SNP as rownames
phenomeans_0 <- as.data.frame(phenomeans_0,col.names = colnames(allele_state_FSHL_lowP))
colnames(phenomeans_0) <- colnames(allele_state_FSHL_lowP)
str(phenomeans_0)

phenomeans_1 <- as.data.frame(phenomeans_1,col.names = colnames(allele_state_FSHL_lowP))
colnames(phenomeans_1) <- colnames(allele_state_FSHL_lowP)

FRQ <- (A1/205) ###FRQ is the frequency of the effect allele divided by the total allele number (205)
FRQ

#the effalleles are the "names" of the effect alleles (0/1) and can be input as a vector, A1 & A2 that were used before are the respective allele counts and are necessary for the FRQ calculation

SNP <- paste(gwasdatafile_FSHL_lowP_MAFfiltered$Chr,gwasdatafile_FSHL_lowP_MAFfiltered$Pos,sep="-")

gwasdatafileFSHL <- cbind.data.frame(SNP,A1=as.vector(eff_alleles),A2=as.vector(noeff_alleles),EFF=gwasdatafile_FSHL_lowP_MAFfiltered$beta,FRQ)
colnames(gwasdatafileFSHL) <- c("SNP","A1","A2","EFF","FRQ")
head(gwasdatafileFSHL)
#write.table(gwasdatafileFSHL, file="/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/FSHL/gwasdatafileFSHL.txt",quote=F,row.names=F)
gwasdatafileFSHL <- read.table(file="/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/FSHL/gwasdatafileFSHL.txt", header=T)



#################################
######2) freqsfile
#################################
#the freqsfile is created using the CreateTraitFile function
source ("/home/bene/Desktop/Bene/Scripts/Publications/Berg&Coop_2014/PolygenicAdaptationCode-master/Scripts/CreateTraitFile.R")	

#the function requires a gwasdatafile (created in 1)) and a gen.data file

#####WRITE A FUNCTION FOR gen.data.file generation
#the gen.data file requires the columns SNP, CLST(all clusters for every SNP) and FRQ with the possibility for further columns for more tests (A1,A2,CHR,POS,IMP,BVAL)
#the clusters are the IDs not regions (as I did before)
#I need all the necessary fields for each ID

#extract the SNPs
SNP <- as.character(gwasdatafileFSHL$SNP)
#in order to call the SNPs from the SNP-matrix there names need to include a space after the chromosome: 1- 92 instead of 1-92:
library(tidyr)
SNP_ <- separate(data.frame(SNP),"SNP",into = c("C","B"), sep="-")
SNP_ <- paste(SNP_$C,SNP_$B,sep="- ")

SNP_freq <- snpmat_beagle_corrected_filtered_lowNA_highMAF[ids,SNP_]
####find which IDs are in which pop, overall create a seperate freqsfile for each pop, merge them in the end and sort by SNP
NE_IDs <- corrected_data_hl[which(corrected_data_hl$population=="Northern_Europe"),"ID.genotype"]
SP_IDs <- corrected_data_hl[which(corrected_data_hl$population=="Spain"),"ID.genotype"]

###??????
#IDs 9580, 9829 & 9931 (already excluded) are not in the X file -> exclude
SP_IDs <- SP_IDs[-c(62,10)]

#for each allele and each ID I need the frequency of 0 and 1
# the gendata file contains SNP, CLST, A1,A2,FRQ

#CLST is just each ID repeated for each SNP
CLST <- rep(c(NE_IDs,SP_IDs),times=length(SNP))
#SNP, needs to be the SNP name repeated for every CLST(=ID)
SNP_in <- rep(SNP,each=length(unique(CLST)))
gendata_FSHL <- data.frame(cbind(SNP=SNP_in,CLST))
head(gendata_FSHL)
#the tricky part is now to get A1, A2 and FRQ; A1 is the name of the effect increasing and A2 of the decreasing allele; I only have allele counts
#and phenotypic data, so I give A1 the name AC_1 is allele state 1 is the increasing one, else it is AC_0 (and vice verse for A2)
#FRQ is the Frequency of A1 over A2 (either 0 or 1 if I take)

#actually A1 and A2, do not change, it is just the allele name, so here I only need to rep the A1/A2 from the GWAS file for each SNP

#so first I just need the allele counts of 1 and 0 for each allele and each ID I am including: essentially it is SNP_freq, but in a different order and as a vector
#I use the eff_alleles from the gwasfile calc
length(rep(gwasdatafileFSHL$A1,each=length(unique(CLST))))#correct
A1 <- rep(gwasdatafileFSHL$A1,each=length(unique(CLST)))
A2 <- rep(gwasdatafileFSHL$A2,each=length(unique(CLST)))

gendata_FSHL <- cbind(gendata_FSHL,A1,A2)
head(gendata_FSHL)
#for FRQ I will loop over each line
#if the allele of the Ind at the position == A1 --> FRQ=1, else FRQ=0
FRQ <- rep(NA,nrow(gendata_FSHL))
for (i in 1:nrow(gendata_FSHL)){
  snp <- data.frame(snp=as.character(gendata_FSHL$SNP[i]))
  snp <- separate(snp,snp,into=c("C","B"))
  snp <- paste(snp$C,snp$B,sep="- ")
  allele <- snpmat_beagle_corrected_filtered_lowNA_highMAF[as.character(gendata_FSHL$CLST[i]),snp]
  if (gendata_FSHL$A1[i] == "AC_0"){
    if (allele==0){
      FRQ[i] <- 1
    } else{
      FRQ[i] <- 0
    }
  } else {
    if (allele==1){
      FRQ[i] <- 1
    } else {
      FRQ[i] <- 0
    }
  }
}

head(FRQ,n=40)
summary(FRQ)
gendata_FSHL <- cbind(gendata_FSHL,FRQ)
head(gendata_FSHL,n=40)#looks good
#write.table(gendata_FSHL, file="/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/FSHL/gendatafileFSHL.txt", quote = F,row.names = F)
gendata_FSHL <- read.table(file="/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/FSHL/gendatafileFSHL.txt", header=T)

setwd("/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/FSHL/")
source ("/home/bene/Desktop/Bene/Scripts/Publications/Berg&Coop_2014/PolygenicAdaptationCode-master/Scripts/CreateTraitFile.R")
CreateTraitFile("gwasdatafileFSHL.txt", "gendatafileFSHL.txt")
#I need to save the files without quotes and rownames + change the SNP identifiers (for a few SNPs I just erased the space between CHR and BP)

#################################
######3 EnvVar files (several files with environmental traits)
#################################
#traits: bioclim variables, latitude,longitude, Hannes seasonal variance data
#Hannes climate data
climate <- read.csv("/home/bene/Desktop/Bene/Phenotype_data/Publications/Dittberner2018/253021-4.csv")
head(climate)
#Bioclim variables:
load("/media/bene/Data/Climate_data/climate1001accessions_19bioclim.Rdata")
head(climate1001)
load("/media/bene/Data/Climate_data/climateall_myaccessions_growingseasondata_Hannes.Rdata")
head(climate_my_accessions)

#create a file with all my accessions
load("/media/bene/Data/Phenotypic_data_output/HLLL_exp_2016/New_Analysis201718/corrected_genotype_means_GWAS_HL_274acc.R")#corrected_data_hl
head(corrected_data_hl)
load("/media/bene/Data/Phenotypic_data_output/HLLL_exp_2016/New_Analysis201718/corrected_genotype_means_GWAS_LL_277acc.R")
head(corrected_data_ll)

#I need all Ids I will have in the dataset




#load full info file
completedata <- read.table("/home/bene/Desktop/Bene/Phenotype_data/Light_Experiment_Traits/CompleteDataHLLL_HighT.txt", header=T)
head(completedata)
completedata$ID.genotype <- as.numeric(as.character(completedata$ID.genotype))
completedata<- completedata[order(completedata$ID.genotype),]
rownames(completedata) <- 1:nrow(completedata)

#only retain accessions from the 1001 genomes project
only1001 <- completedata[1:1198,]
head(only1001)
#extract env parameter
latitude <- unique(only1001[,c("ID.genotype","Latitude","population")])
longitude <- unique(only1001[,c("ID.genotype","Longitude","population")])
rownames(latitude) <- 1:nrow(latitude)
rownames(longitude) <- 1:nrow(longitude)
#change colnames to the right parameters and change population to region identifier
colnames(latitude) <- c("CLST","ENV","REG")
colnames(longitude) <- c("CLST","ENV","REG")
#change the reg identifier: France=1, Spain=2, Sweden=3
latitude$REG <- as.character(latitude$REG)
latitude[which(latitude$REG=="Western_Europe"),"REG"]  <- 1
latitude[which(latitude$REG=="Spain"),"REG"]  <- 2
latitude[which(latitude$REG=="Northern_Europe"),"REG"]  <- 3

longitude$REG <- as.character(longitude$REG)
longitude[which(longitude$REG=="Western_Europe"),"REG"]  <- 1
longitude[which(longitude$REG=="Spain"),"REG"]  <- 2
longitude[which(longitude$REG=="Northern_Europe"),"REG"]  <- 3

head(latitude)
head(longitude)
#####the genotypes 9580,9829,9931 are not in the SNP file --> need to be excluded from the envvar files as well
latitude <- latitude[-c(127,153,204),]
longitude <- longitude[-c(127,153,204),]

#write.table(latitude,"../EnvVar/latitude_206.txt", row.names = F, quote = F)
#write.table(longitude,"../EnvVar/longitude_206.txt", row.names = F, quote = F)









########################################################################################################################
#######RUN THE METHOD
########################################################################################################################
#####test the coop function with this
source ("~/Desktop/Bene/Scripts/Publications/Berg&Coop_2014/PolygenicAdaptationCode-master/Scripts/functions.R")	 

setwd("/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/FSHL/")
###i changed the location, the fulldata is now in 180226
polygen_result <- PolygenicAdaptationFunction ( 
  gwas.data.file ="gwasdatafileFSHL.txt" , 
  freqs.file = "gwasdatafileFSHL.gendatafileFSHL.txt.freqs" , 
  env.var.data.files = list ( "/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/EnvVar/latitude_203.txt" 
  ) , # Note: you can supply as many env.var.data.files concurrently as you want. If only supplying one file it should still be included in a list, e.g. env.var.data.files = list ( "Example/EnvVar/HGDP_LATS_GLOBAL")
  match.pop.file = "/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/Input_t50HL_206/180226/FULLDATA_MAF/matchpop_997.txt" , 
  full.dataset.file = "/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/Input_t50HL_206/180226/FULLDATA_MAF/fulldata_t50HL_merged.txt" , 
  path = "NewResult" , 
  match.categories = c ( "FRQ") ,
  match.bins = list ( c(0,1) ) , 
  cov.SNPs.per.cycle = 5000 , 
  cov.cycles = 4 , 
  null.phenos.per.cycle = 1000 , 
  null.cycles = 2 ,
  load.cov.mat = F ,
  sim.null = T ,
  check.allele.orientation = F
)


#Frequency dataset has data for populations not present in environmental variable datasets.
#load EnvVar and frequency data set
lat <- read.table("/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/EnvVar/latitude_206.txt",header=T)
head(lat)
freq <- read.table("gwasdatafileFSHL.gendatafileFSHL.txt.freqs",header=T)
head(freq)#this one does not have the accession names anyway
gwas <- read.table("gwasdatafileFSHL.txt",header=T)
head(gwas)
#the gendata file I created (and thereby the freqsfile) are not correct, CLST needs to be the ID so they can match to EnvVar

#corrected, run again:
#Environmental variable dataset has data for populations not present frequency dataset.
#I will change the present envvar latitude file for now and create them properly later, now I just want them to have the same pops
latitude <- read.table("/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/EnvVar/latitude_209.txt",header=T)
head(latitude)
nrow(latitude[which(latitude$CLST %in% gendata_FSHL$CLST),])
#only 203 instead of 209--> new latitude
latitude <- latitude[which(latitude$CLST %in% gendata_FSHL$CLST),]
write.table(latitude, file="/media/bene/Data/Output_files/Gwas/GWAS_withinregion_and_Coop/BergandCoop/NEW_2018/EnvVar/latitude_203.txt",row.names = F,quote = F)

#additionally I have to be careful, that the columns of the files are in the correct order!