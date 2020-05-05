#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# parse input arguments
if ((length(args)==0) || (length(args)==1)) {
  stop("Not enough input arguments. Enter two arguments (genotype data and race info for target cohort).", call.=FALSE)
} else if (length(args)>=3) {
    stop("Too many arguments. Enter only two arguments (genotype data and race info for target cohort).", call.=FALSE)
} else if (length(args)==2) {
    target_genome_path <- args[1]
    target_race_path <- args[2]
}

library(tidyverse)
require(readr)
require(dplyr)
require(data.table)


##############################################################
###### set variables

plink <- "/PHShome/bm363/bin/plink2"

# reference genome and labels
ref_genome_path <- "/PHShome/bm363/ancestry_estimation/ref_1kgenome_prunned"
ref_race_path <- "/PHShome/bm363/ancestry_estimation/sam_info_1k.txt"

# the following will be entered as input arguments while running the script
# target genome and ancestry labels
# fc data
#target_genome_path <- "/PHShome/bm363/ancestry_estimation/FC_cleanest_IDs_converted"
#target_race_path <- "/PHShome/bm363/ancestry_estimation/fc_race2.txt"
# sb data
#target_genome_path <- "/PHShome/bm363/ancestry_estimation/sb1"
#target_race_path <- "/PHShome/bm363/ancestry_estimation/sb_race.txt"
# aurora data
#target_genome_path <- "/data/humgen/burook/data_raw/aur/au1"
#target_race_path <- "/PHShome/bm363/ancestry_estimation/aur_race2.txt"

##############################################################
### Merge reference and target datasets

system(paste(plink, " --bfile ", ref_genome_path, " --bmerge ", target_genome_path, ".bed ", target_genome_path, ".bim ", target_genome_path, ".fam --make-bed --out tmpmerged", sep = ""))
# flip problematic snps
system(paste(plink, " --bfile ", target_genome_path, " --flip tmpmerged-merge.missnp --make-bed --out tmpFlip", sep =""))
system(paste(plink, " --bfile ", ref_genome_path, " --bmerge tmpFlip.bed tmpFlip.bim tmpFlip.fam --make-bed --out tmpmerged", sep = ""))
# exclude those that persist
system(paste(plink, " --bfile tmpFlip --exclude tmpmerged-merge.missnp --make-bed --out tmpFlip2", sep =""))
system(paste(plink, " --bfile ", ref_genome_path, " --bmerge tmpFlip2.bed tmpFlip2.bim tmpFlip2.fam --make-bed --out tmpmerged", sep =""))


##############################################################
### QC target dataset
system(paste(plink, " --bfile tmpmerged --maf 0.01  --geno  0.01 --mind 0.1 --hwe 0.0000000001 --out tmpmerg", sep =""))


##############################################################
### Thinning of the merged data

system(paste(plink, " --bfile tmpmerged --maf 0.01  --indep-pairwise  50 10 0.2 --out tmpmerg", sep =""))
system(paste(plink, " --bfile tmpmerged --extract tmpmerg.prune.in --make-bed --out tmp_merged", sep =""))
system(paste("echo ", "Number of markers remaining: ", sep =""))
system(paste("wc -l tmp_merged.bim | cut -f1 -d \" \" ", sep =""))


##############################################################
### Estimate principal components

system(paste(plink, " --bfile tmp_merged --indep-pairwise 50 10 0.2 --out tmpprune2", sep=""))
system(paste(plink, " --bfile tmp_merged --extract tmpprune2.prune.in --genome --out tmpibs2 ", sep=""))
system(paste(plink, " --bfile tmp_merged --read-genome tmpibs2.genome --cluster --pca 5 --out tmppca5", sep=""))



##############################################################
### train and test SVM and tree-based classifier

library(e1071)
library(rpart)
set.seed(17)

eigenvec <- read_delim("tmppca5.eigenvec", col_names = FALSE, delim = " ");
eigenvec <- eigenvec[-1]
colnames(eigenvec) <- c("Sample","P1","P2","P3","P4","P5")

# read labels
ref_race <- read_delim(ref_race_path, col_names = TRUE, delim = "\t")
target_race <- read_delim(target_race_path, col_names = TRUE, delim = "\t")
colnames(target_race) <- c("Sample","Ethnicity_race")
target_race$Sample <- as.character(target_race$Sample)

tmpdata5 <- eigenvec %>% 
  left_join(ref_race, by="Sample") %>%
  left_join(target_race, by="Sample")

# let's save the raw data 
write.table(tmpdata5, "table_raw", row.names = FALSE, col.names = TRUE, quote = FALSE)

# exclude AMR population
tmpdata5 <- tmpdata5 %>% 
  filter(!(SuperPopulation %in% c("AMR")) | (is.na(SuperPopulation)))

# train SVM classifier
tmpsvm1 <- svm(SuperPopulation ~ .,
            data=tmpdata5[(!is.na(tmpdata5$SuperPopulation)),c("SuperPopulation","P1","P2")],
            type="C-classification")

tmpdata5$predsvm <- predict(tmpsvm1, tmpdata5[,c("P1","P2","P3","P4","P5")], decision.values = FALSE, probability = FALSE)


# train tree-based classifier
tmptree1 <- rpart(SuperPopulation ~ .,
            data=tmpdata5[(!is.na(tmpdata5$SuperPopulation)),c("SuperPopulation","P1","P2")],
            method="class", minsplit = 2, minbucket = 1)

tmpdata5$predtree <- predict(tmptree1, tmpdata5[,c("P1","P2")], type="class")


write.table(tmpdata5, "table_pred", row.names = FALSE, col.names = TRUE, quote = FALSE)



# clean up
system("rm tmp*", ignore.stdout=T,ignore.stderr=T)
