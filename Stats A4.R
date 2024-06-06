##BINF6079 Project by Cynthia Du, Jacob Hambly and Robin Zutshi


library(vcfR)
library(factoextra)
library(tidyverse)
library(randomForest)
library(caret)
library(genetics)
library(pROC)
library(plotly)

###### PROCESSING PRDM16 DATASET ------


#*******************************************************************************

# Use 'vcfR' package to read in vcf file
PRDM16_vcf <- read.vcfR("PRDM16(chr1;3069203-3438621).vcf")
PRDM16_vcf

# Have a look at each of the components in vcfR object
PRDM16_vcf@meta
PRDM16_vcf@fix
PRDM16_vcf@gt

# Check dimensions of vcf file
dim(PRDM16_vcf)

# Ensure number of variants is greater than 500
dim(PRDM16_vcf)[1] > 500

#*******************************************************************************

# Extract 'fixed' information from vcfR object into data frame
PRDM16_fix <- as.data.frame(PRDM16_vcf@fix)

# Remove Indels by filtering out observations with more than one allele in the Reference or Alternative column 
PRDM16_fix <- PRDM16_fix %>%
  filter(nchar(REF) == 1 & nchar(ALT) == 1)

#*******************************************************************************

# Extract genotype (GT) information from vcf data
gt_PRDM16 <- extract.gt(PRDM16_vcf, element = "GT", as.numeric = F)

# View genotype calls 
View(gt_PRDM16)

# Check dimensions of genotype data
dim(gt_PRDM16)

# Table of frequencies of genotype call
table(gt_PRDM16)

# Check if there are any NA values in the extracted genotypes
sum(is.na(gt_PRDM16))
which(is.na(gt_PRDM16))

# Remove all NA values in the genotypes
gt_PRDM16 <- as.data.frame(na.omit(gt_PRDM16))

# Encode genotype
gt_PRDM16[gt_PRDM16 == "0|0"] <- 0
gt_PRDM16[gt_PRDM16 == "0|1"] <- 1
gt_PRDM16[gt_PRDM16 == "1|0"] <- 1
gt_PRDM16[gt_PRDM16 == "1|1"] <- 2 

# View genotype calls 
View(gt_PRDM16)

#*******************************************************************************

# Filter genotype calls by keeping only observations of SNP ids remaining in the filtered fix data frame
gt_PRDM16 <- gt_PRDM16 %>%
  filter(rownames(gt_PRDM16) %in% PRDM16_fix$ID)

# Table of frequencies of genotype call
table(as.matrix(gt_PRDM16))

#*******************************************************************************

# Creating empty matrix to store our genotype and allele counts/frequencies
AF_matrix <- matrix(data = NA, nrow = nrow(gt_PRDM16), ncol = 10)

# Renaming columns of empty matrix
colnames(AF_matrix) <- c("SNP_id", "Observed_AA", "Observed_Aa", "Observed_aa", 
                         "p", "q", "Expected_AA", "Expected_Aa", "Expected_aa", 
                         "chi-square p_value")

# For loop to iterate over every genotype call 
for (i in 1:nrow(gt_PRDM16)) {
  # Sum up all count for all genotypes
  Observed_AA <- sum(gt_PRDM16[i,] == 0)
  Observed_Aa <- sum(gt_PRDM16[i,] == 1)
  Observed_aa <- sum(gt_PRDM16[i,] == 2)
  Observed_counts <- c(Observed_AA, Observed_Aa, Observed_aa)
  # Calculate p allele frequency
  p <- (2 * Observed_AA + Observed_Aa) / (2 * (Observed_AA + Observed_Aa + Observed_aa))
  # Calculate q allele frequency
  q <- 1 - p
  # Calculate expected values using p and q
  Expected_AA <- p^2 * nrow(gt_PRDM16)
  Expected_Aa <- 2 * p * q * nrow(gt_PRDM16)
  Expected_aa <- p^2 * nrow(gt_PRDM16)
  Expected_counts <- c(Expected_AA, Expected_Aa, Expected_aa)
  # Combine observed and expected counts into a contingency table
  contingency_table <- rbind(Observed_counts, Expected_counts) 
  # Perform chi-square test
  chi_square_result <- chisq.test(contingency_table)
  # Fill in alele matrix 
  AF_matrix[i,] <- c(rownames(gt_PRDM16)[i], Observed_AA, Observed_Aa, Observed_aa, p, 
                     q, Expected_AA, Expected_Aa, Expected_aa, chi_square_result$p.value)
}

# View the filled in allele frequency matrix
View(AF_matrix)

# Filter out SNPs with MAF less than 0.001 in the populations
AF_matrix <- as.data.frame(AF_matrix) %>%
  filter(q > 0.001)

# Filter out SNPs that violate the HWE assumptions in the populations
AF_matrix <- as.data.frame(AF_matrix) %>%
  filter(`chi-square p_value` > 0.05)

#*******************************************************************************

# Filter genotype calls by keeping only observations of SNP ids remaining in the allele frequency matrix, filtered by MAF > 0.001 and HWE assumptions held
gt_PRDM16 <- gt_PRDM16 %>%
  filter(rownames(gt_PRDM16) %in% AF_matrix$SNP_id)

#*******************************************************************************

# Load in sample metadata file
library(readr)
sample_metadata <- read_tsv("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")[1:4]

# View sample metadata
View(sample_metadata)

# Creating new variable of sample names from extracted genotypes
sample_ids <- colnames(gt_PRDM16)

# Filtering all IGSR sample metadata to only include samples in study
sample_metadata <- sample_metadata %>%
  filter(sample %in% sample_ids)

# Counts by different factors
table(sample_metadata$pop)
table(sample_metadata$super_pop)
table(sample_metadata$gender)

#*******************************************************************************

# Transpose genotype call matrix 
gt_PRDM16 <- t(gt_PRDM16)

# Convert all columns to numeric class, ensuring preservation of row names
gt_PRDM16 <- as.data.frame(lapply(gt_PRDM16, function(x) as.numeric(as.character(x))), 
                           row.names = row.names(gt_PRDM16))

# Write out filtered genotype calls to csv file for PCA analysis
write.csv(gt_PRDM16, file = "prdm16_pca.csv")



###### PROCESSING PAX3 DATASET ------



#loading in gbr vcf
pax3 <- read.vcfR("allpop_pax3_GRCh38.vcf.gz")

#create a pax3 tidy data frame 
pax3_tidy <- vcfR2tidy(pax3, single_frame = TRUE, info_types = TRUE, format_types = TRUE)

#create a data frame for the info columns
pax3_info <- pax3_tidy$dat

#filter the info data frame based on single alleles in the REF and ALT and 'SNP' variant types
pax3_info <- pax3_info %>%
  select(ID, REF, ALT, VT, AF) %>%
  filter(nchar(REF) == 1) %>%
  filter(nchar(ALT) == 1) %>%
  filter(VT == "SNP")

#extracting genotypes
pax3_snps <- extract.gt(pax3, 
                        element = "GT",
                        IDtoRowNames  = T,
                        as.numeric = F,
                        convertNA = T,
                        return.alleles = F)

#making the object containing the genotypes into a dataframe 
pax3_snps <- data.frame(pax3_snps)

#filter the data frame to only contain the filtered SNPs from the info dataframe 
filtered_pax3_snps <- pax3_snps %>%
  filter(row.names(.) %in% pax3_info$ID)

#encode the genotypes into three different levels where Homozygous reference is 0, heterozygous reference-alternate is 1 and homozygous alternate is 2.
filtered_pax3_snps <- as.data.frame(apply(filtered_pax3_snps, c(1,2), function(x) {
  ifelse(x == '0|0', 0, 
         ifelse(x == '1|0' | x == '0|1', 1, 
                ifelse(x == '1|1', 2, x)))
}))


#making a matrix for the minor allele frequency for each SNP
alternate_allele_frequency <- matrix(NA, nrow = nrow(filtered_pax3_snps), ncol = 2)
alternate_allele_frequency <- as.data.frame(alternate_allele_frequency)

#edit the column names and row names
colnames(alternate_allele_frequency) <- c("p_count", "q_count")
rownames(alternate_allele_frequency) <- rownames(filtered_pax3_snps)

#find the p count for each row in filtered_pax3_snps
alternate_allele_frequency$p_count <- apply(filtered_pax3_snps, 1, function(row) {
  sum_0 <- sum(row == 0)
  sum_1 <- sum(row == 1)
  sum_2 <- sum(row == 2)
  p_count <- (2 * sum_0) + sum_1
  q_count <- (2 * sum_2) + sum_1
  return(p_count)
  
}
)

#find the q count for each row in filtered_pax3_snps  
alternate_allele_frequency$q_count <- apply(filtered_pax3_snps, 1, function(row) {
  sum_0 <- sum(row == 0)
  sum_1 <- sum(row == 1)
  sum_2 <- sum(row == 2)
  p_count <- (2 * sum_0) + sum_1
  q_count <- (2 * sum_2) + sum_1
  return(q_count)
  
}
)  

#creating a MAF column
alternate_allele_frequency$MAF <- alternate_allele_frequency$q_count / (alternate_allele_frequency$q_count + alternate_allele_frequency$p_count) 

#creating a new dataframe to filter out any SNPs with a MAF < 0.001
alternate_allele_frequency_filtered <- alternate_allele_frequency %>%
  filter(MAF > 0.001)

#filtering pax3 SNPs based on alternate_allele_frequency_filtered
filtered_pax3_snps_MAF <- filtered_pax3_snps %>%
  filter(row.names(.) %in% row.names(alternate_allele_frequency_filtered))


#creating a matrix to store the genotypes across the samples for each SNP
genotypes <- matrix(NA, nrow = nrow(filtered_pax3_snps_MAF), ncol = 3)

#re-naming the column lables with the genotype labels
colnames(genotypes) <- c("AA", "Aa", "aa")

#adding the record of the number of each of the genotypes for each SNP
for (i in 1:dim(filtered_pax3_snps_MAF)[1]) {
  genotype_counts <- table(factor(filtered_pax3_snps_MAF[i,], levels = 0:2))
  genotypes[i,] <- genotype_counts
}

##calculating the expected HWE for each SNP across the samples and using chi-squared to compare with the observed
hwe_pvalues <- matrix(NA, nrow = nrow(filtered_pax3_snps_MAF), ncol = 1)

colnames(hwe_pvalues) <- "p-value"
rownames(hwe_pvalues) <- rownames(filtered_pax3_snps_MAF)

for (i in 1:dim(filtered_pax3_snps_MAF)[1]) {
  
  observed_counts <- as.numeric(genotypes[i,])
  
  p <- 1 - alternate_allele_frequency_filtered[i,3]
  
  expected_AA <- as.numeric(p^2 * ncol(filtered_pax3_snps_MAF))
  
  expected_Aa <- as.numeric(2*p*ncol(filtered_pax3_snps_MAF) * (1-p))
  
  expected_aa <- as.numeric(((1-p)^2) * ncol(filtered_pax3_snps_MAF))
  
  expected_counts <- c(expected_AA,expected_Aa, expected_aa )
  
  # Combine observed and expected counts into a contingency table
  contingency_table <- rbind(observed_counts, expected_counts)
  
  # Perform chi-square test
  chi_square_result <- chisq.test(contingency_table)
  
  hwe_pvalues[i,] <- chi_square_result$p.value
  
}

#filtering hwe_pvalues
hwe_pvalues <- as.data.frame(hwe_pvalues)

hwe_pvalues_filtered <- hwe_pvalues %>%
  filter(`p-value` >= 0.05)


#filtering pax3 SNPs for those that meet the assumptions of HWE based on observed genotype counts
filtered_pax3_snps_HWE <- filtered_pax3_snps_MAF %>%
  filter(row.names(.) %in% row.names(hwe_pvalues_filtered))


##transforming data so that it is suitable for PCA
t_pax3 <- t(filtered_pax3_snps_HWE)

#exporting t_pax3 as a csv file
write.csv(t_pax3, file = "pax3_pca.csv")




###### PROCESSING TP63 DATASET------
#### LOAD DATA ----

#load raw dataframe
dfRawTP63 <- read.vcfR("TP63 GRCh38 Final.vcf")
dfRawPRDM16 <- read.vcfR("PRDM16(chr1;3069203-3438621).vcf.gz")
dfRawPAX3 <- read.vcfR("allpop_pax3_GRCh38.vcf.gz")


#load sample information

metaSubpop <- read_tsv("https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")[1:2]

#CHECK: check no NAs
metaSubpop %>% is.na() %>% table()

#total number of samples
NSAMPLES <- dim(dfRawTP63)[3]-1 #note: -1 is because of GT column


#### CLEAN DATA ----



#extract information and genotype data and combine both in main dataframe (only SNP VT rows)
dfMain <- cbind(dfRawTP63@fix[, c("CHROM", "POS", "ID", "REF", "ALT", "INFO")], dfRawTP63@gt) %>%
  rbind(., cbind(dfRawPRDM16@fix[, c("CHROM", "POS", "ID", "REF", "ALT", "INFO")], dfRawPRDM16@gt)) %>%
  rbind(., cbind(dfRawPAX3@fix[, c("CHROM", "POS", "ID", "REF", "ALT", "INFO")], dfRawPAX3@gt)) 


#change encoding of all genotypes (0 = homozyous ref, 1 = heterozygous, 2 = homozygous alt)
dfMain <- dfMain %>% as_tibble() %>% filter(nchar(ALT) == 1 & nchar(REF)==1) %>%   #remove multivariant SNPs
  mutate(across (!c(1:7), ~ case_when( .=="0|0" ~ 0,
                                       .=="0|1" ~ 1, 
                                       .=="1|0" ~ 1, 
                                       .=="1|1" ~ 2,
                                       .default = 123456789)))


#Calculating p, q, and filtering based on MAF
dfCleaned <- dfMain %>% 
  mutate(q = rowSums(.[,-c(1:7)])/(NSAMPLES*2), .after="FORMAT") %>%  #This is also the allele frequency
  mutate(p = 1 - q, .after="FORMAT") %>%                              #calculate p
  filter (q > 0.001 & q < 0.999)                                      #filter for MAF > 0.001

#CHECK: check that no q outside expected range
dfCleaned %>% filter(q>1 | q<0) %>% dim()


#Filter based on hardy weinberg (2609 - 2557 = 52 SNPs filtered out)

dfCleaned <- dfCleaned %>% 
  
  #Count actual genotypes
  mutate(act_p2 = rowSums(apply(.[,-c(1:9)], 2, function(x) x == 0)), .after="q") %>%
  mutate(act_2pq = rowSums(apply(.[,-c(1:10)], 2, function(x) x == 1)), .after="q") %>%
  mutate(act_q2 = rowSums(apply(.[,-c(1:11)], 2, function(x) x == 2)), .after="q") %>%
  
  #Calculate HW predicted genotypes (ceiling to reduce divide by 0 errors)
  mutate(hw_p2 = ceiling(p*p*NSAMPLES), .after="q") %>%
  mutate(hw_2pq = ceiling(2*p*q*NSAMPLES), .after="q") %>%
  mutate(hw_q2 = ceiling(q*q*NSAMPLES), .after="q") %>% 
  
  #chi-squared test
  rowwise() %>%
  mutate(chi_squared_p =  chisq.test(rbind(c(act_p2, act_2pq, act_q2), c(hw_p2, hw_2pq, hw_q2)))$p.value, .after="q") 

#Visual inspection of NAN values reveals that they should be assigned a value of 1
dfCleaned %>% filter(is.nan(chi_squared_p)) %>% View()

#Filter out non-hardy weignberg SNPs
dfCleaned <- dfCleaned %>% 
  mutate( chi_squared_p = ifelse(is.na(chi_squared_p), 1.0, chi_squared_p) ) %>%
  filter( chi_squared_p > 0.05 ) %>%
  dplyr::select(!(chi_squared_p:act_p2))



#TESTING
# dfTest <- dfCleaned[, -c(1:2,4:9)] %>% 
#   pivot_longer(-"ID") %>% 
#   pivot_wider(names_from="ID", values_from = value) %>%     
#   inner_join(x=metaSubpop, y=., by = join_by("sample" == "name")) 
# 
# dfP <- aggregate(dfTest[,-c(1:2)], list(dfTest$pop), sum)[,-1]
# dfP <- apply(dfP, c(1, 2), function(x) x/NSAMPLES*2)
# dfP %>% select_if(.>0.001)



#transpose rows and columns & add population group labels
dfTransposed <- dfCleaned[, -c(1:2,4:9)] %>%                          #take only genotype columns
  pivot_longer(-"ID") %>%                                             #swap column and rows
  pivot_wider(names_from="ID", values_from = value) %>%               #add ID column
  inner_join(x=metaSubpop, y=., by = join_by("sample" == "name"))     #add population group labels from meta data


###Rough work
#apply(dfMutated[,-c(1:7)], 1, sum)
#mutate(subPop = metaSubpop$pop[metaSubpop$sample==name], .after="name")
#mutate(subPop = name %in% metaSubpop$sample)
#metaSubpop$pop[metaSubpop$sample=="HG03603"]
#write.csv(dfTransposed, "tp63_pca.csv", row.names = FALSE)

dfFinal=dfTransposed

###### COMBINING CLEANED DATASETS ------

dfTP63<-read.csv("tp63_pca.csv", header = TRUE)
dfPAX3<- read.csv("pax3_pca.csv", header = TRUE)
dfPRDM16<-read.csv("prdm16_pca.csv", header = TRUE)

dfFinal <- inner_join(x=dfTP63, y=dfPAX3, by="sample")
dfFinal <- inner_join(x=dfFinal, y=dfPRDM16, by=join_by("sample"=="X"))

dfFinalNoBEB <- dfFinal %>% filter(pop != "BEB")

#### PCA ----

#compute PCs (uncorrelated linear combination explaining as much variance as possible)
pc <- prcomp(dfFinal[, -c(1:2)], center = TRUE, scale. = FALSE)

#these are the loading vectors (eg. the eigen vectors, explains contribution to the PC) 
#View(pc$rotation)
#check that they are normalized (yes they are, with minor calculation differences)
apply(pc$rotation, 2, function(i) sqrt(sum(i^2))) %>% table()


#PLOT

#Scree Plot (proportion of variance for each PC)
p <- fviz_eig(pc, xlab = "Principal component", addlabels=TRUE, hjust = -0.1, main = "Figure B1. Scree Plot") +
  ylim(0, 10)
print(p)

#Cumulative sum of PC variance plot
plot( 1:30 , summary(pc)$importance[3,c(1:30)], ylim = c(0, 1), 
      main = "Figure B2. Cumulative sum of variance over PC", 
      ylab = "Cumulative variance explained",
      xlab="Principal component", type = "b", xaxt="n")
grid(nx = 30)
axis(1, at = seq(1, 30, by = 1), las=2)

#PCA BIPLOT
fviz_pca_biplot(pc, repel = TRUE,
                select.var = list(contrib=15),
                habillage=dfFinal$pop,
                addEllipses = TRUE,
                label = "var",
                geom = "point", 
                pointsize = 3, 
                title = "Figure B5. PCA biplot by subpopulation with top 15 contributing genes",
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


#poor kmeans clustering
k2 <-kmeans(pc$x[,c(1:2)], centers=3, nstart = 25)
fviz_cluster(k2, data=pc$x[,c(1:2)])


####RANDOM FOREST ----

#split train-test set
set.seed(1000)
testIndex <- sample(2, 282, replace=TRUE, prob= c(0.7, 0.3))

dfFinal$pop <- as.factor(dfFinal$pop)
train <- dfFinal[testIndex==1, ]
test <- dfFinal[testIndex==2, ]


#### ntrees ----

#try random forest with 2 different tree sets
set.seed(1000)
system.time(rf500trees <- randomForest(pop ~., data = train[,-1], proximity=TRUE, ntree=500))   
system.time(rf1500trees <- randomForest(pop ~., data = train[,-1], proximity=TRUE, ntree=1500)) 

#check details
print(rf500trees) 
print(rf1500trees)

#plot
plot(rf500trees)
plot(rf1500trees, main="Random forest error rate with default parameters across 1500 trees")
#NOTE: not much improvement after 400 trees. 500 trees is a good measure
#NOTE: as expected BEB is difficult to converge


#set ntree to use
bestNtree <- 500


### mtrys ----

#tuning mtrys (default is about sqrt(8000) = 80ish )
mtryGrid <- data.frame(x = c(30,seq(from = 65, to = 75, by = 2),150,200)) #num of variables randomly sampled per split

for (i in mtryGrid$x) {
  tempmodel <- randomForest(pop ~., data = train[,-1], mtry=i, ntree = bestNtree, proximity=TRUE)
  mtryGrid$performance[ which(mtryGrid$x == i) ]<- tempmodel$err.rate[nrow(tempmodel$err.rate),1]
}

#convergest at 70
plot(mtryGrid$x, mtryGrid$performance, xlab="mtry", ylab="error rate")
bestMtry <- 70

#save best model
bestmodel <- randomForest(pop ~., data = train[,-1], mtry=bestMtry, ntree = bestNtree, proximity=TRUE)

### actual performance ----

#ACTUAL PERFORMANCE ON TEST DATA
p1 <- predict(bestmodel, test[,-1]) 
confusionMatrix(p1, test[,-1]$pop)
#score of: 95%

#do I need this?
p2 <- predict(bestmodel, newdata=test, type="prob")[,1]
multiclass.roc(test$pop,p2)

#variable importance
varImpPlot(bestmodel)
plot(varImp(bestmodel))



### nodesize ----

#tuning nodesize (default is about sqrt(8000) = 80ish )
nodesizeGrid <- data.frame(x = seq(from = 1, to = 25, by = 5)) 

for (i in nodesizeGrid$x) {
  tempmodel <- randomForest(pop ~., data = train[,-1], mtry=bestMtry, ntree=bestNtree, nodesize=i, proximity=TRUE)
  nodesizeGrid$performance[ which(nodesizeGrid$x == i) ]<- tempmodel$err.rate[nrow(tempmodel$err.rate),1]
}

#convergest at 70
plot(nodesizeGrid$x, nodesizeGrid$performance, xlab="mtry", ylab="error rate")
#bestNodesize <- 

### extra stuff ----

#Tune both

tunegrid <- expand.grid(mtry=seq(50,90,20),nodesize=seq(1,20,5))

#set.seed(88888)
for (i in 1:dim(tunegrid)[1]) {
  set.seed(88888)
  tempmodel <- randomForest(pop ~., data = train[,-1], 
                            mtry=tunegrid$mtry[i], 
                            nodesize=tunegrid$nodesize[i], 
                            ntree = bestNtree, proximity=TRUE)
  tunegrid$train_accuracy[i]<-  1 - tempmodel$err.rate[nrow(tempmodel$err.rate),1] 
  predict_test_rf <- predict(tempmodel, test[,-c(1,2)])
  tunegrid$test_accuracy[i] <- confusionMatrix(predict_test_rf, test$pop)$overall[1]
                                           
  #tunegrid$actual_unseen <- test[,2]
  print(c("model completed ",i))
}

#Save best model
bestIndex <- 9
set.seed(88888)
finalBestModel <- randomForest(pop ~., data = train[,-1], 
                      mtry=tunegrid$mtry[bestIndex], 
                      nodesize=tunegrid$nodesize[bestIndex], 
                      ntree = bestNtree, proximity=TRUE)

finalBestModel


#plot
plot_ly(data = tunegrid, x = ~mtry, y = ~nodesize, z = ~train_accuracy, color = ~test_accuracy) %>%
  layout(title = "3D plot of training for random forest optimization")

#heatmap for best model


#ACTUAL PERFORMANCE ON TEST DATA
p1 <- predict(finalBestModel, test[,-1]) 
confusionMatrix(p1, test[,-1]$pop)
#score of: 95%

#variable importance
top10 <- varImp(finalBestModel) %>% arrange(desc(Overall)) %>% 
  rownames_to_column(., var="snp")
top10 <- top10[c(1:10),]
top10

dfRawPAX3@fix

varImpPlot(finalBestModel)





