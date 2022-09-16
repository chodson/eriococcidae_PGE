#setwd("/Users/christina/projects/cystococcus_snp/")
#library(lme4)
library(nlme)
#library(ggeffects)
library(stargazer)
library(ggplot2)
library(lmerTest)

# Summary of what script does

#we want to know whether males show heterozygous allele expression (i.e. express more than one allele per SNP)
#that would indicate they express alleles inherited from their father, which we wouldn't expect if they have heterochromatic bodies in somatic cells
#we also are looking at the mothers as a comparison since we expect females to express both maternally and paternally inherited SNPs on a whole genome scale
#script below reads in  GATK ASEreadcounter output for each sample which summarizes the reference and alternate allele counts for each SNP
#we apply some filters to the data (keeping only SNPs with a depth of 40 and a low proportion of other bases (i.e. errors))
#then we generate histograms of the distribution of SNPs with uniparental vs. biparental expression for each sample 
#although note we cannot assign whether alleles are paternal or maternal in origin, we look at overall distribution of allelic bias (i.e. one allele expressed vs. two alleles expressed)



### SNP count files for Cystococcus campanidosalis
### Cystococcus campanidorsalis- family LGC3538

#output from ASEreadcounter
ccamp.f.snp<- read.delim("data/transcriptome.vs.3538.f.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953
ccamp.m1.snp<- read.delim("data/transcriptome.vs.3538.m1.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953
ccamp.m2.snp<- read.delim("data/transcriptome.vs.3538.m2.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953


#### Filters applied to count data
#prop.other=(otherBases/(otherBases+totalCount))
#prop.other > 0.05 
#totalCount> 40

ccamp.f.snp$prop.other<-ccamp.f.snp$otherBases/(ccamp.f.snp$otherBases+ccamp.f.snp$totalCount) 
hist(ccamp.f.snp$prop.other, breaks=1000, xlim=c(0,0.2)) ## filter <0.05

ccamp.f.snp<-ccamp.f.snp[ccamp.f.snp$prop.other<0.05,] #31758
ccamp.f.snp$prop.het<-ccamp.f.snp$altCount/ccamp.f.snp$totalCount
hist(ccamp.f.snp$prop.het, breaks = 100)
# het 0.2-0.8
# hom <0.1
hist(ccamp.f.snp$totalCount, breaks=1000, xlim=c(0,300))
ccamp.f.snp<-ccamp.f.snp[ccamp.f.snp$totalCount>40,] #filtering SNPs with totalCount <40, 25398

het.3538.cnp<-ccamp.f.snp[ccamp.f.snp$prop.het>0.2,]
het.3538.cnp<-het.3538.cnp[het.3538.cnp$prop.het<0.8,]#17492 het snps
head(het.3538.cnp)
het.3538.snp.id<-het.3538.cnp$contig

hom.3538.snp<-ccamp.f.snp[ccamp.f.snp$prop.het<0.05,]
hom.3538.snp2<-ccamp.f.snp[ccamp.f.snp$prop.het>0.9,]
hom.3538.snp3<-rbind(hom.3538.snp,hom.3538.snp2)#4961
hom.3538.snp.id<-hom.3538.snp3$contig

hom.3538.snp4<-ccamp.f.snp[ccamp.f.snp$prop.het<0.01,]##2028
hom.3538.snp4.id<-hom.3538.snp$contig

#taking just contig and position IDs for merging

hom.3538.snp<-hom.3538.snp[,c(1,2)]
het.3538.cnp<-het.3538.cnp[,c(1,2)]

#write.table(hom.3538.snp, file='data/hom.3538.f.snp.tsv',sep='\t', quote = F, row.names = F)
#write.table(het.3538.cnp, file='data/het.3538.f.snp.tsv',sep='\t', quote = F, row.names = F)

head(ccamp.m1.snp)
head(ccamp.m2.snp)
ccamp.m1.snp$prop.other<-ccamp.m1.snp$otherBases/(ccamp.m1.snp$otherBases+ccamp.m1.snp$totalCount) 
ccamp.m2.snp$prop.other<-ccamp.m2.snp$otherBases/(ccamp.m2.snp$otherBases+ccamp.m2.snp$totalCount) 
ccamp.m1.snp<-ccamp.m1.snp[ccamp.m1.snp$prop.other<0.05,]
ccamp.m2.snp<-ccamp.m2.snp[ccamp.m2.snp$prop.other<0.05,]
hist(ccamp.m2.snp$totalCount, breaks=1000, xlim=c(0,300))
ccamp.m1.snp<-ccamp.m1.snp[ccamp.m1.snp$totalCount>40,] #filtering SNPs with totalCount <40, #26841
ccamp.m2.snp<-ccamp.m2.snp[ccamp.m2.snp$totalCount>40,] #filtering SNPs with totalCount <40,  #23272
ccamp.m1.snp$prop.het<-ccamp.m1.snp$altCount/ccamp.m1.snp$totalCount
ccamp.m2.snp$prop.het<-ccamp.m2.snp$altCount/ccamp.m2.snp$totalCount


#histograms of prop.het for ccamp.sp
hist(ccamp.m2.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,2300))
hist(ccamp.m1.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,2300))
hist(ccamp.f.snp$prop.het, breaks = 100, col='#BF6B73', ylim=c(0,2300))

#ccamp.fhom.m1<-merge(hom.3538.snp4.id,ccamp.m1.snp)



##################################################################################
### Next family- Cystococcus echiniformis family LGC3571F5
#transcriptome.vs.3571F5.f.csv
#transcriptome.vs.3571F5.m3.csv
#transcriptome.vs.3571F5.m4.csv
# files
cech.3571f5.f.snp<- read.delim("data/transcriptome.vs.3571F5.f.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953
cech.3571f5.m3.snp<- read.delim("data/transcriptome.vs.3571F5.m3.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953
cech.3571f5.m4.snp<- read.delim("data/transcriptome.vs.3571F5.m4.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953

dim(cech.3571f5.f.snp)#17519
dim(cech.3571f5.m3.snp)#18672
dim(cech.3571f5.m4.snp)#17235 

cech.3571f5.f.snp$prop.other<-cech.3571f5.f.snp$otherBases/(cech.3571f5.f.snp$otherBases+cech.3571f5.f.snp$totalCount)
cech.3571f5.m3.snp$prop.other<-cech.3571f5.m3.snp$otherBases/(cech.3571f5.m3.snp$otherBases+cech.3571f5.m3.snp$totalCount)
cech.3571f5.m4.snp$prop.other<-cech.3571f5.m4.snp$otherBases/(cech.3571f5.m4.snp$otherBases+cech.3571f5.m4.snp$totalCount)

hist(cech.3571f5.f.snp$prop.other, breaks=1000, xlim=c(0,0.2))
hist(cech.3571f5.m3.snp$prop.other, breaks=1000, xlim=c(0,0.2))
hist(cech.3571f5.m4.snp$prop.other, breaks=1000, xlim=c(0,0.2))

cech.3571f5.f.snp<-cech.3571f5.f.snp[cech.3571f5.f.snp$prop.other<0.05,]
cech.3571f5.m3.snp<-cech.3571f5.m3.snp[cech.3571f5.m3.snp$prop.other<0.05,]
cech.3571f5.m4.snp<-cech.3571f5.m4.snp[cech.3571f5.m4.snp$prop.other<0.05,]

dim(cech.3571f5.f.snp)#17429 
dim(cech.3571f5.m3.snp)#18511 
dim(cech.3571f5.m4.snp)#17115

cech.3571f5.f.snp<-cech.3571f5.f.snp[cech.3571f5.f.snp$totalCount>40,]
cech.3571f5.m3.snp<-cech.3571f5.m3.snp[cech.3571f5.m3.snp$totalCount>40,]
cech.3571f5.m4.snp<-cech.3571f5.m4.snp[cech.3571f5.m4.snp$totalCount>40,]
dim(cech.3571f5.f.snp)#  12679 
dim(cech.3571f5.m3.snp)#  15191
dim(cech.3571f5.m4.snp)#  13271

cech.3571f5.f.snp$prop.het<-cech.3571f5.f.snp$altCount/cech.3571f5.f.snp$totalCount
cech.3571f5.m3.snp$prop.het<-cech.3571f5.m3.snp$altCount/cech.3571f5.m3.snp$totalCount
cech.3571f5.m4.snp$prop.het<-cech.3571f5.m4.snp$altCount/cech.3571f5.m4.snp$totalCount

hist(cech.3571f5.m4.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,3500))
hist(cech.3571f5.m3.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,3500))
hist(cech.3571f5.f.snp$prop.het, breaks = 100, col='#BF6B73', ylim=c(0,3500))

#hom and het snp list
hom.cech.3571f5.f.snp<-cech.3571f5.f.snp[cech.3571f5.f.snp$prop.het<0.05,]
het.cech.3571f5.f.snp<-cech.3571f5.f.snp[cech.3571f5.f.snp$prop.het>0.2,]
het.cech.3571f5.f.snp<-het.cech.3571f5.f.snp[het.cech.3571f5.f.snp$prop.het<0.8,]

hom.cech.3571f5.f.snp<-hom.cech.3571f5.f.snp[,c(1,2)]
het.cech.3571f5.f.snp<-het.cech.3571f5.f.snp[,c(1,2)]

#write.table(hom.cech.3571f5.f.snp, file='data/hom.cech.3571f5.f.snp',sep='\t', quote = F, row.names = F)
#write.table(het.cech.3571f5.f.snp, file='data/het.cech.3571f5.f.snp',sep='\t', quote = F, row.names = F)





##################################################################################
### Next family- Cystococcus echiniformis family LGC3571F6
#transcriptome.vs.3571F6.f.csv
#transcriptome.vs.3571F6.m.csv
#transcriptome.vs.3571F6.m2.csv
# files
cech.3571F6.f.snp<- read.delim("data/transcriptome.vs.3571F6.f.csv", comment.char = '#', header=T, stringsAsFactors = F)# 
cech.3571F6.m.snp<- read.delim("data/transcriptome.vs.3571F6.m.csv", comment.char = '#', header=T, stringsAsFactors = F)# 
cech.3571F6.m2.snp<- read.delim("data/transcriptome.vs.3571F6.m2.csv", comment.char = '#', header=T, stringsAsFactors = F)# 

dim(cech.3571F6.f.snp) #15870 
dim(cech.3571F6.m.snp)#15510 
dim(cech.3571F6.m2.snp)#12342 

####### Filters applied to raw data
# prop.other >0.05
# total count> 40

cech.3571F6.f.snp$prop.other<-cech.3571F6.f.snp$otherBases/(cech.3571F6.f.snp$otherBases+cech.3571F6.f.snp$totalCount)
cech.3571F6.m.snp$prop.other<-cech.3571F6.m.snp$otherBases/(cech.3571F6.m.snp$otherBases+cech.3571F6.m.snp$totalCount)
cech.3571F6.m2.snp$prop.other<-cech.3571F6.m2.snp$otherBases/(cech.3571F6.m2.snp$otherBases+cech.3571F6.m2.snp$totalCount)

hist(cech.3571F6.f.snp$prop.other, breaks=1000, xlim=c(0,0.2))
hist(cech.3571F6.m.snp$prop.other, breaks=1000, xlim=c(0,0.2))
hist(cech.3571F6.m2.snp$prop.other, breaks=1000, xlim=c(0,0.2))

cech.3571F6.f.snp<-cech.3571F6.f.snp[cech.3571F6.f.snp$prop.other<0.05,]
cech.3571F6.m.snp<-cech.3571F6.m.snp[cech.3571F6.m.snp$prop.other<0.05,]
cech.3571F6.m2.snp<-cech.3571F6.m2.snp[cech.3571F6.m2.snp$prop.other<0.05,]

dim(cech.3571F6.f.snp)#15802
dim(cech.3571F6.m.snp)#15309 
dim(cech.3571F6.m2.snp)#12244 

cech.3571F6.f.snp<-cech.3571F6.f.snp[cech.3571F6.f.snp$totalCount>40,]
cech.3571F6.m.snp<-cech.3571F6.m.snp[cech.3571F6.m.snp$totalCount>40,]
cech.3571F6.m2.snp<-cech.3571F6.m2.snp[cech.3571F6.m2.snp$totalCount>40,]
dim(cech.3571F6.f.snp)#  10681
dim(cech.3571F6.m.snp)# 11624
dim(cech.3571F6.m2.snp)#   8155

cech.3571F6.f.snp$prop.het<-cech.3571F6.f.snp$altCount/cech.3571F6.f.snp$totalCount
cech.3571F6.m.snp$prop.het<-cech.3571F6.m.snp$altCount/cech.3571F6.m.snp$totalCount
cech.3571F6.m2.snp$prop.het<-cech.3571F6.m2.snp$altCount/cech.3571F6.m2.snp$totalCount

hist(cech.3571F6.m2.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,2700))
hist(cech.3571F6.m.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,2700))
hist(cech.3571F6.f.snp$prop.het, breaks = 100, col='#BF6B73', ylim=c(0,2700))

#hom and het snp list
hom.cech.3571F6.f.snp<-cech.3571F6.f.snp[cech.3571F6.f.snp$prop.het<0.05,]
het.cech.3571F6.f.snp<-cech.3571F6.f.snp[cech.3571F6.f.snp$prop.het>0.2,]
het.cech.3571F6.f.snp<-het.cech.3571F6.f.snp[het.cech.3571F6.f.snp$prop.het<0.8,]

hom.cech.3571F6.f.snp<-hom.cech.3571F6.f.snp[,c(1,2)]
het.cech.3571F6.f.snp<-het.cech.3571F6.f.snp[,c(1,2)]


#write.table(hom.cech.3571F6.f.snp, file='data/hom.cech.3571F6.f.snp',sep='\t', quote = F, row.names = F)
#write.table(het.cech.3571F6.f.snp, file='data/het.cech.3571F6.f.snp',sep='\t', quote = F, row.names = F)







##################################################################################
### Next family - Cystococcus echiniformis family LGC3572F4
#transcriptome.vs.3572F4.f.csv
#transcriptome.vs.3572F4.m2.csv
#transcriptome.vs.3572F4.m3.csv


cech.3572F4.f.snp<- read.delim("data/transcriptome.vs.3572F4.f.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953
cech.3572F4.m2.snp<- read.delim("data/transcriptome.vs.3572F4.m2.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953
cech.3572F4.m3.snp<- read.delim("data/transcriptome.vs.3572F4.m3.csv", comment.char = '#', header=T, stringsAsFactors = F)# 31953

dim(cech.3572F4.f.snp)#15981
dim(cech.3572F4.m2.snp)#18101 
dim(cech.3572F4.m3.snp)#17785  

####### Filters applied to raw data
# prop.other >0.05
# total count> 40

cech.3572F4.f.snp$prop.other<-cech.3572F4.f.snp$otherBases/(cech.3572F4.f.snp$otherBases+cech.3572F4.f.snp$totalCount)
cech.3572F4.m2.snp$prop.other<-cech.3572F4.m2.snp$otherBases/(cech.3572F4.m2.snp$otherBases+cech.3572F4.m2.snp$totalCount)
cech.3572F4.m3.snp$prop.other<-cech.3572F4.m3.snp$otherBases/(cech.3572F4.m3.snp$otherBases+cech.3572F4.m3.snp$totalCount)

hist(cech.3572F4.f.snp$prop.other, breaks=1000, xlim=c(0,0.2))
hist(cech.3572F4.m2.snp$prop.other, breaks=1000, xlim=c(0,0.2))
hist(cech.3572F4.m3.snp$prop.other, breaks=1000, xlim=c(0,0.2))

cech.3572F4.f.snp<-cech.3572F4.f.snp[cech.3572F4.f.snp$prop.other<0.05,]
cech.3572F4.m2.snp<-cech.3572F4.m2.snp[cech.3572F4.m2.snp$prop.other<0.05,]
cech.3572F4.m3.snp<-cech.3572F4.m3.snp[cech.3572F4.m3.snp$prop.other<0.05,]

dim(cech.3572F4.f.snp)#15879
dim(cech.3572F4.m2.snp)#17896 
dim(cech.3572F4.m3.snp)#17568

cech.3572F4.f.snp<-cech.3572F4.f.snp[cech.3572F4.f.snp$totalCount>40,]
cech.3572F4.m2.snp<-cech.3572F4.m2.snp[cech.3572F4.m2.snp$totalCount>40,]
cech.3572F4.m3.snp<-cech.3572F4.m3.snp[cech.3572F4.m3.snp$totalCount>40,]
##looks fine, most very small, could filter <0.05
dim(cech.3572F4.f.snp)#  11803
dim(cech.3572F4.m2.snp)#  14312 
dim(cech.3572F4.m3.snp)#  14202 

cech.3572F4.f.snp$prop.het<-cech.3572F4.f.snp$altCount/cech.3572F4.f.snp$totalCount
cech.3572F4.m2.snp$prop.het<-cech.3572F4.m2.snp$altCount/cech.3572F4.m2.snp$totalCount
cech.3572F4.m3.snp$prop.het<-cech.3572F4.m3.snp$altCount/cech.3572F4.m3.snp$totalCount

hist(cech.3572F4.f.snp$prop.het, breaks=100, col='#BF6B73', ylim=c(0,3000))
hist(cech.3572F4.m2.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,3000))
hist(cech.3572F4.m3.snp$prop.het, breaks = 100, col='#41719C', ylim=c(0,3000))

#hom and het snp list
hom.cech.3572F4.f.snp<-cech.3572F4.f.snp[cech.3572F4.f.snp$prop.het<0.05,]
het.cech.3572F4.f.snp<-cech.3572F4.f.snp[cech.3572F4.f.snp$prop.het>0.2,]
het.cech.3572F4.f.snp<-het.cech.3572F4.f.snp[het.cech.3572F4.f.snp$prop.het<0.8,]#17492 het snps



#write.table(hom.cech.3572F4.f.snp, file='data/hom.cech.3572F4.f.snp',sep='\t', quote = F, row.names = F)
#write.table(het.cech.3572F4.f.snp, file='data/het.cech.3572F4.f.snp',sep='\t', quote = F, row.names = F)





### Comparing het and hom proportions for males and females
#I'm considering a allelic bias (prop.het) less than 0.1 and greater than 0.9 as homozygous (as a proxy for uniparental) expression
#for heterozygous (biparental) expression- allelic bias (prop.het) between 0.2-0.8
#I basically want to determine whether sons show evidence of expressing paternally inherited alleles 
#so basically comparing whether expression patterns look similar between mothers and sons- but any heterozygous expression in sons is interesting since very little in mealybugs



#### Family 1 35672F4
#cech.3572F4.f.snp
hom.cech.3572F4.f.snp1<-cech.3572F4.f.snp[cech.3572F4.f.snp$prop.het<0.1,]
hom.cech.3572F4.f.snp2<-cech.3572F4.f.snp[cech.3572F4.f.snp$prop.het>0.9,]
hom.cech.3572F4.f.snp<-rbind(hom.cech.3572F4.f.snp1,hom.cech.3572F4.f.snp2)
hom.cech.3572F4.f.snp.num<-length(hom.cech.3572F4.f.snp$contig)#3939

het.cech.3572F4.f.snp<-cech.3572F4.f.snp[cech.3572F4.f.snp$prop.het>0.2,]
het.cech.3572F4.f.snp<-het.cech.3572F4.f.snp[het.cech.3572F4.f.snp$prop.het<0.8,]#
het.cech.3572F4.f.snp.num<-length(het.cech.3572F4.f.snp$contig)# 6870

#exp props
prop13572F4<-hom.cech.3572F4.f.snp.num/(hom.cech.3572F4.f.snp.num+het.cech.3572F4.f.snp.num)# 0.3644185
prop23572F4<-het.cech.3572F4.f.snp.num/(hom.cech.3572F4.f.snp.num+het.cech.3572F4.f.snp.num)# 0.6355815

#cech.3572F4.m2.snp
hom.cech.3572F4.m2.snp1<-cech.3572F4.m2.snp[cech.3572F4.m2.snp$prop.het<0.1,]
hom.cech.3572F4.m2.snp2<-cech.3572F4.m2.snp[cech.3572F4.m2.snp$prop.het>0.9,]
hom.cech.3572F4.m2.snp<-rbind(hom.cech.3572F4.m2.snp1,hom.cech.3572F4.m2.snp2)
hom.cech.3572F4.m2.snp.num<-length(hom.cech.3572F4.m2.snp$contig)#4792

het.cech.3572F4.m2.snp<-cech.3572F4.m2.snp[cech.3572F4.m2.snp$prop.het>0.2,]
het.cech.3572F4.m2.snp<-het.cech.3572F4.m2.snp[het.cech.3572F4.m2.snp$prop.het<0.8,]#
het.cech.3572F4.m2.snp.num<-length(het.cech.3572F4.m2.snp$contig)# 8046

exp<-c(hom.cech.3572F4.m2.snp.num,het.cech.3572F4.m2.snp.num)
obs<-c(prop13572F4,prop23572F4)
cech3572F4.m2.test<-chisq.test(x=exp,p=obs)

#cech.3572F4.m3.snp
hom.cech.3572F4.m3.snp1<-cech.3572F4.m3.snp[cech.3572F4.m3.snp$prop.het<0.1,]
hom.cech.3572F4.m3.snp2<-cech.3572F4.m3.snp[cech.3572F4.m3.snp$prop.het>0.9,]
hom.cech.3572F4.m3.snp<-rbind(hom.cech.3572F4.m3.snp1,hom.cech.3572F4.m3.snp2)
hom.cech.3572F4.m3.snp.num<-length(hom.cech.3572F4.m3.snp$contig)# 4095

het.cech.3572F4.m3.snp<-cech.3572F4.m3.snp[cech.3572F4.m3.snp$prop.het>0.2,]
het.cech.3572F4.m3.snp<-het.cech.3572F4.m3.snp[het.cech.3572F4.m3.snp$prop.het<0.8,]#
het.cech.3572F4.m3.snp.num<-length(het.cech.3572F4.m3.snp$contig)# 8073

exp<-c(hom.cech.3572F4.m3.snp.num,het.cech.3572F4.m3.snp.num)
obs<-c(prop13572F4,prop23572F4)
cech.3572F4.m3.test<-chisq.test(x=exp,p=obs)

f3572F4count<-c(het.cech.3572F4.f.snp.num, hom.cech.3572F4.f.snp.num)
m2.3572F4count<-c(het.cech.3572F4.m2.snp.num, hom.cech.3572F4.m2.snp.num)
m3.3572F4count<-c(het.cech.3572F4.m3.snp.num, hom.cech.3572F4.m3.snp.num)

f.m2.3572F4.comp<-cbind(f3572F4count,m2.3572F4count)
fisher.test(f.m.3571F6.comp)
f.m3.3572F4.comp<-cbind(f3572F4count,m3.3572F4count)
fisher.test(f.m2.3571F6.comp)



#### Family2   3571F6
#cech.3571F6.f.snp
hom.cech.3571F6.f.snp1<-cech.3571F6.f.snp[cech.3571F6.f.snp$prop.het<0.1,]
hom.cech.3571F6.f.snp2<-cech.3571F6.f.snp[cech.3571F6.f.snp$prop.het>0.9,]
hom.cech.3571F6.f.snp<-rbind(hom.cech.3571F6.f.snp1,hom.cech.3571F6.f.snp2)
hom.cech.3571F6.f.snp.num<-length(hom.cech.3571F6.f.snp$contig)#3596

het.cech.3571F6.f.snp<-cech.3571F6.f.snp[cech.3571F6.f.snp$prop.het>0.2,]
het.cech.3571F6.f.snp<-het.cech.3571F6.f.snp[het.cech.3571F6.f.snp$prop.het<0.8,]#
het.cech.3571F6.f.snp.num<-length(het.cech.3571F6.f.snp$contig)# 6220

#exp props
prop13571F6<-hom.cech.3571F6.f.snp.num/(hom.cech.3571F6.f.snp.num+het.cech.3571F6.f.snp.num)   #0.3663407
prop23571F6<-het.cech.3571F6.f.snp.num/(hom.cech.3571F6.f.snp.num+het.cech.3571F6.f.snp.num)  #0.6336593


#cech.3571F6.m.snp
hom.cech.3571F6.m.snp1<-cech.3571F6.m.snp[cech.3571F6.m.snp$prop.het<0.1,]
hom.cech.3571F6.m.snp2<-cech.3571F6.m.snp[cech.3571F6.m.snp$prop.het>0.9,]
hom.cech.3571F6.m.snp<-rbind(hom.cech.3571F6.m.snp1,hom.cech.3571F6.m.snp2)
hom.cech.3571F6.m.snp.num<-length(hom.cech.3571F6.m.snp$contig)#4159

het.cech.3571F6.m.snp<-cech.3571F6.m.snp[cech.3571F6.m.snp$prop.het>0.2,]
het.cech.3571F6.m.snp<-het.cech.3571F6.m.snp[het.cech.3571F6.m.snp$prop.het<0.8,]
het.cech.3571F6.m.snp.num<-length(het.cech.3571F6.m.snp$contig)#  5989

exp<-c(hom.cech.3571F6.m.snp.num,het.cech.3571F6.m.snp.num)
obs<-c(prop13571F6,prop23571F6)
cech.3571F6.m.test<-chisq.test(x=exp,p=obs)

#cech.3571F6.m2.snp
hom.cech.3571F6.m2.snp1<-cech.3571F6.m2.snp[cech.3571F6.m2.snp$prop.het<0.1,]
hom.cech.3571F6.m2.snp2<-cech.3571F6.m2.snp[cech.3571F6.m2.snp$prop.het>0.9,]
hom.cech.3571F6.m2.snp<-rbind(hom.cech.3571F6.m2.snp1,hom.cech.3571F6.m2.snp2)
hom.cech.3571F6.m2.snp.num<-length(hom.cech.3571F6.m2.snp$contig)#

het.cech.3571F6.m2.snp<-cech.3571F6.m2.snp[cech.3571F6.m2.snp$prop.het>0.2,]
het.cech.3571F6.m2.snp<-het.cech.3571F6.m2.snp[het.cech.3571F6.m2.snp$prop.het<0.8,]#
het.cech.3571F6.m2.snp.num<-length(het.cech.3571F6.m2.snp$contig)# 

exp<-c(hom.cech.3571F6.m2.snp.num,het.cech.3571F6.m2.snp.num)
obs<-c(prop13571F6,prop23571F6)
cech.3571F6.m2.test<-chisq.test(x=exp,p=obs)

f3571F6count<-c(het.cech.3571F6.f.snp.num, hom.cech.3571F6.f.snp.num)
m.3571F6count<-c(het.cech.3571F6.m.snp.num, hom.cech.3571F6.m.snp.num)
m2.3571F6count<-c(het.cech.3571F6.m2.snp.num, hom.cech.3571F6.m2.snp.num)

f.m.3571F6.comp<-cbind(f3571F6count,m.3571F6count)
fisher.test(f.m.3571F6.comp)
f.m2.3571F6.comp<-cbind(f3571F6count,m2.3571F6count)
fisher.test(f.m2.3571F6.comp)


#### Family 3 3571F5
#cech.3571f5.f.snp
hom.cech.3571f5.f.snp1<-cech.3571f5.f.snp[cech.3571f5.f.snp$prop.het<0.1,]
hom.cech.3571f5.f.snp2<-cech.3571f5.f.snp[cech.3571f5.f.snp$prop.het>0.9,]
hom.cech.3571f5.f.snp<-rbind(hom.cech.3571f5.f.snp1,hom.cech.3571f5.f.snp2)
hom.cech.3571f5.f.snp.num<-length(hom.cech.3571f5.f.snp$contig)#4398

het.cech.3571f5.f.snp<-cech.3571f5.f.snp[cech.3571f5.f.snp$prop.het>0.2,]
het.cech.3571f5.f.snp<-het.cech.3571f5.f.snp[het.cech.3571f5.f.snp$prop.het<0.8,]#
het.cech.3571f5.f.snp.num<-length(het.cech.3571f5.f.snp$contig)# 7244

#exp props
prop13571f5<-hom.cech.3571f5.f.snp.num/(hom.cech.3571f5.f.snp.num+het.cech.3571f5.f.snp.num)#
prop23571f5<-het.cech.3571f5.f.snp.num/(hom.cech.3571f5.f.snp.num+het.cech.3571f5.f.snp.num)#

#cech.3571f5.m3.snp
hom.cech.3571f5.m3.snp1<-cech.3571f5.m3.snp[cech.3571f5.m3.snp$prop.het<0.1,]
hom.cech.3571f5.m3.snp2<-cech.3571f5.m3.snp[cech.3571f5.m3.snp$prop.het>0.9,]
hom.cech.3571f5.m3.snp<-rbind(hom.cech.3571f5.m3.snp1,hom.cech.3571f5.m3.snp2)
hom.cech.3571f5.m3.snp.num<-length(hom.cech.3571f5.m3.snp$contig)#3610

het.cech.3571f5.m3.snp<-cech.3571f5.m3.snp[cech.3571f5.m3.snp$prop.het>0.2,]
het.cech.3571f5.m3.snp<-het.cech.3571f5.m3.snp[het.cech.3571f5.m3.snp$prop.het<0.8,]#
het.cech.3571f5.m3.snp.num<-length(het.cech.3571f5.m3.snp$contig)# 21007

exp<-c(hom.cech.3571f5.m3.snp.num,het.cech.3571f5.m3.snp.num)
obs<-c(prop13571f5,prop23571f5)
cech.3571f5.m3.test<-chisq.test(x=exp,p=obs)

#cech.3571f5.m4.snp
hom.cech.3571f5.m4.snp1<-cech.3571f5.m4.snp[cech.3571f5.m4.snp$prop.het<0.1,]
hom.cech.3571f5.m4.snp2<-cech.3571f5.m4.snp[cech.3571f5.m4.snp$prop.het>0.9,]
hom.cech.3571f5.m4.snp<-rbind(hom.cech.3571f5.m4.snp1,hom.cech.3571f5.m4.snp2)
hom.cech.3571f5.m4.snp.num<-length(hom.cech.3571f5.m4.snp$contig)#

het.cech.3571f5.m4.snp<-cech.3571f5.m4.snp[cech.3571f5.m4.snp$prop.het>0.2,]
het.cech.3571f5.m4.snp<-het.cech.3571f5.m4.snp[het.cech.3571f5.m4.snp$prop.het<0.8,]#
het.cech.3571f5.m4.snp.num<-length(het.cech.3571f5.m4.snp$contig)# 

f3571f5count<-c(het.cech.3571f5.f.snp.num, hom.cech.3571f5.f.snp.num)
m3.3571f5count<-c(het.cech.3571f5.m3.snp.num, hom.cech.3571f5.m3.snp.num)
m4.3571f5count<-c(het.cech.3571f5.m4.snp.num, hom.cech.3571f5.m4.snp.num)

f.m3.3571f5.comp<-cbind(f3571f5count,m3.3571f5count)
fisher.test(f.m3.3571f5.comp)
f.m4.3571f5.comp<-cbind(f3571f5count,m4.3571f5count)
fisher.test(f.m4.3571f5.comp)

#### Family 4 3538- C. campanidorsalis
hom.3538.f.snp1<-ccamp.f.snp[ccamp.f.snp$prop.het<0.1,]
hom.3538.f.snp2<-ccamp.f.snp[ccamp.f.snp$prop.het>0.9,]
hom.3538.f.snp<-rbind(hom.3538.snp1,hom.3538.snp2)
hom.3538.f.snp.num<-length(hom.3538.snp$contig)#4961

het.3538.f.snp<-ccamp.f.snp[ccamp.f.snp$prop.het>0.2,]
het.3538.f.snp<-het.3538.f.snp[het.3538.f.snp$prop.het<0.8,]#
het.3538.f.snp.num<-length(het.3538.f.snp$contig)# 17492

#exp props
prop1<-4961/(4961+17492)#0.22
prop2<-17492/(4961+17492)#779

#ccamp.m1.snp
hom.ccamp.m1.snp1<-ccamp.m1.snp[ccamp.m1.snp$prop.het<0.1,]
hom.ccamp.m1.snp2<-ccamp.m1.snp[ccamp.m1.snp$prop.het>0.9,]
hom.ccamp.m1.snp<-rbind(hom.ccamp.m1.snp1,hom.ccamp.m1.snp2)
hom.ccamp.m1.snp.num<-length(hom.ccamp.m1.snp$contig)#3610

het.ccamp.m1.snp<-ccamp.m1.snp[ccamp.m1.snp$prop.het>0.2,]
het.ccamp.m1.snp<-het.ccamp.m1.snp[het.ccamp.m1.snp$prop.het<0.8,]#
het.ccamp.m1.snp.num<-length(het.ccamp.m1.snp$contig)# 21007

exp<-c(hom.ccamp.m1.snp.num,het.ccamp.m1.snp.num)
obs<-c(prop1,prop2)
ccamp.m1.snp.test<-chisq.test(x=exp,p=obs)
# they're different I guess

#ccamp.m2.snp
hom.ccamp.m2.snp1<-ccamp.m2.snp[ccamp.m2.snp$prop.het<0.1,]
hom.ccamp.m2.snp2<-ccamp.m2.snp[ccamp.m2.snp$prop.het>0.9,]
hom.ccamp.m2.snp<-rbind(hom.ccamp.m2.snp1,hom.ccamp.m2.snp2)
hom.ccamp.m2.snp.num<-length(hom.ccamp.m2.snp$contig)#

het.ccamp.m2.snp<-ccamp.m2.snp[ccamp.m2.snp$prop.het>0.2,]
het.ccamp.m2.snp<-het.ccamp.m2.snp[het.ccamp.m2.snp$prop.het<0.8,]#
het.ccamp.m2.snp.num<-length(het.ccamp.m2.snp$contig)# 

exp<-c(hom.ccamp.m2.snp.num,het.ccamp.m2.snp.num)
obs<-c(prop1,prop2)
ccamp.m2.snp.test<-chisq.test(x=exp,p=obs)
#

ccamp.m1.snp.test
ccamp.m2.snp.test


fcount<-c(het.3538.f.snp.num, hom.3538.f.snp.num)
m2.count<-c(het.ccamp.m2.snp.num, hom.ccamp.m2.snp.num)
m1.count<-c(het.ccamp.m1.snp.num, hom.ccamp.m1.snp.num)

f.m1.comp<-cbind(fcount,m1.count)
fisher.test(f.m1.comp)
f.m2.comp<-cbind(fcount,m2.count)
fisher.test(f.m2.comp)


##############################
# Plotting histograms of SNP heterozygosity for each species and sample

#histograms of prop.het for C. campanidorsalis- LGC3538
hist(ccamp.f.snp$prop.het, breaks = 100, col='#BF6B73', ylim=c(0,1000))
hist(ccamp.m1.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,1000))
hist(ccamp.m2.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,1000))

#histograms of prop.het for C. echiniformis- LGC3572F4
hist(cech.3572F4.f.snp$prop.het, breaks=100, col='#BF6B73', ylim=c(0,1000))
hist(cech.3572F4.m2.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,1000))
hist(cech.3572F4.m3.snp$prop.het, breaks = 100, col='#41719C', ylim=c(0,1000))

#histograms of prop.het for C. echiniformis- LGC3571F6
hist(cech.3571F6.m2.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,1000))
hist(cech.3571F6.m.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,1000))
hist(cech.3571F6.f.snp$prop.het, breaks = 100, col='#BF6B73', ylim=c(0,1000))

#histograms of prop.het for C. echiniformis- LGC3571F5
hist(cech.3571f5.m4.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,1000))
hist(cech.3571f5.m3.snp$prop.het, breaks=100, col='#41719C', ylim=c(0,1000))
hist(cech.3571f5.f.snp$prop.het, breaks = 100, col='#BF6B73', ylim=c(0,1000))
