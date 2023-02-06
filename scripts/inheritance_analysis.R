#setwd("/Users//christina//Dropbox//Cystococcus_experiments")
#setwd("/Users/christina/projects/cystococcus_snp/cystococcus_other/")
#setwd("/Users/christina/projects/cystococcus_snp/")
library(ggplot2)
#library(gridExtra)
library(nlme)
#library(ggeffects)
library(stargazer)
library(lmerTest)

# Summary of what this script is:
## basically we want to understand if allele inheritance patterns suggest that Cystococcus species have PGE.
## We have data from galls containing the mother and her offspring for two species of Cystococcus, one of which we thinks has PGE and one that we're not sure about.
## So, below, we plot the data showing whether in the broods of offspring for each family (which we have several of for each species), broods have different inheritance patterns for the paternally vs. maternally derived alleles.
## Then we compare the two species to see if the species that we think has PGE has different inheritance patterns from the one we're not sure about.


alleles<- read.csv("msat_allelebyfamily_copy.csv",header=TRUE, stringsAsFactors=TRUE)

pat.alleles<- alleles[alleles$origin == "pat", ]
mat.alleles<- alleles[alleles$origin == "mat", ]
cech.pat.alleles<- pat.alleles[pat.alleles$Species == "C.echiniformis", ]
ccamp.pat.alleles<- pat.alleles[pat.alleles$Species == "C.campanidorsalis", ]
cech.mat.alleles<- mat.alleles[mat.alleles$Species == "C.echiniformis", ]
ccamp.mat.alleles<- mat.alleles[mat.alleles$Species == "C.campanidorsalis", ]

cechpat<-ggplot(cech.pat.alleles, aes(fill=numalleles, y=numprimer, x=Familyid)) + 
  geom_bar(stat="identity") + theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "Family ID", y="count of primers")
ccamppat<-ggplot(ccamp.pat.alleles, aes(fill=numalleles, y=numprimer, x=Familyid)) + 
  geom_bar(stat="identity") + theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "Family ID", y="count of primers")
cechmat<-ggplot(cech.mat.alleles, aes(fill=numalleles, y=numprimer, x=Familyid)) + 
  geom_bar(stat="identity") + theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "Family ID", y="count of primers")
ccampmat<-ggplot(ccamp.mat.alleles, aes(fill=numalleles, y=numprimer, x=Familyid)) + 
  geom_bar(stat="identity") + theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "Family ID", y="count of primers")

quartz(height=6,width=7)
cechpat1<-ggplot(cech.pat.alleles, aes(fill=numalleles, y=numprimer, x=Familyid)) +scale_fill_manual(values=c("#A8CAE9", "#41719C","grey72")) + 
  geom_bar(stat="identity") + theme_classic(base_size = 18)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "", y="") + theme(legend.position = "none")+ scale_y_continuous(name="", breaks=c(0,2,4,6,8,10),limits=c(0,9))
ccamppat1<-ggplot(ccamp.pat.alleles, aes(fill=numalleles, y=numprimer, x=Familyid)) +scale_fill_manual(values=c("#A8CAE9", "#41719C","grey72")) + 
  geom_bar(stat="identity") + theme_classic(base_size = 18)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "", y="")+ theme(legend.position = "none")+ scale_y_continuous(name="", breaks=c(0,2,4,6,8),limits=c(0,8))
cechmat1<-ggplot(cech.mat.alleles, aes(fill=numalleles, y=numprimer, x=Familyid)) +scale_fill_manual(values=c("#FFC3C7", "#BF6B73","grey72")) + 
  geom_bar(stat="identity") + theme_classic(base_size = 18)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "", y="")+ theme(legend.position = "none")+ scale_y_continuous(name="", breaks=c(0,2,4,6,8,10),limits=c(0,9))
ccampmat1<-ggplot(ccamp.mat.alleles, aes(fill=numalleles, y=numprimer, x=Familyid)) +scale_fill_manual(values=c("#FFC3C7", "#BF6B73","grey72")) + 
  geom_bar(stat="identity") + theme_classic(base_size = 18)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "")+ theme(legend.position = "none")+ scale_y_continuous(name="", breaks=c(0,2,4,6,8),limits=c(0,8))

#+scale_fill_manual(values=c("#B6394F", "#C46F7E","grey72"))




primer<- read.csv("msat_allelebyprimer.csv",header=TRUE, stringsAsFactors=TRUE)

pat.primer<- primer[primer$origin == "pat", ]
mat.primer<- primer[primer$origin == "mat", ]
cech.pat.primer<- pat.primer[pat.primer$species == "C.echiniformis", ]
ccamp.pat.primer<- pat.primer[pat.primer$species == "C.campanidorsalis", ]
cech.mat.primer<- mat.primer[mat.primer$species == "C.echiniformis", ]
ccamp.mat.primer<- mat.primer[mat.primer$species == "C.campanidorsalis", ]

cechpatprim<-ggplot(cech.pat.primer, aes(fill=numalleles, y=count, x=Primer)) + 
  geom_bar(stat="identity") + theme_classic()
cechmatprim<-ggplot(cech.mat.primer, aes(fill=numalleles, y=count, x=Primer)) + 
  geom_bar(stat="identity") + theme_classic()
ccamppatprim<-ggplot(ccamp.pat.primer, aes(fill=numalleles, y=count, x=Primer)) + 
  geom_bar(stat="identity") + theme_classic()
ccampmatprim<-ggplot(ccamp.mat.primer, aes(fill=numalleles, y=count, x=Primer)) + 
  geom_bar(stat="identity") + theme_classic()
]
#+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x = "Family ID", y="count of primers")

###all works, does the same thing as above code but summarizing the number of alleles that are paternally or maternally derived for each primer

##things left to do:
###plot comparing results for each species (with some sort of confidence interval)
alleles<- read.csv("msat_allelebyfamily.csv",header=TRUE, stringsAsFactors=TRUE)
pat.alleles<- alleles[alleles$origin == "pat", ]
mat.alleles<- alleles[alleles$origin == "mat", ]
cech.pat.alleles<- pat.alleles[pat.alleles$Species == "C.echiniformis", ]
ccamp.pat.alleles<- pat.alleles[pat.alleles$Species == "C.campanidorsalis", ]
cech.mat.alleles<- mat.alleles[mat.alleles$Species == "C.echiniformis", ]
ccamp.mat.alleles<- mat.alleles[mat.alleles$Species == "C.campanidorsalis", ]

#getting means for each species and sex and allele num

Cechpat<-aggregate(cech.pat.alleles, by=list(cech.pat.alleles$numalleles), FUN=mean)
Cechpat<- Cechpat[c(1,6)]
Cechpatlen<-aggregate(cech.pat.alleles, by=list(cech.pat.alleles$numalleles), FUN=length)
Cechpatlen<-Cechpatlen[c(1,6)]
Cechpatsd<-aggregate(cech.pat.alleles, by=list(cech.pat.alleles$numalleles), FUN=sd)
Cechpatsd<-Cechpatsd[c(1,6)]
CI<-1.96*Cechpatsd$numprimer/sqrt(Cechpatlen$numprimer)
Species1<-c('Cech','Cech','Cech')
parentorg<-c('pat','pat','pat')
Cechpat<-cbind(Cechpat,CI,Species1,parentorg)
#this worked

Cechmat<-aggregate(cech.mat.alleles, by=list(cech.mat.alleles$numalleles), FUN=mean)
Cechmat<- Cechmat[c(1,6)]
Cechmatsd<-aggregate(cech.mat.alleles, by=list(cech.mat.alleles$numalleles), FUN=sd)
Cechmatsd<- Cechmatsd[c(1,6)]
Cechmatlen<-aggregate(cech.mat.alleles, by=list(cech.mat.alleles$numalleles), FUN=length)
Cechmatlen<- Cechmatlen[c(1,6)]
CI<-1.96*Cechmatsd$numprimer/sqrt(Cechmatlen$numprimer)
Species1<-c('Cech','Cech','Cech')
parentorg<-c('mat','mat','mat')
Cechmat<-cbind(Cechmat,CI,Species1,parentorg)

Ccamppat<-aggregate(ccamp.pat.alleles, by=list(ccamp.pat.alleles$numalleles), FUN=mean)
Ccamppat<- Ccamppat[c(1,6)]
Ccamppatsd<-aggregate(ccamp.pat.alleles, by=list(ccamp.pat.alleles$numalleles), FUN=sd)
Ccamppatsd<- Ccamppatsd[c(1,6)]
Ccamppatlen<-aggregate(ccamp.pat.alleles, by=list(ccamp.pat.alleles$numalleles), FUN=length)
Ccamppatlen<- Ccamppatlen[c(1,6)]
CI<-1.96*Ccamppatsd$numprimer/sqrt(Ccamppatlen$numprimer)
Species1<-c('Ccamp','Ccamp','Ccamp')
parentorg<-c('pat','pat','pat')
Ccamppat<-cbind(Ccamppat,CI,Species1,parentorg)

Ccampmat<-aggregate(ccamp.mat.alleles, by=list(ccamp.mat.alleles$numalleles), FUN=mean)
Ccampmat<-Ccampmat[c(1,6)]
Ccampmatsd<-aggregate(ccamp.mat.alleles, by=list(ccamp.mat.alleles$numalleles), FUN=sd)
Ccampmatsd<-Ccampmatsd[c(1,6)]
Ccampmatlen<-aggregate(ccamp.mat.alleles, by=list(ccamp.mat.alleles$numalleles), FUN=length)
Ccampmatlen<-Ccampmatlen[c(1,6)]
CI<-1.96*Ccampmatsd$numprimer/sqrt(Ccampmatlen$numprimer)
Species1<-c('Ccamp','Ccamp','Ccamp')
parentorg<-c('mat','mat','mat')
Ccampmat<-cbind(Ccampmat,CI,Species1,parentorg)

speciesbyparent<-rbind(Cechpat,Cechmat,Ccampmat,Ccamppat)
#worked, dataframe of means of number of alleles coming from each parent for each species

cechbyparent<-rbind(Cechpat,Cechmat)
ccampbyparent<-rbind(Ccampmat,Ccamppat)

allelecech<-ggplot(cechbyparent, aes(fill=Group.1, y=numprimer, x=parentorg)) + 
  geom_bar(stat="identity") + theme_classic()
allelesccamp<-ggplot(ccampbyparent, aes(fill=Group.1, y=numprimer, x=parentorg)) + 
  geom_bar(stat="identity") + theme_classic()

allelecech
allelesccamp
#plotted the cumulative amount not by parent of origin, need to get means and separate by parent


###plotting genotypes for multiple allele families
multmat<- read.csv("mult_mating_msat.csv",header=TRUE, stringsAsFactors=TRUE)

mult1<- multmat[multmat$Family == "LGC_2450", ]
mult2<- multmat[multmat$Family == "LGC_0628", ]
mult3<- multmat[multmat$Family == "LGC_1267", ]
mult4<- multmat[multmat$Family == "LGC_2525", ]

mult1_plot<-ggplot(mult1, aes(y=Offspring, x=ComboID)) + 
  geom_bar(stat="identity", fill=c("#0072B2")) + theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title="LGC_2450 C. campanidorsalis",x = "Allele", y="Offspring Count")
mult2_plot<-ggplot(mult2, aes(y=Offspring, x=ComboID)) + 
  geom_bar(stat="identity",fill=c("#D55E00")) + theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title="LGC_0628 C. echiniformis",x = "Allele", y="Offspring Count")
mult3_plot<-ggplot(mult3, aes(y=Offspring, x=ComboID)) + 
  geom_bar(stat="identity",fill=c("#D55E00")) + theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title="LGC_1267 C. echiniformis",x = "Allele", y="Offspring Count")
mult4_plot<-ggplot(mult4, aes(y=Offspring, x=ComboID)) + 
  geom_bar(stat="identity",fill=c("#D55E00")) + theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title="LGC_2525 C. echiniformis",x = "Allele", y="Offspring Count")
##above worked well enough, basically they just first import the data showing the allele combinations in potential multiple mating cases, then it subsets the data by family, then plots the number of offspring for each allele combination


###code from prev
#p<-ggplot(FAweek,aes(Week,mean, colour=Type))+geom_point()+scale_color_manual(values=c( "#9999CC", "#66CC99","#CC6666"))+ylim(0,8)
#p<-p+theme_classic()+geom_errorbar(aes(ymin=mean-CI,ymax=mean+CI,width=0.2))+aes(xlab="",ylab="")


## Test of whether number of alleles from mother and father differs between Cech and Ccamp
#data sets
msat_data<- read.delim("Msat_inheritance.tsv",header=T, stringsAsFactors = F)
msat_data$observation<-seq(1:38)
model_msat<-glmer(cbind(one_allelle,two_alleles)~species+parent+ (1|observation), data=msat_data, family=binomial)
summary(model_msat)

png('out1.png')
stargazer(model_msat, type = "text", digits = 3, star.cutoffs = c(0.05, 0.01, 0.001), digit.separator = "")
dev.off()
# number of primers with one vs two alleles inherited differs between parent but not between species



# Chi square of whether two pat alleles in broods is likely due to multiple mating
multmat<- read.csv("mult_mating_msat.csv",header=TRUE, stringsAsFactors=TRUE)

mult1<- multmat[multmat$Family == "LGC_2450", ]
mult2<- multmat[multmat$Family == "LGC_0628", ]
mult3<- multmat[multmat$Family == "LGC_1267", ]
mult4<- multmat[multmat$Family == "LGC_2525", ]

mult1_plot<-ggplot(mult1, aes(y=Offspring, x=ComboID)) + 
  geom_bar(stat="identity", fill=c("#0072B2")) + theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title="LGC_2450 C. campanidorsalis",x = "Allele", y="Offspring Count")
mult2_plot<-ggplot(mult2, aes(y=Offspring, x=ComboID)) + 
  geom_bar(stat="identity",fill=c("#D55E00")) + theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title="LGC_0628 C. echiniformis",x = "Allele", y="Offspring Count")
mult3_plot<-ggplot(mult3, aes(y=Offspring, x=ComboID)) + 
  geom_bar(stat="identity",fill=c("#D55E00")) + theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title="LGC_1267 C. echiniformis",x = "Allele", y="Offspring Count")
mult4_plot<-ggplot(mult4, aes(y=Offspring, x=ComboID)) + 
  geom_bar(stat="identity",fill=c("#D55E00")) + theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(title="LGC_2525 C. echiniformis",x = "Allele", y="Offspring Count")

#example code
#observed =      # no. observed
#expected =     # expected proportions

#chisq.test(x = observed,
#           p = expected) 

# Family data
# family 2450
obs1=c(0,7,8,1)
exp1=c(0.25,0.25,0.25,0.25)
chisq.test(x=obs1, p=exp1)

# family 628
obs2=c(0,11,9,0)
exp2=c(0.25,0.25,0.25,0.25)
chisq.test(x=obs2, p=exp2)

# family 1267
obs3=c(7,0,0,0,0,0,0,13)
exp3=c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
chisq.test(x=obs3, p=exp3)

# family 2525
obs4=c(0,0,0,0,0.5,0.5,12,0,0,5,0,0,0,0,0,0)
exp4=c(0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625,0.0625)
chisq.test(x=obs4, p=exp4)
