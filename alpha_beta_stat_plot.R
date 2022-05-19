library(qiime2R)
library(phyloseq)
library(vegan)
library(microbiome)
library(dunn.test)
library(btools) #to calculate Faith's PD
library(ggplot2)
library(plyr)
require(tidyverse)
require(dplyr)
require(lubridate)
require(ggpubr)
require(openxlsx)
require(cluster)
require(ecodist)

#Generate physeq
physeq <- qza_to_phyloseq(features = "microbiome/table2.qza",
                          tree="microbiome/rooted-tree.qza",
                          taxonomy="microbiome/taxonomy.qza",
                          metadata = "microbiome/metadata.txt")
##metadata to dataframe
meta<- sample_data(physeq)%>% data.frame
meta$SampleID <- rownames(meta)

#Plotting for rarefaction curve
rareplot = ranacapa::ggrare(physeq, se= FALSE)+theme_classic()

#Alpha diversity indices
diversity = estimate_richness(physeq,measures = c("Observed","Shannon"))
pd = estimate_pd(physeq) #Faith's PD 
#Make new file for alpha diversity data
diversity = merge(diversity,pd,by="row.names")
diversity = left_join(diversity,meta,by=c("Row.names"="SampleID")) %>% rename("SampleID"="Row.names")

#Statistics & plots for alpha diversity
##Observed richness
dunn.test (diversity$Observed, g=diversity$Season, method="bonferroni", kw=TRUE, label=TRUE,
           wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

p1<-ggplot(data=diversity,aes(x=Season, y= Observed,color=Season))+
  geom_boxplot(lwd=1.5)+
  ggtitle("(a) Observed richness")+
  theme_classic()+
  scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))+
  theme(axis.text.x = element_blank())+labs(x= "")
p1

##Shannon 
dunn.test (diversity$Shannon, g=diversity$Season, method="bonferroni", kw=TRUE, label=TRUE,
           wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

p2<-ggplot(data=diversity,aes(x=Season, y= Shannon,color=Season))+
  geom_boxplot(lwd=1.5)+
  ggtitle("(b) Shannon's index")+
  theme_classic()+
  scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))+
  theme(axis.text.x = element_blank())+labs(x= "")
p2

##Faith's PD
dunn.test (diversity$PD, g=diversity$Season, method="bonferroni", kw=TRUE, label=TRUE,
           wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

p3<-ggplot(data=diversity,aes(x=Season, y= PD,color=Season))+
  geom_boxplot(lwd=1.5)+
  ggtitle("(c) Faith's PD")+
  theme_classic()+
  scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p3

ggarrange(p1,p2,p3,ncol=1,legend = "none")


#Plots for beta diversity
##NMDS by Bray-curtis
bray = ordinate(physeq,method="NMDS",distance="bray")
p1.1<- plot_ordination(physeq,bray,color="Season",shape="Indiv")+
  theme_classic()+geom_point(size=3.5)+
  ggtitle("(a) Bray-Curtis")+theme(text = element_text(size = 15))+
  scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))+ 
  stat_ellipse(type = "t",level=0.95,aes(color=Season,group=Season)) 
p1.2<- plot_ordination(physeq,bray,color="Indiv",shape="Season")+
  theme_classic()+geom_point(size=3.5)+ggtitle("(a) Bray-Curtis")+
  theme(text = element_text(size = 15))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  stat_ellipse(type = "t",level=0.95,aes(color=Indiv,group=Indiv)) 

##PCoA plot by unweighted Unifrac
out.uf <- ordinate(physeq, "PCoA", "unifrac", weighted=FALSE)
p2.1<-plot_ordination(physeq,out.uf,color="Season",shape="Indiv")+theme_classic()+
  geom_point(size=3.5)+ggtitle("(b) Unweighted UniFrac")+
  theme(text = element_text(size = 15))+
  scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))+ 
  stat_ellipse(type = "t",level=0.95,aes(color=Season,group=Season)) 
p2.2<- plot_ordination(physeq,out.uf,color="Indiv",shape="Season")+
  theme_classic()+geom_point(size=3.5)+ggtitle("(b) Unweighted UniFrac")+
  theme(text = element_text(size = 15))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  stat_ellipse(type = "t",level=0.95,aes(color=Indiv,group=Indiv)) 


##PCoA plot by Weighted Unifrac
out.wuf <- ordinate(physeq, "PCoA", "unifrac", weighted=TRUE) 
p3.1<- plot_ordination(physeq,out.wuf,color="Season",shape="Indiv")+
  theme_classic()+geom_point(size=3.5)+ggtitle("(c) Weighted UniFrac")+
  theme(text = element_text(size = 15))+
  scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))+ 
  stat_ellipse(type = "t",level=0.95,aes(color=Season,group=Season)) 
p3.2<- plot_ordination(physeq,out.wuf,color="Indiv",shape="Season")+
  theme_classic()+geom_point(size=3.5)+ggtitle("(c) Weighted UniFrac")+
  theme(text = element_text(size = 15))+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  stat_ellipse(type = "t",level=0.95,aes(color=Indiv,group=Indiv)) 


###Beta diversity plots color by Season
ggarrange(p1.1,p2.1,p3.1, common.legend = TRUE,ncol=1,legend="right")

###Beta diversity plots color by Individual
ggarrange(p1.2,p2.2,p3.2, common.legend = TRUE,ncol=1,legend="right")

#Dispersion plots
##Dispersion plot of Bray curtis
d.bray = phyloseq::distance(physeq,"bray") #distance bray curtis
disp=betadisper(d.bray, sample_data(physeq)$Indiv)
boxplot(disp,border=c("#E69F00","#0072B2","#009E73"),lwd=4,main="SD_Bray",xlab=" ")

##Dispersion plot of unweighted UniFrac
unifrac_dist = phyloseq::distance(physeq, "unifrac")
disp=betadisper(unifrac_dist, sample_data(physeq)$Season)
boxplot(disp,border=c("#E69F00","#0072B2","#009E73"),lwd=4,main="SD_unweighted UniFrac",xlab=" ")

##Dispersion plot of weighted UniFrac
wunifrac_dist = phyloseq::distance(physeq, "wunifrac")
disp=betadisper(wunifrac_dist, sample_data(physeq)$Season)
boxplot(disp,border=c("#E69F00","#0072B2","#009E73"),lwd=4,main="SD_weighted UniFrac",xlab=" ")

#PERMANOVA for beta diversity
adonis2(unifrac_dist ~ sample_data(physeq)$Indiv*sample_data(physeq)$Seas) 
adonis2(wunifrac_dist~ sample_data(physeq)$Indiv*sample_data(physeq)$Seas) 
adonis2(d.bray~ sample_data(physeq)$Indiv*sample_data(physeq)$Seas) 

