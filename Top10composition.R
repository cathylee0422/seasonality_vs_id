#Composition barplots at different taxonomic rank (phylum, family and genus level)
library(ggplot2)
library(plyr)
require(tidyverse)
require(dplyr)
library(phyloseq)

# Set color for the plot
mycolors   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                "#44AA99","#332288", "#000000")

#Composition for Phylum
##Create data frame used for plotting (top 10 phylum + other)
tax_table(physeq)[tax_table(physeq)[, "Phylum"] == "", "Phylum"] <- "Unclassified phylum"
pseq.phylum <- physeq %>% transform_sample_counts(function(x) x/sum(x)) %>% tax_glom(taxrank="Phylum")
dat <- psmelt(pseq.phylum)
dat$Phylum <- as.character(dat$Phylum)
medians <- ddply(dat, ~Phylum, function(x) c(median=median(x$Abundance)))
medians = medians %>% arrange(desc(median))
top10.phylum <- c(medians$Phylum[1:10],"Other")

other.phylum <- medians$Phylum[11:length(medians$Phylum)]
dat[dat$Phylum %in% other.phylum,]$Phylum <- 'Other'
dat$Phylum <- factor(dat$Phylum,levels=c(top10.phylum)) #order phylum by relative abundance (most to fewest)

dat <- dat %>% group_by(Sample) %>%arrange(ymd(Date)) #sort samples by date
samp.ord <-dat$Sample
dat$Sample = factor(dat$Sample, unique(samp.ord))

##Plots by individual
###Relative abundance
p1<-dat %>% filter(Indiv=="Bo")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),
        axis.ticks = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Bo")+
  guides(fill=guide_legend(ncol=1))

###Add color coded bar annotation for season
leg <- dat %>% filter(Indiv=="Bo")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p1<-p1+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)

p1 


p2<-dat %>% filter(Indiv=="Em")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.ticks = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Em")

leg <- dat %>% filter(Indiv=="Em")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p2<-p2+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)

p2

p3<-dat %>% filter(Indiv=="Fl")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),
        axis.ticks = element_blank(),axis.title.y = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Fl")

leg <- dat %>% filter(Indiv=="Fl")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p3<-p3+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)

phylum<- ggarrange(p1,p2,p3,nrow=1,common.legend = TRUE,legend = "right")

#ggsave("barplot_phylum.pdf",width=8,height=11)


#Composition for genus
tax_table(physeq)[tax_table(physeq)[, "Genus"] == "", "Genus"] <- "Unclassified genus"
pseq.genus <- physeq %>% transform_sample_counts(function(x) x/sum(x)) %>% tax_glom(taxrank="Genus")
dat <- psmelt(pseq.genus)
dat$Genus<- as.character(dat$Genus)
medians <- ddply(dat, ~Genus, function(x) c(median=median(x$Abundance)))
medians = medians %>% arrange(desc(median))
top10.genus <- c(medians$Genus[1:10],"Other")

other.genus <- medians$Genus[11:length(medians$Genus)]
dat[dat$Genus %in% other.genus,]$Genus <- 'Other'
dat$Genus <- factor(dat$Genus,levels=c(top10.genus))

dat <- dat %>% group_by(Sample) %>%arrange(ymd(Date))
samp.ord <-dat$Sample
dat$Sample = factor(dat$Sample, unique(samp.ord))

#plot
p1<-dat %>% filter(Indiv=="Bo")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),
        axis.ticks = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Bo")+
  guides(fill=guide_legend(ncol=1))

leg <- dat %>% filter(Indiv=="Bo")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p1<-p1+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)
p1 

p2<-dat %>% filter(Indiv=="Em")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.ticks = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Em")

leg <- dat %>% filter(Indiv=="Em")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p2<-p2+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)

p2

p3<-dat %>% filter(Indiv=="Fl")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),
        axis.ticks = element_blank(),axis.title.y = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Fl")

leg <- dat %>% filter(Indiv=="Fl")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p3<-p3+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)

genus<-ggarrange(p1,p2,p3,nrow=1,common.legend = TRUE,legend = "right")

#ggsave("barplot_genus.pdf",width=8,height=11)

# Composition for family
tax_table(physeq)[tax_table(physeq)[, "Family"] == "", "Family"] <- "Unclassified family"
pseq.family <- physeq %>% transform_sample_counts(function(x) x/sum(x)) %>% tax_glom(taxrank="Family")
dat <- psmelt(pseq.family)
dat$Family <- as.character(dat$Family)
medians <- ddply(dat, ~Family, function(x) c(median=median(x$Abundance)))
medians = medians %>% arrange(desc(median))
top10.family <- c(medians$Family[1:10],"Other")

other.family <- medians$Family[11:length(medians$Family)]
dat[dat$Family %in% other.family,]$Family <- 'Other'

dat <- dat %>% group_by(Sample) %>%arrange(ymd(Date))
samp.ord <-dat$Sample
dat$Sample = factor(dat$Sample, unique(samp.ord))

#plot
p1<-dat %>% filter(Indiv=="Bo")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Family)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),
        axis.ticks = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Bo")+
  guides(fill=guide_legend(ncol=1))

leg <- dat %>% filter(Indiv=="Bo")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p1<-p1+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)

p1 

p2<-dat %>% filter(Indiv=="Em")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Family)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),axis.title.y = element_blank(),
        axis.ticks = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Em")

leg <- dat %>% filter(Indiv=="Em")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p2<-p2+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)

p2

p3<-dat %>% filter(Indiv=="Fl")%>%
  ggplot(aes(x=Sample, y=Abundance, fill=Family)) + 
  geom_bar(stat="identity",position="fill") +  
  scale_fill_manual(values=mycolors)+ 
  scale_color_manual(values=mycolors)+ theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1),axis.title.x = element_blank(),
        axis.ticks = element_blank(),axis.title.y = element_blank(),plot.margin = unit(c(1,1,1,1), "lines"))+ggtitle("Fl")

leg <- dat %>% filter(Indiv=="Fl")%>%ggplot(aes(x = Sample, y = 0)) + geom_point(aes(color = Season), shape = 15, size = 4, show.legend = F) + 
  theme_classic() + 
  theme(axis.title = element_blank(), axis.line = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(), 
        plot.margin = unit(c(0,0,0,0), "cm")) + scale_color_manual(values=c("#E69F00","#0072B2","#009E73"))
p3<-p3+ annotation_custom(ggplotGrob(leg),  ymin = -.2, ymax = 0, xmin = -Inf, xmax = Inf)

family <- ggarrange(p1,p2,p3,nrow=1,common.legend = TRUE,legend = "right")

#ggsave("barplot_family.pdf",width=8,height=11)

ggarrange(phylum,family,genus,nrow=3)

