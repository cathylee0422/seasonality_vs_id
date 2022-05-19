#Venn diagram on core microbiome between seasons
### code source https://microbiome.github.io/tutorials/core_venn.html

library(phyloseq)
library(microbiome)
library(tidyverse)
library(dplyr)

season <- unique(sample_data(physeq)$Season)
pseq.rel <- transform_sample_counts(physeq,function(x) x / sum(x))

list_core <- c() # an empty object to store information

for (n in season){ # for each variable n in season
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Season == n) # Choose sample from season by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 75% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each season.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

mycols <- adjustcolor(c('Fruit/Seed'="#E69F00", 'Leaf'="#009E73", 'Invertebrate'="#0072B2"),alpha.f=0.5) 
plot(eulerr::venn(list_core),
     fills = mycols,cex=10.5)

#List of core ASV and respective taxonomy in each season
print(list_core)

lf.core <- data.frame(tax_table(pseq.rel)) %>% 
  rownames_to_column("ASV")  %>% 
  filter(ASV %in% list_core$Leaf) %>% mutate(Season="Leaf")

fs.core <- data.frame(tax_table(pseq.rel)) %>% 
  rownames_to_column("ASV")  %>% 
  filter(ASV %in% list_core$`Fruit/Seed`) %>% mutate(Season="Fruit/Seed")

ot.core <- data.frame(tax_table(pseq.rel)) %>% 
  rownames_to_column("ASV")  %>% 
  filter(ASV %in% list_core$Invertebrate) %>% mutate(Season="Invertebrate")
