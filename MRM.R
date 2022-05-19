#MRM analysis
require(phyloseq)
require(tidyverse)
require(dplyr)
require(cluster)
require(ecodist)
require(openxlsx)

#calculate gower distance for individual
ind_df = meta %>%
  dplyr::select(Indiv) %>% 
  mutate(Indiv = as.factor(Indiv)) %>%
  as.data.frame 

ind_d = cluster::daisy(ind_df, metric='gower')
options(repr.plot.height=5, repr.plot.width=14)
plot(hclust(ind_d), cex=0.4)

#calculate euclidean distance using 4 major diet category (euclidean)
diet <- read.xlsx("monthlydwintake_prop_4categories.xlsx",detectDates=TRUE)
meta <- sample_data(physeq) %>% data.frame %>% left_join(diet,by = c("ym"="yearmonth")) #match by year and month
row.names(meta)<-meta$SampleID

#calculate euclidean distance based % dry weight intake of 4 major food items
diet_df = meta %>%
  dplyr::select(fs,lf,am,ot) %>%
  as.data.frame 

diet_d = cluster::daisy(diet_df, metric='euclidean')
options(repr.plot.height=5, repr.plot.width=14)
plot(hclust(diet_d), cex=0.4)


#calculate respecitive gower distance based on intakes of major food items
##fruit&seed
fs_df = meta %>%
  dplyr::select(fs) %>%
  as.data.frame 
fs_d = cluster::daisy(fs_df, metric='gower')

##leaf
lf_df = meta %>%
  dplyr::select(lf) %>%
  as.data.frame 
lf_d= cluster::daisy(lf_df, metric='gower')

##invertebrate
am_df = meta %>%
  dplyr::select(am) %>%
  as.data.frame 
am_d= cluster::daisy(am_df, metric='gower')

##other
ot_df = meta %>%
  dplyr::select(ot) %>%
  as.data.frame 
ot_d= cluster::daisy(ot_df, metric='gower')

#Caculate euclidean distance based on alpha diversity
row.names(diversity) <- diversity$SampleID
#Observed
ob_df <- diversity %>% 
  dplyr::select(Observed) %>% mutate(Observed=Observed %>% as.numeric) %>%
  as.data.frame()

ob_d = cluster::daisy(ob_df, metric='euclidean')
options(repr.plot.height=5, repr.plot.width=14)
plot(hclust(ob_d), cex=0.4)

#Shannon
sh_df <- diversity %>% 
  dplyr::select(Shannon) %>% mutate(Shannon=Shannon %>% as.numeric) %>%
  as.data.frame()

sh_d = cluster::daisy(sh_df, metric='euclidean')
options(repr.plot.height=5, repr.plot.width=14)
plot(hclust(sh_d), cex=0.4)

#PD
pd_df <- diversity %>% 
  dplyr::select(PD) %>% mutate(PD=PD %>% as.numeric) %>%
  as.data.frame()

pd_d = cluster::daisy(pd_df, metric='euclidean')
options(repr.plot.height=5, repr.plot.width=14)
plot(hclust(pd_d), cex=0.4)


#beta div were previously calculated and save as follows
d.wuf = phyloseq::distance(physeq, "wunifrac") #weighted UniFrac
d.uf = phyloseq::distance(physeq,"uunifrac") #unweighted UniFrac
d.bray = phyloseq::distance(physeq,"bray") #bray curtis

## checking if everything is overlapping, it should return "character(0)"
setdiff(ind_d %>% labels, d.wuf %>% labels) %>% print
setdiff(d.wuf  %>% labels, ind_d %>% labels) %>% print
setdiff(sample_names(physeq), d.wuf %>% labels) %>% print
setdiff(d.wuf %>% labels, sample_names(physeq)) %>% print

# Ordering the matrices in the same way , find function dist_mtx_order.R in github
X = labels(ob_d)
d.wuf_o = dist_mtx_order(d.wuf , X)
d.uf_o = dist_mtx_order(d.uf , X)
d.bray_o = dist_mtx_order(d.bray , X)
sh_d_o = dist_mtx_order(sh_d , X)
pd_d_o = dist_mtx_order(pd_d, X)
ind_d_o = dist_mtx_order(ind_d , X)
diet_d_o = dist_mtx_order(diet_d , X)
fs_d_o = dist_mtx_order(fs_d , X)
lf_d_o = dist_mtx_order(lf_d , X)
am_d_o = dist_mtx_order(am_d , X)
ot_d_o = dist_mtx_order(ot_d , X)

#MRM
##Observed richness
MRM(ob_d ~ ind_d_o + diet_d_o,mrank=TRUE,nperm=1000)
MRM(ob_d ~ ind_d_o + fs_d_o + lf_d_o + am_d_o + ot_d_o,mrank=TRUE,nperm=1000)

##Shannon
MRM(sh_d_o ~ ind_d_o + diet_d_o,mrank=TRUE,nperm=1000)
MRM(sh_d_o ~ ind_d_o + fs_d_o + lf_d_o + am_d_o + ot_d_o,mrank=TRUE,nperm=1000)

##Faith's PD
MRM(pd_d_o ~ ind_d_o + diet_d_o,mrank=TRUE,nperm=1000)
MRM(pd_d_o ~ ind_d_o + fs_d_o + lf_d_o + am_d_o + ot_d_o,mrank=TRUE,nperm=1000)

##Weighted UniFrac
MRM(d.wuf_o ~ ind_d_o + diet_d_o,mrank=TRUE,nperm=1000)
MRM(d.wuf_o ~ ind_d_o + fs_d_o + lf_d_o + am_d_o + ot_d_o,mrank=TRUE,nperm=1000)

##Unweighted UniFrac
MRM(d.uf_o ~ ind_d_o + diet_d_o,mrank=TRUE,nperm=1000)
MRM(d.uf_o ~ ind_d_o + fs_d_o + lf_d_o + am_d_o + ot_d_o,mrank=TRUE,nperm=1000)

##Bray
MRM(d.bray_o ~ ind_d_o + diet_d_o,mrank=TRUE,nperm=1000)
MRM(d.bray_o ~ ind_d_o + fs_d_o + lf_d_o + am_d_o + ot_d_o,mrank=TRUE,nperm=1000)

