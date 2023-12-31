
#title: "Hyena Microbiome"
#author: "Dominik Schmid"
#date: "24/11/2023"


####Introduction

#This scipt the analysis of 301 gut microbial samples from 12 hyenas
#The original data comes from Rojas et al. 2023 Msystem doi:



### This skript seeks to address whether gut microbial rhythms are present in hyena ###

##################################

##### clean up 
rm(list=ls())


### load packages 
library(vctrs)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(devtools)
library(CoDaSeq)
library(zCompositions)
library(tibble)
library(ggplot2)
library(dplyr)
library(vegan)
library(picante)
library(gridExtra)
library(pairwiseAdonis)
library(lme4)
library(lmerTest)
library(microbiome)
library(DT)
library(eulerr)
library(microbiomeutilities)
library(ggrepel)
library(ranacapa)
library(sjPlot)
library(sjmisc)
library(MuMIn)
library(jtools)
library(car)
library(plyr)

# Finally, set a random seed for reproducibility
set.seed(777)

###Upload ASV, taxonomy and meta data and compile into phyloseq object thats saved as RDS file (based on previous filtering)

ps = readRDS("~/Downloads/hyena_phyloseq.rds")

### lets remind ourselves of some basics about our phyloseq object 

microbiome::summarize_phyloseq(ps)

mean(sample_sums(ps)) ## on average we still have more than 30 000 reads from our samples 
sd(sample_sums(ps))

### turn row names into sample ids
sample_data(ps)$sample.ID<-row.names(sample_data(ps))

# remove two samples at midday and the NA

ps  <- ps %>% subset_samples(!sample.ID %in% c("P180"#, #NA
                                               #"P181", "P291" #afternoon
))

#format the time variable
#sample_data(ps)$sample_time<-as.POSIXct(sample_data(ps)$sample_time, format ="%H:%M:%S")
sample_data(ps)$sample_time<-as.POSIXct(paste(sample_data(ps)$sample_date, sample_data(ps)$sample_time), format ="%Y-%m-%d %H:%M:%S")

# make time category
sample_data(ps)$day_time<-ifelse(sample_data(ps)$sample_time <= as.POSIXct('2023-11-24 12:00'), "dawn", "dusk")

#time since sunrise
meta <- as(sample_data(ps), "data.frame")
light<-suncalc::getSunlightTimes(meta$sample_date, lat=1.41, lon=35.51, tz = "CET")
sample_data(ps)$sunrise<-light$sunrise
sample_data(ps)$goldenHour<-light$goldenHour
sample_data(ps)$hours_after_sunrise<-sample_data(ps)$sample_time-sample_data(ps)$sunrise
sample_data(ps)$hours_after_sunrise<-as.numeric(round(sample_data(ps)$hours_after_sunrise/60, 2))
sample_data(ps)$day_time<-ifelse(sample_data(ps)$hours_after_sunrise <= 8, "dawn", "dusk")
sample_data(ps)$time_frames<-ifelse(sample_data(ps)$hours_after_sunrise<=2, "2-hours", 
                                    ifelse(sample_data(ps)$hours_after_sunrise>2 & sample_data(ps)$hours_after_sunrise<=4, "4-hours", 
                                           ifelse(sample_data(ps)$hours_after_sunrise>4 & sample_data(ps)$hours_after_sunrise<=6, "6-hours", 
                                                  ifelse(sample_data(ps)$hours_after_sunrise>6 & sample_data(ps)$hours_after_sunrise<=12, "12-hours", 
                                                         ifelse(sample_data(ps)$hours_after_sunrise>12 & sample_data(ps)$hours_after_sunrise<=14, "14-hours", 
                                                                ifelse(sample_data(ps)$hours_after_sunrise>14 & sample_data(ps)$hours_after_sunrise<=16, "16-hours", "NA"))))))

#format matriline 
sample_data(ps)$matriline<-as.factor(sample_data(ps)$matriline)


#meta data
meta <- as(sample_data(ps), "data.frame")
str(meta)

######## now for the microbiome ########

ps_class<-tax_glom(ps, taxrank = "Class")

sample_data(ps_class)<-sample_data(ps_class)[,c("sample.ID")]
ps_class <- microbiome::transform(ps_class, "compositional") #### ps was not previously transformed - i.e. this is the first time - just to keep in mind 
class_per_sample<-psmelt(ps_class)

ddply(class_per_sample, c("Class"), summarise, mean=round(mean(Abundance),4))


ps_genus<-tax_glom(ps, taxrank = "Genus")

sample_data(ps_genus)<-sample_data(ps_genus)[,c("sample.ID")]
ps_genus <- microbiome::transform(ps_genus, "compositional") #### ps was not previously transformed - i.e. this is the first time - just to keep in mind 
genus_per_sample<-psmelt(ps_genus)

ddply(genus_per_sample, c("Genus"), summarise, mean=round(mean(Abundance),4))


#ps_family<-tax_glom(ps, taxrank = "Family")

#sample_data(ps_family)<-sample_data(ps_family)[,c("sample.ID")]
#ps_family <- microbiome::transform(ps_family, "compositional") #### ps was not previously transformed - i.e. this is the first time - just to keep in mind 
#family_per_sample<-psmelt(ps_family)

#ddply(family_per_sample, c("Family"), summarise, mean=round(mean(Abundance),4))

##### find most common classes

top<-as_tibble(get_group_abundances(ps_class, level="Class", group="Class", transform = "compositional") %>%arrange (-mean_abundance))
print(top,n=15)

##### classes that are more prominent than 1%
top_class<-c("Clostridia","Bacilli",  "Actinobacteria", "Bacteroidia", "Fusobacteriia", 
             "Erysipelotrichia", "Coriobacteriia","Gammaproteobacteria")

class_per_sample$class_plot<-as.character(ifelse(class_per_sample$Class %in% top_class, class_per_sample$Class, "Other"))

class_per_sample$class_plot<-factor(class_per_sample$class_plot, 
                                    levels = c("Clostridia","Bacilli",  "Actinobacteria", "Bacteroidia", "Fusobacteriia", 
                                               "Erysipelotrichia", "Coriobacteriia","Gammaproteobacteria", "Other"))

plot.meta<-meta[,c("sample.ID", "hours_after_sunrise")]

class_per_sample<-merge(class_per_sample, plot.meta, by="sample.ID")

class_plot_per_sample<-ggplot(class_per_sample, aes(y=Abundance, x=reorder(hours_after_sunrise, sample.ID), fill=class_plot))+
  geom_col(position = "fill", size=0, width=1)+
  scale_fill_manual(values=c("#DDAAD0","#D376BA","#CA3FA5","#E58890","#F2AC86", "#FFD07B", "#FEC457", "#FDB833", "grey65"))+
  ylab("Relative Abundance")+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank(), legend.position="bottom", axis.title.y = element_text(colour="black",size=12,face="bold"))+
  scale_y_continuous(expand = c(0,0))+ labs(fill ="Bacterial Class")

##### per genus
ps_genus<-tax_glom(ps, taxrank = "Genus")
core<-microbiome::core(ps_genus, detection = 1.5, prevalence = 0.85)
core.abundance <- sample_sums(core(ps_genus, detection = 2, prevalence = .85))/sample_sums(ps_genus)
mean(core.abundance)
sd(core.abundance)

top<-as_tibble(get_group_abundances(core, level="Genus", group="Genus", transform = "compositional") %>%arrange (-mean_abundance))
print(top)
#1 Lachnospiraceae          Lachnospiraceae                 0.247        0.173  
#2 Clostridium_ss1          Clostridium_ss1                 0.121        0.125  
#3 Clostridiales_unclass    Clostridiales_unclass           0.101        0.140  
#4 Peptoniphilus            Peptoniphilus                   0.0713       0.120  
#5 Peptoclostridium         Peptoclostridium                0.0636       0.0840 
#6 Fusobacterium            Fusobacterium                   0.0550       0.0953 
#7 Bacteroides              Bacteroides                     0.0450       0.0895 
#8 Erysipelotrichaceae      Erysipelotrichaceae             0.0414       0.0456 
#9 Paeniclostridium         Paeniclostridium                0.0404       0.0665 
#10 Clostridiaceae_1         Clostridiaceae_1                0.0354       0.0607 
#11 Alloprevotella           Alloprevotella                  0.0318       0.0713 
#12 Clostridium_ss7          Clostridium_ss7                 0.0308       0.0534 
#13 Streptococcus            Streptococcus                   0.0237       0.0581 
#14 Clostridiales_XIII       Clostridiales_XIII              0.0216       0.0243 
#15 Coriobacteriales_unclass Coriobacteriales_unclass        0.0198       0.0318 
#16 Slackia                  Slackia                         0.0166       0.0294 
#17 Peptostreptococcus       Peptostreptococcus              0.0123       0.0238 
#18 Peptococcus              Peptococcus                     0.00944      0.0159 
#19 Coprococcus_3            Coprococcus_3                   0.00894      0.0137 
#20 Eggerthellaceae          Eggerthellaceae                 0.00351      0.00492


top_genus<-c("Lachnospiraceae", "Clostridium_ss1", "Clostridiales_unclass", "Clostridium_ss7", "Clostridiaceae_1",  
             "Clostridiales_XIII", "Paeniclostridium", "Peptoclostridium", "Peptoniphilus",  "Peptostreptococcus",
             "Peptococcus", "Coprococcus_3", "Streptococcus", ### all Clostridia
             "Fusobacterium", 
             "Bacteroides", "Alloprevotella",
              "Erysipelotrichaceae",  
             "Coriobacteriales_unclass", "Slackia") ###"Eggerthellaceae", 

genus_per_sample$genus_plot<-as.character(ifelse(genus_per_sample$Genus %in% top_genus, genus_per_sample$Genus, "Other"))

genus_per_sample$genus_plot<-factor(genus_per_sample$genus_plot, 
                                    levels = c("Lachnospiraceae", "Clostridium_ss1", "Clostridiales_unclass", "Clostridium_ss7", "Clostridiaceae_1",  
                                               "Clostridiales_XIII", "Paeniclostridium", "Peptoclostridium", "Peptoniphilus",  "Peptostreptococcus",
                                               "Peptococcus", "Coprococcus_3", "Streptococcus", ### all Clostridia
                                               "Fusobacterium", 
                                               "Bacteroides", "Alloprevotella",
                                               "Erysipelotrichaceae",  
                                               "Coriobacteriales_unclass", "Slackia", #"Erysipelotrichaceae", #removed for vis because legend too big
                                               "Other"))

plot.meta<-meta[,c("sample.ID", "hours_after_sunrise")]

genus_per_sample<-merge(genus_per_sample, plot.meta, by="sample.ID")

genus_plot_per_sample<-ggplot(genus_per_sample, aes(y=Abundance, x=reorder(hours_after_sunrise, sample.ID), fill=genus_plot))+
  geom_col(position = "fill", size=0, width=1)+
  scale_fill_manual(values=c("#2F242C","#402350","#502274","#7A27A5","#A42CD6", "#D050FF","#D270E9","#D490D3","#D897B6","#EEB4B3","#F3C7C7","#E9D4BF", "#F0E849",
                             "#CECA82",
                             "#9DC488", "#85C18C",
                             "#6CBE8F",
                             "#3AB795", "#3D7564",
                             "grey65"))+
  ylab("Relative Abundance")+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), 
        strip.text = element_blank(), legend.position = "bottom", axis.title.y = element_text(colour="black",size=12,face="bold"))+
  scale_y_continuous(expand = c(0,0))+ labs(fill ="Bacterial Genus")


ggarrange(
  class_plot_per_sample,
  genus_plot_per_sample, 
  labels = c("A", "B"),
  nrow=2, ncol=1, align = c("h")
)


#####  alpha diversity metrices and some stats ######
#####################################################

# create alpha diversity metrices (choose choa1 over observed because almost identical)
# Explore phylogenetic tree and create PD diversity index

#order sample based on OTUs
#create tree plot and a merged tree plot
#calculate PD diversity index (library: picante)

#*this needs to be redone when we have actual tree infromation


#first we rarefy 

#total = 250
minimum = min(sample_sums(ps)) #find median sampling depth

standf = function(x, t=minimum) round(t * (x / sum(x))) #standardise by this sampling depth 

ps.rare = transform_sample_counts(ps, standf) #transform data accordingly 


otu.table <- as.data.frame(otu_table(ps.rare))
meta <- as(sample_data(ps.rare), "data.frame")
df.pd <- pd(t(otu.table), phy_tree(ps.rare), include.root=T) #what is SR? species richness?
meta$Phyogenetic_diversity <- df.pd$PD

alpha = estimate_richness(ps.rare, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")) #generate df with various diversity metices 

#add indicies to meta data (like we did with PD)
meta$Shannon <- alpha$Shannon 
#meta$Simpson <- alpha$Simpson 
#meta$InvSimpson <- alpha$InvSimpson
#meta$Chao1 <- alpha$Chao1 
meta$Observed <- alpha$Observed

##### modelling using gams
observed_gam <- mgcv::gam(log(Observed)~
                             s(age_yrs, bs="cr", k=10)+
                             matriline+
                             s(hours_after_sunrise, bs = "cr", k=6) + # specifically interested in this variable
                            #s(sample_year, bs="cr")+ 
                            s(hyenaID2, bs = "re"),
                           data=meta,
                           family = gaussian)

observed_gam <- mgcv::gam(log(Observed)~
                            s(age_yrs, bs="cr", k=10)+
                            matriline+
                            time_frames + # specifically interested in this variable
                            #s(sample_year, bs="cr")+ 
                            s(hyenaID2, bs = "re"),
                          data=meta,
                          family = gaussian)

print(summary(observed_gam)) # no treatment effect (great)
mgcv::gam.check(observed_gam)
plot(observed_gam)

library(mgcViz)
observed_plot<-ggplot(meta, aes(hours_after_sunrise, Observed, colour=hours_after_sunrise))+
  geom_point(size=2.5)+
  stat_smooth(method="gam", lty="dotted", colour="black")+
  xlab("hours since sunrise") + ylab("Observed ASV richness")+theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position = "bottom") +
  scale_colour_gradient(low="#C2E6AA",high="#E76E8E")+labs(fill = "")

observed_plot<-ggplot(meta, aes(fct_relevel(time_frames, "2-hours", "4-hours", "6-hours", "12-hours", "14-hours", "16-hours"), 
                                Observed, fill=time_frames))+
  geom_boxplot()+
  stat_smooth(method="gam", lty="dotted", colour="black")+
  xlab("hours since sunrise") + ylab("Observed ASV richness")+theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position = "bottom") +
  scale_fill_manual(values=c("#DA9B99","#DE8C95","#E76E8E","#C2E6AA","#C8D4A6","#CEC1A2"))+labs(fill = "")

##### modelling using gams
PD_gam <- mgcv::gam(log(Phyogenetic_diversity)~
                            s(age_yrs, bs="cr", k=4)+
                            matriline+
                            s(hours_after_sunrise, bs = "cr", k=6) + # specifically interested in this variable
                            #s(sample_year, bs="cr")+
                            s(hyenaID2, bs = "re"),
                          data=meta,
                          family = gaussian)

print(summary(PD_gam)) # no treatment effect (great)
mgcv::gam.check(PD_gam)
plot(PD_gam)


##### modelling using gams
shannon_gam <- mgcv::gam(log(Shannon)~
                           s(age_yrs, bs="cr", k=4)+
                           matriline+
                           s(hours_after_sunrise, bs = "cr", k=6) + # specifically interested in this variable
                           #s(sample_year, bs="cr")+
                           s(hyenaID2, bs = "re"),
                         data=meta,
                         family = gaussian)

print(summary(shannon_gam)) # no treatment effect (great)
mgcv::gam.check(shannon_gam)
plot(shannon_gam)


##### same on unrarefied data 

meta <- as(sample_data(ps), "data.frame")
df.pd <- pd(t(otu.table), phy_tree(ps), include.root=T) #what is SR? species richness?
meta$Phyogenetic_diversity <- df.pd$PD

alpha = estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")) #generate df with various diversity metices 

#add indicies to meta data (like we did with PD)
meta$Shannon <- alpha$Shannon 
#meta$Simpson <- alpha$Simpson 
#meta$InvSimpson <- alpha$InvSimpson
#meta$Chao1 <- alpha$Chao1 
meta$Observed <- alpha$Observed

##### modelling using gams
observed_gam <- mgcv::gam(log(Observed)~
                            s(age_yrs, bs="cr", k=4)+
                            matriline+
                            s(hours_after_sunrise, bs = "cr", k=6) + # specifically interested in this variable
                            #s(sample_year, bs="cr")+ 
                            s(hyenaID2, bs = "re"),
                          data=meta,
                          family = gaussian)

print(summary(observed_gam)) # no treatment effect (great)
mgcv::gam.check(observed_gam)
plot(observed_gam)

observed_gam <- mgcv::gam(log(Observed)~
                            s(age_yrs, bs="cr", k=4)+
                            matriline+
                            time_frames + # specifically interested in this variable
                            #s(sample_year, bs="cr")+ 
                            s(hyenaID2, bs = "re"),
                          data=meta,
                          family = gaussian)

print(summary(observed_gam)) # no treatment effect (great)
mgcv::gam.check(observed_gam)
plot(observed_gam)

######## beta diversity index ######
#ps.beta  <- ps.rare %>% #### unrarefied data produces the same outcome
#  subset_samples(!prey_abundance_avg_month %in% NA)

ps.beta<-ps.rare
meta.beta<- as(sample_data(ps.beta), "data.frame")
# calculate Bray-Curtis distance using the phyloseq package
#dist_bray = phyloseq::distance(ps, method="bray") ### abundance matters
# calculate Jaccard distance using the phyloseq package
#dist_jacc = phyloseq::distance(ps, method="jaccard") ### presence/absence
# calculate weighted Unifrag distance using the phyloseq package
dist_wUni = phyloseq::distance(ps.beta, method="weighted unifrac") ### abundance matters + phylogeny matters 
# calculate unweighted distance using the phyloseq package
dist_Uni = phyloseq::distance(ps.beta, method="unifrac") ### presence/absence + phylogeny matters


#### PERMANOVA (since in rojas et al. 2023 year and prey abundance had no effect we did not include it here)
adonis2(dist_Uni ~ hours_after_sunrise + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
        method="unifrac", by="margin")
# rarified results

#day_time    1    0.755 0.01424 4.3184  0.001 ***
#  matriline   3    1.564 0.02949 2.9808  0.321    
#age_yrs     1    0.327 0.00618 1.8728  0.592    
#Residual  288   50.356 0.94985                  
#Total     293   53.015 1.00000

#hours_after_sunrise   1    0.729 0.01359 4.1802  0.003 **
#  matriline             3    1.631 0.03041 3.1171  0.315   
#age_yrs               1    0.336 0.00626 1.9261  0.598   
#Residual            292   50.925 0.94965                 
#Total               297   53.625 1.00000 

#time_frames   5    1.493 0.02767 1.7145  0.003 **
#  matriline     3    1.593 0.02952 3.0494  0.247   
#age_yrs       1    0.359 0.00666 2.0639  0.482   
#Residual    290   50.499 0.93593                 
#Total       299   53.956 1.00000 

adonis2(dist_wUni ~ time_frames + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
        method="weighted unifrac", by="margin")
# rarified results

#day_time    1   0.1230 0.00989 2.9921  0.052 .
#matriline   3   0.3346 0.02690 2.7138  0.330  
#age_yrs     1   0.1408 0.01132 3.4272  0.254  
#Residual  288  11.8358 0.95153                
#Total     293  12.4387 1.00000 

#hours_after_sunrise   1   0.1194 0.00951 2.9219  0.063 .
#matriline             3   0.3564 0.02839 2.9084  0.438  
#age_yrs               1   0.1464 0.01167 3.5849  0.246  
#Residual            292  11.9283 0.95019                
#Total               297  12.5536 1.00000 

#time_frames   5   0.2505 0.01981 1.2223  0.321
#matriline     3   0.3457 0.02734 2.8121  0.266
#age_yrs       1   0.1618 0.01279 3.9471  0.138
#Residual    290  11.8842 0.93986              
#Total       299  12.6446 1.00000 

##### now for unrarefied data
#ps.beta  <- ps %>% #### unrarefied data produces the same outcome
#  subset_samples(!prey_abundance_avg_month %in% NA)

ps.beta.unrarefied<-ps
meta.beta.unrarefied<- as(sample_data(ps.beta.unrarefied), "data.frame")
# calculate Bray-Curtis distance using the phyloseq package
#dist_bray = phyloseq::distance(ps, method="bray") ### abundance matters
# calculate Jaccard distance using the phyloseq package
#dist_jacc = phyloseq::distance(ps, method="jaccard") ### presence/absence
# calculate weighted Unifrag distance using the phyloseq package
dist_wUni = phyloseq::distance(ps.beta.unrarefied, method="weighted unifrac") ### abundance matters + phylogeny matters 
# calculate unweighted distance using the phyloseq package
dist_Uni = phyloseq::distance(ps.beta.unrarefied, method="unifrac") ### presence/absence + phylogeny matters

#### PERMANOVA (since in rojas et al. 2023 year and prey abundance had no effect we did not include it here)
adonis2(dist_Uni ~ hours_after_sunrise + matriline + age_yrs, data = meta.beta.unrarefied, strata=meta.beta.unrarefied$hyenaID2, 
        method="unifrac", by="margin")
# UN-rarified results

#day_time    1    0.539 0.01102 3.3098  0.001 ***
#matriline   3    1.081 0.02210 2.2123  0.155    
#age_yrs     1    0.364 0.00743 2.2319  0.116    
#Residual  288   46.927 0.95895                  
#Total     293   48.936 1.00000

#hours_after_sunrise   1    0.532 0.01073 3.2679  0.002 **
#  matriline             3    1.117 0.02251 2.2853  0.081 . 
#age_yrs               1    0.381 0.00769 2.3412  0.094 . 
#Residual            292   47.560 0.95887                 
#Total               297   49.600 1.00000

#time_frames   5    1.280 0.02566 1.5766  0.001 ***
#  matriline     3    1.087 0.02180 2.2320  0.129    
#age_yrs       1    0.411 0.00825 2.5341  0.041 *  
#  Residual    290   47.090 0.94400                  
#Total       299   49.883 1.00000

adonis2(dist_wUni ~ hours_after_sunrise + matriline + age_yrs, data = meta.beta.unrarefied, strata=meta.beta.unrarefied$hyenaID2, 
        method="weighted unifrac", by="margin")
# UN-rarified results

#day_time    1   0.1213 0.00994 3.0085  0.054 .
#matriline   3   0.3305 0.02706 2.7318  0.311  
#age_yrs     1   0.1413 0.01157 3.5039  0.228  
#Residual  288  11.6154 0.95107                
#Total     293  12.2130 1.00000

#hours_after_sunrise   1   0.1178 0.00956 2.9381  0.060 .
#matriline             3   0.3520 0.02856 2.9268  0.383  
#age_yrs               1   0.1469 0.01192 3.6641  0.215  
#Residual            292  11.7070 0.94974                
#Total               297  12.3266 1.00000

#time_frames   5   0.2474 0.01992 1.2300  0.308
#matriline     3   0.3412 0.02748 2.8280  0.200
#age_yrs       1   0.1623 0.01307 4.0361  0.122
#Residual    290  11.6637 0.93933              
#Total       299  12.4170 1.00000


######## is it the same if we just consider core taxa?
#ps.beta  <- ps %>% #### unrarefied data produces the same outcome
#  subset_samples(!prey_abundance_avg_month %in% NA)
ps.beta.core<-ps
meta.beta.core<- as(sample_data(ps.beta.core), "data.frame")
ps_genus<-tax_glom(ps.beta.core, taxrank = "Genus") 
core<-microbiome::core(ps_genus, detection = 1, prevalence = 0.85)

top<-as_tibble(get_group_abundances(core, level="Genus", group="Genus", transform = "compositional") %>%arrange (-mean_abundance))
print(top,n=26)

#core = subset_taxa(core, Genus!="Peptococcus" & Genus!="Coprococcus_3" & Genus!="Eggerthellaceae")

dist_wUni = phyloseq::distance(core, method="weighted unifrac") ### abundance matters + phylogeny matters 
dist_Uni = phyloseq::distance(core, method="unifrac") ### presence/absence + phylogeny matters

#meta.beta.core$Seq_depth<-as.numeric(sample_sums(ps_genus)) #has no impact whatsoever
#meta.beta.core$Seq_depth_scaled <- as.numeric(scales::rescale(meta.beta.core$Seq_depth, to = c(-1,1)))

unweighted_core_PERM<-adonis2(dist_Uni ~ time_frames + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
        method="unifrac", by="margin")

#hours_after_sunrise   1  0.01994 0.01468 4.4209  0.013 *
#  matriline             3  0.01654 0.01217 1.2226  0.026 *
#  age_yrs               1  0.00756 0.00556 1.6752  0.137  
#Residual            292  1.31694 0.96928                
#Total               297  1.35868 1.00000

#time_frames   5  0.02725 0.01549 0.9325  0.487  
#matriline     3  0.02226 0.01265 1.2694  0.116  
#age_yrs       1  0.01514 0.00861 2.5903  0.049 *
#  Residual    290  1.69511 0.96359                
#Total       299  1.75916 1.00000

weighted_core_PERM<-adonis2(dist_wUni ~ time_frames + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
        method="weighted unifrac", by="margin")

#hours_after_sunrise   1   0.2057 0.01242 3.8349  0.031 *
#  matriline             3   0.5481 0.03310 3.4066  0.425  
#age_yrs               1   0.1639 0.00990 3.0565  0.456  
#Residual            292  15.6595 0.94559                
#Total               297  16.5606 1.00000  

#time_frames   5   0.4687 0.02364 1.4678  0.097 .
#matriline     3   0.6192 0.03124 3.2320  0.200  
#age_yrs       1   0.2456 0.01239 3.8460  0.252  
#Residual    290  18.5211 0.93430                
#Total       299  19.8234 1.00000 

####
ordination_wUni<-ordinate(core, method="PCoA", distance=dist_wUni)

weighted <- plot_ordination(core, ordination_wUni) + 
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = day_time)) + # put ellipses around group centroids
  geom_point(aes(
    #fill = day_time, 
    col = hours_after_sunrise), # fill color by infection status
    size = 3) +# make points size 4
    #pch = 21, # Make points circular with a border 
    #col = hours_after_sunrise) + # the border colour should be black
  #labs(subtitle = "a)") + # insert subtitle
  #ggtitle("Unweighted Unifrac") + # insert title
  theme_bw() + scale_fill_manual(values=c("#C2E6AA", "#E76E8E")) +
  labs(fill = "time of day", colour="hours since sunrise") + scale_colour_gradient(low="#C2E6AA",high="#E76E8E")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold"))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))# center align the plot title
weighted



ggarrange(
  observed_plot,
  weighted, 
  labels = c("A", "B"),
  nrow=1, ncol=2, align = c("h"), legend ="bottom", common.legend = T
)






##### GAMS 
#continuous variables typically s(var, ...)
#categories need to be factors 
#interactions coded as s(var, by = interaction_term)
#k = number of basic functions fitted to a line 
#REML or use sp = # to adjust smoothing 
#bs = 

#edf is effective degree of freedoms
#edf= 1 is a linear fit 
#edf= 2 is a quadratic term 
#edf> 2 the more wiggly the line becomes

#a significant smooth term is one where you cannot draw a horizontal line through the 95% confidence interval.

### prep for GAMMs

### we look specifically at the core taxa 
clr_scaled <-microbiome::transform(core, transform = "clr")

#### 
sample_data(clr_scaled)<-sample_data(clr_scaled)[,c( "sample.ID", "hyenaID2", "age_yrs", "matriline", "hours_after_sunrise", "time_frames", "day_time")]
core_genera_melt<-psmelt(clr_scaled)
head(core_genera_melt)
meta.gam<-as(sample_data(clr_scaled), "data.frame")
meta.gam$Seq_depth<-as.numeric(sample_sums(core))
meta.gam$Seq_depth_scaled <- as.numeric(scales::rescale(meta.gam$Seq_depth, to = c(-1,1)))
core_genera_gam<- merge(core_genera_melt, meta.gam, by.x = "Sample", by.y = "sample.ID")
names(core_genera_gam)
head(core_genera_gam)
str(core_genera_gam)
core_genera_gam$time_frames.x<-fct_relevel(core_genera_gam$time_frames.x, "2-hours", "4-hours", "6-hours", "12-hours", "14-hours", "16-hours")

#### GAMMS
library(mgcv)
library(tidygam)
library(tidymv)

######### Lachnospiraceae #########
Lachnospiraceae<-subset(core_genera_gam, Genus == "Lachnospiraceae")

Lachnospiraceae_gam <- mgcv::gam(Abundance~
                      s(age_yrs.x, bs="cr", k=4)+
                      matriline.x+
                      s(hours_after_sunrise.x, bs = "cr", k=6) + # specifically interested in this variable
                      #s(Seq_depth_scaled)+
                      s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                      data=Lachnospiraceae,
                      family = gaussian)

print(summary(Lachnospiraceae_gam)) # no treatment effect (great)
gam.check(Lachnospiraceae_gam)
plot(Lachnospiraceae_gam)

#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    2.7811     0.1361  20.434   <2e-16 ***
#  matriline.x2   0.3846     0.1902   2.022   0.0441 *  
#  matriline.x3   0.4728     0.1845   2.563   0.0109 *  
#  matriline.x4   0.1174     0.1919   0.612   0.5410    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
#  edf Ref.df     F p-value  
#s(age_yrs.x)                        1.000  1.000 0.021  0.8852  
#s(hours_after_sunrise.x)            1.509  1.864 0.380  0.6257  
#s(hyenaID2.x,hours_after_sunrise.x) 4.030 11.000 0.605  0.0938 

####### Clostridium_ss1 ########
Clostridium_ss1<-subset(core_genera_gam, Genus == "Clostridium_ss1")

Clostridium_ss1_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr", k=4)+
                                   matriline.x+
                                   s(hours_after_sunrise.x, bs = "cr", k=6) + # specifically interested in this variable
                                   #s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                 data=Clostridium_ss1,
                                 family = gaussian)

print(summary(Clostridium_ss1_gam)) #strong treatment effect (fewer in Staphylococcus in treatment)
gam.check(Clostridium_ss1_gam) # plus day effect (random)
plot(Clostridium_ss1_gam)

#matriline.x2 -0.23137    0.31654  -0.731    0.465    
#matriline.x3 -0.01733    0.29787  -0.058    0.954    
#matriline.x4  0.10522    0.31757   0.331    0.741    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
#  edf Ref.df     F  p-value    
#s(age_yrs.x)                        1.000  1.000 2.986   0.0851 .  
#s(hours_after_sunrise.x)            1.680  2.098 1.058   0.3194    
#s(hyenaID2.x,hours_after_sunrise.x) 7.483 11.000 3.439 1.21e-06

###########  Clostridiales_unclass    ##############
Clostridiales_unclass<-subset(core_genera_gam, Genus == "Clostridiales_unclass")

Clostridiales_unclass_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr", k=4)+
                                   matriline.x+
                                   s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                   #s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                 data=Clostridiales_unclass,
                                 family = gaussian)

print(summary(Clostridiales_unclass_gam)) # no effect of treatment but day (random effect)
gam.check(Clostridiales_unclass_gam)
plot(Clostridiales_unclass_gam)

###########  Clostridiaceae_1    ##############
Clostridiaceae_1<-subset(core_genera_gam, Genus == "Clostridiaceae_1")

Clostridiaceae_1_gam <- mgcv::gam(Abundance~
                                         s(age_yrs.x, bs="cr", k=4)+
                                         matriline.x+
                                         s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                         #s(Seq_depth_scaled)+
                                         s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                       data=Clostridiaceae_1,
                                       family = gaussian)

print(summary(Clostridiaceae_1_gam)) # some effect of treatment (more Corynebacterium in treatment)
gam.check(Clostridiaceae_1_gam)
plot(Clostridiaceae_1_gam)

#Estimate Std. Error t value Pr(>|t|)   
#(Intercept)    0.5356     0.3415   1.568  0.11797   
#matriline.x2  -1.2057     0.4561  -2.643  0.00866 **
#  matriline.x3  -0.5100     0.4353  -1.172  0.24230   
#matriline.x4  -0.1003     0.4620  -0.217  0.82830   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Approximate significance of smooth terms:
#  edf Ref.df     F  p-value    
#s(age_yrs.x)                        2.605  2.884 3.348   0.0115 *  
#  s(hours_after_sunrise.x)            1.000  1.000 3.322   0.0694 .  
#s(hyenaID2.x,hours_after_sunrise.x) 7.085 11.000 2.406 8.72e-05 ***

Clostridiaceae_1_gam <- mgcv::gamm(Abundance~
                                     s(age_yrs.x, bs="cr", k=4)+
                                     matriline.x+
                                     s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                     #s(Seq_depth_scaled)+
                                     s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                   data=Clostridiaceae_1)

sjPlot::tab_model(Clostridiaceae_1_gam$gam,  show.stat = TRUE, show.se = TRUE, title = "Table S3: GAMM predicting observed ASV richness. Note that the estimate for smooth terms is actually EDF (estimated degrees of freedom), it is not an 'estimate'. EDF reflects the degree of non-linearity of a curve (1 = linear)",
                  dv.labels = "Observed ASV richness",
                  pred.labels = c("Intercept", "matriline.x", "age_yrs.x (non-linear)", "hours_after_sunrise.x (non-linear)"))


###########  Clostridium_ss7    ##############
Clostridium_ss7<-subset(core_genera_gam, Genus == "Clostridium_ss7")

Clostridium_ss7_gam <- mgcv::gam(Abundance~
                                    s(age_yrs.x, bs="cr", k=4)+
                                    matriline.x+
                                    s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                    #s(Seq_depth_scaled)+
                                    s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                  data=Clostridium_ss7,
                                  family = gaussian)

print(summary(Clostridium_ss7_gam)) # no effect of treatment (good)
gam.check(Clostridium_ss7_gam)
plot(Clostridium_ss7_gam)


###########  Clostridiales_XIII    ##############
Clostridiales_XIII<-subset(core_genera_gam, Genus == "Clostridiales_XIII")

Clostridiales_XIII_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr", k=4)+
                                   matriline.x+
                                   s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                   #s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                 data=Clostridiales_XIII,
                                 family = gaussian)

print(summary(Clostridiales_XIII_gam)) ## no effect of treatment on abundance (good)
gam.check(Clostridiales_XIII_gam)
plot(Clostridiales_XIII_gam)

#s(age_yrs.x)                        2.095   2.49 2.812  0.0408 *  
#  s(hours_after_sunrise.x)            1.000   1.00 4.556  0.0337 *  
#  s(hyenaID2.x,hours_after_sunrise.x) 5.762  11.00 1.984 9.6e-05 ***

###########  Paeniclostridium    ##############
Paeniclostridium<-subset(core_genera_gam, Genus == "Paeniclostridium")

Paeniclostridium_gam <- mgcv::gam(Abundance~
                                      s(age_yrs.x, bs="cr", k=4)+
                                      matriline.x+
                                      s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                      #s(Seq_depth_scaled)+
                                      s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                    data=Paeniclostridium,
                                    family = gaussian)

print(summary(Paeniclostridium_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Paeniclostridium_gam) # and day effect (random)
plot(Paeniclostridium_gam)

############ Peptoclostridium ###########
Peptoclostridium<-subset(core_genera_gam, Genus == "Peptoclostridium")

Peptoclostridium_gam <- mgcv::gam(Abundance~
                                    s(age_yrs.x, bs="cr", k=6)+
                                    matriline.x+
                                    s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                    #s(Seq_depth_scaled)+
                                    s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                  data=Peptoclostridium,
                                  family = gaussian)

print(summary(Peptoclostridium_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Peptoclostridium_gam) # and day effect (random)
plot(Peptoclostridium_gam)

#s(age_yrs.x)                        4.06  4.592 5.033 0.000358 ***
#  s(hours_after_sunrise.x)            1.66  2.076 2.912 0.054415 .  
#s(hyenaID2.x,hours_after_sunrise.x) 7.93 11.000 3.099 9.36e-06 ***


Peptoclostridium_gam <- mgcv::gamm(Abundance~
                                    s(age_yrs.x, bs="cr", k=4)+
                                    matriline.x+
                                    s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                    #s(Seq_depth_scaled)+
                                    s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                  data=Peptoclostridium)

sjPlot::tab_model(Peptoclostridium_gam$gam,  show.stat = TRUE, show.se = TRUE, title = "Table S3: GAMM predicting observed ASV richness. Note that the estimate for smooth terms is actually EDF (estimated degrees of freedom), it is not an 'estimate'. EDF reflects the degree of non-linearity of a curve (1 = linear)",
                  dv.labels = "Observed ASV richness",
                  pred.labels = c("Intercept", "matriline.x", "age_yrs.x (non-linear)", "hours_after_sunrise.x (non-linear)"))


############ Peptoniphilus ###########
Peptoniphilus<-subset(core_genera_gam, Genus == "Peptoniphilus")

Peptoniphilus_gam <- mgcv::gam(Abundance~
                                    s(age_yrs.x, bs="cr", k=4)+
                                    matriline.x+
                                    s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                    #s(Seq_depth_scaled)+
                                    s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                  data=Peptoniphilus,
                                  family = gaussian)

print(summary(Peptoniphilus_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Peptoniphilus_gam) # and day effect (random)
plot(Peptoniphilus_gam)

############ Peptostreptococcus ###########
Peptostreptococcus<-subset(core_genera_gam, Genus == "Peptostreptococcus")

Peptostreptococcus_gam <- mgcv::gam(Abundance~
                                 s(age_yrs.x, bs="cr", k=4)+
                                 matriline.x+
                                 s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                 #s(Seq_depth_scaled)+
                                 s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                               data=Peptostreptococcus,
                               family = gaussian)

print(summary(Peptostreptococcus_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Peptostreptococcus_gam) # and day effect (random)
plot(Peptostreptococcus_gam)

############ Peptococcus ###########
Peptococcus<-subset(core_genera_gam, Genus == "Peptococcus")

Peptococcus_gam <- mgcv::gam(Abundance~
                                      s(age_yrs.x, bs="cr", k=4)+
                                      matriline.x+
                                      s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                      #s(Seq_depth_scaled)+
                                      s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                    data=Peptococcus,
                                    family = gaussian)

print(summary(Peptococcus_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Peptococcus_gam) # and day effect (random)
plot(Peptococcus_gam)

#s(age_yrs.x)                        2.351  2.715 4.227 0.00469 ** 
#  s(hours_after_sunrise.x)            1.459  1.792 4.104 0.01408 *  
#  s(hyenaID2.x,hours_after_sunrise.x) 8.435 11.000 6.301 < 2e-16

############ Coprococcus_3 ###########
Coprococcus_3<-subset(core_genera_gam, Genus == "Coprococcus_3")

Coprococcus_3_gam <- mgcv::gam(Abundance~
                               s(age_yrs.x, bs="cr", k=4)+
                               matriline.x+
                               s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                               #s(Seq_depth_scaled)+
                               s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                             data=Coprococcus_3,
                             family = gaussian)

print(summary(Coprococcus_3_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Coprococcus_3_gam) # and day effect (random)
plot(Coprococcus_3_gam)

#s(age_yrs.x)                        2.809   2.97  6.235 0.000423 ***
#  s(hours_after_sunrise.x)            1.000   1.00 11.028 0.001014 ** 
#  s(hyenaID2.x,hours_after_sunrise.x) 6.938  11.00  2.967 3.19e-06 ***
  
############ Streptococcus ###########
Streptococcus<-subset(core_genera_gam, Genus == "Streptococcus")

Streptococcus_3_gam <- mgcv::gam(Abundance~
                                 s(age_yrs.x, bs="cr", k=4)+
                                 matriline.x+
                                 s(hours_after_sunrise.x, bs = "cr", k=8) + # specifically interested in this variable
                                 #s(Seq_depth_scaled)+
                                 s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                               data=Streptococcus,
                               family = gaussian)

print(summary(Streptococcus_3_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Streptococcus_3_gam) # and day effect (random)
plot(Streptococcus_3_gam)

#s(age_yrs.x)                        2.695  2.927 4.121 0.00488 **
#  s(hours_after_sunrise.x)            5.812  6.501 2.770 0.00976 **
#  s(hyenaID2.x,hours_after_sunrise.x) 4.017 11.000 0.684 0.04999 *
  

Streptococcus<-subset(core_genera_gam, Genus == "Streptococcus")

Streptococcus_3_gam <- mgcv::gamm(Abundance~
                                   s(age_yrs.x, bs="cr", k=4)+
                                   matriline.x+
                                   s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                                   #s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                 data=Streptococcus)

sjPlot::tab_model(Streptococcus_3_gam$gam,  show.stat = TRUE, show.se = TRUE, title = "Table S3: GAMM predicting observed ASV richness. Note that the estimate for smooth terms is actually EDF (estimated degrees of freedom), it is not an 'estimate'. EDF reflects the degree of non-linearity of a curve (1 = linear)",
                  dv.labels = "Observed ASV richness",
                  pred.labels = c("Intercept", "matriline.x", "age_yrs.x (non-linear)", "hours_after_sunrise.x (non-linear)"))

############ Fusobacterium ###########
Fusobacterium<-subset(core_genera_gam, Genus == "Fusobacterium")

Fusobacterium_3_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr", k=4)+
                                   matriline.x+
                                   s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                                   #s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                 data=Fusobacterium,
                                 family = gaussian)

print(summary(Fusobacterium_3_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Fusobacterium_3_gam) # and day effect (random)
plot(Fusobacterium_3_gam)

############ Bacteroides ###########
Bacteroides<-subset(core_genera_gam, Genus == "Bacteroides")

Bacteroides_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr", k=4)+
                                   matriline.x+
                                   s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                                   #s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                 data=Bacteroides,
                                 family = gaussian)

print(summary(Bacteroides_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Bacteroides_gam) # and day effect (random)
plot(Bacteroides_gam)

############ Alloprevotella ###########
Alloprevotella<-subset(core_genera_gam, Genus == "Alloprevotella")

Alloprevotella_gam <- mgcv::gam(Abundance~
                               s(age_yrs.x, bs="cr", k=4)+
                               matriline.x+
                               s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                               #s(Seq_depth_scaled)+
                               s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                             data=Alloprevotella,
                             family = gaussian)

print(summary(Alloprevotella_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Alloprevotella_gam) # and day effect (random)
plot(Alloprevotella_gam)

############ Erysipelotrichaceae ###########
Erysipelotrichaceae<-subset(core_genera_gam, Genus == "Erysipelotrichaceae")

Erysipelotrichaceae_gam <- mgcv::gam(Abundance~
                                  s(age_yrs.x, bs="cr", k=4)+
                                  matriline.x+
                                  s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                                  #s(Seq_depth_scaled)+
                                  s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                data=Erysipelotrichaceae,
                                family = gaussian)

print(summary(Erysipelotrichaceae_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Erysipelotrichaceae_gam) # and day effect (random)
plot(Erysipelotrichaceae_gam)

############ Coriobacteriales_unclass ###########
Coriobacteriales_unclass<-subset(core_genera_gam, Genus == "Coriobacteriales_unclass")

Coriobacteriales_unclass_gam <- mgcv::gam(Abundance~
                                       s(age_yrs.x, bs="cr", k=4)+
                                       matriline.x+
                                       s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                                       #s(Seq_depth_scaled)+
                                       s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                     data=Coriobacteriales_unclass,
                                     family = gaussian)

print(summary(Coriobacteriales_unclass_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Coriobacteriales_unclass_gam) # and day effect (random)
plot(Coriobacteriales_unclass_gam)

############ Slackia ###########
Slackia<-subset(core_genera_gam, Genus == "Slackia")

Slackia_gam <- mgcv::gam(Abundance~
                                            s(age_yrs.x, bs="cr", k=4)+
                                            matriline.x+
                                            s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                                            #s(Seq_depth_scaled)+
                                            s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                          data=Slackia,
                                          family = gaussian)

print(summary(Slackia_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Slackia_gam) # and day effect (random)
plot(Slackia_gam)

############ Erysipelotrichaceae ###########
Erysipelotrichaceae<-subset(core_genera_gam, Genus == "Erysipelotrichaceae")

Erysipelotrichaceae_gam <- mgcv::gam(Abundance~
                           s(age_yrs.x, bs="cr", k=4)+
                           matriline.x+
                           s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                           #s(Seq_depth_scaled)+
                           s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                         data=Erysipelotrichaceae,
                         family = gaussian)

print(summary(Erysipelotrichaceae_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Erysipelotrichaceae_gam) # and day effect (random)
plot(Erysipelotrichaceae_gam)

#######  ´Éxtra Enterococcus #####
Enterococcus<-subset(core_genera_gam, Genus == "Enterococcus")

Enterococcus_gam <- mgcv::gam(Abundance~
                                  s(age_yrs.x, bs="cr", k=4)+
                                  matriline.x+
                                  s(hours_after_sunrise.x, bs = "cc", k=8) + # specifically interested in this variable
                                  #s(Seq_depth_scaled)+
                                  s(hyenaID2.x, hours_after_sunrise.x, bs = "re"),
                                data=Enterococcus,
                                family = gaussian)

print(summary(Enterococcus_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Enterococcus_gam) # and day effect (random)
plot(Enterococcus_gam)

#s(age_yrs.x)                        2.807  2.970 20.467  <2e-16 ***
#  s(hours_after_sunrise.x)            1.552  1.928  3.923  0.0172 *  
#  s(hyenaID2.x,hours_after_sunrise.x) 8.901 11.000  8.911  <2e-16 ***




##########
sig_genera<-rbind(Clostridiaceae_1, Streptococcus, Coprococcus_3, Peptococcus, Peptoclostridium, Clostridiales_XIII, Clostridiaceae_1)

genus_plot<-ggplot(sig_genera, aes(hours_after_sunrise.x, Abundance, colour=Genus, fill=Genus))+
  stat_smooth(method="gam", alpha=0.2, size=1.5) + 
  scale_colour_manual(values=c("#A42CD6", "#D050FF", "#E9D4BF", "#D490D3", "#F3C7C7", "#F0E849"))+
  scale_fill_manual(values=c("#A42CD6", "#D050FF", "#E9D4BF", "#D490D3", "#F3C7C7", "#F0E849"))+
  theme_bw()+ xlab("hours after sunrise") +
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))# center align the plot title


age_genera<-subset(core_genera_gam, Genus!="Peptostreptococcus" & Genus!="Lachnospiraceae" & Genus!= "Clostridium_ss1")

age_plot<-ggplot(age_genera, aes(age_yrs.x, Abundance, colour=Genus, fill=Genus))+
  stat_smooth(method="gam", alpha=0.1, size=1.5, aes(fill=Genus)) +
  scale_colour_manual(values=c("#85C18C", #Alloprevotella
                               "#3AB795", #Bacteroides
                               "#A42CD6", #Clostridiaceae_1
                               "#502274", #Clostridiales_unclass
                               "#D050FF", #Clostridiales_XIII
                               "#7A27A5", #Clostridium_ss7
                               "#E9D4BF", #Coprococcus_3
                               "#3AB795", #Coriobacteriales_unclass
                               "#327021", #Eggerthellaceae
                               "#6CBE8F", #Erysipelotrichaceae
                               "#CECA82", #Fusobacterium
                               "#D270E9", #Paeniclostridium
                               "#D490D3", #Peptoclostridium
                               "#EEB4B3", #Peptococcus
                               "#D897B6", #Peptoniphilus
                               "#3D7564", #"Slackia"
                               "#F0E849")) + #Streptococcus)+
                               scale_fill_manual(values=c("#85C18C", #Alloprevotella
                                                            "#3AB795", #Bacteroides
                                                            "#A42CD6", #Clostridiaceae_1
                                                            "#502274", #Clostridiales_unclass
                                                            "#D050FF", #Clostridiales_XIII
                                                            "#7A27A5", #Clostridium_ss7
                                                            "#E9D4BF", #Coprococcus_3
                                                            "#3AB795", #Coriobacteriales_unclass
                                                            "#327021", #Eggerthellaceae
                                                            "#6CBE8F", #Erysipelotrichaceae
                                                            "#CECA82", #Fusobacterium
                                                            "#D270E9", #Paeniclostridium
                                                            "#D490D3", #Peptoclostridium
                                                            "#EEB4B3", #Peptococcus
                                                            "#D897B6", #Peptoniphilus
                                                            "#3D7564", #"Slackia"
                                                            "#F0E849"))+ #Streptococcus)+
  theme_bw()+xlab("age in years")+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))# center align the plot title


ggarrange(
  genus_plot,
  age_plot, 
  labels = c("A", "B"),
  nrow=1, ncol=2, align = c("h")
)



