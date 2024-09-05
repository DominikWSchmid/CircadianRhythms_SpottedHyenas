
#title: "Hyena Microbiome"
#author: "Dominik Melville"
#date: "26/07/2023"


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
#library(DT)
#library(eulerr)
library(microbiomeutilities)
library(ggrepel)
#library(sjPlot)
#library(sjmisc)
library(MuMIn)
#library(jtools)
library(car)
library(plyr)

# Finally, set a random seed for reproducibility
set.seed(777)

###Upload ASV, taxonomy and meta data and compile into phyloseq object thats saved as RDS file (based on previous filtering)

ps = readRDS("/Users/dom/Dropbox/Hyena/hyena_phyloseq.rds")

### lets remind ourselves of some basics about our phyloseq object 

microbiome::summarize_phyloseq(ps)

### remove singeltons 
ps <- filter_taxa(ps, function (x) {sum(x > 0) >1}, prune=TRUE)

#mean(sample_sums(ps)) 
#sd(sample_sums(ps))

### turn row names into sample ids
sample_data(ps)$sample.ID<-row.names(sample_data(ps))

# remove one samples with NA
ps  <- ps %>% subset_samples(!sample.ID %in% c("P180"#, #no time specified
))

#format the time variable
sample_data(ps)$sample_time<-as.POSIXct(paste(sample_data(ps)$sample_date, sample_data(ps)$sample_time), format ="%Y-%m-%d %H:%M:%S", tz="Africa/Nairobi")

#calcualte time since sunrise
meta <- as(sample_data(ps), "data.frame")
light<-suncalc::getSunlightTimes(meta$sample_date, lat=1.482, lon=35.130, tz = "Africa/Nairobi")
sample_data(ps)$sunrise<-light$sunrise
sample_data(ps)$hours_after_sunrise<-sample_data(ps)$sample_time-sample_data(ps)$sunrise
sample_data(ps)$hours_after_sunrise<-as.numeric(round((sample_data(ps)$hours_after_sunrise/60)/60, 2))

#make time category
sample_data(ps)$day_time<-ifelse(sample_data(ps)$hours_after_sunrise <= 5.00, "morning", "afternoon")

#format matriline 
sample_data(ps)$matriline<-as.factor(sample_data(ps)$matriline)
sample_data(ps)$seq_depth<-as.numeric(sample_sums(ps))
sample_data(ps)$hours_after_sunrise_rescaled<-as.numeric(scales::rescale(sample_data(ps)$hours_after_sunrise, to =c(-1,1)))
sample_data(ps)$prey_abundance_avg_month_rescaled<-as.numeric(scales::rescale(sample_data(ps)$prey_abundance_avg_month, to =c(-1,1)))

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

library(ggplot2)

class_plot_per_sample<-ggplot(class_per_sample, aes(y=Abundance, x=reorder(hours_after_sunrise, sample.ID), fill=class_plot))+
  geom_col(position = "fill", size=0, width=1)+
  scale_fill_manual(values=c("#7082BC","#E9D752", "#A59A86", "#A4BDF1", "#675524", "#D1C9A3", "#EFEBF5", "#A6983E", "grey75"))+
  ylab("Relative Abundance")+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank(), legend.position="bottom", axis.title.y = element_text(colour="black",size=12,face="bold"))+
  scale_y_continuous(expand = c(0,0))+ labs(fill ="Bacterial Class")+geom_vline(xintercept = 103.5, size=1.5, linetype="dashed", colour="black")+
  ggtitle("morning                                                                                                                                           afternoon")



##### per genus
ps_genus<-tax_glom(ps, taxrank = "Genus")
core<-microbiome::core(ps_genus, detection = 0, prevalence = 0.85)
core.abundance <- sample_sums(core(ps_genus, detection =0, prevalence = .85))/sample_sums(ps_genus)
mean(core.abundance)
sd(core.abundance)

top<-as_tibble(get_group_abundances(core, level="Genus", group="Genus", transform = "compositional") %>%arrange (-mean_abundance))
print(top,n=26)
#1 Lachnospiraceae          Lachnospiraceae                 0.188        0.133  
#2 Micrococcaceae           Micrococcaceae                  0.0941       0.198  
#3 Clostridium_ss1          Clostridium_ss1                 0.0940       0.103  
#4 Clostridiales_unclass    Clostridiales_unclass           0.0806       0.119  
#5 Peptoniphilus            Peptoniphilus                   0.0591       0.107  
#6 Peptoclostridium         Peptoclostridium                0.0562       0.0820 
#7 Enterococcus             Enterococcus                    0.0553       0.0880 
#8 Fusobacterium            Fusobacterium                   0.0508       0.0905 
#9 Bacteroides              Bacteroides                     0.0416       0.0847 
#10 Erysipelotrichaceae      Erysipelotrichaceae             0.0349       0.0419 
#11 Paeniclostridium         Paeniclostridium                0.0307       0.0508 
#12 Alloprevotella           Alloprevotella                  0.0298       0.0684 
#13 Clostridiaceae_1         Clostridiaceae_1                0.0267       0.0457 
#14 Clostridium_ss7          Clostridium_ss7                 0.0235       0.0441 
#15 Streptococcus            Streptococcus                   0.0187       0.0508 
#16 Clostridiales_XIII       Clostridiales_XIII              0.0174       0.0217 
#17 Faecalitalea             Faecalitalea                    0.0159       0.0330 
#18 Coriobacteriales_unclass Coriobacteriales_unclass        0.0152       0.0255 
#19 Hathewaya                Hathewaya                       0.0131       0.0243 
#20 Slackia                  Slackia                         0.0125       0.0216 
#21 Peptostreptococcus       Peptostreptococcus              0.0101       0.0195 
#22 Collinsella              Collinsella                     0.00895      0.0257 
#23 Peptococcus              Peptococcus                     0.00806      0.0143 
#24 Coprococcus_3            Coprococcus_3                   0.00736      0.0124 
#25 Ruminococcaceae          Ruminococcaceae                 0.00514      0.0115 
#26 Eggerthellaceae          Eggerthellaceae                 0.00260      0.00346


top_genus<-c("Lachnospiraceae", "Clostridium_ss1", "Clostridiales_unclass", "Clostridium_ss7", "Clostridiaceae_1",  
             "Clostridiales_XIII", "Paeniclostridium", "Peptoclostridium", "Peptoniphilus",  "Peptostreptococcus",
             "Peptococcus", "Coprococcus_3", "Streptococcus", ### all Clostridia
             "Fusobacterium", "Enterococcus", "Faecalitalea", "Collinsella", 
             "Bacteroides", "Alloprevotella", "Micrococcaceae",
              "Erysipelotrichaceae",  "Hathewaya",
             "Coriobacteriales_unclass", "Slackia","Eggerthellaceae", "Ruminococcaceae") 

genus_per_sample$genus_plot<-as.character(ifelse(genus_per_sample$Genus %in% top_genus, genus_per_sample$Genus, "Other"))

genus_per_sample$genus_plot<-factor(genus_per_sample$genus_plot, 
                                    levels = c("Lachnospiraceae", "Clostridium_ss1", "Clostridiales_unclass", "Clostridium_ss7", "Clostridiaceae_1",  
                                               "Clostridiales_XIII", "Paeniclostridium", "Peptoclostridium", "Peptoniphilus",  "Peptostreptococcus",
                                               "Peptococcus", "Coprococcus_3", "Hathewaya", "Ruminococcaceae", ### all Clostridia
                                               "Streptococcus", "Enterococcus", ### Bacilli
                                               "Micrococcaceae", "Collinsella", #Actinobacteria
                                               "Bacteroides", "Alloprevotella", #Bacteroidia
                                               "Fusobacterium", #Fusobacteria
                                               "Faecalitalea", "Erysipelotrichaceae", #Erysipelotrichia
                                               "Coriobacteriales_unclass", "Slackia", "Eggerthellaceae", #Coriobacteriia
                                               "Other"))

plot.meta<-meta[,c("sample.ID", "hours_after_sunrise")]

#str_sort(plot.meta$hours_after_sunrise, numeric = TRUE)

genus_per_sample<-merge(genus_per_sample, plot.meta, by="sample.ID")

genus_plot_per_sample<-ggplot(genus_per_sample, aes(y=Abundance, x=reorder(hours_after_sunrise, sample.ID), fill=genus_plot))+
  geom_col(position = "fill", size=0, width=1)+
  scale_fill_manual(values=c("#2E417F", "#334684", "#384B89", "#425593", "#5265A4","#5669A7","#5A6DAA","#6174B0", 
                             "#697BB6","#7082BC","#7B8DC5","#8697CE","#91A1D7","#9BABDF", 
                             "#F2E15F","#D2C558",
                             "#988E7B", "#A7A185", 
                             "#AEC4F3", "#94B0EB", 
                             "#58491D", 
                             "#B3A982", "#D9D2B8", 
                             "#E8E3E1", "#EFEBF5", "#EAE3F7",
                             "grey75"))+
  ylab("Relative Abundance")+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), 
        strip.text = element_blank(), legend.position = "bottom", axis.title.y = element_text(colour="black",size=12,face="bold"))+
  scale_y_continuous(expand = c(0,0))+ labs(fill ="Bacterial Genus")+geom_vline(xintercept = 103.5, size=1.5, linetype="dashed", colour="black")+
  ggtitle("morning                                                                                                                                           afternoon")


ggarrange(
  class_plot_per_sample,
  genus_plot_per_sample, 
  labels = c("A", "B"),
  nrow=2, ncol=1, align = c("h")
)


#####  alpha diversity metrices and some stats ######
#####################################################

#### prey abundance averaged for the last month #####
# This proxy was used in the previous publication by Rojas et al. 2023 #
# However, it is a very coarse proxy since there is no immediate #
# connection to prey capture success and it, unsurprisingly, explained #
# very little variation (if any). For the purpose of this analysis #
# we decided to abandon the variable, which is also not available for all #
# hyenas. Yet, we provide code here to test its influence yourself #
# simply un-comment some currently commented-out lines #
# for further questions get in touch with: dominikwerner.schmid@uni-ulm.de #

#first we rarefy 

ps.rarefied <- rarefy_even_depth(ps, sample.size = 2900, rngseed = 123)


otu.table <- as.data.frame(otu_table(ps.rarefied))
meta <- as(sample_data(ps.rarefied), "data.frame")
df.pd <- pd(t(otu.table), phy_tree(ps.rarefied), include.root=T) #what is SR? species richness?
meta$Phyogenetic_diversity <- df.pd$PD

alpha = estimate_richness(ps.rarefied, measures = c("Observed", "Shannon")) #generate df with various diversity metices 

#add indicies to meta data (like we did with PD)
meta$Shannon <- alpha$Shannon 
meta$Observed <- alpha$Observed

#rescale or refit parameters
meta$hyenaID2<-as.numeric(meta$hyenaID2)
meta$hours_after_sunrise_rescaled<-as.numeric(scales::rescale(meta$hours_after_sunrise, to =c(-1,1)))
meta$prey_abundance_avg_month_rescaled<-as.numeric(scales::rescale(meta$prey_abundance_avg_month, to =c(-1,1)))

##### modelling using gams
observed_gam <- mgcv::gam(log(Observed)~
                             s(age_yrs, bs="cr")+
                             matriline+
                             s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                             #s(prey_abundance_avg_month_rescaled, bs = "cr")+
                             s(hyenaID2, bs = "re"),
                           data=meta,
                           family = gaussian)

observed_gam <- mgcv::gam(log(Observed)~
                            s(age_yrs, bs="cr")+
                            matriline+
                            day_time + # specifically interested in this variable
                            #s(prey_abundance_avg_month_rescaled, bs = "cr")+
                            s(hyenaID2, bs = "re"),
                          data=meta,
                          family = gaussian)

print(summary(observed_gam)) # no treatment effect (great)
mgcv::gam.check(observed_gam)
plot(observed_gam)

library(mgcViz)
observed_plot<-ggplot(meta, aes(hours_after_sunrise, Observed, colour=hours_after_sunrise))+
  geom_point(size=3)+
  stat_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=5), lty="dotted", colour="black")+
  xlab("hours since sunrise") + ylab("Observed ASV richness")+theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position = "none") +
  scale_colour_gradient(low="#F7ED5D",high="#1538A0")+labs(fill = "")

sample_size_plot<-ggplot(meta, aes(hours_after_sunrise, colour=hours_after_sunrise))+
  geom_histogram(colour="black", fill="grey75", bins=15)+theme_bw() + ylab("Sample Sizes")+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank())

Alpha.SampleSize_plot<-ggarrange(
  sample_size_plot,
  observed_plot,
  labels = c("A", "B"),
  nrow=2, ncol=1, align = c("v"), heights=c(1,3)
)

##### modelling using gams
PD_gam <- mgcv::gam((Phyogenetic_diversity)~ #same for Shannon
                            s(age_yrs, bs="cr")+
                            matriline+
                            s(hours_after_sunrise_rescaled, bs = "cr") + 
                            #s(prey_abundance_avg_month_rescaled, bs = "cr") +
                            s(hyenaID2, bs = "re"),
                          data=meta,
                          family = gaussian)

print(summary(PD_gam)) 
mgcv::gam.check(PD_gam)
plot(PD_gam)


##### modelling using gams
shannon_gam <- mgcv::gam((Phyogenetic_diversity)~ #same for Shannon
                           s(age_yrs, bs="cr")+
                           matriline+
                           day_time + 
                           #s(prey_abundance_avg_month_rescaled, bs = "cr")+
                           s(hyenaID2, bs = "re"),
                         data=meta,
                         family = gaussian)

print(summary(shannon_gam)) 
mgcv::gam.check(shannon_gam)
plot(shannon_gam)



otu.table <- as.data.frame(otu_table(ps))
meta_unrarefied <- as(sample_data(ps), "data.frame")
alpha = estimate_richness(ps, measures = c("Observed", "Shannon")) #generate df with various diversity metices 
meta_unrarefied$Observed_unrarefied<-alpha$Observed
meta_unrarefied<-subset(meta_unrarefied, seq_depth>=2900)
meta$Observed_unrarefied<-meta_unrarefied$Observed_unrarefied
ggplot(meta, aes(Observed, Observed_unrarefied))+geom_point()
cor.test(meta$Observed, meta$Observed_unrarefied)
### rarefying makes sense ### 
### would anyway, given Weiss 2023 ###



######## beta diversity index ######

ps.beta<-ps.rarefied
meta.beta<- as(sample_data(ps.beta), "data.frame")
dist_wUni = phyloseq::distance(ps.beta, method="weighted unifrac") ### abundance matters + phylogeny matters 
dist_Uni = phyloseq::distance(ps.beta, method="unifrac") ### presence/absence + phylogeny matters
dist_bray = phyloseq::distance(ps.beta, method="bray") ### abundance matters 
dist_aitchison = phyloseq::distance(ps.beta, method="euclidean") ### presence/absence 

#### PERMANOVA (since in rojas et al. 2023 year and prey abundance had no effect we did not include it here)
adonis2(dist_wUni ~ hours_after_sunrise_rescaled + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
        method="wunifrac", by="margin")
adonis2(dist_wUni ~ day_time + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
        method="wunifrac", by="margin")

adonis2(dist_Uni ~ hours_after_sunrise_rescaled + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
        method="unifrac", by="margin")
adonis2(dist_Uni ~ day_time + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
        method="unifrac", by="margin")

anova(betadisper(dist_Uni, meta.beta$hours_after_sunrise_rescaled))
permutest(betadisper(dist_wUni, meta.beta$hours_after_sunrise_rescaled))

anova(betadisper(dist_wUni, meta.beta$day_time))
permutest(betadisper(dist_wUni, meta.beta$day_time))


#adonis2(dist_bray ~ hours_after_sunrise_rescaled + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
#        method="wunifrac", by="margin")
#adonis2(dist_bray ~ day_time + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
#        method="wunifrac", by="margin")

#adonis2(dist_aitchison ~ hours_after_sunrise_rescaled + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
#        method="unifrac", by="margin")
#adonis2(dist_aitchison ~ day_time + matriline + age_yrs, data = meta.beta, strata=meta.beta$hyenaID2, 
#        method="unifrac", by="margin")

#anova(betadisper(dist_Uni, meta.beta$hours_after_sunrise_rescaled))
#permutest(betadisper(dist_wUni, meta.beta$hours_after_sunrise_rescaled))

#anova(betadisper(dist_wUni, meta.beta$day_time))
#permutest(betadisper(dist_wUni, meta.beta$day_time))

#### but for completeness we have included prey abundance here
#make phyloseq removing the those samples with prey abundance = NA
#ps.beta_prey  <- ps.rarefied %>% 
#  subset_samples(!prey_abundance_avg_month %in% NA)

#meta.beta_prey<- as(sample_data(ps.beta_prey), "data.frame")
#prey_dist_wUni = phyloseq::distance(ps.beta_prey, method="weighted unifrac") ### abundance matters + phylogeny matters 
#prey_dist_Uni = phyloseq::distance(ps.beta_prey, method="unifrac") ### presence/absence + phylogeny matters


#adonis2(prey_dist_Uni ~ hours_after_sunrise_rescaled + matriline + age_yrs +prey_abundance_avg_month_rescaled, data = meta.beta_prey, strata=meta.beta_prey$hyenaID2, 
#        method="unifrac", by="margin")
# rarified results

#adonis2(prey_dist_wUni ~ hours_after_sunrise_rescaled + matriline + age_yrs + prey_abundance_avg_month_rescaled, data = meta.beta_prey, strata=meta.beta_prey$hyenaID2, 
#        method="weighted unifrac", by="margin")
# rarified results


######## is it the same if we just consider core taxa?
ps.beta.core<-ps
meta.beta.core<- as(sample_data(ps.beta.core), "data.frame")
ps_genus<-tax_glom(ps.beta.core, taxrank = "Genus") 
core<-microbiome::core(ps_genus, detection = 0, prevalence = 0.85)

top<-as_tibble(get_group_abundances(core, level="Genus", group="Genus", transform = "compositional") %>%arrange (-mean_abundance))
print(top,n=26)

dist_wUni = phyloseq::distance(core, method="weighted unifrac") ### abundance matters + phylogeny matters 
dist_Uni = phyloseq::distance(core, method="unifrac") ### presence/absence + phylogeny matters
dist_bray = phyloseq::distance(core, method="bray") ### abundance matters  
dist_aitchison = phyloseq::distance(core, method="euclidean") ### presence/absence 

adonis2(dist_Uni ~ hours_after_sunrise_rescaled + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
        method="unifrac", by="margin")
adonis2(dist_Uni ~ day_time + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
        method="unifrac", by="margin")

adonis2(dist_wUni ~ hours_after_sunrise_rescaled + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
        method="weighted unifrac", by="margin")
adonis2(dist_wUni ~ day_time + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
        method="weighted unifrac", by="margin")

anova(betadisper(dist_wUni, meta.beta.core$hours_after_sunrise_rescaled))
permutest(betadisper(dist_wUni, meta.beta.core$hours_after_sunrise_rescaled))

anova(betadisper(dist_wUni, meta.beta.core$day_time))
permutest(betadisper(dist_wUni, meta.beta.core$day_time))


#adonis2(dist_bray ~ hours_after_sunrise_rescaled + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
#        method="weighted unifrac", by="margin")
#adonis2(dist_bray ~ day_time + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
#        method="weighted unifrac", by="margin")

#adonis2(dist_aitchison ~ hours_after_sunrise_rescaled + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
#        method="weighted unifrac", by="margin")
#adonis2(dist_aitchison ~ day_time + matriline + age_yrs, data = meta.beta.core, strata=meta.beta.core$hyenaID2, 
#        method="weighted unifrac", by="margin")

#anova(betadisper(dist_aitchison, meta.beta.core$hours_after_sunrise_rescaled))
#permutest(betadisper(dist_aitchison, meta.beta.core$hours_after_sunrise_rescaled))

#anova(betadisper(dist_aitchison, meta.beta.core$day_time))
#permutest(betadisper(dist_aitchison, meta.beta.core$day_time))

###
#### but for completeness we have included prey abundance here
#make phyloseq removing the those samples with prey abundance = NA
#ps.beta.core_prey  <- ps.beta.core %>% 
#  subset_samples(!prey_abundance_avg_month %in% NA)

#meta.beta.core_prey<- as(sample_data(ps.beta.core_prey), "data.frame")
#prey_dist_wUni = phyloseq::distance(ps.beta.core_prey, method="weighted unifrac") ### abundance matters + phylogeny matters 
#prey_dist_Uni = phyloseq::distance(ps.beta.core_prey, method="unifrac") ### presence/absence + phylogeny matters

#adonis2(prey_dist_wUni ~ day_time + matriline + age_yrs +prey_abundance_avg_month_rescaled, data = meta.beta.core_prey, strata=meta.beta.core_prey$hyenaID2, 
#        method="weighted unifrac", by="margin")
# rarified results

#adonis2(prey_dist_wUni ~ hours_after_sunrise_rescaled + matriline + age_yrs + prey_abundance_avg_month_rescaled, data = meta.beta.core_prey, strata=meta.beta.core_prey$hyenaID2, 
#        method="weighted unifrac", by="margin")


####
ordination<-ordinate(core, method="DCA", distance="unweighted Unifrac")

ord_plot <- plot_ordination(core, ordination) + 
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = day_time)) + # put ellipses around group centroids
  geom_point(aes(
    #fill = day_time, 
    col = hours_after_sunrise), # fill color by infection status
    size = 3) +# make points size 4
    #pch = 21, # Make points circular with a border 
    #col = hours_after_sunrise) + # the border colour should be black
  #labs(subtitle = "a)") + # insert subtitle
  #ggtitle("Unweighted Unifrac") + # insert title
  theme_bw() + scale_fill_manual(values=c("#1538A0", "#F7ED5D")) +
  labs(fill = "time of day", colour="hours since sunrise") + scale_colour_gradient(low="#F7ED5D",high="#1538A0")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold"))+
  ylab("Axis 2 [25.2%]")+ xlab("Axis 1 [41.6%]")+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))# center align the plot title
ord_plot


ordination<-ordinate(ps.rarefied, method="DCA", distance="bray")

ord_plot_rare <- plot_ordination(ps.rarefied, ordination) + 
  stat_ellipse(geom = "polygon", alpha = 0.3, aes(fill = day_time)) + # put ellipses around group centroids
  geom_point(aes(
    #fill = day_time, 
    col = hours_after_sunrise), # fill color by infection status
    size = 3) +# make points size 4
  #pch = 21, # Make points circular with a border 
  #col = hours_after_sunrise) + # the border colour should be black
  #labs(subtitle = "a)") + # insert subtitle
  #ggtitle("Unweighted Unifrac") + # insert title
  theme_bw() + scale_fill_manual(values=c("#1538A0", "#F7ED5D")) +
  labs(fill = "time of day", colour="hours since sunrise") + scale_colour_gradient(low="#F7ED5D",high="#1538A0")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold"))+
  ylab("Axis 2 [25.2%]")+ xlab("Axis 1 [41.6%]")+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))# center align the plot title
ord_plot_rare


ggarrange(
  Alpha.SampleSize_plot,
  ord_plot, 
  labels = c("",  "C"),
  nrow=1, ncol=2, legend ="bottom", common.legend = T
)



###### joint species distribution modelling #####
library(gllvm)

##### transform data to account for its compositionality
ps.gllvm.clr <-microbiome::transform((otu_table(core)), transform = "clr")

genus<-data.frame(tax_table(core)[,"Genus"])
taxa_names(core)<-genus$Genus

ASV<-data.frame(t(otu_table(core)))

core.clr <- data.frame(sample_data(core))
core.clr$seq_depth<-as.numeric(scales::rescale(sample_sums(core), to =c(-1,1)))
core.clr<-core.clr[,c(1,3,8,15,16,18,19)]

gllvm1<-gllvm(ASV, 
              core.clr, #otu table and meta data as data frame 
              formula = ~ 
                hours_after_sunrise_rescaled + age_yrs + matriline + seq_depth, 
              family = "negative.binomial", 
              row.eff = ~(1|hyenaID), #random effect
              num.lv = 3)

par(mfrow = c(1, 2))
plot(gllvm1, which = 1:2, var.colors = 1, n.plot = 20) #looks good
summary(gllvm1)

##### building a graph ##### 

df <- coef(gllvm1)
est_df <- data.frame(df$Intercept)
est_df2 <- data.frame(df$Xcoef)
est_df3 <- merge(est_df, est_df2, by=0)

#reorder genera

row.names(est_df3)<-est_df3$Row.names
est_df3<-est_df3[colnames(ASV),]

#turn df into long format 

names(est_df3)[1]<-"Genus"
names(est_df3)[2]<-"Intercept"
estimates_df <- gather(est_df3, Treatment, Estimate, names(est_df3)[2]:names(est_df3)[ncol(est_df3)], factor_key=TRUE)

# extract CIs

confint_df <- data.frame(confint(gllvm1))
confint_df <- rbind(confint_df[grep("^Xcoef", rownames(confint_df)), ], confint_df[grep("^Intercept", rownames(confint_df)), ]) 

# adding the variable levels as a column

variables <- colnames(est_df3)[3:ncol(est_df3)]
variables <- c(variables, "Intercept")
variables1 <- rep(variables, nrow(est_df3))
variables2 <- variables1[order(match(variables1, variables))]

confint_df$Treatment<-variables2

# column with taxa names

confint_df$Genus <- rep(colnames(ASV), length(unique(confint_df$Treatment)))

# now have estimates and confidence intervals as seperate data frames, but they are in different formats. Need to massage them into one dataframe for plotting.

merged <- merge(estimates_df, confint_df, by =c("Treatment", "Genus"))

names(merged)[4]<-"CI_lower"
names(merged)[5]<-"CI_upper"

unique(merged$Treatment)


merged$Sig <- !data.table::between(0, merged$CI_lower, merged$CI_upper) 
merged$Genus<-as.factor(merged$Genus)

write.csv(merged, file="/Users/dominikschmid/Dropbox (Personal)/Hyena/gllvm_results_29.07.24.csv")

subset_merged<-subset(merged, Treatment != "Intercept" & 
                        Treatment != "SeqDepth")

variable_names <- c('age_yrs' = "age in years", 
                    'hours_after_sunrise_rescaled' = "hours after sunrise", 
                    'matriline2' = "matriline 2",
                    'matriline3' = "matriline 3", 
                    'matriline4' = "matriline 4")


ggplot(subset_merged, aes(y = Genus, x=Estimate, fill= Sig))+
  facet_wrap(~Treatment, scales="free_x", nrow=1, labeller = as_labeller(variable_names))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey") + 
  geom_errorbar(aes(xmin=CI_lower, xmax = CI_upper),width=0.5, size=0.5, colour="black")+
  geom_point(size=2.5, pch=21)+ 
  theme_classic(base_size = 12) +
  scale_fill_manual(values=c("grey85", "cornflowerblue"))+
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(face="bold"), axis.title.y = element_blank())+
  theme(axis.text.y = element_text(colour=subset_merged$Colour)) +
  scale_y_discrete(label = function(x) abbreviate(x, minlength=7, TRUE))


subset_merged2<-subset(merged, Treatment == "hours_after_sunrise_rescaled" )

gllvm.plot<-ggplot(subset_merged2, aes(y = Genus, x=Estimate, fill= Sig))+
  facet_wrap(~Treatment, scales="free_x", nrow=1, labeller = as_labeller(variable_names))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey") + 
  geom_errorbar(aes(xmin=CI_lower, xmax = CI_upper),width=0.5, size=0.5, colour="black")+
  geom_point(size=2.5, pch=21)+ 
  theme_classic(base_size = 12) +
  scale_fill_manual(values=c("grey85", "cornflowerblue"))+
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(face="bold"), axis.title.y = element_blank())+
  theme(axis.text.y = element_text(colour=subset_merged$Colour)) +
  scale_y_discrete(label = function(x) abbreviate(x, minlength=7, TRUE))

###### same when including prey abundance ######

#core_prey  <- core %>% 
#  subset_samples(!prey_abundance_avg_month %in% NA)


##### transform data to account for its compositionality
#ps.prey.gllvm.clr <-microbiome::transform((otu_table(core_prey)), transform = "clr")

#genus<-data.frame(tax_table(core_prey)[,"Genus"])
#taxa_names(core_prey)<-genus$Genus

#ASV<-data.frame(t(otu_table(core_prey)))

#core.prey.clr <- data.frame(sample_data(core_prey))
#core.prey.clr$seq_depth<-as.numeric(scales::rescale(sample_sums(core_prey), to =c(-1,1)))
#core.prey.clr<-core.prey.clr[,c(1,3,8,15,16,18,19,20)]

#gllvm1<-gllvm(ASV, 
#              core.prey.clr, #otu table and meta data as data frame 
#              formula = ~ 
#                hours_after_sunrise_rescaled + age_yrs + matriline + seq_depth + prey_abundance_avg_month_rescaled, 
#              family = "negative.binomial", 
#              row.eff = ~(1|hyenaID), #random effect
#              num.lv = 3)

#par(mfrow = c(1, 2))
#plot(gllvm1, which = 1:2, var.colors = 1, n.plot = 20) #looks good
#summary(gllvm1)

##### building a graph ##### 

#df <- coef(gllvm1)
#est_df <- data.frame(df$Intercept)
#est_df2 <- data.frame(df$Xcoef)
#est_df3 <- merge(est_df, est_df2, by=0)

#reorder genera

#row.names(est_df3)<-est_df3$Row.names
#est_df3<-est_df3[colnames(ASV),]

#turn df into long format 

#names(est_df3)[1]<-"Genus"
#names(est_df3)[2]<-"Intercept"
#estimates_df <- gather(est_df3, Treatment, Estimate, names(est_df3)[2]:names(est_df3)[ncol(est_df3)], factor_key=TRUE)

# extract CIs

#confint_df <- data.frame(confint(gllvm1))
#confint_df <- rbind(confint_df[grep("^Xcoef", rownames(confint_df)), ], confint_df[grep("^Intercept", rownames(confint_df)), ]) 

# adding the variable levels as a column

#variables <- colnames(est_df3)[3:ncol(est_df3)]
#variables <- c(variables, "Intercept")
#variables1 <- rep(variables, nrow(est_df3))
#variables2 <- variables1[order(match(variables1, variables))]

#confint_df$Treatment<-variables2

# column with taxa names

#confint_df$Genus <- rep(colnames(ASV), length(unique(confint_df$Treatment)))

# now have estimates and confidence intervals as seperate data frames, but they are in different formats. Need to massage them into one dataframe for plotting.

#merged <- merge(estimates_df, confint_df, by =c("Treatment", "Genus"))

#names(merged)[4]<-"CI_lower"
#names(merged)[5]<-"CI_upper"

#unique(merged$Treatment)


#merged$Sig <- !data.table::between(0, merged$CI_lower, merged$CI_upper) 
#merged$Genus<-as.factor(merged$Genus)

#write.csv(merged, file="/Users/dom/Dropbox/Hyena/gllvm_results_29.07.24.csv")

#subset_merged<-subset(merged, Treatment != "Intercept" & 
#                        Treatment != "SeqDepth")

#variable_names <- c('age_yrs' = "age in years", 
#                    'hours_after_sunrise_rescaled' = "hours after sunrise", 
#                    'matriline2' = "matriline 2",
#                    'matriline3' = "matriline 3", 
#                    'matriline4' = "matriline 4", 
#                    'seq_depth' = "sequencing depth",
#                    'prey_abundance_avg_month_rescaled' = "average prey abundance")


#ggplot(subset_merged, aes(y = Genus, x=Estimate, fill= Sig))+
#  facet_wrap(~Treatment, scales="free_x", nrow=1, labeller = as_labeller(variable_names))+
#  geom_vline(xintercept=0, linetype="dashed", colour="grey") + 
#  geom_errorbar(aes(xmin=CI_lower, xmax = CI_upper),width=0.5, size=0.5, colour="black")+
#  geom_point(size=2.5, pch=21)+ 
#  theme_classic(base_size = 12) +
#  scale_fill_manual(values=c("grey85", "cornflowerblue"))+
#  theme(legend.position = "none") +
#  theme(axis.text.y = element_text(face="bold"), axis.title.y = element_blank())+
#  theme(axis.text.y = element_text(colour=subset_merged$Colour)) +
#  scale_y_discrete(label = function(x) abbreviate(x, minlength=7, TRUE))

#subset_merged3<-subset(merged, Treatment == "prey_abundance_avg_month_rescaled" )

#prey.plot<-ggplot(subset_merged3, aes(y = Genus, x=Estimate, fill= Sig))+
#  facet_wrap(~Treatment, scales="free_x", nrow=1, labeller = as_labeller(variable_names))+
#  geom_vline(xintercept=0, linetype="dashed", colour="grey") + 
#  geom_errorbar(aes(xmin=CI_lower, xmax = CI_upper),width=0.5, size=0.5, colour="black")+
#  geom_point(size=2.5, pch=21)+ 
#  theme_classic(base_size = 12) +
#  scale_fill_manual(values=c("grey85", "cornflowerblue"))+
#  theme(legend.position = "none") +
#  theme(axis.text.y = element_text(face="bold"), axis.title.y = element_blank())+
#  theme(axis.text.y = element_text(colour=subset_merged$Colour)) +
#  scale_y_discrete(label = function(x) abbreviate(x, minlength=7, TRUE))



##### GAMS ######
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
sample_data(clr_scaled)<-sample_data(clr_scaled)[,c( "sample.ID", "hyenaID2", "age_yrs", "matriline", "hours_after_sunrise", "day_time")]
core_genera_melt<-psmelt(clr_scaled)
head(core_genera_melt)
meta.gam<-as(sample_data(clr_scaled), "data.frame")
meta.gam$Seq_depth<-as.numeric(sample_sums(core))
meta.gam$Seq_depth_scaled <- as.numeric(scales::rescale(meta.gam$Seq_depth, to = c(-1,1)))
meta.gam$hours_after_sunrise_rescaled <- as.numeric(scales::rescale(meta.gam$hours_after_sunrise, to = c(-1,1)))
core_genera_gam<- merge(core_genera_melt, meta.gam, by.x = "Sample", by.y = "sample.ID")
names(core_genera_gam)
head(core_genera_gam)
str(core_genera_gam)

#### GAMMS
library(mgcv)
library(tidygam)
library(tidymv)

######### Faecalitalea #########
Faecalitalea<-subset(core_genera_gam, Genus == "Faecalitalea")

Faecalitalea_gam <- mgcv::gam(Abundance~
                      s(age_yrs.x, bs="cr")+
                      matriline.x+
                      s(hours_after_sunrise_rescaled, bs = "cr") + 
                      s(Seq_depth_scaled)+
                      s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                      data=Faecalitalea,
                      family = gaussian)

print(summary(Faecalitalea_gam)) # no treatment effect (great)
gam.check(Faecalitalea_gam)
plot(Faecalitalea_gam)


####### Collinsella ########
Collinsella<-subset(core_genera_gam, Genus == "Collinsella")

Collinsella_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr")+
                                   matriline.x+
                                   s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                                   s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                                 data=Collinsella,
                                 family = gaussian)

print(summary(Collinsella_gam)) #strong treatment effect (fewer in Staphylococcus in treatment)
gam.check(Collinsella_gam) # plus day effect (random)
plot(Collinsella_gam)


###########  Enterococcus    ##############
Enterococcus<-subset(core_genera_gam, Genus == "Enterococcus")

Enterococcus_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr")+
                                   matriline.x+
                                   s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                                   s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                                 data=Enterococcus,
                                 family = gaussian)

print(summary(Enterococcus_gam)) # no effect of treatment but day (random effect)
gam.check(Enterococcus_gam)
plot(Enterococcus_gam)

###########  Ruminococcaceae    ##############
Ruminococcaceae<-subset(core_genera_gam, Genus == "Ruminococcaceae")

Ruminococcaceae_gam <- mgcv::gam(Abundance~
                                         s(age_yrs.x, bs="cr")+
                                         matriline.x+
                                         s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                                         s(Seq_depth_scaled)+
                                         s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                                       data=Ruminococcaceae,
                                       family = gaussian)

print(summary(Ruminococcaceae_gam)) 
gam.check(Ruminococcaceae_gam)
plot(Ruminococcaceae_gam)

###########  Clostridiales_XIII    ##############
Clostridiales_XIII<-subset(core_genera_gam, Genus == "Clostridiales_XIII")

Clostridiales_XIII_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr")+
                                   matriline.x+
                                   s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                                   s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                                 data=Clostridiales_XIII,
                                 family = gaussian)

print(summary(Clostridiales_XIII_gam)) ## no effect of treatment on abundance (good)
gam.check(Clostridiales_XIII_gam)
plot(Clostridiales_XIII_gam)


############ Peptoclostridium ###########
Peptoclostridium<-subset(core_genera_gam, Genus == "Peptoclostridium")

Peptoclostridium_gam <- mgcv::gam(Abundance~
                                    s(age_yrs.x, bs="cr")+
                                    matriline.x+
                                    s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                                    s(Seq_depth_scaled)+
                                    s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                                  data=Peptoclostridium,
                                  family = gaussian)

print(summary(Peptoclostridium_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Peptoclostridium_gam) # and day effect (random)
plot(Peptoclostridium_gam)


############ Peptococcus ###########
Peptococcus<-subset(core_genera_gam, Genus == "Peptococcus")

Peptococcus_gam <- mgcv::gam(Abundance~
                                      s(age_yrs.x, bs="cr")+
                                      matriline.x+
                                      s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                                      s(Seq_depth_scaled)+
                                      s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                                    data=Peptococcus,
                                    family = gaussian)

print(summary(Peptococcus_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Peptococcus_gam) # and day effect (random)
plot(Peptococcus_gam)


############ Coprococcus_3 ###########
Coprococcus_3<-subset(core_genera_gam, Genus == "Coprococcus_3")

Coprococcus_3_gam <- mgcv::gam(Abundance~
                               s(age_yrs.x, bs="cr")+
                               matriline.x+
                               s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                               s(Seq_depth_scaled)+
                               s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                             data=Coprococcus_3,
                             family = gaussian)

print(summary(Coprococcus_3_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Coprococcus_3_gam) # and day effect (random)
plot(Coprococcus_3_gam)

############ Streptococcus ###########
Streptococcus<-subset(core_genera_gam, Genus == "Streptococcus")

Streptococcus_3_gam <- mgcv::gam(Abundance~
                                 s(age_yrs.x, bs="cr")+
                                 matriline.x+
                                 s(hours_after_sunrise_rescaled, bs = "cr") + # specifically interested in this variable
                                 s(Seq_depth_scaled)+
                                 s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                               data=Streptococcus,
                               family = gaussian)

print(summary(Streptococcus_3_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Streptococcus_3_gam) # and day effect (random)
plot(Streptococcus_3_gam)

############ Fusobacterium ###########
Fusobacterium<-subset(core_genera_gam, Genus == "Fusobacterium")

Fusobacterium_3_gam <- mgcv::gam(Abundance~
                                   s(age_yrs.x, bs="cr")+
                                   matriline.x+
                                   s(hours_after_sunrise_rescaled, bs = "cc") + # specifically interested in this variable
                                   s(Seq_depth_scaled)+
                                   s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                                 data=Fusobacterium,
                                 family = gaussian)

print(summary(Fusobacterium_3_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Fusobacterium_3_gam) # and day effect (random)
plot(Fusobacterium_3_gam)


############ Erysipelotrichaceae ###########
Erysipelotrichaceae<-subset(core_genera_gam, Genus == "Erysipelotrichaceae")

Erysipelotrichaceae_gam <- mgcv::gam(Abundance~
                                  s(age_yrs.x, bs="cr")+
                                  matriline.x+
                                  s(hours_after_sunrise_rescaled, bs = "cc") + # specifically interested in this variable
                                  s(Seq_depth_scaled)+
                                  s(hyenaID2.x, hours_after_sunrise_rescaled, bs = "re"),
                                data=Erysipelotrichaceae,
                                family = gaussian)

print(summary(Erysipelotrichaceae_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Erysipelotrichaceae_gam) # and day effect (random)
plot(Erysipelotrichaceae_gam)


##########
sig_genera<-rbind(Clostridiales_XIII, Coprococcus_3, Enterococcus, Ruminococcaceae, Faecalitalea, Peptoclostridium, Peptococcus, Streptococcus)

genus_plot<-ggplot(sig_genera, aes(hours_after_sunrise.x, Abundance, colour=Genus, fill=Genus))+
  stat_smooth(method="gam", alpha=0.2, size=1.5) + 
  scale_colour_manual(values=c("#5669A7", "#8697CE", "#D2C558", "#B3A982", "#6174B0", "#7B8DC5", "#9BABDF","#F2E15F"))+
  scale_fill_manual(values=c("#5669A7", "#8697CE", "#D2C558", "#B3A982", "#6174B0", "#7B8DC5", "#9BABDF","#F2E15F"))+
  theme_bw()+ xlab("hours after sunrise") + 
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))# center align the plot title

  ggarrange(
    gllvm.plot,
    genus_plot, 
    labels = c("A", "B"),
    nrow=1, ncol=2, widths=c(1,2), align = c("h")
  )
  