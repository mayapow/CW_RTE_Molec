#Maya Powell
#ITS2 processing
#Curacao samples from RTE

setwd("~/Documents/Castillo_Lab/CW_2022/CW_RTE_ITS2")
#load libraries
library(phyloseq)
library('ggplot2')
library('Rmisc')
library(cowplot)
library("ggpubr")
library("vegan")
library(dplyr)
#remotes::install_github("KarstensLab/microshades")
#library("microshades")
#remotes::install_github("mikemc/speedyseq")
library("speedyseq")
library("microViz")

#COMING BACK TO THIS SCRIPT? READ IN PS OBJECTS AT LINE #80

#Check for low read samples
counts <- read.csv('symportal_profilecounts_rte.csv',header=TRUE,row.names=1,check.names=FALSE)
plot(rowSums(counts)) 
counts$sum <- rowSums(counts)
print(counts$sum == "0")
#1 is zero
counts.no0 <- subset(counts, sum != 0) 
print(counts.no0$sum == "0")
counts.no0 <- counts.no0 %>% select(-sum)
plot(rowSums(counts.no0))
#ignore any samples in negative control and just remove for ITS2 profiles and post-med sequences

#BEFORE READING IN SAMPLE DATAFRAME I JUST MANUALLY DELETED THE WATER SAMPLE
#BECAUSE IT HAS 0 READS, and then read in below
#sample dataframe
samdf<-read.csv("symportal_metadata_rte.csv")
head(samdf)
#use this line if it accidentally adds in the X variable as the first column
#samdf <- samdf %>% select(-X)
rownames(samdf) <- samdf$sample_name

#Make all phyloseq objects
#ran this once and then read them in later
#import taxa info
taxa <- read.csv("symportal_taxa.csv",header=TRUE)
rownames(taxa) <- as.factor(taxa$ITS2.type.profile)
mtaxa <- as.matrix(taxa)
# import counts (absolute abundance from its2 type profiles)
mcounts <- as.matrix(counts.no0)

# Construct phyloseq object 
#this one is just the ITS2 profiles
ps <- phyloseq(otu_table(mcounts, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(mtaxa))
ps
#write how many taxa and samples you have here
#33 taxa and 193 samples
#save file to read in later
saveRDS(ps,file="ps.its2.RDS")

#now we are using the all counts table (post-med)
counts.all <- read.csv("symportal_allcounts_rte.csv",header=T,row.names=1)
plot(rowSums(counts.all)) 
#1 zero, removing
counts.all$sum <- rowSums(counts.all)
print(counts.all$sum == "0")
counts.all.no0 <- subset(counts.all, sum != 0) 
print(counts.all.no0$sum == "0")
counts.all.no0 <- counts.all.no0 %>% select(-sum)
plot(rowSums(counts.all.no0)) 

# import counts (absolute abundance from its2 type profiles)
mcounts.all <- as.matrix(counts.all.no0)
# Construct phyloseq object 
ps.all <- phyloseq(otu_table(mcounts.all, taxa_are_rows=FALSE),
                   sample_data(samdf))
ps.all
#9236 taxa and 194 samples
saveRDS(ps.all,"ps.all.its2.RDS")

###Read in phyloseq objects
#ps = ITS2 type profiles
#ps.all = post-med sequences
#XXX = pre-med sequences
ps <- readRDS("ps.its2.RDS")
#ps <- ps %>% ps_mutate(transplant_time = paste(transplant_group, timepoint))
#saveRDS(ps, "ps.its2.RDS")
ps.all <- readRDS("ps.all.its2.RDS")

#make relative abundance dataframe
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))

##Bar plots
#absolute abundance by sample (just to look at them, not need to save)
plot_bar(ps, "full_id", fill="ITS2.type.profile")

#absolute abundance by sample faceted by site and coral species
plot_bar(ps,"ITS2.type.profile", fill="ITS2.type.profile",
         facet_grid=coral_species~reef_bay)

##Type profile by sample##
#all samples
#plot_bar(ps,"sample_name")
gg.bar <- plot_bar(ps.rel,"sample_name",fill="Majority.ITS2.sequence")+
  geom_bar(stat="identity")+
  theme_classic()+
  facet_grid(coral_species~transplant_group, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Site")+
  ylab("Relative Abundance") 
  #theme(axis.text.x=element_blank(),axis.ticks.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar
ggsave(gg.bar,file="sym.barplot.all.png",h=8,w=35)    

#now do this for all datasets (relative) split by Coral.Species
ps.ss.rel <- subset_samples(ps.rel, coral_species=="Siderastrea siderea")
ps.pp.rel <- subset_samples(ps.rel, coral_species=="Porites sp.")

#siderastrea siderea
gg.bar.ss <- plot_bar(ps.ss.rel,"full_id",fill="Majority.ITS2.sequence")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(timepoint~transplant_group, scales = "free")+
  ylab("Relative Abundance")+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.ss
ggsave(gg.bar.ss,file="sym.barplot.ss.transplant.timepoint.png",h=15,w=15)

#porites spp
gg.bar.pp <- plot_bar(ps.pp.rel,"full_id",fill="Majority.ITS2.sequence")+
  geom_bar(stat="identity")+
  theme_classic(base_size=22)+
  facet_wrap(timepoint~transplant_group, scales = "free")+
  ylab("Relative Abundance")+ 
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank()) #for some reason this won't just add to the original plot
gg.bar.pp
ggsave(gg.bar.pp,file="sym.barplot.pp.transplant.timepoint.png",h=15,w=25)

#put these graphs together
#gg.panel.spp <- ggarrange(gg.bar.ss,gg.bar.pp,nrow=1,ncol=3,labels="AUTO")
#gg.panel.spp
#ggsave(gg.panel.spp,file="bac.div.site.png",height=8)

###MP December 2nd, 2024

##########################################################
#### Summary of Dominant ITS2 Majority Types and DIVs ####
##########################################################
#based off code by Hannah Aichelman 
#https://github.com/hannahaichelman/CrypticCorals/tree/main/ITS2_Symbiodiniaceae
#Ander please look through this code at this link

#just majority ITS2 types
ps <- readRDS("ps.its2.RDS")
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
seq <- data.frame(ps.rel@otu_table)
seq.test <- seq 
#making the test just to check to make sure it doesn't change the identities of the taxa bc I'm ~paranoid~ hehe
samdf <- data.frame(ps.rel@sam_data)

# change column names to majority its2 sequence
sym_taxa = read.csv(file = "symportal_taxa.csv", header = TRUE) 

sym_taxa$Majority.ITS2.sequence #you need to do it this way because you need to add the 1s etc to not have duplicates
colnames(seq) =  c("A4.1","A4.2","A4.3","A1.1","A1.2","A1.3","A4.4",
                   "B19.1","B19.2","B19.3","B19.4","C3.1","C1.1",
                   "C47a","C3.2","C3.3","C1.2","C1.3","C45a.1",
                   "C42.2/C45a","C1.4","C1.5","C1.6","C45a.2","C3.4",
                   "C1.7","C3.5","C1.8","D1.1","D1.2","D1.3","D1.4","D1.5")

# make new data frame and sum columns with the same majority its2 sequence
seqtab.rel.sums.rte <- seq %>%
  mutate(A1_sum = rowSums(select(., starts_with("A1")))) %>%
  mutate(A4_sum = rowSums(select(., starts_with("A4")))) %>%
  mutate(B19_sum = rowSums(select(., starts_with("B19")))) %>%
  mutate(C1_sum = rowSums(select(., starts_with("C1")))) %>%
  mutate(C3_sum = rowSums(select(., starts_with("C3")))) %>%
  mutate(C42.2_sum = rowSums(select(., starts_with("C42")))) %>%
  mutate(C45a_sum = rowSums(select(., starts_with("C45a")))) %>%
  mutate(C47a_sum = C47a) %>%
  mutate(D1_sum = rowSums(select(., starts_with("D1")))) %>%
  rownames_to_column(var = "id") %>%
  select(id, contains("_sum"))

its2.sums.rel.rte <- left_join(seqtab.rel.sums.rte, samdf.16s, by = "id")

#use seqtab rel sum to make graphs
#"A4","B19","B2","C46","C3","C47a","C1","C45","C42","D1"
taxa.maj <- data.frame(colnames(seqtab.rel.sums.rte))
taxa.maj <- data.frame(taxa.maj[-c(1), ])
colnames(taxa.maj)[1] <- 'maj_its2'
taxa.maj$genus <- c("A","A","B","C","C","C","C","C","D")
rownames(taxa.maj) <- taxa.maj$maj_its2

taxa.maj$maj_its2 = as.factor(taxa.maj$maj_its2)
taxa.maj$genus = as.factor(taxa.maj$genus)

#seqtab rel sums
rownames(seqtab.rel.sums.rte) <- seqtab.rel.sums.rte$id
seqtab.rel.sums.rte.maj <- seqtab.rel.sums.rte %>% select(-id)
#col names of seqtab.rel.sums.rte.maj match rownames of taxa.maj
#rownames of seqtab.rel.sums.rte.maj match rownames of samdf

#make ps object of majority type sums
ps.its2.rte.maj <- phyloseq(sample_data(samdf),
                        otu_table(seqtab.rel.sums.rte.maj,taxa_are_rows=FALSE),
                        tax_table(as.matrix(taxa.maj)))
ps.its2.rte.maj
saveRDS(ps.its2.rte.maj,file = "ps.its2.rte.maj")
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 9 taxa and 193 samples ]:
#   sample_data() Sample Data:        [ 193 samples by 20 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 9 taxa by 2 taxonomic ranks ]:
#   taxa are columns
ps.its2.rte.maj <- readRDS("ps.its2.rte.maj")

plot_bar(ps.its2.rte.maj, x="id",fill="maj_its2")+
  theme_classic()

its2_colors = c("A1_sum" = "#222255","A4_sum" = "#ffaabb","B19_sum" =  "#aaaa00",
                "C1_sum" = "#eedd88","C3_sum" = "#77aadd",
                "C42.2_sum" = "#225555","C45a_sum" = "#44bb99",
                "C47a_sum" = "#99ddff","D1_sum" = "#994455")

ps.ss.maj <- subset_samples(ps.its2.rte.maj, coral_species=="Siderastrea siderea")
ps.pp.maj <- subset_samples(ps.its2.rte.maj, coral_species=="Porites sp.")

###SIDERASTREA SIDEREA###
ss.site <- plot_bar(ps.ss.maj, x="parent", fill="maj_its2") +
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~transplant_time, scales = "free")+
  scale_fill_manual(name = "Majority ITS2", values = its2_colors) +
  theme(axis.text.x=element_text(angle = 90))+
  #theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle("Siderastrea siderea")
ss.site
ggsave(ss.site,file="its2.ss.site.time.pdf",h=15,w=25)
ggsave(ss.site,file="its2.ss.site.time.labels.pdf",h=20,w=25) 

###BRANCHING PORITES SP###
pp.site <- plot_bar(ps.pp.maj, x="parent", fill="maj_its2") +
  theme_classic(base_size = 30)+
  ylab("Relative Abundance")+
  facet_wrap(~transplant_time, scales = "free")+
  scale_fill_manual(name = "Majority ITS2", values = its2_colors) +
  #theme(axis.text.x=element_text(angle = 90))+
  theme(axis.text.x=element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank())+
  ggtitle("Branching Porites sp.")
pp.site
ggsave(pp.site,file="its2.pp.site.time.pdf",h=15,w=25) 
ggsave(pp.site,file="its2.pp.site.time.labels.pdf",h=20,w=25) 

#convert to factors and numeric as needed for different variables
its2.sums.rel.rte = its2.sums.rel.rte %>%
  mutate_at(c(2:10), as.numeric)

its2.sums.rel.rte = its2.sums.rel.rte %>%
  mutate_at(c(11:26), as.factor)

# add in dominant and minor distinctions to use in other plots
its2.dom.rel.rte <- its2.sums.rel.rte %>%
  mutate(dominant_type = case_when(A1_sum >= 0.7 ~ "A1",
                                   A4_sum >= 0.7 ~ "A4",
                                   B19_sum >= 0.7 ~ "B19",
                                   C1_sum >= 0.7 ~ "C1",
                                   C3_sum >= 0.7 ~ "C3",
                                   C42.2_sum >= 0.7 ~ "C42.2",
                                   C45a_sum >= 0.7 ~ "C45a",
                                   C47a_sum >= 0.7 ~ "C47a",
                                   D1_sum >= 0.7 ~ "D1")) %>%
  mutate(minor_type = case_when(A1_sum < 0.7 & A1_sum > 0.0 ~ "A1",
                                A4_sum < 0.7 & A4_sum > 0.0 ~ "A4",
                                B19_sum < 0.7 & B19_sum > 0.0 ~ "B19",
                                C1_sum < 0.7 & C1_sum > 0.0 ~ "C1",
                                C3_sum < 0.7 & C3_sum > 0.0 ~ "C3",
                                C42.2_sum < 0.7 & C42.2_sum > 0.0 ~ "C42.2",
                                C45a_sum < 0.7 & C45a_sum > 0.0 ~ "C45a",
                                C47a_sum < 0.7 & C47a_sum > 0.0 ~ "C47a",
                                D1_sum < 0.7 & D1_sum > 0.0 ~ "D1"))
#make dominant type dataframe
write.csv(its2.dom.rel.rte, file = "ITS2.dominanttype.CW.RTE.csv", row.names = FALSE)
its2.dom.rel.rte <- read.csv("ITS2.dominanttype.CW.RTE.csv")
its2.dom.rel.rte <- its2.dom.rel.rte %>% 
  mutate(treatment = case_when(transplant_group == "Bay Initial" ~ "native",
                               transplant_group == "Reef Initial" ~ "native",
                               transplant_group == "Bay Native" ~ "native",
                               transplant_group == "Reef Native" ~ "native",
                               transplant_group == "Bay Transplant" ~ "transplant",
                               transplant_group == "Reef Transplant" ~ "transplant")
)

ss.its2.dom.rte <- filter(its2.dom.rel.rte, coral_species=="Siderastrea siderea")
pp.its2.dom.rte <- filter(its2.dom.rel.rte, coral_species=="Porites sp.")

# summarize proportion of individuals with more than 70% of the sym types
maj_types_spp <- its2.dom.rel.rte %>%
  group_by(coral_species) %>%
  dplyr::summarize(n = n(),
            #n_maj_A1 = sum(A1_sum >= 0.7),
            n_maj_A4 = sum(A4_sum >= 0.7),
            #n_maj_B19 = sum(B19_sum >= 0.7),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C42.2 = sum(C42.2_sum >= 0.7),
            n_maj_C45a = sum(C45a_sum >= 0.7),
            n_maj_C47a = sum(C47a_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            #p_maj_A1 = n_maj_A1/n,
            p_maj_A4 = n_maj_A4/n,
            #p_maj_B19 = n_maj_B19/n,
            p_maj_C1 = n_maj_C1/n,
            p_maj_C3 = n_maj_C3/n,
            p_maj_C42.2 = n_maj_C42.2/n,
            p_maj_C45a = n_maj_C45a/n,
            p_maj_C47a = n_maj_C47a/n,
            p_maj_D1 = n_maj_D1/n)
#based on this - removing samples from the species that do not have any majority of certain sym types
#pp keep = A4, C3, C42.2, C45a, C47a
#ss keep = C1, C3, D1
#all remove = A1, B19

#####CORAL SPECIES#####
View(maj_types_spp)

####SSID######
#transplant_initial_site (aka native site)
#transplant_final_site (aka transplant site)
#timepoint
#transplant_group_sampled (where the coral was sampled)
#transplant_time (transplant and timepoint)

#native site
ss.its2.dom.rte %>%
  group_by(transplant_initial_site) %>%
  dplyr::summarize(n = n(),
            n_maj_C1 = sum(C1_sum >= 0.7),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_D1 = sum(D1_sum >= 0.7),
            p_maj_C1 = n_maj_C1/n,
            p_maj_C3 = n_maj_C3/n,
            p_maj_D1 = n_maj_D1/n)

####BRANCHING PORITES SP######
#A4, C3, C42.2, C45a, C47a
#transplant_initial_site (aka native site)
#transplant_final_site (aka transplant)
#timepoint
#transplant_group_sampled (where the coral was sampled)
#transplant_time (transplant timepoint)
pp.its2.dom.rte %>%
  group_by(transplant_time) %>%
  summarize(n = n(),
            n_maj_A4 = sum(A4_sum >= 0.7),
            n_maj_C3 = sum(C3_sum >= 0.7),
            n_maj_C42.2 = sum(C42.2_sum >= 0.7),
            n_maj_C45a = sum(C45a_sum >= 0.7),
            n_maj_C47a = sum(C47a_sum >= 0.7),
            p_maj_A4 = n_maj_A4/n,
            p_maj_C3 = n_maj_C3/n,
            p_maj_C42.2 = n_maj_C42.2/n,
            p_maj_C45a = n_maj_C45a/n,
            p_maj_C47a = n_maj_C47a/n)

#####STATS#####

###Trying multinomial models using nnet

library(nnet,package_version())
library(nnet); packageVersion("nnet")
library(car); packageVersion("car")
library(emmeans); packageVersion("emmeans")

#coral species
spp_model <-multinom(maj_ITS2 ~ coral_species, data = its2.dom.rel.rte)
Anova(spp_model)

spp_em <- emmeans(spp_model, ~ coral_species | maj_ITS2)
pairs(spp_em)

#siderastrea siderea
ss_null <- multinom(maj_ITS2 ~ 1, data = ss.its2.dom.rte)
ss_origin <- multinom(maj_ITS2 ~ native_site, data = ss.its2.dom.rte)
ss_treat <- multinom(maj_ITS2 ~ transplant_site, data = ss.its2.dom.rte)
ss_time <- multinom(maj_ITS2 ~ timepoint, data = ss.its2.dom.rte)
ss_add <- multinom(maj_ITS2 ~ native_site + transplant_site, data = ss.its2.dom.rte)
ss_ot <- multinom(maj_ITS2 ~ native_site * transplant_site, data = ss.its2.dom.rte)
ss_model <- multinom(maj_ITS2 ~ timepoint + native_site + transplant_site, data = ss.its2.dom.rte)
ss_group <- multinom(maj_ITS2 ~ timepoint + native_site * transplant_site, data = ss.its2.dom.rte)
ss_model_int <- multinom(maj_ITS2 ~ timepoint * native_site * transplant_site, data = ss.its2.dom.rte)
#interactive model

#summary(ss_model)
anova(ss_null, ss_origin, ss_treat, ss_time, ss_add, ss_ot, ss_model, ss_group, ss_model_int)
AIC(ss_null, ss_origin, ss_treat, ss_time, ss_add, ss_ot, ss_model, ss_group, ss_model_int)

Anova(ss_ot, type = 3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: maj_ITS2
#                       LR Chisq Df Pr(>Chisq)    
# native_site             38.192  2  5.091e-09 ***
# treatment               32.725  2  7.831e-08 ***
# native_site:treatment   25.493  2  2.913e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

ss_t_em <- emmeans(ss_ot, ~ transplant_site | maj_ITS2)
pairs(ss_t_em)

#branching porites sp.
pp_null <- multinom(maj_ITS2 ~ 1, data = pp.its2.dom.rte)
pp_origin <- multinom(maj_ITS2 ~ native_site, data = pp.its2.dom.rte)
pp_treat <- multinom(maj_ITS2 ~ transplant_site, data = pp.its2.dom.rte)
pp_time <- multinom(maj_ITS2 ~ timepoint, data = pp.its2.dom.rte)
pp_add <- multinom(maj_ITS2 ~ native_site + transplant_site, data = pp.its2.dom.rte)
pp_ot <- multinom(maj_ITS2 ~ native_site * transplant_site, data = pp.its2.dom.rte)
pp_model <- multinom(maj_ITS2 ~ timepoint + native_site + transplant_site, data = pp.its2.dom.rte)
pp_group <- multinom(maj_ITS2 ~ timepoint + native_site * transplant_site, data = pp.its2.dom.rte)
pp_model_int <- multinom(maj_ITS2 ~ timepoint * native_site * transplant_site, data = pp.its2.dom.rte)
#interactive model

#summary(pp_model)
anova(pp_null, pp_origin, pp_treat, pp_time, pp_add, pp_ot, pp_model, pp_group, pp_model_int)
AIC(pp_null, pp_origin, pp_treat, pp_time, pp_add, pp_ot, pp_model, pp_group, pp_model_int)

Anova(pp_ot)
# Response: maj_ITS2
#                       LR Chisq Df Pr(>Chisq)    
# native_site             73.669  4  3.809e-15 ***
# treatment               68.216  4  5.400e-14 ***
# native_site:treatment    0.000  4          1    

pp_t_em <- emmeans(pp_ot, ~ native_site:transplant_site | maj_ITS2)
pairs(pp_t_em)


####STATS FROM CAS - trying CLMMs###

library("ordinal")
library(MASS)

#notes from meeting with Cas Grupstra
#ITS2 profiles - assigned to sym type, no cutoff, just did analysis on C15
#UID = sym type
#clmm = mixed model for ordinal data
#clmm UID ~ lineage*type + (1|site), data = ordinal)
#Did full to null models
#ordinal::anova.clmm()
#differences in counts
#microsatellite data - just don't include in the analysis for now
#majority 0 could define as the larger proportion
#do it both ways (majority & combined) and put one in supplement
#separate analysis - which ones have the background abdunaces/are more mixed
#dominated vs mixed analysis as well? - run it and if it's interesting then dig in
#Raw ITS-2
#reads were processed by Symportal (Hume et al. 2019) to
#produce defining intragenomic sequence variant (DIV) profiles
#for each coral colony (final n = 73; Table S1). All colonies
#were dominated by C15 DIVs. Cumulative link mixed models
#(clmm) were used to test for differences in dominant C15
#DIV associations. A null model was generated that included
#only site as a random effect, and a full model included an interaction
#between lineage and reef type, as well as site as a
#random effect. Two additional models included reef type or
#lineage as fixed effects, with site as a random effect. The function
#anova.clmm in ordinal v2022.11.16 was used for model
#selection, and Anova.clmm in rcompanion v.2.4.21 assessed
#effect sizes.

its2.dom.rel.rte <- read.csv("ITS2.dominanttype.CW.RTE.csv")

ss.its2.dom.rte <- filter(its2.dom.rel.rte, coral_species=="Siderastrea siderea")
pp.its2.dom.rte <- filter(its2.dom.rel.rte, coral_species=="Porites sp.")
write.csv(ss.its2.dom.rte, file = "ss.ITS2.dominanttype.CW.RTE.csv", row.names = FALSE)
write.csv(pp.its2.dom.rte, file = "pp.ITS2.dominanttype.CW.RTE.csv", row.names = FALSE)

#columns to use:
#timepoint
#native_site
#transplant_site
#parent
#maj_its2

#For my models do:
#cumulative link mixed models (LMM) were used to test for the effects of origin (bay or reef), treatment (native or transplant), and timepoint (T0, T4, T12)
#origin, treatment, and timepoint
#Random effect of genotype (nested within origin)

ftable(xtabs(~ maj_ITS2 + native_site + transplant_site + timepoint, data = ss.its2.dom.rte))

#siderastrea siderea
a <- clmm(factor(maj_ITS2) ~(1|parent), data = ss.its2.dom.rte)
b <- clmm(factor(maj_ITS2) ~native_site + (1|parent), data = ss.its2.dom.rte)
c <- clmm(factor(maj_ITS2) ~transplant_site+(1|parent), data = ss.its2.dom.rte)
d <- clmm(factor(maj_ITS2) ~timepoint+(1|parent), data = ss.its2.dom.rte)
e <- clmm(factor(maj_ITS2) ~transplant_site+native_site+ (1|parent), data = ss.its2.dom.rte)
f <- clmm(factor(maj_ITS2) ~timepoint+native_site+(1|parent), data = ss.its2.dom.rte)
g <- clmm(factor(maj_ITS2) ~timepoint+transplant_site+(1|parent), data = ss.its2.dom.rte)
h <- clmm(factor(maj_ITS2) ~timepoint+transplant_site+native_site+(1|parent), data = ss.its2.dom.rte)
i <- clmm(factor(maj_ITS2) ~timepoint*transplant_site+native_site+(1|parent), data = ss.its2.dom.rte)
j <- clmm(factor(maj_ITS2) ~timepoint+transplant_site*native_site+(1|parent), data = ss.its2.dom.rte)
k <- clmm(factor(maj_ITS2) ~transplant_site+timepoint*native_site+(1|parent), data = ss.its2.dom.rte)
l <- clmm(factor(maj_ITS2) ~transplant_site*native_site*timepoint+(1|parent), data = ss.its2.dom.rte)

anova(a,h,j,l)

# Likelihood ratio tests of cumulative link models:
#  
#   formula:                                                                    link: threshold:
# a factor(maj_ITS2) ~ (1 | parent)                                             logit flexible  
# b factor(maj_ITS2) ~ native_site + (1 | parent)                               logit flexible  
# d factor(maj_ITS2) ~ timepoint + (1 | parent)                                 logit flexible  
# f factor(maj_ITS2) ~ timepoint + native_site + (1 | parent)                   logit flexible  
# g factor(maj_ITS2) ~ timepoint + transplant_site + (1 | parent)               logit flexible  
# h factor(maj_ITS2) ~ timepoint + transplant_site + native_site + (1 | parent) logit flexible  
# j factor(maj_ITS2) ~ timepoint + transplant_site * native_site + (1 | parent) logit flexible  
# i factor(maj_ITS2) ~ timepoint * transplant_site + native_site + (1 | parent) logit flexible  
# k factor(maj_ITS2) ~ transplant_site + timepoint * native_site + (1 | parent) logit flexible  
# l factor(maj_ITS2) ~ transplant_site * native_site * timepoint + (1 | parent) logit flexible  
# 
#   no.par    AIC  logLik  LR.stat df Pr(>Chisq)    
# a      3 197.76 -95.882                           
# b      4 199.28 -95.639   0.4878  1    0.48490    
# d      5 195.22 -92.612   6.0532  1    0.01388 *  
# f      6 196.90 -92.449   0.3252  1    0.56848    
# g      6 180.57 -84.287  16.3250  0               
# h      7 180.68 -83.341   1.8907  1    0.16912    
# j      8 159.13 -71.564  23.5543  1  1.214e-06 ***
# i      9 184.61 -83.307 -23.4853  1    1.00000    
# k      9 184.63 -83.314  -0.0148  0               
# l     12 167.02 -71.508  23.6126  3  3.009e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#model selection
anova(a,h)
#   no.par    AIC  logLik LR.stat df Pr(>Chisq)    
# a      3 197.76 -95.882                          
# h      7 180.68 -83.341  25.082  4  4.843e-05 ***

#effect of random
anova(h,hnr)
hnr <- clm(factor(maj_ITS2) ~timepoint+transplant_site+native_site, data = ss.its2.dom.rte)
#     no.par    AIC  logLik LR.stat df Pr(>Chisq)    
# hnr      6 191.90 -89.949                          
# h        7 180.68 -83.341  13.214  1  0.0002778 ***
  
Anova(h, type = 2)
summary(h)

# Analysis of Deviance Table (Type III tests)
# 
# Response: factor(maj_ITS2)
#                 Df   Chisq Pr(>Chisq)    
# timepoint        2  8.9007  0.0116746 *  
# transplant_site  1 14.6096  0.0001322 ***
# native_site      1  1.8075  0.1788039    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

hnr <- clm(factor(maj_ITS2) ~timepoint+transplant_site+native_site, data = ss.its2.dom.rte)
h_pred <- clmm2(factor(maj_ITS2) ~ timepoint+transplant_site+native_site, random = parent, Hess = TRUE, data = ss.its2.dom.rte)
summary(CLMM2_End_S2_min)

marginal = emmeans(glm.3e.1, ~ m_y*reef_bay*area)

summary(g)
Anova(g, type = 3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: factor(maj_ITS2)
#                       Df   Chisq Pr(>Chisq)    
# timepoint              2  8.9393  0.0114513 *  
# transplant_final_site  1 13.2422  0.0002737 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#branching porites sp
#getting convergence errors due to low sample sizes because of C
#need to collapse all C to C

pp.its2.dom.rte <- pp.its2.dom.rte %>% 
  mutate(maj_ITS2_por= case_when(maj_ITS2 == "A4" ~ "A4",
                                 maj_ITS2 == "C3" ~ "C",
                                 maj_ITS2 == "C42.2" ~ "C",
                                 maj_ITS2 == "C45a" ~ "C",
                                 maj_ITS2 == "C47a" ~ "C",
                                 ))


ftable(xtabs(~ maj_ITS2 + native_site + transplant_site + timepoint, data = pp.its2.dom.rte))

z <- clmm(factor(maj_ITS2) ~(1|parent), data = pp.its2.dom.rte)
y <- clmm(factor(maj_ITS2) ~native_site + (1|parent), data = pp.its2.dom.rte)
x <- clmm(factor(maj_ITS2) ~transplant_site+(1|parent), data = pp.its2.dom.rte)
w <- clmm(factor(maj_ITS2) ~timepoint+(1|parent), data = pp.its2.dom.rte)
v <- clmm(factor(maj_ITS2) ~transplant_site+native_site+ (1|parent), data = pp.its2.dom.rte)
u <- clmm(factor(maj_ITS2) ~timepoint+native_site+(1|parent), data = pp.its2.dom.rte)
t <- clmm(factor(maj_ITS2) ~timepoint+transplant_site+(1|parent), data = pp.its2.dom.rte)
s <- clmm(factor(maj_ITS2) ~timepoint+transplant_site+native_site+(1|parent), data = pp.its2.dom.rte)
r <- clmm(factor(maj_ITS2) ~timepoint*transplant_site+native_site+(1|parent), data = pp.its2.dom.rte)
q <- clmm(factor(maj_ITS2) ~timepoint+transplant_site*native_site+(1|parent), data = pp.its2.dom.rte)
p <- clmm(factor(maj_ITS2) ~transplant_site+timepoint*native_site+(1|parent), data = pp.its2.dom.rte)
o <- clmm(factor(maj_ITS2) ~transplant_site*native_site*timepoint+(1|parent), data = pp.its2.dom.rte)

z <- clmm(factor(maj_ITS2) ~(1|parent), data = pp.its2.dom.rte)
s <- clmm(factor(maj_ITS2) ~timepoint+transplant_site+native_site+(1|parent), data = pp.its2.dom.rte)

anova(a,b,d,f,g,h,i,j,k,l)
anova(h)
#same here for porites for AIC - g is best, then i and j

Anova(h, type = 3)

#Cas code:
Sym_ordinal  <- read.csv("Sym_table_ordinal.csv")

ftable(xtabs(~ UID + Type + Lineage, data = Sym_ordinal))

?clm()
#type here is reef type (extreme vs normal)
x <- clmm(factor(UID) ~Lineage*Type+(1|Site), data = Sym_ordinal)
y <- clmm(factor(UID) ~Lineage+Type+(1|Site), data = Sym_ordinal)
z <- clmm(factor(UID) ~Type+(1|Site), data = Sym_ordinal)
aa <- clmm(factor(UID) ~Lineage+(1|Site), data = Sym_ordinal)
ab <- clmm(factor(UID) ~(1|Site), data = Sym_ordinal)

