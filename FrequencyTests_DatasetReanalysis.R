#-----------------------------------------------#
# Re-analysis of Orozco-Terwengel 2012 Dataset  #
#                                               #
# Author: R. Axel W. Wiberg                     #
# Created: 24.01.2016                           #
# Last Modified: 06.02.2017                     #
#-----------------------------------------------#
# Clean environment
rm(list=ls(all = TRUE))
# Set maximum number of decimal positions to a
# high number.
options(scipen = 999)

# Load libraries
# ####
library(ggplot2)
#library(vcd) # Another version of Woolf-test is available from here.
library(dplyr)
library(plyr)
library(reshape)
library(DescTools) # BreslowDayTest from here
library(pscl)
library(qvalue)

#-----------------------------------#
# SOURCE THE CUSTOM FUNCTIONS:      #
# From FrequencyTests_Functions.R   #
#-----------------------------------#
#
source("~/Desktop/Data/RData/RScripts/FrequencyTests_Functions.R")

#------------------------#
# ORIGINAL DATA ANALYSIS #
#------------------------# 
setwd("~/Desktop/Data/RData/allele_frequency_tests/")

# Load original Data (Using the CMH-Test)
# ####
OtW_data<-read.table("BF37.sync",header = FALSE)
colnames(OtW_data)<-c("Chr","pos","ref",
                      "BR1","BR3","BR2",
                      "F15R4","F15R5","F23R1",
                      "F27R5","F37R4","F37R5",
                      "F37R1","cmh_p")

OtW_data$data <- rep("CMH-test",nrow(OtW_data))
# Order by chromosome and position
OtW_data<-OtW_data[order(OtW_data$Chr,OtW_data$pos),]
nrow(OtW_data)
#How many are "HET"
nrow(OtW_data[grep("Het",OtW_data$Chr),])
# Remove the "HET" chromosomes
OtW_data_nh <- OtW_data[grep("Het",OtW_data$Chr,invert=TRUE),]
nrow(OtW_data_nh)
# Remoge the Het label
OtW_data$Chr <- gsub("Het","",OtW_data$Chr)
nrow(OtW_data)
OtW_data_nh$SNP <- paste(OtW_data_nh$Chr,"_",OtW_data_nh$pos,sep = "")


head(OtW_data)
str(OtW_data)
nrow(OtW_data)
# Load table of inclusion SNPs
# SNPs on: X, 2L and 3R were removed for various reasons.
incl_snps <- read.table("OtW_all_included_snps.tab", header = FALSE)
colnames(incl_snps)<-c("Chr","pos","ref","alt","incl")
#order by chromosome and position
incl_snps <- incl_snps[order(incl_snps$Chr,incl_snps$pos),]
head(incl_snps)
str(incl_snps)
nrow(incl_snps)
#How many are "HET"
nrow(incl_snps[grep("Het",incl_snps$Chr),])
#Remove the "Het" tag
incl_snps$Chr <- gsub("Het","",incl_snps$Chr)
incl_snps$SNP <- paste(incl_snps$Chr,"_",incl_snps$pos,sep = "")
length(unique(incl_snps$SNP))

#how many SNPs are excluded
length(OtW_data$SNP[OtW_data$SNP %in% incl_snps$SNP])


# print original top 2000 SNPs
OtW_cmh_top2000<-head(OtW_data_nh[
  order(OtW_data_nh$cmh_p),], n = 2000)

write.table(OtW_cmh_top2000$SNP,file = "~/Desktop/Data/allele_frequency_analysis_project/orozco-terwengel_2012/BF37_cmh_top2000_all.list",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# In format for GOwinda
# Top SNPs
head(OtW_cmh_top2000)
write.table(OtW_cmh_top2000[,c(1,2)],file = "~/Desktop/Data/allele_frequency_analysis_project/orozco-terwengel_2012/BF37_cmh_top2000_all_GO.tab",
            quote = FALSE, col.names = FALSE, row.names = FALSE,
            sep = "\t")


# All SNPs
write.table(OtW_data_nh[,c(1,2)],file = "~/Desktop/Data/allele_frequency_analysis_project/orozco-terwengel_2012/BF37_cmh_allsnps_GO.tab",
            quote = FALSE, col.names = FALSE, row.names = FALSE,
            sep = "\t")
# ####

#
# Original Manhattan Plot: Figure 3 A)
# ####
# Make nice midpoint vectors
OtW_data_nh$absPos <- seq(1:nrow(OtW_data_nh))
OtW_data_nh$Chr <- factor(OtW_data_nh$Chr)
chrNum=length(unique(OtW_data_nh$Chr))
c <- 1
midvec <- vector(length = chrNum)
for (i in unique(OtW_data_nh$Chr)){
  ndx <- which(OtW_data_nh$Chr == i)
  SubPos <- OtW_data_nh$absPos[ndx]
  midvec[c] <-((max(SubPos) - min(SubPos))/2) + min(SubPos)
  c <- c + 1
}
# Get the top 2000 SNPs
OtW_data_top2000<-head(OtW_data_nh[order(OtW_data_nh$cmh_p),],n = 2000)
nrow(OtW_cmh_top2000)
head(OtW_data_nh)
nrow(OtW_data_nh)

manhplot <- ggplot() +
  geom_point(data = OtW_data_nh[
    !(OtW_data_nh$SNP %in% OtW_data_top2000$SNP),], 
             aes(absPos,-log10(cmh_p), 
                 colour = Chr),alpha = 1/5) +
  geom_point(data = OtW_data_top2000, 
             aes(absPos,-log10(cmh_p)),
             colour = "darkblue",alpha = 1/5) +
  scale_colour_manual(limits = as.character(unique(OtW_data_nh$Chr)),
                      values = c(rep(c(
                        "grey50","black"),3),"grey50"), 
                      guide = FALSE) +
  scale_x_continuous(labels=unique(OtW_data_nh$Chr), 
                     breaks = midvec) +
  geom_text(aes(x=500,y=40,label="A)"),size=8)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45)) +
  xlab("") +
  ylab("-log10(CMH p-value)")
manhplot + my.theme  +
  theme(axis.text.x  = element_text(size = 12,face = "bold", 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        plot.title = element_text(hjust=0))

# ####

#----------------------------------------------------------------#
# Print a subset of the sync file for the re-analysis:           #
# So that only the final generation (37) and the base generation #
# are present in the file.                                       #
#----------------------------------------------------------------#
# ####
re_dat <- OtW_data[,c(1,2,3,4,11,5,12,6,13)]
head(re_dat)
nrow(re_dat)
write.table(x = re_dat,file = "BF37_reduced.sync",sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
# ####

#---------------------#
# RE-ANALYSIS OF DATA #
#---------------------# 
#
# Load re-analysis (QB-GLM)
# ####
# t_p = Treatment effect p-value
###
# Raw counts
###
OtW_glm_datar<-read.table("BF37_qbglm_unp_cov10-500.rout",
                          header = FALSE)
colnames(OtW_glm_datar)<-c("Chr","pos","ref",
                          "BR1","F37R4",
                          "BR3","F37R5",
                          "BR2","F37R1",
                          "t_p")

OtW_glm_datar$SNP <- paste(OtW_glm_datar$Chr,"_",OtW_glm_datar$pos,sep = "")
OtW_glm_datar$data <- rep("Raw Counts",nrow(OtW_glm_datar))

# What is the Bonferroni threshhold?
bonf_raw <- 0.05/nrow(OtW_glm_datar)
# ####

###
# Rescaled to neff
# ####
OtW_glm_dataneff<-read.table("BF37_qbglm_unp_neff_cov10-500.rout",
                             header = FALSE)
colnames(OtW_glm_dataneff)<-c("Chr","pos","ref",
                              "BR1","F37R4",
                              "BR3","F37R5",
                              "BR2","F37R1",
                              "t_p")

OtW_glm_dataneff$SNP <- paste(OtW_glm_dataneff$Chr,
                              "_",OtW_glm_dataneff$pos,sep = "")
OtW_glm_dataneff$data <- rep("Scaled (neff)",nrow(OtW_glm_dataneff))

# What is the Bonferroni threshhold?
bonf_neff <- 0.05/nrow(OtW_glm_dataneff)

# ####






###
# Rescaled to 1000
# ####
OtW_glm_data_s1000<-read.table("BF37_qbglm_unp_scale1000_cov10-500.rout",
                         header = FALSE)
colnames(OtW_glm_data_s1000)<-c("Chr","pos","ref",
                          "BR1","F37R4",
                          "BR3","F37R5",
                          "BR2","F37R1",
                          "t_p")

OtW_glm_data_s1000$SNP <- paste(OtW_glm_data_s1000$Chr,
                                "_",OtW_glm_data_s1000$pos,sep = "")
OtW_glm_data_s1000$data <- rep("Scaled (1000)",
                               nrow(OtW_glm_data_s1000))

# What is the Bonferroni threshhold?
bonf_s1000 <- 0.05/nrow(OtW_glm_data_s1000)

# ####

# Rescaled to 100
# ####
OtW_glm_data_s100<-read.table("BF37_qbglm_unp_scale100_cov10-500.rout",
                              header = FALSE)
colnames(OtW_glm_data_s100)<-c("Chr","pos","ref",
                               "BR1","F37R4",
                               "BR3","F37R5",
                               "BR2","F37R1",
                               "t_p")
OtW_glm_data_s100$SNP <- paste(OtW_glm_data_s100$Chr,
                               "_",OtW_glm_data_s100$pos,sep = "")
OtW_glm_data_s100$data <- rep("Scaled (100)",
                              nrow(OtW_glm_data_s100))

# What is the Bonferroni threshhold?
bonf_s100 <- 0.05/nrow(OtW_glm_data_s100)

# ####

#
# Manhattan plot (GLM: all re-analyses combined)
# ####
# Make nice midpoint vectors

# Raw Counts
head(OtW_glm_datar)
OtW_glm_datar$absPos <- seq(1:nrow(OtW_glm_datar))
OtW_glm_datar$Chr <- factor(OtW_glm_datar$Chr)
unique(OtW_glm_datar$Chr)
unique(OtW_data_nh$Chr)
chrNum=length(unique(OtW_glm_datar$Chr))
c <- 1
midvec <- vector(length = chrNum)
for (i in unique(OtW_glm_datar$Chr)){
  ndx <- which(OtW_glm_datar$Chr == i)
  SubPos <- OtW_glm_datar$absPos[ndx]
  midvec[c] <-((max(SubPos) - min(SubPos))/2) + min(SubPos)
  c <- c + 1
}
# Get top 2000 SNPs
OtW_glm_datar_top2000<-head(OtW_glm_datar[
  order(OtW_glm_datar$t_p),], n = 2000)
nrow(OtW_glm_datar_top2000)
# Get SNPs that pass bonferroni threshhold
bonf_raw
OtW_glm_datar_bonf<-OtW_glm_datar[
  OtW_glm_datar$t_p < bonf_raw,]
nrow(OtW_glm_datar_bonf)

# Get SNPs that pass B-H threshhold
# Order tests
OtW_glm_datar_bh<-OtW_glm_datar[order(OtW_glm_datar$t_p),]
# Make B-H column: 
OtW_glm_datar$BH <- (0.05*seq(1,nrow(OtW_glm_datar_bh)))/
  nrow(OtW_glm_datar_bh)
# Which SNPs have t_p <= BH
OtW_glm_datar_bhsig<-OtW_glm_datar_bh[
  OtW_glm_datar_bh$t_p <= OtW_glm_datar_bh$BH,]
nrow(OtW_glm_datar_bhsig)
tail(OtW_glm_datar_bhsig)

# Get SNPs that have q-values < 0.05
# Save t_p as list
write.table(OtW_glm_datar$t_p,"OtW_glm_datar_t_p.list",
            row.names = FALSE,col.names = FALSE,quote=FALSE)
# Read the p-values
pvals<-scan("OtW_glm_datar_t_p.list")
# Convert to q-values with default settings
qobj<-qvalue(pvals)
qplot(qobj)
# Attach q-values to data
OtW_glm_datar$t_q<-qobj$qvalues
# Which SNPs have t_q < 0.05
OtW_glm_datar_q<-OtW_glm_datar[
  OtW_glm_datar$t_q < 0.05,]
nrow(OtW_glm_datar_q)
tail(OtW_glm_datar_q)


# Scale Counts (neff)
OtW_glm_dataneff$absPos <- seq(1:nrow(OtW_glm_dataneff))
OtW_glm_dataneff$Chr <- factor(OtW_glm_dataneff$Chr)
chrNum=length(unique(OtW_glm_dataneff$Chr))
c <- 1
midvec <- vector(length = chrNum)
for (i in unique(OtW_glm_dataneff$Chr)){
  ndx <- which(OtW_glm_dataneff$Chr == i)
  SubPos <- OtW_glm_dataneff$absPos[ndx]
  midvec[c] <-((max(SubPos) - min(SubPos))/2) + min(SubPos)
  c <- c + 1
}
# Get top 2000 SNPs
OtW_glm_dataneff_top2000<-head(OtW_glm_dataneff[
  order(OtW_glm_dataneff$t_p),], n = 2000)
nrow(OtW_glm_dataneff_top2000)
# Get SNPs that pass bonferroni threshhold
bonf_neff
OtW_glm_dataneff_bonf<-OtW_glm_dataneff[
  OtW_glm_dataneff$t_p < bonf_neff,]
nrow(OtW_glm_dataneff_bonf)

# Get SNPs that pass B-H threshhold
# Order tests
OtW_glm_dataneff_bh<-OtW_glm_dataneff[order(OtW_glm_dataneff$t_p),]
# Make B-H column: 
OtW_glm_dataneff_bh$BH <- (0.05*seq(1,nrow(OtW_glm_dataneff_bh)))/
  nrow(OtW_glm_dataneff_bh)
# Which SNPs have t_p <= BH
OtW_glm_dataneff_bh_sig<-OtW_glm_dataneff_bh[
  OtW_glm_dataneff_bh$t_p <= OtW_glm_dataneff_bh$BH,]
nrow(OtW_glm_dataneff_bh_sig)


# Get SNPs that have q-values < 0.05
# Save t_p as list
write.table(OtW_glm_dataneff$t_p,"OtW_glm_dataneff_t_p.list",
            row.names = FALSE,col.names = FALSE,quote=FALSE)
# Read the p-values
pvals<-scan("OtW_glm_dataneff_t_p.list")
# Convert to q-values with default settings
qobj<-qvalue(pvals)
qplot(qobj)
# Attach q-values to data
OtW_glm_dataneff$t_q<-qobj$qvalues
# Which SNPs have t_q < 0.05
OtW_glm_dataneff_q<-OtW_glm_dataneff[
  OtW_glm_dataneff$t_q < 0.05,]
nrow(OtW_glm_dataneff_q)
tail(OtW_glm_dataneff_q)










# Scale Counts (1000)
OtW_glm_data_s1000$absPos <- seq(1:nrow(OtW_glm_data_s1000))
OtW_glm_data_s1000$Chr <- factor(OtW_glm_data_s1000$Chr)
chrNum=length(unique(OtW_glm_data_s1000$Chr))
c <- 1
midvec <- vector(length = chrNum)
for (i in unique(OtW_glm_data_s1000$Chr)){
  ndx <- which(OtW_glm_data_s1000$Chr == i)
  SubPos <- OtW_glm_data_s1000$absPos[ndx]
  midvec[c] <-((max(SubPos) - min(SubPos))/2) + min(SubPos)
  c <- c + 1
}
# Get top 2000 SNPs
OtW_glm_data_s1000_top2000<-head(OtW_glm_data_s1000[
  order(OtW_glm_data_s1000$t_p),], n = 2000)

# Get SNPs that pass bonferroni threshhold
bonf_s1000
OtW_glm_data_s1000_bonf<-OtW_glm_data_s1000[
  OtW_glm_data_s1000$t_p < bonf_s1000,]
nrow(OtW_glm_data_s1000_bonf)
# There are only 33 that pass bonferroni 
# plot the actual allele frequency difference for all of them
# Save the syncfile
write.table(OtW_glm_data_s1000_bonf[,seq(1,9)],
            "OtW_glm_data_s1000_bonf.sync",
            row.names=FALSE,col.names=FALSE,quote=FALSE,
            sep = "\t")

# Get SNPs that pass B-H threshhold
# Order tests
OtW_glm_data_s1000_bh<-OtW_glm_data_s1000[order(OtW_glm_data_s1000$t_p),]
# Make B-H column: 
OtW_glm_data_s1000_bh$BH <- (0.05*seq(1,nrow(OtW_glm_data_s1000_bh)))/
  nrow(OtW_glm_data_s1000_bh)
nrow(OtW_glm_data_s1000_bonf)
head(OtW_glm_data_s1000)
# Which SNPs have t_p <= BH
OtW_glm_data_s1000_bh_sig<-OtW_glm_data_s1000_bh[
  OtW_glm_data_s1000_bh$t_p <= OtW_glm_data_s1000_bh$BH,]
nrow(OtW_glm_data_s1000_bh_sig)
tail(OtW_glm_data_s1000_bh_sig)

# Get SNPs that have q-values < 0.05
# Save t_p as list
write.table(OtW_glm_data_s1000$t_p,"OtW_glm_data_s1000_t_p.list",
            row.names = FALSE,col.names = FALSE,quote=FALSE)
# Read the p-values
pvals<-scan("OtW_glm_data_s1000_t_p.list")
# Convert to q-values with default settings
qobj<-qvalue(pvals)
qplot(qobj)
# Attach q-values to data
OtW_glm_data_s1000$t_q<-qobj$qvalues
# Which SNPs have t_q < 0.05
OtW_glm_data_s1000_q<-OtW_glm_data_s1000[
  OtW_glm_data_s1000$t_q < 0.05,]
nrow(OtW_glm_data_s1000_q)
tail(OtW_glm_data_s1000_q)


# Scale Counts (100)
OtW_glm_data_s100$absPos <- seq(1:nrow(OtW_glm_data_s100))
OtW_glm_data_s100$Chr <- factor(OtW_glm_data_s100$Chr)
chrNum=length(unique(OtW_glm_data_s100$Chr))
c <- 1
midvec <- vector(length = chrNum)
for (i in unique(OtW_glm_data_s100$Chr)){
  ndx <- which(OtW_glm_data_s100$Chr == i)
  SubPos <- OtW_glm_data_s100$absPos[ndx]
  midvec[c] <-((max(SubPos) - min(SubPos))/2) + min(SubPos)
  c <- c + 1
}
# Get top 2000 SNPs
OtW_glm_data_s100_top2000<-head(OtW_glm_data_s100[
  order(OtW_glm_data_s100$t_p),], n = 2000)
# Get SNPs that pass bonferroni threshhold
bonf_s100
OtW_glm_data_s100_bonf<-OtW_glm_data_s100[
  OtW_glm_data_s100$t_p < bonf_s100,]
nrow(OtW_glm_data_s100_bonf)
# Sample X SNPs randomly from the ones that pass 
# bonferroni to plot the actual allele frequency difference
OtW_glm_data_s100_bonf_sub<-OtW_glm_data_s100_bonf[
  sample(seq(1,nrow(OtW_glm_data_s100_bonf)),100),seq(1,9)]
# Save the syncfile
write.table(OtW_glm_data_s100_bonf_sub,
            "OtW_glm_data_s100_bonf_sub.sync",
            row.names=FALSE,col.names=FALSE,quote=FALSE,
            sep = "\t")

# Get SNPs that pass B-H threshhold
# Order tests
OtW_glm_data_s100_bh<-OtW_glm_data_s100[order(OtW_glm_data_s100$t_p),]
# Make B-H column: 
OtW_glm_data_s100_bh$BH <- (seq(1,nrow(OtW_glm_data_s100_bh))/
  nrow(OtW_glm_data_s100_bh))*0.05
# Which SNPs have t_p <= BH
OtW_glm_data_s100_bh_sig<-OtW_glm_data_s100_bh[
  OtW_glm_data_s100_bh$t_p <= OtW_glm_data_s100_bh$BH,]
nrow(OtW_glm_data_s100_bh_sig)
tail(OtW_glm_data_s100_bh_sig)

# Get SNPs that have q-values < 0.05
# Save t_p as list
write.table(OtW_glm_data_s100$t_p,"OtW_glm_data_s100_t_p.list",
            row.names = FALSE,col.names = FALSE,quote=FALSE)
# Read the p-values
pvals<-scan("OtW_glm_data_s100_t_p.list")
# Convert to q-values with default settings
qobj<-qvalue(pvals)
qplot(qobj)
# Attach q-values to data
OtW_glm_data_s100$t_q<-qobj$qvalues
# Which SNPs have t_q < 0.05
OtW_glm_data_s100_q<-OtW_glm_data_s100[
  OtW_glm_data_s100$t_q < 0.05,]
nrow(OtW_glm_data_s100_q)
tail(OtW_glm_data_s100_q)










# Plot data: Figure 3 B) and C)
dat <- rbind(OtW_glm_datar,
             OtW_glm_dataneff)


dat_top2000 <- rbind(OtW_glm_datar_top2000,
                     OtW_glm_dataneff_top2000)

dat_bonf <-rbind(OtW_glm_datar_bonf,
                 OtW_glm_dataneff_bonf)


plot_lab<-data.frame(x=c(500,500),y=c(7,7),
                     lab=c("B)","C)"),
                     data=c("Raw Counts","Scaled (neff)"))

# Draw the plot ####
manhplot_glm <- ggplot() +
  geom_point(data = dat, 
             aes(absPos,-log10(t_p),
                 colour = Chr),alpha = 1/5) +
  geom_point(data = dat_top2000, 
             aes(absPos,-log10(t_p)),
             colour = "darkblue",alpha = 1/5) +
  geom_point(data = dat_bonf, 
             aes(absPos,-log10(t_p)),
             shape=21,colour = "red",size=3) +
  scale_colour_manual(limits = as.character(unique(dat$Chr)),
                      values = c(rep(c(
                        "grey50","black"),3)), 
                      guide = FALSE) +
  scale_x_continuous(labels=unique(dat$Chr), 
                     breaks = midvec) +
  geom_text(data=plot_lab,aes(x=x,y=y,label=lab),size=8)+
  xlab("") +
  ylab("-log10(T p-value)") +
  facet_grid(data~.,scales = "free_y")

manhplot_glm + my.theme +
  theme(axis.text.x  = element_text(size = 12,face = "bold", 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        plot.title = element_text(hjust=0))


# Figure S5
dat <- rbind(OtW_glm_data_s1000,
             OtW_glm_data_s100)

dat_bonf<-rbind(OtW_glm_data_s1000_bonf,
                OtW_glm_data_s100_bonf)

dat_top2000<-rbind(OtW_glm_data_s1000_top2000,
                   OtW_glm_data_s100_top2000)

manhplot_glm <- ggplot() +
  geom_point(data = dat, 
             aes(absPos,-log10(t_p),
                 colour = Chr),alpha = 1/5) +
  geom_point(data = dat_top2000, 
             aes(absPos,-log10(t_p)),
             colour = "darkblue",alpha = 1/5) +
  geom_point(data = dat_bonf, 
             aes(absPos,-log10(t_p)),
             shape=21,colour = "red") +
  scale_colour_manual(limits = as.character(unique(dat$Chr)),
                      values = c(rep(c(
                        "grey50","black"),3)), 
                      guide = FALSE) +
  scale_x_continuous(labels=unique(dat$Chr), 
                     breaks = midvec) +
  xlab("") +
  ylab("-log10(T p-value)") +
  facet_grid(data~.,scales = "free_y")

manhplot_glm + my.theme +
  theme(axis.text.x  = element_text(size = 12,face = "bold", 
                                    colour = "black",
                                    vjust = 0.6,
                                    hjust = 0.5),
        plot.title = element_text(hjust=0))



# Look at some example SNPs
# NEFF: those that pass bonferroni correction
#3R 10295611
matrix<-array(c(25,25,25,9,39,35,39,13,21,33,21,12),dim=c(2,2,3))
data1<-get_dat(matrix,zeroes=1)
data1$snp<-rep("3R_10295611",nrow(data1))
res<-glm(cbind(A_Cnt,Tot_Cnt-A_Cnt)~tr_l,
         family="quasibinomial",data=data1)
n_rows<-nrow(summary(res)$coefficients)
#3R 17760815
matrix<-array(c(41,31,12,0,48,31,14,0,19,31,5,0),dim=c(2,2,3))
data2<-get_dat(matrix,zeroes=1)
data2$snp<-rep("3R_17760815",nrow(data2))
#3R 20122954
matrix<-array(c(36,43,13,0,64,42,23,0,23,43,8,0),dim=c(2,2,3))
data3<-get_dat(matrix,zeroes=1)
data3$snp<-rep("3R_20122954",nrow(data3))

alldata<-rbind(data1,data2,data3)

ggplot()+
  geom_point(data=alldata,aes(tr_l,A_Cnt/Tot_Cnt, colour=rep),
             size = 3,position="jitter")+
  ylim(0,1)+
  xlab("Treatment")+
  ylab("Allele Frequency")+
  facet_grid(snp~.,scales="free_y")+
  my.theme +theme(text=element_text(size=10),
                  axis.text = element_text(size=10),
                  strip.text = element_text(size=10))
# ####