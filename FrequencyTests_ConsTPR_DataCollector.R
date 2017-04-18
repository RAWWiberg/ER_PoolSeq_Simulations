#!/usr/bin/env Rscript
#-------------------------------------------------------#
# R Script for collecting data from simulation results. #
#                                                       #
# Author: R. Axel W. Wiberg                             #
# Created: Apr 2017                                     #
# Last Modified: 05.04.2017                             #
#                                                       # 
#-------------------------------------------------------#

# Clean environment
#rm(list = ls(all = TRUE))

#---------------#
# Load libraries
#---------------#
#library(ggplot2)
#install.packages("vcd")
#library(vcd) # Another version of Woolf-test is available from here.
#install.packages("dplyr")
library(dplyr)
#install.packages("plyr")
library(plyr)
#install.packages("reshape")
library(reshape) # melt() function from here.
#install.packages("DescTools")
library(DescTools) # BreslowDayTest from here
#install.packages("pscl")
library(pscl)
# Set FST and mcov
args <- commandArgs(trailingOnly = TRUE)
fst=as.numeric(args[1])#0.1
mcov=as.numeric(args[2])#200
dir=args[3]
K<-c(4,10)
#-----------------------------------#
# SOURCE THE FUNCTIONS:
# From FrequencyTests_Functions.R
#-----------------------------------#
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
scriptdir<-dirname(script.name)
source(paste(scriptdir,"/FrequencyTests_Functions.R",sep=""))

# Set the working directory
cat("Looking for data in: ",dir,"\n")
setwd(dir)

#----------------------------------------#
# Consistency of True Positive Rate (TPR)#
#----------------------------------------#

# Initialise a dataframe to store data
tprdatcons<-data.frame(
  npops=vector(length=28),
  cmh_tp_rates=vector(length=28),
  qbinglm_unp_tp_rates=vector(length=28),
  lm_unp_tp_rates=vector(length=28),
  sim=vector(length=28),
  snps=vector(length=28))

# Load the data, perform calculations and store the results
# results files are quite big so they are loaded one at a time.
i <- 1
for(k in K){
  ktxt<-paste("k=",k,sep="")
  k <- as.character(k)
  for(snps in c(10000,1000000)){
    if(snps == 10000){
      print(ktxt)
      snptxt<-paste("SNPs=",snps,sep="")
      print(snptxt)
      for(sim in seq(1,10)){
        simtxt<-paste("sim=",sim,sep="")
        print(simtxt)
        sim_dat<-read.table(paste(
          "k=",k,
          "_fst=",fst,
          "_N=100_mcov=",mcov,
          "_res=0_SNPs=10000_p_tp=0.01_sel_diff=0.2_SIM",
          sim,"/FrequencyTest_Simulations_k=",k,
          "_fst=",fst,
          "_N=100_mcov=",mcov,
          "_res=0_SNPs=10000_p_tp=0.01_sel_diff=0.2_SIM",
          sim,".csv",
          sep=""),
          sep = "\t",header = TRUE)
          # TPR = Proportion of true positives that are in the top 1%.
        # CMH-test
        tprdatcons$cmh_tp_rates[i]<-nrow(sim_dat[
          sim_dat$tp == "1" &
            sim_dat$cmh_pval < quantile(
              sim_dat$cmh_pval,0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
        # QB-GLM
        tprdatcons$qbinglm_unp_tp_rates[i]<-nrow(sim_dat[
          sim_dat$tp == "1" &
            sim_dat$qbinglm_unp_l_pval < quantile(
              sim_dat$qbinglm_unp_l_pval,0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
        #LM (t-test)
        tprdatcons$lm_unp_tp_rates[i]<-nrow(sim_dat[
          sim_dat$tp == "1" &
            sim_dat$lm_unp_pval < quantile(
              sim_dat$lm_unp_pval,0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
        
        tprdatcons$npops[i]=ktxt
        
        tprdatcons$sim[i]=simtxt
        
        tprdatcons$snps[i]=snptxt
        
        i<-i+1
      }
    }else if(snps == 1000000){
      snptxt<-paste("SNPs=",snps,sep="")
      print(snptxt)
      for(sim in seq(1,4)){
        simtxt<-paste("sim=",sim,sep="")
        print(simtxt)
        sim_dat<-read.table(paste(
          "k=",k,
          "_fst=",fst,
          "_N=100_mcov=",mcov,
          "_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIM",
          sim,"/FrequencyTest_Simulations_k=",k,
          "_fst=",fst,
          "_N=100_mcov=",mcov,
          "_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIM",
          sim,".csv",
          sep=""),
          sep = "\t",header = TRUE)
        # TPR = Proportion of true positives that are in the top 1%.
        
        # CMH-test
        tprdatcons$cmh_tp_rates[i]<-nrow(sim_dat[
          sim_dat$tp == "1" &
            sim_dat$cmh_pval < quantile(
              sim_dat$cmh_pval,0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
        # QB-GLM
        tprdatcons$qbinglm_unp_tp_rates[i]<-nrow(sim_dat[
          sim_dat$tp == "1" &
            sim_dat$qbinglm_unp_l_pval < quantile(
              sim_dat$qbinglm_unp_l_pval,0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
        # LM (t-test)
        tprdatcons$lm_unp_tp_rates[i]<-nrow(sim_dat[
          sim_dat$tp == "1" &
            sim_dat$lm_unp_pval < quantile(
              sim_dat$lm_unp_pval,0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
        
        tprdatcons$npops[i]=ktxt
        
        tprdatcons$sim[i]=simtxt
        
        tprdatcons$snps[i]=snptxt
        i<-i+1
      }
    }
  }
}
rm(sim_dat)


# Re-order the factor levels for the npops (Nr. Replicates) column
# this makes the plots more intuitive.
tprdatcons$npops <- factor(tprdatcons$npops, 
                           levels = c("k=4","k=10"))
# "melt" the data for easier plotting.
tprdatcons_m <- melt(tprdatcons,
                 measure.vars = c(
                   "cmh_tp_rates",
                   "qbinglm_unp_tp_rates",
                   "lm_unp_tp_rates"))
colnames(tprdatcons_m)<-c("npops","sim","snps","test","tpr")

tprdatcons_m

# Save the data.frame as an R object to avoid having to perform
# calculations again.
write.table(tprdatcons_m,paste("FST=",fst,
			   "_mcov=",mcov,
			   "_tprdatcons_melted.tab"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(tprdatcons,paste("FST=",fst,
			   "_mcov=",mcov,
			   "_tprdatcons.tab"),quote=FALSE,row.names=FALSE,sep="\t")

#save(list = ls(all=TRUE), file = "tprcons.RData",envir=.GlobalEnv)

