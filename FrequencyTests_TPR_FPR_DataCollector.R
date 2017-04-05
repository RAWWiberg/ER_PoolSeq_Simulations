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
library(ggplot2)
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
#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)
# Set FST and mcov
args <- commandArgs(trailingOnly = TRUE)
fst=args[1]#0.1
mcov=args[2]#200
#-----------------------------------#
# SOURCE THE FUNCTIONS:
# From FrequencyTests_Functions.R
#-----------------------------------#
source("~/Desktop/Data/RData/RScripts/FrequencyTests_Functions.R")
# Set the working directory
setwd("~/Desktop/Data/allele_frequency_analysis_project/data/")

#--------------------------#
# False Positive Rate (FPR)#
#--------------------------#
# Set alpha levels for which FPR should be computed
alpha<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5)

# Initialise a dataframe to store data
fprdat<-data.frame(alpha=vector(),
                   fpr_glm=vector(),
                   fpr_cmh=vector(),
                   fpr_cmh_woo=vector(),
                   fpr_qglm_unp=vector(),
                   fpr_gtest=vector(),
                   fpr_ttest_unp=vector(),
                   scale=vector(),
                   npops=vector())
# Load the data, perform calculations and store the results
# results files are quite big so they are loaded one at a time.
for(k in c(2,3,4,10)){
  ktxt<-paste("k=",k,sep="")
  print(ktxt)
  for(res in c(1,0)){
    print(res)
    if(res == 0){
      k <- as.character(k)
      sim_dat<-read.table(paste(
        "k=",k,
        "_fst=",fst,
	"0.2_N=100_mcov=",mcov,
	"_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
        "_fst=",fst,
	"_N=100_mcov",mcov,
	"_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
        sep=""),
        sep = "\t",header = TRUE)
      sim_dat<-sim_dat[sim_dat$tp == "0",]
      scale<-unique(sim_dat$scale)
      
      # Initialise a temporary internal dataframe for the results
      fpr<-data.frame(alpha=rep(alpha,length(unique(sim_dat$scale))),
                      fpr_glm=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_cmh=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_cmh_woo=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_qglm_unp=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_gtest=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                    
                      fpr_ttest_unp=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      scale=rep(unique(sim_dat$scale),each=length(alpha)),
                      npops=rep(unique(sim_dat$npops),each=length(alpha)))

      # Compute the FDR at different levels of alpha
      for(a in alpha){
        fpr$fpr_cmh[
          fpr$alpha==a & 
            fpr$scale==scale & 
            fpr$npops==ktxt]<-table(sim_dat$cmh_pval<a)["TRUE"]/nrow(sim_dat)
        fpr$fpr_cmh_woo[
          fpr$alpha==a & 
            fpr$scale==scale & 
            fpr$npops==ktxt]<-table(sim_dat$cmh_pval[
              sim_dat$woo_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
                sim_dat$woo_pval > 0.05,])
        fpr$fpr_glm[
          fpr$alpha==a & 
            fpr$scale==scale & 
            fpr$npops==ktxt]<-table(sim_dat$glm_l_pval[
              sim_dat$glm_lxp_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
                sim_dat$glm_lxp_pval > 0.05,])
        fpr$fpr_qglm_unp[
          fpr$alpha==a & 
           fpr$scale==scale & 
           fpr$npops==ktxt]<-table(sim_dat$qglm_unp_l_pval<a)["TRUE"]/nrow(sim_dat)
        fpr$fpr_gtest[
          fpr$alpha==a & 
           fpr$scale==scale & 
           fpr$npops==ktxt]<-table(sim_dat$g_lxc_pval[
             sim_dat$g_lxcxpx_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
               sim_dat$g_lxcxpx_pval > 0.05,])
        fpr$fpr_ttest_unp[
          fpr$alpha==a & 
           fpr$scale==scale & 
           fpr$npops==ktxt]<-table(sim_dat$lm_unp_pval<a)["TRUE"]/nrow(sim_dat)
      }
      fprdat<-rbind(fprdat,fpr)
    }else if(res == 1){
      for(sc in c(100,200,"neff")){
        k <- as.character(k)
        sim_dat<-read.table(paste(
          "k=",k,
          "_fst=",fst,
	  "_N=100_mcov=",mcov,
	  "_res=1_scale=",sc,
          "_SNPs=100000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
          "_fst=",fst,
	  "_N=100_mcov=",mcov,
	  "_res=1_scale=",sc,"_SNPs=100000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
          sep=""),
          sep = "\t",header = TRUE)
        sim_dat<-sim_dat[sim_dat$tp == "0",]
        scale<-unique(sim_dat$scale)

        # Initialise a temporary internal dataframe for the results
        fpr<-data.frame(alpha=rep(alpha,length(unique(sim_dat$scale))),
                        fpr_glm=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_cmh=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_cmh_woo=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_qglm_unp=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_gtest=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_ttest_unp=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        scale=rep(unique(sim_dat$scale),each=length(alpha)),
                        npops=rep(unique(sim_dat$npops),each=length(alpha)))
        # Compute the FDR at different levels of alpha
        for(a in alpha){
          fpr$fpr_cmh[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$cmh_pval<a)["TRUE"]/nrow(sim_dat)
          fpr$fpr_cmh_woo[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$cmh_pval[
                sim_dat$woo_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
                  sim_dat$woo_pval > 0.05,])
          fpr$fpr_glm[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$glm_l_pval[
                sim_dat$glm_lxp_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
                  sim_dat$glm_lxp_pval > 0.05,])
          fpr$fpr_qglm_unp[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$qglm_unp_l_pval<a)["TRUE"]/nrow(sim_dat)
          fpr$fpr_gtest[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$g_lxc_pval[
                sim_dat$g_lxcxpx_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
                  sim_dat$g_lxcxpx_pval > 0.05,])
          fpr$fpr_ttest_unp[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$lm_unp_pval<a)["TRUE"]/nrow(sim_dat)
          }
        fprdat<-rbind(fprdat,fpr)
      }
    }
  }
  
}
rm(sim_dat)

fprdat_m<-melt(data=fprdat,
               measure.vars = c("fpr_cmh",
                                "fpr_cmh_woo",
                                "fpr_glm",
                                "fpr_qglm_unp",
                                "fpr_gtest",
                                "fpr_ttest_unp"))
colnames(fprdat_m)<-c("alpha","scale","npops","test","fpr")
head(fprdat_m)

# Save the data.frame as an R object to avoid having to perform
# calculations again.
#save(list = ls(all=TRUE), file = "fpr.RData",envir=.GlobalEnv)
write.table(fprdat_m,paste("FST="fst,
			   "_mcov=",mcov,
			   "fprdat_melted.tab"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(fprdat,paste("FST="fst,
			   "_mcov=",mcov,
			   "fprdat.tab"),quote=FALSE,row.names=FALSE,sep="\t")
rm(sim_dat,fprdat,fprdat_m)

#-------------------------#
# True Positive Rate (TPR)#
#-------------------------#
# Initialise a dataframe to store data
tprdat<-data.frame(
  npops=vector(length=4),
  scale=vector(length=4),
  cmh_tp_rates=vector(length=4),
  cmh_woo_tp_rates=vector(length=4),
  glm_l_tp_rates=vector(length=4),
  qglm_unp_tp_rates=vector(length=4),
  gtest_tp_rates=vector(length=4),
  ttest_unp_tp_rates=vector(length=4))

# Load the data, perform calculations and store the results
# results files are quite big so they are loaded one at a time.
i <- 1
for(k in c(2,3,4,10)){
  ktxt<-paste("k=",k,sep="")
  print(ktxt)
  for(res in c(1,0)){
    print(res)
    if(res == 0){
      k <- as.character(k)
      sim_dat<-read.table(paste(
        "k=",k,
        "_fst=",fst,
	"_N=100_mcov=",mcov,
	"_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
        "_fst=",fst,
        "_N=100_mcov=",mcov,
	"_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
        sep=""),
        sep = "\t",header = TRUE)
      scale<-unique(sim_dat$scale)
      # TPR = Proportion of true positives that are in the top 1%.
      # CMH-test
      tprdat$cmh_tp_rates[i]<-nrow(sim_dat[
        sim_dat$tp == "1" & 
          sim_dat$cmh_pval < quantile(
            sim_dat$cmh_pval,0.01,na.rm=TRUE),])/nrow(
              sim_dat[sim_dat$tp == "1",])
      
      tprdat$cmh_woo_tp_rates[i]<-nrow(sim_dat[
        sim_dat$tp == "1" & 
          sim_dat$woo_pval > 0.05 & 
          sim_dat$cmh_pval < quantile(
            sim_dat$cmh_pval[
              sim_dat$woo_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
      
      #Bin GLM
      tprdat$glm_l_tp_rates[i]<-nrow(sim_dat[
        sim_dat$glm_lxp_pval > 0.05 & 
          sim_dat$tp == "1" & 
          sim_dat$glm_l_pval < quantile(
            sim_dat$glm_l_pval[
              sim_dat$glm_lxp_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
      
      #Qbin-GLM
      tprdat$qglm_unp_tp_rates[i]<-nrow(sim_dat[
        sim_dat$tp == "1" &
          sim_dat$qglm_unp_l_pval < quantile(
            sim_dat$qglm_unp_l_pval,0.01,na.rm=TRUE),])/nrow(
              sim_dat[sim_dat$tp == "1",])
      
      #G-test
      tprdat$gtest_tp_rates[i]<-nrow(sim_dat[
        sim_dat$g_lxcxpx_pval > 0.05 &
          sim_dat$tp == "1" &
            sim_dat$g_lxc_pval < quantile(
              sim_dat$g_lxc_pval[
                sim_dat$g_lxcxpx_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                  sim_dat[sim_dat$tp == "1",])
      
      #T-test
      tprdat$ttest_unp_tp_rates[i]<-nrow(sim_dat[
        sim_dat$tp == "1" &
          sim_dat$lm_unp_pval < quantile(
            sim_dat$lm_unp_pval,0.01,na.rm=TRUE),])/nrow(
              sim_dat[sim_dat$tp == "1",])
      tprdat$npops[i]=ktxt
      tprdat$scale[i]="CT=var."
      i<-i+1
      } else if(res == 1){
        for (sc in c(100,200,"neff")){
          k <- as.character(k)
          sim_dat<-read.table(paste(
            "k=",k,
            "_fst=",fst,
            "_N=100_mcov=",mcov,
            "_res=1_scale=",sc,"_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
            "_fst=",fst,
            "_N=100_mcov=",mcov,
            "_res=1_scale=",sc,"_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
            sep=""),
            sep = "\t",header = TRUE)
          # TPR = Proportion of true positives that are in the top 1%.
          # CMH-test
          tprdat$cmh_tp_rates[i]<-nrow(sim_dat[
            sim_dat$tp == "1" & 
              sim_dat$cmh_pval < quantile(
                sim_dat$cmh_pval,0.01,na.rm=TRUE),])/nrow(
                  sim_dat[sim_dat$tp == "1",])
        
          tprdat$cmh_woo_tp_rates[i]<-nrow(sim_dat[
            sim_dat$tp == "1" & 
              sim_dat$woo_pval > 0.05 & 
              sim_dat$cmh_pval < quantile(
                sim_dat$cmh_pval[
                  sim_dat$woo_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                    sim_dat[sim_dat$tp == "1",])
        
          tprdat$glm_l_tp_rates[i]<-nrow(sim_dat[
            sim_dat$glm_lxp_pval > 0.05 & 
              sim_dat$tp == "1" & 
              sim_dat$glm_l_pval < quantile(
                sim_dat$glm_l_pval[
                  sim_dat$glm_lxp_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                    sim_dat[sim_dat$tp == "1",])
        
          tprdat$qglm_unp_tp_rates[i]<-nrow(sim_dat[
            sim_dat$tp == "1" &
              sim_dat$qglm_unp_l_pval < quantile(
                sim_dat$qglm_unp_l_pval,0.01,na.rm=TRUE),])/nrow(
                  sim_dat[sim_dat$tp == "1",])
        
          tprdat$gtest_tp_rates[i]<-nrow(sim_dat[
            sim_dat$g_lxcxpx_pval > 0.05 &
              sim_dat$tp == "1" &
                sim_dat$g_lxc_pval < quantile(
                  sim_dat$g_lxc_pval[
                    sim_dat$g_lxcxpx_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                      sim_dat[sim_dat$tp == "1",])
        
          tprdat$ttest_unp_tp_rates[i]<-nrow(sim_dat[
            sim_dat$tp == "1" &
              sim_dat$lm_unp_pval < quantile(
                sim_dat$lm_unp_pval,0.01,na.rm=TRUE),])/nrow(
                  sim_dat[sim_dat$tp == "1",])
          tprdat$npops[i]=ktxt
          tprdat$scale[i]=paste("CT=",sc)
          i<-i+1
        }
      }
  }
}
rm(sim_dat)

# "melt" the data for easier plotting.
tprdat_m <- melt(tprdat,
                 measure.vars = c(
                   "cmh_tp_rates",
                   "cmh_woo_tp_rates",
                   "glm_l_tp_rates",
                   "qglm_unp_tp_rates",
                   "gtest_tp_rates",
                   "ttest_unp_tp_rates")
                 )
colnames(tprdat_m)<-c("npops","scale","test","tpr")

# Save the data.frame as an R object to avoid having to perform
# calculations again.
write.table(tprdat_m,paste("FST="fst,
			   "_mcov=",mcov,
			   "tprdat_melted.tab"),quote=FALSE,row.names=FALSE,sep="\t")
write.table(tprdat,paste("FST="fst,
			   "_mcov=",mcov,
			   "tprdat.tab"),quote=FALSE,row.names=FALSE,sep="\t")
#save(list = ls(all=TRUE), file = "tpr.RData",envir=.GlobalEnv)

