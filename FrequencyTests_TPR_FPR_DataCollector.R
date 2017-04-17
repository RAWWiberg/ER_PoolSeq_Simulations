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
K<-c(2,3,4,10)
scales<-c(100,200,"neff")
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

#--------------------------#
# False Positive Rate (FPR)#
#--------------------------#
# Set alpha levels for which FPR should be computed
alpha<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5)

# Initialise a dataframe to store data
fprdat<-data.frame(alpha=vector(),
                   fpr_binglm_ni=vector(),
                   fpr_binglm=vector(),
                   fpr_cmh=vector(),
                   fpr_cmh_woo=vector(),
                   fpr_qbinglm_unp=vector(),
                   fpr_gtest=vector(),
                   fpr_lm_unp=vector(),
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
	"_N=100_mcov=",mcov,
	"_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
        "_fst=",fst,
	"_N=100_mcov=",mcov,
	"_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
        sep=""),
        sep = "\t",header = TRUE)
      sim_dat<-sim_dat[sim_dat$tp == "0",]
      scale<-unique(sim_dat$scale)
      
      # Initialise a temporary internal dataframe for the results
      fpr<-data.frame(alpha=rep(alpha,length(unique(sim_dat$scale))),
                      fpr_binglm=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_binglm_ni=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_cmh=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_cmh_woo=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_qbinglm_unp=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                      
                      fpr_gtest=vector(
                        length=length(alpha)*length(unique(sim_dat$scale))),
                    
                      fpr_lm_unp=vector(
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
        
        fpr$fpr_binglm[
          fpr$alpha==a & 
            fpr$scale==scale & 
            fpr$npops==ktxt]<-table(sim_dat$binglm_l_pval[
              sim_dat$binglm_lxp_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
                sim_dat$binglm_lxp_pval > 0.05,])
        
        fpr$fpr_binglm_ni[
          fpr$alpha==a & 
            fpr$scale==scale & 
            fpr$npops==ktxt]<-table(sim_dat$binglm_l_p_pval_ni<a)["TRUE"]/nrow(sim_dat)
        
        fpr$fpr_qbinglm_unp[
          fpr$alpha==a & 
           fpr$scale==scale & 
           fpr$npops==ktxt]<-table(sim_dat$qbinglm_unp_l_pval<a)["TRUE"]/nrow(sim_dat)
        
        fpr$fpr_gtest[
          fpr$alpha==a & 
           fpr$scale==scale & 
           fpr$npops==ktxt]<-table(sim_dat$g_lxc_pval[
             sim_dat$g_lxcxpx_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
               sim_dat$g_lxcxpx_pval > 0.05,])
       
         fpr$fpr_lm_unp[
          fpr$alpha==a & 
           fpr$scale==scale & 
           fpr$npops==ktxt]<-table(sim_dat$lm_unp_pval<a)["TRUE"]/nrow(sim_dat)
      }
      fprdat<-rbind(fprdat,fpr)
    }else if(res == 1){
      for(sc in scales){
        k <- as.character(k)
        sim_dat<-read.table(paste(
          "k=",k,
          "_fst=",fst,
	  "_N=100_mcov=",mcov,
	  "_res=1_scale=",sc,
          "_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
          "_fst=",fst,
	  "_N=100_mcov=",mcov,
	  "_res=1_scale=",sc,"_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
          sep=""),
          sep = "\t",header = TRUE)
        sim_dat<-sim_dat[sim_dat$tp == "0",]
        scale<-unique(sim_dat$scale)

        # Initialise a temporary internal dataframe for the results
        fpr<-data.frame(alpha=rep(alpha,length(unique(sim_dat$scale))),
                        
                        fpr_binglm=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_binglm_ni=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_cmh=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_cmh_woo=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_qbinglm_unp=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_gtest=vector(
                          length=length(alpha)*length(unique(sim_dat$scale))),
                        
                        fpr_lm_unp=vector(
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
          
          fpr$fpr_binglm[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$binglm_l_pval[
                sim_dat$binglm_lxp_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
                  sim_dat$binglm_lxp_pval > 0.05,])
          
          fpr$fpr_binglm_ni[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$binglm_l_p_pval_ni<a)["TRUE"]/nrow(sim_dat)
          
          fpr$fpr_qbinglm_unp[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$qbinglm_unp_l_pval<a)["TRUE"]/nrow(sim_dat)
          
          fpr$fpr_gtest[
            fpr$alpha==a & 
              fpr$scale==scale & 
              fpr$npops==ktxt]<-table(sim_dat$g_lxc_pval[
                sim_dat$g_lxcxpx_pval > 0.05]<a)["TRUE"]/nrow(sim_dat[
                  sim_dat$g_lxcxpx_pval > 0.05,])
          
          fpr$fpr_lm_unp[
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
                                "fpr_binglm",
                                "fpr_binglm_ni",
                                "fpr_qbinglm_unp",
                                "fpr_gtest",
                                "fpr_lm_unp"))
colnames(fprdat_m)<-c("alpha","scale","npops","test","fpr")
head(fprdat_m)

# Save the data.frame as an R object to avoid having to perform
# calculations again.
#save(list = ls(all=TRUE), file = "fpr.RData",envir=.GlobalEnv)
write.table(fprdat_m,paste("FST=",fst,
			   "_mcov=",mcov,
			   "_fprdat_melted.tab",sep=""),quote=FALSE,row.names=FALSE,sep="\t")
write.table(fprdat,paste("FST=",fst,
			   "_mcov=",mcov,
			   "_fprdat.tab",sep=""),quote=FALSE,row.names=FALSE,sep="\t")
rm(sim_dat,fprdat,fprdat_m)

#-------------------------#
# True Positive Rate (TPR)#
#-------------------------#
# Calculated as the proportion of true positives that have a 
# p-value < the 1st percentile of the p-value distribution
# i.e. the proportion of true positives that lie in the tail of
# the p-value distribution.

# Initialise a dataframe to store data
tprdat<-data.frame(
  npops=vector(length=16),
  scale=vector(length=16),
  cmh_tp_rates=vector(length=16),
  cmh_woo_tp_rates=vector(length=16),
  binglm_l_tp_rates=vector(length=16),
  binglm_l_ni_tp_rates=vector(length=16),
  qbinglm_unp_tp_rates=vector(length=16),
  gtest_tp_rates=vector(length=16),
  lm_unp_tp_rates=vector(length=16))

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
      tprdat$binglm_l_tp_rates[i]<-nrow(sim_dat[
        sim_dat$binglm_lxp_pval > 0.05 & 
          sim_dat$tp == "1" & 
          sim_dat$binglm_l_pval < quantile(
            sim_dat$binglm_l_pval[
              sim_dat$binglm_lxp_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])

      tprdat$binglm_l_ni_tp_rates[i]<-nrow(sim_dat[
        sim_dat$tp == "1" &
          sim_dat$binglm_l_p_pval_ni < quantile(
            sim_dat$binglm_l_p_pval_ni,0.01,na.rm=TRUE),])/nrow(
              sim_dat[sim_dat$tp == "1",])
      
      #Qbin-GLM
      tprdat$qbinglm_unp_tp_rates[i]<-nrow(sim_dat[
        sim_dat$tp == "1" &
          sim_dat$qbinglm_unp_l_pval < quantile(
            sim_dat$qbinglm_unp_l_pval,0.01,na.rm=TRUE),])/nrow(
              sim_dat[sim_dat$tp == "1",])
      
      #G-test
      tprdat$gtest_tp_rates[i]<-nrow(sim_dat[
        sim_dat$g_lxcxpx_pval > 0.05 &
          sim_dat$tp == "1" &
            sim_dat$g_lxc_pval < quantile(
              sim_dat$g_lxc_pval[
                sim_dat$g_lxcxpx_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                  sim_dat[sim_dat$tp == "1",])
      
      #T-test (LM)
      tprdat$lm_unp_tp_rates[i]<-nrow(sim_dat[
        sim_dat$tp == "1" &
          sim_dat$lm_unp_pval < quantile(
            sim_dat$lm_unp_pval,0.01,na.rm=TRUE),])/nrow(
              sim_dat[sim_dat$tp == "1",])
      
      tprdat$npops[i]=ktxt
      
      tprdat$scale[i]="CT=var."
      i<-i+1
      } else if(res == 1){
        for (sc in scales){
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
          #Bin GLM
          tprdat$binglm_l_tp_rates[i]<-nrow(sim_dat[
            sim_dat$binglm_lxp_pval > 0.05 & 
              sim_dat$tp == "1" & 
              sim_dat$binglm_l_pval < quantile(
                sim_dat$binglm_l_pval[
                  sim_dat$binglm_lxp_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                    sim_dat[sim_dat$tp == "1",])
          
          tprdat$binglm_l_ni_tp_rates[i]<-nrow(sim_dat[
            sim_dat$tp == "1" &
              sim_dat$binglm_l_p_pval_ni < quantile(
                sim_dat$binglm_l_p_pval_ni,0.01,na.rm=TRUE),])/nrow(
                  sim_dat[sim_dat$tp == "1",])
          
          #QBin GLM
          tprdat$qglm_unp_tp_rates[i]<-nrow(sim_dat[
            sim_dat$tp == "1" &
              sim_dat$qbinglm_unp_l_pval < quantile(
                sim_dat$qbinglm_unp_l_pval,0.01,na.rm=TRUE),])/nrow(
                  sim_dat[sim_dat$tp == "1",])
          #G-test
          tprdat$gtest_tp_rates[i]<-nrow(sim_dat[
            sim_dat$g_lxcxpx_pval > 0.05 &
              sim_dat$tp == "1" &
                sim_dat$g_lxc_pval < quantile(
                  sim_dat$g_lxc_pval[
                    sim_dat$g_lxcxpx_pval > 0.05],0.01,na.rm=TRUE),])/nrow(
                      sim_dat[sim_dat$tp == "1",])
          #T-test (LM)
          tprdat$lm_unp_tp_rates[i]<-nrow(sim_dat[
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
                   "binglm_l_tp_rates",
                   "binglm_l_ni_tp_rates",
                   "qbinglm_unp_tp_rates",
                   "gtest_tp_rates",
                   "lm_unp_tp_rates")
                 )
colnames(tprdat_m)<-c("npops","scale","test","tpr")

# Save the data.frame as an R object to avoid having to perform
# calculations again.
write.table(tprdat_m,paste("FST=",fst,
			   "_mcov=",mcov,
			   "_tprdat_melted.tab",sep=""),quote=FALSE,row.names=FALSE,sep="\t")
write.table(tprdat,paste("FST=",fst,
			   "_mcov=",mcov,
			   "_tprdat.tab",sep=""),quote=FALSE,row.names=FALSE,sep="\t")
#save(list = ls(all=TRUE), file = "tpr.RData",envir=.GlobalEnv)