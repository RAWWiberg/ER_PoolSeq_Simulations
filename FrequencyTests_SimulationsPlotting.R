#-------------------------------------------------------#
# R Script for producing plots from simulation results. #
#                                                       #
# Author: R. Axel W. Wiberg                             #
# Created: Jan 2017                                     #
# Last Modified: 06.02.2017                             #
#                                                       # 
#-------------------------------------------------------#

# Clean environment
rm(list = ls(all = TRUE))

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
        "_fst=0.2_N=100_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
        "_fst=0.2_N=100_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
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
          "_fst=0.2_N=100_res=1_scale=",sc,"_SNPs=100000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
          "_fst=0.2_N=100_res=1_scale=",sc,"_SNPs=100000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
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
        "_fst=0.2_N=100_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
        "_fst=0.2_N=100_res=0_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
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
            "_fst=0.2_N=100_res=1_scale=",sc,"_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA/FrequencyTest_Simulations_k=",k,
            "_fst=0.2_N=100_res=1_scale=",sc,"_SNPs=1000000_p_tp=0.01_sel_diff=0.2_SIMDATA.csv",
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
#save(list = ls(all=TRUE), file = "tpr.RData",envir=.GlobalEnv)



#-------------------------#
# True Positive Rate (TPR)#
#-------------------------#

# Initialise a dataframe to store data
tprdatcons<-data.frame(
  npops=vector(length=28),
  cmh_tp_rates=vector(length=28),
  qglm_unp_tp_rates=vector(length=28),
  ttest_unp_tp_rates=vector(length=28),
  sim=vector(length=28),
  snps=vector(length=28))

# Load the data, perform calculations and store the results
# results files are quite big so they are loaded one at a time.
i <- 1
for(k in c(4,10)){
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
          "_fst=0.2_N=100_res=1_scale=100_SNPs=10000_p_tp=0.01_sel_diff=0.2_",
          sim,"/FrequencyTest_Simulations_k=",k,
          "_fst=0.2_N=100_res=1_scale=100_SNPs=10000_p_tp=0.01_sel_diff=0.2_",
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
        tprdatcons$qglm_unp_tp_rates[i]<-nrow(sim_dat[
          sim_dat$tp == "1" &
            sim_dat$qglm_unp_l_pval < quantile(
              sim_dat$qglm_unp_l_pval,0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
        tprdatcons$ttest_unp_tp_rates[i]<-nrow(sim_dat[
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
          "_fst=0.2_N=100_res=1_scale=100_SNPs=1000000_p_tp=0.01_sel_diff=0.2_",
          sim,"/FrequencyTest_Simulations_k=",k,
          "_fst=0.2_N=100_res=1_scale=100_SNPs=1000000_p_tp=0.01_sel_diff=0.2_",
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
        tprdatcons$qglm_unp_tp_rates[i]<-nrow(sim_dat[
          sim_dat$tp == "1" &
            sim_dat$qglm_unp_l_pval < quantile(
              sim_dat$qglm_unp_l_pval,0.01,na.rm=TRUE),])/nrow(
                sim_dat[sim_dat$tp == "1",])
        tprdatcons$ttest_unp_tp_rates[i]<-nrow(sim_dat[
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
                   "qglm_unp_tp_rates",
                   "ttest_unp_tp_rates"))
colnames(tprdatcons_m)<-c("npops","sim","snps","test","tpr")

tprdatcons_m

# Save the data.frame as an R object to avoid having to perform
# calculations again.
#save(list = ls(all=TRUE), file = "tprcons.RData",envir=.GlobalEnv)



#---------#
# Plotting#
#---------#
source("~/Desktop/Data/RData/RScripts/FrequencyTests_Functions.R")
setwd("~/Desktop/Data/allele_frequency_analysis_project/data/")
# Load saved R objects
load("fpr.RData")
load("tpr.RData")
load("tprcons.RData")
#----------#
# Plot FPR #
#----------#
fpr_v_alp<-ggplot(data=fprdat_m)+
  xlab(expression(paste(alpha,"-threshold")))+
  ylab("False Positive Rate")+
  scale_x_continuous(limits=c(min(alpha),max(alpha)),
                     breaks=c(0.0001,0.001,0.01,0.1,0.5),trans="log10")+
  scale_y_continuous(limits=c(0.0001,1),
                     breaks=c(0.0001,0.001,0.01,0.1,1),trans="log10")+
  geom_abline(intercept=0,slope=1,colour="grey")+
  geom_abline(intercept=log10(5),slope=1,colour="grey",linetype="dashed")+
  geom_text(aes(0.00025,0.00025*5,label="5"),size = 3)+
  geom_abline(intercept=log10(2),slope=1,colour="grey",linetype="dashed")+
  geom_text(aes(0.00025,0.00025*2,label="2"),size = 3)+
  geom_abline(intercept=1,slope=1,colour="grey",linetype="dashed")+
  geom_text(aes(0.00025,0.00025*10,label="10"),size = 3)+
  geom_point(data=fprdat_m,aes(x=alpha,y=fpr,shape=test,colour=test),
             size = 3)+
  scale_shape_manual("",labels=c("CMH","CMH+Woolf","Binomial\nGLM",
                                   "Quasibinomial\nGLM","G-test","LM"),
                     values=c(16,17,15,3,7,10))+
  scale_colour_manual("",labels = c("CMH","CMH+Woolf","Binomial\nGLM",
                                "Quasibinomial\nGLM","G-test","LM"),
                  values=c('#1b9e77','#d95f02','#7570b3',
                           '#e7298a','#66a61e','#e6ab02'))+
  facet_grid(scale~npops)
fpr_v_alp +annotation_logticks() + my.theme + 
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(7,"mm"))

#----------#
# Plot FPR #
#----------#
# Re-order the factor levels for the npops (Nr. Replicates) column
# this makes the plots more intuitive.
tprdat_m$npops<-factor(tprdat_m$npops,levels=c("k=2","k=3","k=4","k=10"))

tp_plot<-ggplot(data=tprdat_m)+
  geom_point(aes(x=npops,y=tpr,shape=test,colour=test),
             colour = "black",
             size = 3,
             position=position_jitter(width = 0.05))+
  xlab("Nr. Replicated Treatment Lines")+
  ylab("True Positive Rate")+
  scale_shape_manual("",labels=c("CMH","CMH+Woolf","Binomial\nGLM",
                                 "Quasibinomial\nGLM","G-test","LM"),
                     values=c(16,17,15,3,7,10))+
  scale_colour_manual("",labels = c("CMH","CMH+Woolf","Binomial\nGLM",
                                    "Quasibinomial\nGLM","G-test","LM"),
                      values=c('#1b9e77','#d95f02','#7570b3',
                               '#e7298a','#66a61e','#e6ab02'))+
  facet_grid(scale~.)

tp_plot + my.theme + theme(
  title = element_text(size = 10,face="bold"),
  axis.text.x = element_text(size = 10,face = "bold"),
  axis.text.y = element_text(size = 10,face = "bold"),
  axis.ticks.x = element_blank(),
  strip.text = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.key.size = unit(7,"mm"))

#-------------------------#
# Plot consistency of FPR #
#-------------------------#
tpconsplot<-ggplot()+
  geom_boxplot(data=tprdatcons_m,aes(npops,tpr,fill=test))+
  xlab("Nr Replicated Treatment Lines")+
  ylab("True Positive Rate")+
  scale_fill_manual("",
                      values=c('#1b9e77','#e7298a','#e6ab02'),
                      labels=c("CMH-test","Quasibinomial\nGLM",
                               "LM"))+
  facet_wrap(~snps)
tpconsplot+my.theme + 
  theme(axis.text = element_text(size = 12,face ="bold"),
        strip.text=element_text(size = 12),
        title = element_text(size=12,face="bold"),
        legend.text = element_text(size = 12),
        legend.key.size = unit(7,"mm"))

#----------------------------------------------------------------#
# Plot Distributions of Quasibinomial GLMS and CMH-test p-values #
#----------------------------------------------------------------#
# Load the data: QBGLM p-values
qbglmdat<-read.table("qbglm_tp0.tab",
                     sep="\t",
                     header=FALSE)
colnames(qbglmdat)<-c("pval","npops","scale")

# Plot the histogram
qbglm_his<-ggplot()+
  geom_histogram(data=qbglmdat,aes(pval),binwidth=0.05)+
  xlab("p-value")+
  ylab("Count")+
  facet_grid(scale~npops)

qbglm_his + my.theme + 
  theme(axis.text.x = element_text(size = 12,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12))

# Load the data: CMH-test p-values
cmhdat<-read.table("cmh_tp0.tab",
                   sep="\t",
                   header=FALSE)
colnames(cmhdat)<-c("pval","npops","scale")
# Plot the histogram
cmh_his<-ggplot()+
  geom_histogram(data=cmhdat,
                 aes(pval),binwidth=0.05)+
  xlab("p-value")+
  ylab("Count")+
  facet_grid(scale~npops)
cmh_his + my.theme + 
  theme(axis.text.x = element_text(size = 12,
                                   angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12))

#-------------------------------------------------------------------------#
# Plot Mean vs. SD of allele frequency differences across all simulations #
#-------------------------------------------------------------------------#
# Load the data
data<-read.table("mean_sd_tp0.tab", header = FALSE)
colnames(data)<-c("sd_diffs","mean_diffs","npops","scale","tp")
# Plot distributions of mean_diffs and sd_diffs

# Plot mean_diffs vs sd_diffs and correlation coefficients (Figure S3)
sd_v_mean_plot <- ggplot() +
  geom_point(data=data, aes(mean_diffs,sd_diffs),
             size = 3,
             alpha = 1/10) +
  xlab("Mean Allele Frequency Difference") +
  ylab("SD of Allele Frequency Differences") +
  facet_grid(npops~scale,scales="free_x")

# Get correlation coefficients and p-values as text for plots
# ####
text_dat<-data.frame(p=vector(length=length(unique(data$scale))),
                     r=vector(length=length(unique(data$scale))),
                     npops = rep(unique(data$npops),
                                 length(unique(data$scale))),
                     scale = rep(unique(data$scale),each=
                                   length(unique(data$npops))))
cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=2" &
                                      data$scale=="CT=200"],
                y = data$sd_diffs[data$npops =="k=2" &
                                    data$scale=="CT=200"],
                method="spearman")

text_dat$p[
  text_dat$npops == "k=2" & 
    text_dat$scale == "CT=200"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=2" & 
    text_dat$scale == "CT=200"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=3" &
                                      data$scale=="CT=200"],
                y = data$sd_diffs[data$npops =="k=3" &
                                    data$scale=="CT=200"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=3" & 
    text_dat$scale == "CT=200"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=3" & 
    text_dat$scale == "CT=200"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=4" &
                                      data$scale=="CT=200"],
                y = data$sd_diffs[data$npops =="k=4" &
                                    data$scale=="CT=200"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=4" & 
    text_dat$scale == "CT=200"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=4" & 
    text_dat$scale == "CT=200"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=10" &
                                      data$scale=="CT=200"],
                y = data$sd_diffs[data$npops =="k=10" &
                                    data$scale=="CT=200"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=10" & 
    text_dat$scale == "CT=200"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=10" & 
    text_dat$scale == "CT=200"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=2" &
                                      data$scale=="CT=100"],
                y = data$sd_diffs[data$npops =="k=2" &
                                    data$scale=="CT=100"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=2" & 
    text_dat$scale == "CT=100"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=2" & 
    text_dat$scale == "CT=100"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=3" &
                                      data$scale=="CT=100"],
                y = data$sd_diffs[data$npops =="k=3" &
                                    data$scale=="CT=100"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=3" & 
    text_dat$scale == "CT=100"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=3" & 
    text_dat$scale == "CT=100"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=4" &
                                      data$scale=="CT=100"],
                y = data$sd_diffs[data$npops =="k=4" &
                                    data$scale=="CT=100"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=4" & 
    text_dat$scale == "CT=100"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=4" & 
    text_dat$scale == "CT=100"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=10" &
                                      data$scale=="CT=100"],
                y = data$sd_diffs[data$npops =="k=10" &
                                    data$scale=="CT=100"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=10" & 
    text_dat$scale == "CT=100"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=10" & 
    text_dat$scale == "CT=100"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=2" &
                                      data$scale=="CT=var."],
                y = data$sd_diffs[data$npops =="k=2" &
                                    data$scale=="CT=var."],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=2" & 
    text_dat$scale == "CT=var."]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=2" & 
    text_dat$scale == "CT=var."]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=3" &
                                      data$scale=="CT=var."],
                y = data$sd_diffs[data$npops =="k=3" &
                                    data$scale=="CT=var."],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=3" & 
    text_dat$scale == "CT=var."]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=3" & 
    text_dat$scale == "CT=var."]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=4" &
                                      data$scale=="CT=var."],
                y = data$sd_diffs[data$npops =="k=4" &
                                    data$scale=="CT=var."],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=4" & 
    text_dat$scale == "CT=var."]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=4" & 
    text_dat$scale == "CT=var."]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=10" &
                                      data$scale=="CT=var."],
                y = data$sd_diffs[data$npops =="k=10" &
                                    data$scale=="CT=var."],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=10" & 
    text_dat$scale == "CT=var."]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=10" & 
    text_dat$scale == "CT=var."]<- as.character(round(cor.t$estimate,2))


cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=2" &
                                      data$scale=="CT=neff"],
                y = data$sd_diffs[data$npops =="k=2" &
                                    data$scale=="CT=neff"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=2" & 
    text_dat$scale == "CT=neff"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=2" & 
    text_dat$scale == "CT=neff"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=3" &
                                      data$scale=="CT=neff"],
                y = data$sd_diffs[data$npops =="k=3" &
                                    data$scale=="CT=neff"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=3" & 
    text_dat$scale == "CT=neff"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=3" & 
    text_dat$scale == "CT=neff"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=4" &
                                      data$scale=="CT=neff"],
                y = data$sd_diffs[data$npops =="k=4" &
                                    data$scale=="CT=neff"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=4" & 
    text_dat$scale == "CT=neff"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=4" & 
    text_dat$scale == "CT=neff"]<- as.character(round(cor.t$estimate,2))

cor.t<-cor.test(x = data$mean_diffs[data$npops =="k=10" &
                                      data$scale=="CT=neff"],
                y = data$sd_diffs[data$npops =="k=10" &
                                    data$scale=="CT=neff"],
                method="spearman")
text_dat$p[
  text_dat$npops == "k=10" & 
    text_dat$scale == "CT=neff"]<- as.character(round(cor.t$p.value,2))
text_dat$r[
  text_dat$npops == "k=10" & 
    text_dat$scale == "CT=neff"]<- as.character(round(cor.t$estimate,2))

print(text_dat)
# ####

# Plot Figure S3 with text labels 
sd_v_mean_plot + my.theme+ theme(axis.text = element_text(size=12))+
  geom_text(data=text_dat,aes(x=0.5,y=0.9,
                              label=paste("p = ",p,sep="")),size=4) +
  geom_text(data=text_dat,aes(x=0.5,y=1,
                              label=paste("rho = ",r,sep="")),size=4)+
  scale_x_continuous(breaks = c(-0.5,0,0.5))

#------------------------------------------------------------------------#
# Plot the histogram of the mean allele frequency differences (Figure S1)#
#------------------------------------------------------------------------#
mean_hist_plot <- ggplot() +
  geom_histogram(data=data, aes(mean_diffs)) +
  xlab("Mean Allele Frequency Difference") +
  ylab("Count") +
  facet_grid(npops~scale,scales="free_x")

mean_hist_plot + my.theme+ theme(axis.text = element_text(size=12))

#-------------------------------------------------------------------------#
# Plot the histogram of the SD of allele frequency differences (Figure S2)#
#-------------------------------------------------------------------------#
sd_hist_plot <- ggplot() +
  geom_histogram(data=data, aes(sd_diffs)) +
  xlab("SD Allele Frequency Difference") +
  ylab("Count") +
  facet_grid(npops~scale,scales="free_x")

sd_hist_plot + my.theme+ theme(axis.text = element_text(size=12))


