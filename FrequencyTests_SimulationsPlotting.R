#-------------------------------------------------------#
# R Script for producing plots from simulation results. #
#                                                       #
# Author: R. Axel W. Wiberg                             #
# Created: Jan 2017                                     #
# Last Modified: 18.05.2017                             #
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
source("~/packages/ER_PoolSeq_Simulations/FrequencyTests_Functions.R")
source("~/RData/RScripts/ggplot_theme.R")
# Set the working directory
setwd("~/PhD/allele_frequency_tests/data/")
alpha<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5)
#----------#
# Plotting #
#----------#
# Load data
fst<-0.1
fprdat_m<-read.table(paste("FST=",fst,"_mcov=200_fprdat_melted.tab",sep=""),
                     header=TRUE,sep="\t")
fprdat<-read.table(paste("FST=",fst,"_mcov=200_fprdat.tab",sep=""),
                   header=TRUE,sep="\t")
head(fprdat_m)
# Re-order the npops column
levels(fprdat_m$npops)
fprdat_m$npops<-factor(fprdat_m$npops,levels=c("k=2","k=3","k=4","k=10"))

tprdat_m<-read.table(paste("FST=",fst,"_mcov=200_tprdat_melted.tab",sep=""),
                     header=TRUE,sep="\t")
tprdat<-read.table(paste("FST=",fst,"_mcov=200_tprdat.tab",sep=""),
                   header=TRUE,sep="\t")
# Re-order the npops column
levels(tprdat_m$npops)
tprdat_m$npops<-factor(tprdat_m$npops,levels=c("k=2","k=3","k=4","k=10"))

tprconsdat<-read.table("tprcons.RData")
#----------#
# Plot FPR #
#----------#
levels(fprdat_m$test)
fpr_v_alp<-ggplot(data=fprdat_m)+
  xlab(expression(paste(alpha,"-threshold")))+
  ylab("False Positive Rate")+
  scale_x_continuous(limits=c(min(alpha),max(alpha)),
                     breaks=c(0.0001,0.001,0.01,0.1,0.5),trans="log10")+
  scale_y_continuous(limits=c(0.00001,1),
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
  scale_shape_manual("",labels = c("Binomial\nGLM (1)","Binomial\nGLM (2)",
                                   "CMH-test","CMH-test+Woolf-test",
                                   "G-test","LM","Quasibinomial\nGLM"),
                     values=c(16,17,15,3,7,10,23))+
  scale_colour_manual("",labels = c("Binomial\nGLM (1)","Binomial\nGLM (2)",
                                    "CMH-test","CMH-test+Woolf-test",
                                    "G-test","LM","Quasibinomial\nGLM"),
                      values=c('#1b9e77','#d95f02','#7570b3',
                               '#e7298a','#66a61e','#e6ab02',
                               'black'))+
  facet_grid(scale~npops)
fpr_v_alp +annotation_logticks() + my.theme + 
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=1,vjust=1),
#        panel.grid.major = element_line(colour = "grey"),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.key.size = unit(7,"mm"),
        legend.position = "right")


#----------#
# Plot TPR #
#----------#
head(tprdat_m)
levels(tprdat_m$test)
tprdat_m$test<-factor(tprdat_m$test,levels=c("binglm_l_tp_rates",
                                             "binglm_l_ni_tp_rates",
                                             "cmh_tp_rates",
                                             "cmh_woo_tp_rates",
                                             "gtest_tp_rates",
                                             "lm_unp_tp_rates",
                                             "qbinglm_unp_tp_rates"))
tp_plot<-ggplot(data=tprdat_m)+
  geom_point(aes(x=npops,y=tpr,shape=test,colour=test),
             size = 3,
             position=position_jitter(width = 0.05))+
  xlab("Nr. Replicated Treatment Lines")+
  ylab("True Positive Rate")+
  scale_shape_manual("",labels = c("Binomial\nGLM (1)","Binomial\nGLM (2)",
                                   "CMH-test","CMH-test+Woolf-test",
                                   "G-test","LM","Quasibinomial\nGLM"),
                     values=c(16,17,15,3,7,10,23))+
  scale_colour_manual("",labels = c("Binomial\nGLM (1)","Binomial\nGLM (2)",
                                    "CMH-test","CMH-test+Woolf-test",
                                    "G-test","LM","Quasibinomial\nGLM"),
                      values=c('#1b9e77','#d95f02','#7570b3',
                               '#e7298a','#66a61e','#e6ab02',
                               'black'))+
  facet_grid(scale~.)

tp_plot + my.theme + theme(
  title = element_text(size = 10,face="bold"),
  axis.text.x = element_text(size = 10,face = "bold"),
  axis.text.y = element_text(size = 10,face = "bold"),
  axis.ticks.x = element_blank(),
  strip.text = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.key.size = unit(7,"mm"),
  legend.position = "right")

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
# reorder npops column
qbglmdat$npops<-factor(qbglmdat$npops,levels=c("k=2","k=3","k=4","k=10"))

# Plot the histogram (FIGURE S5)
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
rm(qbglm_dat)
# Load the data: CMH-test p-values
cmhdat<-read.table("cmh_tp0.tab",
                   sep="\t",
                   header=FALSE)
colnames(cmhdat)<-c("pval","npops","scale")
# reorder npops column
cmhdat$npops<-factor(cmhdat$npops,levels=c("k=2","k=3","k=4","k=10"))
# Plot the histogram (FIGURE S4)
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
rm(cmhdat)
#-------------------------------------------------------------------------#
# Plot Mean vs. SD of allele frequency differences across all simulations #
#-------------------------------------------------------------------------#
# Load the data
data<-read.table("mean_sd_tp0.tab", header = FALSE)
colnames(data)<-c("sd_diffs","mean_diffs","npops","scale")
# reorder npops column
data$npops<-factor(data$npops,levels=c("k=2","k=3","k=4","k=10"))

# Plot distributions of mean_diffs and sd_diffs
#----------------------------------------------------------------------#
# Plot mean_diffs vs sd_diffs and correlation coefficients (Figure S3) #
#----------------------------------------------------------------------#
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


