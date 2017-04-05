#!/usr/bin/env Rscript
#----------------------------------------------------------------------#
# This is a BASH shell script that will run the simulations with user  #
# defined parameters (see the suppelementary "pipeline.txt" file)      #
#                                                                      #
# Author: R. Axel W. Wiberg                                            #
# Created: Jan 2017                                                    #
# Last Modified: 06.02.2017                                            #
#----------------------------------------------------------------------#
#
# Load libraries
library("optparse")
library("methods")
# Set maximum number of decimal positions to a
# high number.
options(scipen = 999)

# Define the command-line options list
option_list = list(
  make_option(c("-f","--fst"), type="numeric", default=0, 
              help="FST for the simulation [default = 0.2]", 
              metavar="numeric"),
  
  make_option(c("-N","--Npool"), type="integer", default=0, 
              help="Sample size of each pool [default= 100]", 
              metavar="integer"),
  
  make_option(c("-c","--scale"), type="character", default="Null", 
              help="what scale to re-scale counts [default = neff]", 
              metavar="numeric"),
  
  make_option(c("-r","--rescale"), type="integer", default=2, 
              help="BOOL (0/1): whether to rescale counts [default = 1]", 
              metavar="numeric"),
  
  make_option(c("-k","--Nreplicates"), type="integer", default=0, 
              help="number of replicate lines [default = 2]", 
              metavar="integer"),

  make_option(c("-m","--mcov"), type="integer", default=40, 
              help="mean of covereage distribution [default = 40]", 
              metavar="integer"),

  make_option(c("-n","--NSNPs"), type="integer", default=0, 
              help="number of simulations to run [default = 100]", 
              metavar="integer"),

  make_option(c("-t","--prop_tp"), type="numeric", default=0, 
              help="Proportion of SNPs that are 'True Positives' 
              [default = 0.01 (1%)]", 
              metavar="integer"),

  make_option(c("-s","--selection_diff"), type="numeric", default=0, 
              help="The difference between line 1 and line 2 due to selection
              [default = 0.2]", 
              metavar="integer"),

  make_option(c("-d","--handle"), type="character",default="", 
              help="", 
              metavar="character"),
  
  make_option(c("-o","--outdir"),default="Null", 
               help="The target output directory 
               (default is current directory)", 
               metavar="character")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print_help(opt_parser)

if(opt$fst==0){
  opt$fst<-0.2
  cat("\n\nUsing default: FST = 0.2\n")
  cat("If you want to change this, see the help above.\n\n")
  
}
if(opt$Npool==0){
  opt$Npool<-100
  cat("Using default: Npool = 100\n")
  cat("If you want to change this, see the help above.\n\n")
  
}
if(opt$NSNPs==0){
  opt$NSNPs<-100
  cat("Using default: NSNPs = 100\n")
  cat("If you want to change this, see the help above.\n\n")
  
}
if(opt$Nreplicates==0){
  opt$Nreplicates<-2
  cat("Using default: Nreplicates = 2\n")
  cat("If you want to change this, see the help above.\n\n")
  
}
if(opt$rescale==2){
  opt$rescale<-1
  cat("Using default: rescale = TRUE\n")
  cat("If you want to change this, see the help above.\n\n")
  
}
if(opt$scale == "Null"){
  if(opt$rescale == 1){
    opt$scale<-"neff"
    cat("Using default: Scale = neff\n")
    cat("If you want to change this, see the help above.\n\n")
  }
}
if(opt$prop_tp==0){
  opt$prop_tp<-0.01
  cat("Using default: proporton true positive = 0.01\n")
  cat("If you want to change this, see the help above.\n\n")
  
}
if(opt$selection_diff==0){
  opt$selection_diff<-0.2
  cat("Using default: selection differential = 0.2\n")
  cat("If you want to change this, see the help above.\n\n")
  
}
if(opt$outdir=="Null"){
  opt$outdir<-getwd()
  cat("Using default: current working directory as output",getwd(),"\n")
  cat("If you want to change this, see the help above.\n\n")
  
}

#print(opt)
#--------------------#
# Load the functions
#--------------------#
# Find the location of the current script and source: 
# FequencyTests_Functions.R
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
scriptdir<-dirname(script.name)
cat("Loading functions from: ",
    paste(scriptdir,"/FequencyTests_Functions.R\n",sep=""))
source(paste(scriptdir,"/FrequencyTests_Functions.R",sep=""))
cat("Loading functions: DONE\n\n")

#------------------------------#
# Define simulation parameters #
#------------------------------#
# Simulation parameters:
sim_fst<-opt$fst
sim_N<-opt$Npool
sim_nsims<-opt$NSNPs
sim_k<-opt$Nreplicates
sim_rscale<-opt$rescale
mcov<-opt$mcov
sim_scale<-opt$scale
p_tp<-opt$prop_tp
sel_diff<-opt$selection_diff
handle<-opt$handle

# Make sure scale option is numeric
if((sim_scale != "neff") & (sim_scale != "Null")){
  sim_scale <- as.integer(sim_scale)
}

# Make a directory for the results:
if(sim_rscale == 1){
  dataset<-paste(paste("k=",sim_k,sep=""),
                 paste("_fst=",sim_fst,sep=""),
                 paste("_N=",sim_N,sep=""),
                 paste("_mcov=",mcov,sep=""),
                 paste("_res=",sim_rscale,sep=""),
                 paste("_scale=",sim_scale,sep=""),
                 paste("_SNPs=",as.integer(sim_nsims),sep=""),
                 paste("_p_tp=",p_tp,sep=""),
                 paste("_sel_diff=",sel_diff,sep=""),
                 paste("_",handle,sep=""),
                 sep = "")
}else{
  dataset<-paste(paste("k=",sim_k,sep=""),
                 paste("_fst=",sim_fst,sep=""),
                 paste("_N=",sim_N,sep=""),
                 paste("_mcov=",mcov,sep=""),
                 paste("_res=",sim_rscale,sep=""),
                 paste("_SNPs=",as.integer(sim_nsims),sep=""),
                 paste("_p_tp=",p_tp,sep=""),
                 paste("_sel_diff=",sel_diff,sep=""),
                 paste("_",handle,sep=""),
                 sep = "")  
}
subdir<-dataset
# Set the directory for the results.
dir<-opt$outdir
cat("Writing output to directory:\n",paste(dir,"/",subdir,sep=""),"\n\n")
outdir<-paste(dir,"/",subdir,"/",sep="")
dir.create(paste(dir,"/",subdir,sep=""))
setwd(paste(dir,"/",subdir,sep=""))

#-----------------------#
# Run Simulations
#-----------------------#
cat("\nSimulation Parameters Are:\n",dataset,"\n\n")
e<-new.env()
assign(dataset,data_gen(als=sim_nsims,
         k=sim_k, l=2,
         mcov=mcov,
         cons=0,
         res=sim_rscale,scale=sim_scale,
         N=sim_N,a=0.2,b=0.2,fst=sim_fst,
         prop_tp=p_tp,selection_diff = sel_diff
         ),envir=e)

data<-e[[dataset]]

#--------------------------#
# Write the data as output
#--------------------------#
write.table(e[[dataset]],paste(outdir,
                          "FrequencyTest_Simulations_",
                          dataset,".csv",sep=""),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
# END