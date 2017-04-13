#-------------------------------------------------------------------#
# Functions used in simulating allele frequency differences and the #
# assessment of performance of different tests of allele frequency  #
# differences.                                                      #
#                                                                   #
# Also includes some test data used to test the functions           #
# Also includes a custom ggplot2 theme                              #
#                                                                   #
# Author: R. Axel W. Wiberg                                         #
# Created: 05.08.2016                                               #
# Last Modified: 02.04.2017                                         #
#-------------------------------------------------------------------#

#-----------------------------------#
# LOAD LIBRARIES NEEDED:
#-----------------------------------#
library(ggplot2)
library(truncnorm)
#-----------------------------------#
# SOME EXAMPLE DATA:
#-----------------------------------#
# Example data to test scripts #### 
# matrixk2: k = 2, non-consistent effect, zero counts present
# matrixk2: k = 3, non-consistent effect, zero counts present
# matrixk2: k = 4, non-consistent effect, zero counts present
matrixk2 <- array(c(66,90,0,1,72,60,0,0),dim=c(2,2,2),
                  dimnames=list(c("H","L"),c("A","T"),
                                c("Pop1","Pop2")))

matrixk3 <- array(c(66,90,0,1,72,60,0,0,69,0,21,72),dim=c(2,2,3),
                  dimnames=list(c("H","L"),c("A","T"),
                                c("Pop1","Pop2","Pop3")))

matrixk4 <- array(c(66,90,0,1,72,60,0,0,69,0,21,72,69,0,21,72),dim=c(2,2,4),
                  dimnames=list(c("H","L"),c("A","T"),
                                c("Pop1","Pop2","Pop3","Pop4")))
# ####

#-----------------------------------#
# GGPLOT2 THEME:                    #
#-----------------------------------#
my.theme <- theme(panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(NULL), 
                  #element_line(colour = "grey"),
                  panel.grid.minor = element_line(NULL),
                  panel.border = element_rect(colour = "grey50", 
                                              fill = NA),
                  strip.background = element_rect(colour = "black", 
                                                  fill = "white"), 
                  strip.text = element_text(size = 20),
                  text = element_text(size = 20,
                                      colour = "black"),
                  axis.text = element_text(size = 20, 
                                           colour = "black"))

#-----------------------------------#
# DEFINE FUNCTIONS:                 #
#-----------------------------------#

#
# FUNCTION: calculate SE
se <- function(x){
  sd(x)/sqrt(length(x))}
#
# FUNCTION: logit and logistic functions
logit <- function(x){log(x/(1-x))}

logistic <- function(x){1/(1+exp(-x))}

#
# FUNCTION: Woolf-test
# The script comes from the help page for the mantelhaen.test()
#?mantelhaen.test()
woolf.test <- function(x) {
  x <- x + 1 / 2
  k <- dim(x)[3]
  or <- apply(x, 3, function(x) (x[1,1]*x[2,2])/(x[1,2]*x[2,1]))
  w <-  apply(x, 3, function(x) 1 / sum(1 / x))
  woolf <- sum(w * (log(or) - weighted.mean(log(or), w)) ^ 2)
  df <- k-1
  p <- 1 - pchisq(woolf, df)
  dat <- c(woolf,df,p)
  names(dat) <- c("Woolf", "df", "p-value")
  dat
}

#
# FUNCTION: The G-test
# G-test direct G computations
# Reference; Sokal and Rohlf (1969, 1981, 2015) Biometry, 
# W.H.Frieman and Co. San Francisco
G_test<-function(array,add=FALSE,correction="none"){
  if(add == TRUE & any(array == 0)){
    array<-array+1
  }
  if(correction=="cont"){
    f_hat<-G_test_fhat(array)$f_hat
    for(i in 1:length(array)){
      if(array[i] < f_hat[i]){
        array[i] <- array[i] + 0.5
      }
      if(array[i] > f_hat[i]){
        array[i] <- array[i] - 0.5    
      }
    }
  }
  # L = line
  # P = population
  # C = allele counts
  #number of levels for each table
  l <- dim(array)[1]
  c <- dim(array)[2]
  p <- dim(array)[3]
  grand <- sum(array)*log(sum(array))
  #cat(c("grand",grand,"\n")) # Test line
  # Get PxL
  PxLa<-array(dim=c(p,l))
  for(i in seq(1,p,1)){
    for(j in seq(1,l,1)){
      PxLa[i,j] <- sum(array[j,,i])
    }
  }
  PXLv <- vector(length=(l*p))
  for(i in seq(1,(l*p))){
    #print(array[i]*log(array[i]))
    PXLv[i] <- PxLa[i]*log(PxLa[i])
  }
  PxL <- sum(PXLv)
  #cat(c("PxL",PxL,"\n")) # Test line
  # Get PxC
  PxCa<-array(dim=c(p,c))
  for(i in seq(1,p,1)){
    for(j in seq(1,c,1)){
      PxCa[i,j] <- sum(array[,j,i])
    }
  }
  PXCv <- vector(length=(p*c))
  for(i in seq(1,(p*c))){
    #print(array[i]*log(array[i]))
    PXCv[i] <- PxCa[i]*log(PxCa[i])
  }
  PxC <- sum(PXCv)
  #cat(c("PxC",PxC,"\n")) # Test line
  # Get LxC
  LXCa<-array(dim=c(l,c))
  for(i in seq(1,l,1)){
    for(j in seq(1,c,1)){
      LXCa[i,j]<-sum(array[i,j,])  
    }
  }
  LXCv <- vector(length=(l*c))
  for(i in seq(1,(l*c))){
    #print(array[i]*log(array[i]))
    LXCv[i] <- LXCa[i]*log(LXCa[i])
  }
  LxC <- sum(LXCv)
  #cat(c("LxC",LxC,"\n")) # Test line
  # Get PxLxC
  PXLXCv <- vector(length=(l*c*p))
  for(i in seq(1,(l*c*p))){
    #print(array[i]*log(array[i]))
    PXLXCv[i] <- array[i]*log(array[i])
  }
  PxLxC <- sum(PXLXCv)
  #cat(c("PxLxC",PxLxC,"\n")) # Test line
  Pv <- vector(length=p)
  for(i in seq(1,p,1)){
    Pv[i] <- sum(array[,,i])*log(sum(array[,,i]))
  }
  P <- sum(Pv)
  #cat(c("P",P,"\n")) # Test line
  Lv <- vector(length=l)
  for(i in seq(1,l,1)){
    Lv[i] <- sum(array[i,,])*log(sum(array[i,,]))
  }
  L <- sum(Lv)
  #cat(c("L",L,"\n")) # Test line
  Cv <- vector(length=c)
  for(i in seq(1,c,1)){
    Cv[i] <- sum(array[,i,])*log(sum(array[,i,]))
  }
  C <- sum(Cv)
  #cat(c("C",C,"\n")) # Test line
  PxLi_G <- 2*(PxL-(P+L)+grand)
  PxCi_G <- 2*(PxC-(P+C)+grand)
  LxCi_G <- 2*(LxC-(L+C)+grand)
  PxLxCi_G <- 2*(PxLxC-(P+L+C)+(2*grand))
  PxLxCx_G <- 2*(PxLxC-(PxL+PxC+LxC)+(P+L+C)-grand)
  PxLi_df <- (p*l)-(p-1)-(l-1)-1 
  PxCi_df <- (p*c)-(p-1)-(c-1)-1 
  LxCi_df <- (l*c)-(l-1)-(c-1)-1
  PxLxCi_df <- (p*l*c)-(p-1)-(l-1)-(c-1)-1
  PxLxCx_df <- (p-1)*(l-1)*(c-1)
  PxLi_p <- pchisq(PxLi_G,df=PxLi_df,lower.tail=FALSE) 
  PxCi_p <- pchisq(PxCi_G,df=PxCi_df,lower.tail=FALSE)
  LxCi_p <- pchisq(LxCi_G,df=LxCi_df,lower.tail=FALSE)
  PxLxCi_p <- pchisq(PxLxCi_G,df=PxLxCi_df,lower.tail=FALSE)
  PxLxCx_p <- pchisq(PxLxCx_G,df=PxLxCx_df,lower.tail=FALSE)
  
  printarray <- matrix(c(PxLi_G,PxCi_G,LxCi_G,PxLxCi_G,PxLxCx_G,
                         PxLi_df,PxCi_df,LxCi_df,PxLxCi_df,PxLxCx_df,
                         PxLi_p,PxCi_p,LxCi_p,PxLxCi_p,PxLxCx_p),nrow=5,
                       dimnames=list(c("PxL independence: ","PxC independence: ",
                                       "LxC independence: ","PxLxC independence: ","PxLxC interaction: "),
                                     c("G","df","p-value")))
  return(printarray)
}

#
# FUNCTION: Get GLM results from an array ####
#get GLM data-set from k-way table
get_glm_dat <- function(array,zeroes=1){
  if(zeroes == 1){
    if(any(array == 0)){
      array <- array+1
    }
  }
  A_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  Tot_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  tr_l <- vector(length = dim(array)[3]*dim(array)[1])
  rep <- vector(length = dim(array)[3]*dim(array)[1])
  j <-1
  for(k in seq(1,dim(array)[3],1)){
    for(i in seq(1,dim(array)[1],1)){
      #      print(c(i,j,k))
      A_Cnt[j]<-array[i,1,k]
      Tot_Cnt[j]<-sum(array[i,,k])
      tr_l[j] <- as.character(i)
      rep[j] <- as.character(k)
      j <- j + 1
    }
  }
  d<-data.frame("A_Cnt"=A_Cnt,"Tot_Cnt"=Tot_Cnt,"tr_l"=tr_l,"rep"=rep)
  mod <- anova(glm(
    cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$rep+d$tr_l+d$tr_l:d$rep,
                   family = "binomial"),test="LRT")
  return(mod)
}
#
# FUNCTION: convert array to data.frame
get_dat <- function(array,zeroes=1){
  if(zeroes == 1){
    if(any(array == 0)){
      array<-array+1
    }
  }
  A_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  Tot_Cnt <- vector(length = dim(array)[3]*dim(array)[1])
  tr_l <- vector(length = dim(array)[3]*dim(array)[1])
  rep <- vector(length = dim(array)[3]*dim(array)[1])
  j <-1
  for(k in seq(1,dim(array)[3],1)){
    for(i in seq(1,dim(array)[1],1)){
      #      print(c(i,j,k))
      A_Cnt[j]<-array[i,1,k]
      Tot_Cnt[j]<-sum(array[i,,k])
      tr_l[j] <- as.character(i)
      rep[j] <- as.character(k)
      j <- j + 1
    }
  }
  d<-data.frame("A_Cnt"=A_Cnt,"Tot_Cnt"=Tot_Cnt,"tr_l"=tr_l,"rep"=rep)
  return(d)
}

# Effective sample size formula #
# Reference; Feder et al., 2012. Plos ONE. 7:e48588
# n = nr. chromosomes in pool
# m = read coverage
neff <- function(n,m){round(((m*n)-1)/(m+n))}

# FUNCTION: get mean and sd allele frequencies from .sync format data frame
#testdat<-OtW_data_top2000[1,c(4,11,5,12,6,13)]
maj_als<-function(snp,als,dat){
  totcount <- rep(0,4)
  for(r in 1:ncol(dat)){
    rep <- as.numeric(unlist(strsplit(
      as.character(dat[snp,r]),split=":")))[1:4]
    #print(rep)
    for(a in 1:length(rep)){
      totcount[a] <- totcount[a]+rep[a]
    }
  }
  mals_indx <- order(totcount,decreasing = TRUE)[c(1,2)]
  #print(mals_indx)
  #mals <- als[mals_indx]
  return(mals_indx)
}
als <- c("A","T","C","G")

# FUNCTION: Make arrays from .sync file


#
# FUNCTION: Data generator
# Main Simulation Function
# als = nr. simulations (SNPs)
# k = nr. replicates
# l = nr. treatment lines per k
# mcov = mean coverage (mean of a poisson distribution) 
# cons = (bool) 1: force consistent allele frequency differences
#               0: allow some inconsistency
# resc = (bool) 1: rescale the allele counts
#               0: don't rescale the allele counts
# scale = (if resc = 1) the total count (CT) to which allele f
#         requencies should be scaled. 
#         If scale = "neff" use the "effective sample size
#         formula of Kolaczkowski et al., (2011) Genetics. 187:245-260
# N = size of pool samples taken from treatment lines.
#     number of chromosomes/alleles is then n = 2*N
# a = alpha parameter of beta distribution for ancestral population 
#     allele frequencies
# b = beta parameter of beta distribution for ancestral population 
#     allele frequencies
# fst = fst value between ancestral population and the
#       treatment lines
# selection_diff = the average allele frequency difference between 
#                  treat_line 1 and treat_line 2 due to selection
data_gen <- function(als=100,k=2,l=2,
                     mcov=40,cons=1,resc="neff",scale=1,
                     N=100,a=0.2,b=0.2,fst = 0.2,prop_tp=0.01,
                     selection_diff = 0.2){
  # Create some empty vectors to store results
  count_dat<-vector(mode = "character", length = als)
  af_dat<-vector(mode = "character", length = als)  

  # CMH-test results
  cmh_pval<-vector(mode = "numeric", length = als)
  cmh_cOR<-vector(mode = "numeric", length = als)
  cmh_chi<-vector(mode = "numeric", length = als)
  
  # Woolf-test results
  woo_pval<-vector(mode = "numeric", length = als)
  woo_W<-vector(mode = "numeric", length = als)
  
  # G-test results
  g_lxcxpx_pval<-vector(mode = "numeric",length=als)
  g_lxc_pval<-vector(mode="numeric",length=als)
  g_lxcxpx_g<-vector(mode = "numeric",length=als)
  g_lxc_g<-vector(mode="numeric",length=als)
  
  # Binomial GLM results
  binglm_lxpx_pval<-vector(mode="numeric",length=als)
  binglm_l_pval<-vector(mode="numeric",length=als)
  binglm_lxpx_dev<-vector(mode="numeric",length=als)
  binglm_l_dev<-vector(mode="numeric",length=als)
  binglm_w<-vector(mode="numeric",length=als)
  # (without interactions)
  binglm_l_p_pval_ni<-vector(mode="numeric",length=als)
  
  # Quasibinomial GLM results
  qbinglm_l_p_pval<-vector(mode="numeric",length=als)
  qbinglm_l_unp_pval<-vector(mode="numeric",length=als)
  
  # LM (t-test) results: [unpaired]
  m_lm_unp_pval<-vector(mode="numeric",length=als)
  m_lm_unp_w<-vector(mode="numeric",length=als)
  # LM (t-test) results: [paired]
  m_lm_p_pval<-vector(mode="numeric",length=als)
  m_lm_p_w<-vector(mode="numeric",length=als)

  # Other Statistics
  sd_diffs <- vector(mode = "numeric",length = als)
  mean_diffs <- vector(mode = "numeric",length = als)
  diff_s <- vector(mode = "numeric",length = als)
  
  # Simulation parameters
  npops_v<-rep(paste("k=",k,sep=""),length=als)
  if(resc == 0){
    scale<-"var."
  }
  fst_v<-rep(fst,als)
  scale_v<-rep(paste("CT=",scale,sep=""),als)
  sample_v<-paste(npops_v,scale_v,sep="-")
  sample_v<-vector(mode = "numeric",length = als)
  tp_v<-vector(mode="numeric",length=als)
  # Print the simulation parameters
  # if this feature is desired uncomment the (3) following lines.
  #cat(paste("Nsnps: ", als,"Nreps: ",k," NTreatLines: ",l,
  #            " Rescale: ",resc,"Scale: ",scale,"Consistent: ",cons,
  #            "FST: ",fst, "N: ",N,"\n\n"))

  # Initialise progress bar
  pb<-txtProgressBar(min=1,max=als,style=3) 

  # How many independent chromosomes/alleles in the pool (n)
  n=2*N
  i=1
  # Degrees of freedom (for the G-test)
  df<-k-1
  
  # Produce a beta distribution which excludes the limits 0 and 1
  # 0 < x < 1, n = 10000 deviates
  beta_dist <- rbeta(n=10000,shape1=a,shape2=b)
  # Remove and 0s and 1s, we want no fixations in the base
  # population.
  beta_dist <- beta_dist[beta_dist < 1 & beta_dist > 0]
  
  # Produce a coverage distribution which excludes coverage < 16
  # either a negative binomial
  cov_dist<-rnbinom(n=10000,size=2,mu=mcov)
  # or a poisson
  #cov_dist<-rpois(n=10000,lambda = mcov)
  # Remove coverage values < 10
  cov_dist<-cov_dist[cov_dist>=10]
  # If consistency is not forced
  if(cons == 0){
    # "als" is the total number of SNPs
    # "tp_als" is the total number of "true_positives"
    tp_als <- round(prop_tp*als)
    for(tp_snp in seq(1,tp_als)){
      
    }

    # FOR EACH SNP
    while(i <= als){
      if (i <= tp_als){
        tp <- "1"
      }else{
        tp <- "0"
      }
      # START POPULATING CELLS FOR AN ARRAY 
      # Initialise an array with k partial lx2 tables
      m<-array(dim=c(l,2,k))
      allele_counts<-vector()
      # Initialise a vector of allele frequencies
      allele_fs <- vector()
      # Sample a frequency of A allele in base
      # population.
      pA <- sample(beta_dist,size=1,replace=TRUE)
      # Start filling the array cells
      
      # FOR EACH REPLICATE
      # (there are two treatment lines per replicate)
      for(repl in seq(1,k,1)){
        # Vector for "A" allele counts
        Acs<-vector(length=l)
        # Vector for "a" allele counts
        acs<-vector(length=l)
        # Vector for "A" allele frequencies
        Afss<-vector(length=l)
        # Vector for "a" allele frequencies
        afss<-vector(length=l)

        # FOR EACH TREATMENT LINE
        for(treat_line in seq(1,l,1)){
          # Sample an allele frequency from the
          # truncated normal distribution
          fA<-rtruncnorm(n=1,
                         mean=pA,
                         sd=sqrt(fst*pA*(1-pA)),
                         a=0,b=1)
          # if the SNP is a true positive, add the differential or
          # take to fixation.
          if (tp == "1" & treat_line == 2){
            fA <- min(1,fA + selection_diff)
          }
          
          # Print the allele frequencies in the base population and
          # in the population after drift.
          # If this feature is desired uncomment the following (1) lines
          #cat(c("pA: ",pA," fA: ",fA))
          
          # Set the total row counts for each 
          # partial table.
          # CT is either fixed, scaled to neff or sampled from a 
          # poisson distribution with mean mcov.
          if(resc==1){
            if(is.numeric(scale)){
            totc<-scale
            }else if (scale == "neff"){
              totc<-sample(cov_dist,size = 1,replace = TRUE)
            }
          }else if(resc==0){
            totc<-sample(cov_dist,size = 1,replace = TRUE)
          }
          
          
          # FIRST ROUND OF BINOMIAL SAMPLING:
          # ALLELE FREQUENCY IN THE POOL
          # Sample the "A" allele count in the sample from 
          # the population
          # Sample n chromosomes from the end generation of the experimental
          # evolution study using a binomial distribution with probability fA
          endAc<-rbinom(n=1,size=n,prob=fA)
          endf<-endAc/n
          # SECOND ROUND OF BINOMIAL SAMPLING:
          # ALLELE FREQUENCY AMONG THE READS
          # Sample the "A" allele count among the reads using totc and endf
          poolAc<-rbinom(n=1,size=totc,prob=endf)
          poolAf <- poolAc/totc
          poolac<-totc-poolAc
          # if scale == "neff", re-compute totc using the neff function
          if(scale=="neff"){
            totc<-neff(n,totc)
            # Re-compute the "A" allele count
            poolAc<-round(poolAf*totc)
            # Re-compute the "a" allele count
            poolac<-totc-poolAc
          }
          # Save the allele data
          Acs[treat_line]<-poolAc
          acs[treat_line]<-poolac
          Afss[treat_line]<-poolAf
          afss[treat_line]<-poolac/totc
          m[treat_line,1,repl]<-poolAc
          m[treat_line,2,repl]<-poolac
        }
        allele_fs<-c(allele_fs,c(Afss,afss))
        allele_counts<-c(allele_counts,c(Acs,acs))
      }
      #FINISHED POPULATING "CELLS" WITH COUNTS FOR SNP
      
      # Run some checks on the allele frequencies being generated.
      # Some checks are redundant.
      # This is because it does not make sense to try to analyse SNPs
      # where one allele is fixed in every line for example. This would
      # not be done with real data.
      
      
      # Check if *all* the treatment lines are fixed 
      # for one or the other allele. If they are, throw them out and 
      # continue with the loop.
      # This will remove SNPs where alleles are fixed in all populations.
      if(all(allele_fs==1|allele_fs==0)){
        next
      # Check if any SINGLE treatment line has minor allele frequency < 0.05.
      # if TRUE, throw them out and continue with the loop.
      }else if(any(allele_fs < 0.05)){
        next

      # Check if *all* the treatment lines have allele frequencies 
      # < 0.05 (5%) or > 0.95 (95%).
      # if TRUE, throw them out and continue with the loop.
      # This should remove stuff where alleles are very close to fixation
      # in all populations.
      }else if(all(allele_fs<0.05|allele_fs>0.95)){
        next
      }

      # RUN THE STATISTICAL TESTS
      # Calculate differences in A allele between lines.
      Adiffs <- vector()
      for(repl in seq(1,k,1)){
        Adiffs <- c(Adiffs,diff(m[,1,repl]/(m[,1,repl]+m[,2,repl])))
      }
      #print(m) # this is a script tester line.

      # Run the CMH- and Wool-test for the SNP
      # If there are any 0s, add 1 to all cells in m
      # before running cmh-test and woolf-test
      if(any(0 %in% m)){
        cmh_dat <- mantelhaen.test(m+1)
        woo_dat <- woolf.test(m+1)
        
      }else{
        cmh_dat <- mantelhaen.test(m)
        woo_dat <- woolf.test(m)
      }

      # Run the G-test for the SNP
      g_dat <- G_test(m,add=TRUE)
      
      # Create data.frame for the SNP
      d <- get_dat(m,zeroes=1)

      # Run GLM with binomial error distribution: [unpaired]
      # (without interaction)
      binglm_res_ni<-glm(cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$tr_l,
                        family = "binomial")
      
      # (with interaction)
      binglm_res<- tryCatch(anova(glm(
        cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$rep+d$tr_l+d$tr_l:d$rep,
        family = "binomial"),test="LRT"),
        warning=function(w){suppressWarnings(anova(glm(
          cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$rep+d$tr_l+d$tr_l:d$rep,
          family = "binomial"),test="LRT"))})

      # FOR QUASIBINOMIAL GLMS TESTING FOR AN INTERACTION IS NOT POSSIBLE: 
      # WE RUN OUT OF DEGREES OF FREEDOM FOR THE RESIDUALS. 
      # IT WORKS FOR THE BINOMIAL GLMs ONLY BECAUSE DISPERSION IS TAKEN TO BE 1
      # SO WE CAN TEST ON A NORMAL DISTRIBUTION.
      
      # Run GLM with quasibinomial error distribution: [paired]
      qbinglm_res_p<-glm(cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$rep+d$tr_l,
        family = "quasibinomial")

      # Run GLM with quasibinomial error distribution: [unpaired]
      qubinglm_res_unp<-glm(cbind(d$A_Cnt,d$Tot_Cnt-d$A_Cnt)~d$tr_l,
        family = "quasibinomial")

      # Run an lm ("t-test"): [unpaired]
      m_lm_unp <- tryCatch(summary(lm(
        I(A_Cnt/Tot_Cnt)~tr_l,data=d))$coefficients,
        warning=function(w){suppressWarnings(summary(lm(
          I(A_Cnt/Tot_Cnt)~tr_l,data=d))$coefficients)})

      # Run an lm ("t-test"): [paired]
      m_lm_p <- tryCatch(summary(lm(
        I(A_Cnt/Tot_Cnt)~rep+tr_l,data=d))$coefficients,
        warning=function(w){suppressWarnings(summary(lm(
          I(A_Cnt/Tot_Cnt)~rep+tr_l,data=d))$coefficients)})
      
      # Save all the results
      count_dat[i] <- paste(allele_counts,collapse=",")
      af_dat[i] <- paste(allele_fs,collapse=",")

      # CMH-Test results
      woo_pval[i] <- woo_dat["p-value"]
      woo_W[i] <- woo_dat["Woolf"]
      cmh_pval[i] <- cmh_dat$p.value
      cmh_cOR[i] <- cmh_dat$estimate
      cmh_chi[i] <- cmh_dat$statistic

      # G-test results
      g_lxcxpx_pval[i] <- g_dat[5,3]
      g_lxc_pval[i] <- g_dat[3,3]
      g_lxcxpx_g[i] <- g_dat[5,1]
      g_lxc_g[i] <- g_dat[3,1]
      
      # Binomial GLM results: [unpaired]
      # (without interaction)
      row<-nrow(summary(binglm_res_ni)$coefficients)
      binglm_l_p_pval_ni[i]<-summary(binglm_res_ni)$coefficients[row,4]
      
      # (with interaction)
      # Check if the Binomial GLM has produced any warnings
      if(is(binom_res[[2]],"warning")){
        glm_d <- binglm_res[[1]]
        binglm_lxpx_pval[i] <- glm_d[[5]][4]
        binglm_l_pval[i] <- glm_d[[5]][3]
        binglm_lxpx_dev[i] <- glm_d[[2]][4]
        binglm_l_dev[i] <- glm_d[[2]][3]
        binglm_w[i] <- "1"
      }else{
        binglm_d <- binglm_res
        binglm_lxpx_pval[i] <- glm_d[[5]][4]
        binglm_l_pval[i] <- glm_d[[5]][3]
        binglm_lxpx_dev[i] <- glm_d[[2]][4]
        binglm_l_dev[i] <- glm_d[[2]][3]
        binglm_w[i] <- "0"
      }

      # Quasibinomial GLM results: [paired]
      # "treatment" effect will be last p-value
      row<-nrow(summary(quasbinom_res_p)$coefficients)
      qbinglm_l_p_pval[i]<-summary(quasbinom_res_p)$coefficients[row,4]
     
      # Quasibinomial GLM results: [unpaired]
      # "treatment" effect will be last p-value
      row<-nrow(summary(quasbinom_res_unp)$coefficients)
      qbinglm_l_unp_pval[i]<-summary(quasbinom_res_unp)$coefficients[
        row,4]
     
      # LM (t-test) results: [unpaired]
      # Check if the LM (t-test) has produced any warnings
      if(is(m_lm_unp[[2]],"warning")){
        row<-nrow(m_lm_unp[[1]])
        m_lm_unp_pval[i] <- m_lm_unp[[1]][row,4]
        m_lm_unp_w[i] <- "1"
      }else{
        row<-nrow(m_lm_unp)
        m_lm_unp_pval[i] <- m_lm_unp[row,4]
        m_lm_unp_w[i] <- "0"
      }
     
      # LM (t-test) results: [paired]
      # Check if the LM (t-test) has produced any warnings
      if(is(m_lm_p[[2]],"warning")){
        row<-nrow(m_lm_p[[1]])
        m_lm_p_pval[i] <- m_lm_p[[1]][row,4]
        m_lm_p_w[i] <- "1"
      }else{
        row<-nrow(m_lm_p)
        m_lm_p_pval[i] <- m_lm_p[row,4]
        m_lm_p_w[i] <- "0"
      }
     
      sd_diffs[i] <- sd(Adiffs)
      mean_diffs[i] <- mean(Adiffs)

      diff_s[i] <- sd(mean(Adiffs)-Adiffs)#/length(Adiffs)
      tp_v[i] <- tp
      setTxtProgressBar(pb,i)
      
     i=i+1
    }
    snp_sim_dat <- data.frame("allele_counts" = count_dat,
                         "allele_freqs" = af_dat,
                         "df" = rep(df,als),
                         "woo_pval" = woo_pval,
                         "woo_W" = woo_W,
                         "cmh_pval" = cmh_pval,
                         "cmh_cOR" = cmh_cOR,
                         "cmh_chi" = cmh_chi,
                         "g_lxcxpx_pval"=g_lxcxpx_pval,
                         "g_lxc_pval"=g_lxc_pval,
                         "g_lxcxpx_g"=g_lxcxpx_g,
                         "g_lxc_g"=g_lxc_g,
                         "binglm_lxp_pval"=glm_lxpx_pval,
                         "binglm_l_pval"=glm_l_pval,
                         "binglm_l_dev"=glm_l_dev,
                         "binglm_lxp_dev"=glm_lxpx_dev,
                         "binglm_w"=glm_w,
                         "binglm_l_p_pval_ni"=binom_l_p_pval,
                         "lm_unp_pval"=m_lm_unp_pval,
                         "lm_p_pval"=m_lm_p_pval,
                         "qbinglm_p_l_pval"=quasbinom_l_p_pval,
                         "qbinglm_unp_l_pval"=quasbinom_l_unp_pval,
                         "sd_diffs"=sd_diffs,
                         "mean_diffs"=mean_diffs,
                         "diff_s"=diff_s,
                         "npops"=npops_v,
                         "scale"=scale_v,
                         "sample"=sample_v,
                         "fst"=fst_v,
                         "tp"=tp_v)
    snp_sim_dat
  }else{
    # If consistency is forced
    # THIS PART OF SIMULATION FUNCTION WAS NEVER USED
    # BUT EXISTS TO PROVIDE THE OPTION.
    #
    while(i <= als){
      counts <- vector()
      allele_fs <- vector()
      # 
      if(resc == 1){
        totc1 <- scale
        totc2 <- scale
      }else{
        totc1 <- sample(cts,1)
        totc2 <- sample(cts,1)    
      }
      # Sample "A" allele frequency from possible frequencies
      # for line 1
      Af1 <- sample(fs,1)
      # Compute "A" allele count
      Ac1 <- round(Af1*totc1)
      # Compute "a" allele count
      ac1 <- totc1-Ac1
      # Sample "A" allele frequency from possible frequencies
      # for line 2  
      Af2 <- sample(fs,1)
      # Compute "A" allele count
      Ac2 <- round(Af2*totc2)
      # Compute "a" allele count
      ac2 <- totc2-Ac2
      allele_counts <- rep(c(Ac1,Ac2,ac1,ac2),k)
      allele_fs <- c(Af1,(ac1/totc1),Af2,(ac2/totc2))
      m <- array(counts,dim=c(l,2,k))
      # Calculate differences in A allele between lines.
      Adiffs <- vector()
      for(pop in seq(1,k,1)){
        Adiffs <- c(Adiffs,diff(m[,1,pop]/(m[,1,pop]+m[,2,pop])))
      }
      # Run the CMH-test for the SNP
      cmh_dat <- mantelhaen.test(m)
      # Run the Woolf-test for the SNP
      woo_dat <- woolf.test(m)
      # Run the G-test for the SNP
      g_dat <- G_test(m,add=TRUE)
      # Run the Binomial GLM
      res<-tryCatch(d<-get_glm_dat(m,zeroes=1),
                    warning=function(w){
                      suppressWarnings(d<-get_glm_dat(m,zeroes=1))
                                        return(list(d,w))})
      count_dat[i] <- paste(allele_counts,collapse=",")
      af_dat[i] <- paste(allele_fs,collapse=",")
      #bd_p[i] <- bd_dat$p.value
      #bd_chi[i] <- bd_dat$statistic
      woo_p[i] <- woo_dat["p-value"]
      woo_W[i] <- woo_dat["Woolf"]
      cmh_p[i] <- cmh_dat$p.value
      cmh_cOR[i] <- cmh_dat$estimate
      cmh_chi[i] <- cmh_dat$statistic
      g_lxcxpx_p[i] <- g_dat[5,3]
      g_lxc_p[i] <- g_dat[3,3]
      g_lxcxpx_g[i] <- g_dat[5,1]
      g_lxc_g[i] <- g_dat[3,1]
      if(is(res[[2]],"warning")){
        glm_d <- res[[1]]
        glm_lxpx_p[i] <- glm_d[[5]][4]
        glm_l_p[i] <- glm_d[[5]][2]
        glm_lxpx_dev[i] <- glm_d[[2]][4]
        glm_l_dev[i] <- glm_d[[2]][2]
        glm_w[i] <- "1"
      }else{
        glm_d <- res
        glm_lxpx_p[i] <- glm_d[[5]][4]
        glm_l_p[i] <- glm_d[[5]][2]
        glm_lxpx_dev[i] <- glm_d[[2]][4]
        glm_l_dev[i] <- glm_d[[2]][2]
        glm_w[i] <- "0"
      }
      sd_diffs[i] <- sd(Adiffs)
      mean_diffs[i] <- mean(Adiffs)
      sd_abs_diffs[i] <- sd(abs(Adiffs))
      mean_abs_diffs[i] <- mean(abs(Adiffs))
      diff_s[i] <- sd(mean(Adiffs)-Adiffs)#/length(Adiffs))
      setTxtProgressBar(pb, i)
      i=i+1
    }
    snp_sim_dat <- data.frame("allele_counts" = count_dat,
                         "allele_freqs" = af_dat,
                         "df" = rep(df,als),
                         #"bd_p" = bd_p,
                         #"bd_chi" = bd_chi,
                         "woo_p" = woo_p,
                         "woo_W" = woo_W,
                         "cmh_p" = cmh_p,
                         "cmh_cOR" = cmh_cOR,
                         "cmh_chi" = cmh_chi,
                         "g_lxcxpx_p"=g_lxcxpx,
                         "g_lxc_p"=g_lxc,
                         "g_lxcxpx_g"=g_lxcxpx,
                         "g_lxc_g"=g_lxc,
                         "glm_lxp_p"=glm_lxpx_p,
                         "glm_l_p"=glm_l_p,
                         "glm_l_dev"=glm_l_dev,
                         "glm_lxp_dev"=glm_lxpx_dev,
                         "glm_w"=glm_w,
                         "sd_diffs"=sd_diffs,
                         "mean_diffs"=mean_diffs,
                         "sd_abs_diffs"=sd_abs_diffs,
                         "mean_abs_diffs"=mean_abs_diffs,
                         "diff_s"=diff_s)
    snp_sim_dat
    
  }
}

# NOTES: To recreate the contingency table array from the "counts" column:
# copy and paste an entry in the "counts" column to the array() function.
# e.g. m <- array(c([paste_counts_element_here]),dim = c(l,2,k))

# Produce a test dataset
#testsim <- data_gen(als=100,k=4,l=2, 
#                    mincov=16,maxcov=400,N=100,
#                    a=0.2,b=0.2,fst=0.2,cons=0,res=1,scale=100)
#head(testsim,n=1)
#m <- array(c(),dim = c(2,2,4))
#Adiffs <- vector()
#k <- 4
#for(pop in seq(1,k,1)){
#  Adiffs <- c(Adiffs,diff(m[,1,pop]/(m[,1,pop]+m[,2,pop])))
#}
#Adiffs
#m
#mean(Adiffs)
#sd(Adiffs)
#get_glm_dat(m)
#get_dat(m)
