#!/usr/bin/env Rscript

#### ABOUT ####
# Input to this script should be a multisample 
# VCF containing only biallelic sites. The output
# will be a hard-filtered VCF containing putative somatic variants

#### INPUT CHECK ####
suppressMessages(
  if (!require("optparse",quietly = T)){
    install.packages("optparse",quiet = T,repos = "https://cloud.r-project.org")
    suppressMessages(library("optparse",quietly = T,verbose = F))
})

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
    help="REQUIRED: Input VCF containing only biallelic sites [default = %default]", metavar="VCF"),
  make_option(c("-o", "--output"), type="character", default="out.vcf.gz", 
    help="OPTIONAL: Output VCF (gzipped) [default = %default]", metavar="VCF"),
  make_option(c("-c","--cores"), type="integer", default=1,
    help="OPTIONAL: Cores [default = %default]",metavar="Integer")
  )
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Missing input (--input)", call.=FALSE)
}

if (file.exists(opt$output)){
  stop(paste(opt$output,"already exists"), call.=FALSE)
}

#### LOAD PACKAGES ####
cat(">>> Configuring environment <<<\n")

suppressMessages(
  if (!require("tidyverse",quietly = T)){
    install.packages("tidyverse",quiet = T,repos = "https://cloud.r-project.org")
    library("tidyverse",quietly = T,verbose = F)
})
suppressMessages(
  if (!require("stringr",quietly = T)){
    install.packages("stringr",quiet = T,repos = "https://cloud.r-project.org")
    library("stringr",quietly = T,verbose = F)
})
suppressMessages(
  if (!require("vcfR",quietly = T)){
    install.packages("vcfR",quiet = T,repos = "https://cloud.r-project.org")
    library("vcfR",quietly = T,verbose = F)
})
suppressMessages(
  if (!require("MASS",quietly = T)){
    install.packages("MASS",quiet = T,repos = "https://cloud.r-project.org")
    library("MASS",quietly = T,verbose = F)
})
suppressMessages(
  if (!require("pbmcapply",quietly = T)){
    install.packages("pbmcapply",quiet = T,repos = "https://cloud.r-project.org")
    library("pbmcapply",quietly = T,verbose = F)
})
suppressMessages(
  if (!require("VGAM",quietly = T)){
    install.packages("VGAM",quiet = T,repos = "https://cloud.r-project.org")
    library("VGAM",quietly = T,verbose = F)
})

#### FILTERING FUNCTIONS ####
flag_missing_data<-function(tsv,bad_threshold=0.2,name="missing_data",soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Print info
  cat(paste0("Flagging sites where missing data ≥ ",bad_threshold),"\n")
  
  #Create vector of loci where NA genotype calls ≥ threshold
  flagged_sites<-tsv %>% 
    group_by(locus) %>%
    summarise(missing_pct=sum(is.na(GT))/n()) %>%
    filter(missing_pct >= bad_threshold) %>%
    pull(locus) %>%
    unique()
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    return(tsv)
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  } 
}
flag_DP<-function(tsv,low_percentile_cutoff=0.1,cutoff_low=NA,high_percentile_cutoff=0.99,cutoff_high=NA,name="DP",soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Estimate poisson-distributed mean
  site.df<-tsv %>%
    summarise(meanDP=mean(DP,na.rm = T),
              .by = locus)
  
  if (is.na(cutoff_low)){
    pois_mu<-site.df %>%
      pull(meanDP) %>%
      round() %>%
      as.integer() %>%
      fitdistr(.,densfun="poisson") %>%
      .$estimate
    cutoff_low<-qpois(p=low_percentile_cutoff,lambda=pois_mu)
  }
  if (is.na(cutoff_high)){
    pois_mu<-site.df %>%
      pull(meanDP) %>%
      round() %>%
      as.integer() %>%
      fitdistr(.,densfun="poisson") %>%
      .$estimate
    cutoff_high<-qpois(p=high_percentile_cutoff,lambda=pois_mu)
  }
  
  cutoff_low<-max(cutoff_low,8) # DP=8 is the lowest value DP where AD=2 can fail lower-tail binomial test
  cutoff_high<-min(cutoff_high,40) #DP=40 is the default high in Mitchell et al. 2022
  
  #Print info
  cat("Flagging sites for abberrant depth (using Poisson distributed read counts):\n")
  cat(paste0("|--- 10th percentile low cutoff: ",low_percentile_cutoff),"\n")
  cat(paste0("|--- 1st percentile high cutoff: ",high_percentile_cutoff),"\n")
  cat(paste0("|--- Low count cutoff: ",cutoff_low),"\n")
  cat(paste0("|--- High count cutoff: ",cutoff_high),"\n")
  
  #Create vector of loci where mean HET DP < cutoff_low
  site.df<-tsv %>%
    filter(GT == 1) %>%
    summarise(meanDP=mean(DP,na.rm = T),
              .by = locus)
  
  flagged_sites<-site.df %>% 
    filter(meanDP < cutoff_low | meanDP > cutoff_high,
           .by = locus) %>%
    pull(locus) %>% 
    unique()
  
  #Append to vector of loci where mean DP < cutoff_low
  site.df<-tsv %>%
    summarise(meanDP=mean(DP,na.rm = T),
              .by = locus)
  
  
  flagged_sites<-c(flagged_sites,
                   site.df %>%
                     filter(meanDP < cutoff_low | meanDP > cutoff_high,
                            .by = locus) %>%
                     pull(locus) %>%
                     unique()) %>%
    unique()
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    return(tsv)
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  } 
}
flag_het_AB_threshold<-function(tsv,bad_thresh=0.2,ADmin=1,name="het_AB_threshold",soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Print info
  cat(paste0("Flagging sites with AB ≤ ",bad_thresh, " (Minimum variant AD ≥ ",ADmin,")"),"\n")
  
  #Create vector of loci AB < bad_thresh (using only hets with altAD ≥ ADmin)
  site.df<-tsv %>%
    mutate(altAD = as.integer(round(AB*DP))) %>%
    filter(GT==1 & altAD >= ADmin) %>%
    summarise(AD=sum(altAD,na.rm = T),
              DP=sum(DP,na.rm = T),
              .by = locus) %>%
    mutate(AB = AD/DP)
  
  flagged_sites<-site.df %>% 
    filter(AB <= bad_thresh) %>%
    pull(locus) %>%
    unique()
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    return(tsv)
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  } 
}
flag_potential_germline<-function(tsv,cores=1,padj_thresh=0.05,p_correction="fdr",name="potential_germline",soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Print info
  cat("Flagging sites as with evidence of germline heterozygosity:\n")
  cat(paste0("|--- Aggregated AB < 0.5 (lower-tail binomial p-value ≥ ",padj_thresh,", ",p_correction," corrected)"),"\n")
  
  #Create vector of loci where mutant read count among heterozygous samples fails binomial test for 50% allele balance
  site.df<-tsv %>%
    mutate(altAD = as.integer(round(AB*DP))) %>%
    summarise(AD=sum(altAD,na.rm = T),
              DP=sum(DP,na.rm = T),
              .by = locus)
  
  site.df<-pbmclapply(1:nrow(site.df), function(x) {
    df<-site.df[x,]
    #print(df)
    binom<-binom.test(x=df$AD,
                      n=df$DP,
                      alternative="less",
                      p=0.5)
    df$binom.p<-binom$p.value
    return(df)
  },mc.cores = as.integer(cores)) %>% 
    bind_rows()
  
  site.df$padj<-p.adjust(site.df$binom.p,method = p_correction)
  
  flagged_sites<-site.df %>% 
    filter(padj >= padj_thresh) %>%
    pull(locus) %>%
    unique()
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    return(tsv)
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  } 
}
flag_low_HET_variance<-function(tsv,name="low_het_var",soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Print info
  cat("Flagging sites where all samples are HET or none are HET\n")
  
  #Create vector of loci where all called genotypes are het or none are het
  flagged_sites<-tsv %>% 
    filter(sum(GT == 1,na.rm = T) == 0,
           .by = locus) %>%
    pull(locus) %>% 
    unique()
  flagged_sites<-c(flagged_sites,
                   tsv %>%
                     filter((sum(GT == 1,na.rm = T)/sum(!is.na(GT))) == 1,
                            .by=locus) %>%
                     pull(locus) %>%
                     unique())
  flagged_sites<-unique(flagged_sites)
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    return(tsv)
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  } 
}
flag_aberrant_AB<-function(tsv,cores=1,ADmin=3,padj_cutoff=0.05,p_correction="fdr",name="aberrant_AB",soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Print info
  cat(paste0("Flagging sites with abberrant AB among HET samples (with variant AD ≥ ",ADmin,")"),"\n")
  cat(paste0("|--- Lower-tail binomial p-value ≤ ",padj_cutoff,", ",p_correction," corrected)"),"\n")
  
  #Create vector of loci where VAF fails lower-tail binomial test (using only samples with altAD >= ADmin)
  site.df<-tsv %>%
    mutate(altAD = as.integer(round(AB*DP))) %>%
    filter(altAD >= ADmin) %>%
    #filter(GT==1) %>%
    summarise(AD=sum(altAD,na.rm = T),
              DP=sum(DP,na.rm = T),
              .by = locus)
  #return(site.df)
  
  site.df<-pbmclapply(1:nrow(site.df), function(x) {
    df<-site.df[x,]
    binom<-binom.test(x=df$AD,
                      n=df$DP,
                      alternative="less",
                      p=0.5)
    df$binom.p<-binom$p.value
    return(df)
  },mc.cores = as.integer(cores)) %>%
    bind_rows()
  
  site.df$padj<-p.adjust(site.df$binom.p,method = p_correction)
  
  flagged_sites<-site.df %>%
    filter(padj <= padj_cutoff) %>%
    pull(locus) %>%
    unique()
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    return(tsv)
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  } 
}
estimateRho_gridml = function(AD_vec,DP_vec,soft=T) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(AD_vec)/sum(DP_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=AD_vec, size=DP_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
} #Function call for the beta-binomial artifact filter
flag_artifact<-function(tsv,cores=1,bad_threshold=0.1,name="potential_artifact",return_rho=F,soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Print info
  cat(paste0("Flagging sites with underdispersed AB (ML beta-binomial rho ≤ ",bad_threshold,")"),"\n")
  
  #Remove sites if ML estimate for beta-binomial rho is <= 0.1
  ##Create base dataframe that recodes DP=0 as DP=1
  df<-tsv %>%
    mutate(altAD = as.integer(round(AB*DP)))
  df[df$DP == 0 | is.na(df$DP),]$DP<-1
  
  ##Prepare two n,j matrices, where n = samples, j = sites, and cells = AD (or DP)
  df.AD<- df %>% 
    dplyr::select(locus,sample,altAD) %>%
    pivot_wider(names_from = sample,values_from = altAD) %>%
    as.data.frame()
  rownames(df.AD)<-df.AD$locus
  df.AD<-df.AD %>% dplyr::select(!locus)
  
  df.DP<- df %>% 
    dplyr::select(locus,sample,DP) %>%
    pivot_wider(names_from = sample,values_from = DP) %>%
    as.data.frame()
  rownames(df.DP)<-df.DP$locus
  df.DP<-df.DP %>% dplyr::select(!locus)
  
  ##Prepare results df
  df<-data.frame(locus=rownames(df.DP),rhoest=NA)
  
  ##Compute best-fit rho
  df<-pbmclapply(1:nrow(df), function(x) {
    var_row<-df[x,]
    rho<-estimateRho_gridml(AD_vec = as.numeric(df.AD[x,]),
                            DP_vec = as.numeric(df.DP[x,]))
    var_row$rhoest<-rho
    return(var_row)
  },mc.cores = as.integer(cores)) %>%
    bind_rows()
  
  flagged_sites<-df %>% 
    filter(rhoest <= bad_threshold) %>%
    pull(locus) %>%
    unique()
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    if ( return_rho ){
      tsv<-left_join(tsv,
                     df %>% dplyr::select(locus,rhoest)) %>%
        rename("MLrho"="rhoest")
      return(tsv)
    } else {
      return(tsv)
    }
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  }
}
flag_high_alt_in_ref<-function(tsv,cores=1,padj_thresh=0.05,p_correction="fdr",name="high_alt_in_ref",soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Print info
  cat("Flagging sites with high variant allele count among homREF samples\n")
  cat(paste0("|--- Upper-tail poisson test p-value ≤ ",padj_thresh,", ",p_correction," corrected)"),"\n")
  
  #Create vector of loci where mutant read among homozygous ref samples fails poisson test where number of events = AD(alt) and rate = 0.01*depth
  total.df<-tsv %>% 
    summarise(DP=sum(DP,na.rm = T),
              .by = locus)
  ref.df<-tsv %>% 
    mutate(altAD = as.integer(round(AB*DP))) %>%
    filter(GT == 0) %>%
    summarise(homref_altAD=sum(altAD,na.rm = T),
              .by = locus)
  site.df<-left_join(ref.df,total.df,by="locus") %>% filter(!homref_altAD == 0)
  
  site.df<-pbmclapply(1:nrow(site.df), function(x) {
    df<-site.df[x,]
    #print(df)
    pois<-poisson.test(x=df$homref_altAD,
                       r=0.01*df$DP,
                       alternative="greater",
                       conf.level=0.90)
    df$pois.p<-pois$p.value
    return(df)
  },mc.cores = as.integer(cores)) %>% 
    bind_rows()
  
  #site.df$q<-qvalue(site.df$pois.p,fdr.level = fdr)$qvalues
  site.df$padj<-p.adjust(site.df$pois.p,method = p_correction)
  
  flagged_sites<-site.df %>% 
    filter(padj <= padj_thresh) %>%
    pull(locus) %>%
    unique()
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    return(tsv)
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  } 
}
flag_lowGQ<-function(tsv,bad_GQthresh=30,name="lowGQ",soft=T){
  #If first filtering step, initiate FILTER column annotation
  if ((! "FILTER" %in% colnames(tsv))){
    tsv$FILTER<-"PASS"
  }
  
  #Print info
  cat(paste0("Flagging sites where mean GQ < ",bad_GQthresh),"\n")
  
  #Pull sites where mean GQ <= bad_GQthresh 
  site.df<-tsv %>%
    filter(DP > 0) %>% 
    summarise(meanGQ=mean(GQ),
              .by = locus)
  
  flagged_sites<-site.df %>% 
    filter(meanGQ < bad_GQthresh) %>%
    pull(locus) %>%
    unique()
  
  #Pull sites where mean GQ(het) <= bad_GQthresh
  site.df<-tsv %>%
    filter(GQ==1) %>%
    summarise(meanGQ=mean(GQ),
              .by = locus)
  
  flagged_sites<-c(flagged_sites,
                   site.df %>%
                     filter(meanGQ < bad_GQthresh) %>%
                     pull(locus) %>%
                     unique()) %>%
    unique()
  
  if (( length(flagged_sites) == 0 )){
    return(tsv)
  }
  
  #Annotate FILTER field
  if ( soft == T ){
    tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER<-paste(tsv[tsv$locus %in% flagged_sites & !tsv$FILTER == "PASS",]$FILTER,
                                                                             name,
                                                                             sep = ",")
    tsv[tsv$locus %in% flagged_sites & tsv$FILTER == "PASS",]$FILTER<-name
    return(tsv)
  }
  #Hard filter
  if ( soft == F ){
    return(tsv[!tsv$locus %in% flagged_sites,] %>% filter(FILTER == "PASS"))
  } 
}

#### READ INPUT ####
cat(">>> Reading input <<<\n")
in_vcf<-read.vcfR(opt$input,verbose=T)
nc<-as.integer(opt$cores)

#### CONVERT INPUT TO FILTERABLE TABLE ####
#Load loci from FIX matrix
cat(">>> Preparing VCF for filtering <<<\n")
loci<-in_vcf@fix %>%
  .[,c("CHROM","POS")] %>%
  apply(.,1,function(x) paste(x[1],x[2],sep=":"))

#Load FORMAT data
vcf.df<-extract_gt_tidy(in_vcf,
                format_fields = c("GT","AD","DP","GQ"),
                gt_column_prepend = "",
                alleles = F) %>%
  arrange(Key) %>%
  mutate(locus = loci[Key]) %>%
  dplyr::select(locus,Indiv,GT,AD,DP,GQ)

#Recode genotypes
vcf.df<-vcf.df %>% 
  mutate(GT=case_when(GT %in% c("0/0","0|0") ~ as.factor(0),
                      GT %in% c("0/1","0|1","1/0","1|0") ~ as.factor(1),
                      GT %in% c("1/1","1|1") ~ as.factor(2),
                      GT == "./." ~ NA))

#Compute VAF as AB
vcf.df$ALT<-str_split(vcf.df$AD,pattern = ",") %>%
  purrr::map(.,2) %>%
  unlist() %>%
  as.integer()
vcf.df$AB<-as.integer(vcf.df$ALT)/as.integer(vcf.df$DP)
vcf.df[is.nan(vcf.df$AB),]$AB<-NA

#Subset to relevant fields, group, rename
vcf.df<-vcf.df %>% 
  dplyr::select(locus,Indiv,GT,AB,DP,GQ)
vcf.df<-vcf.df %>% rename("sample"="Indiv")
  
#### APPLY FILTERS ####
# Order of filters prioritizes functions that share information across sites (`flag_DP`) or require fdr correction (`flag_potential_germline`)
# Following these filters, hard-filtering is applied to reduce compute time
cat(">>> Filtering VCF <<<\n")
vcf.df.filtered<-flag_DP(vcf.df,
    cutoff_high = 40,
    soft = F) %>%
  flag_potential_germline(cores = nc,
   p_correction = "none",
   padj_thresh = 0.001, 
   soft = F) %>%
  flag_aberrant_AB(cores = nc,
   p_correction = "none",
   ADmin = 1,
   soft = F) %>%
  flag_high_alt_in_ref(cores = nc,
   p_correction = "none",
   soft = F) %>%
  flag_missing_data(bad_threshold = 0.2,
    soft = F) %>%
  flag_het_AB_threshold(ADmin = 1,
    bad_thresh = 1/3,
    soft = F) %>%
  flag_low_HET_variance(soft = F) %>%
  flag_lowGQ(bad_GQthresh = 30, 
   soft = F) %>%
  flag_artifact(cores = nc,
    soft = F)

PASS<-unique(vcf.df.filtered$locus)
out_vcf<-in_vcf[loci %in% PASS]

#### PREPARE OUTPUT ####
cat(paste0(">>> Writing output to ",opt$output," <<<\n"))
out_vcf@meta[length(out_vcf@meta)+1]<-paste0("##SomaticFilterAmalgam.R hard-filtered; Date=",date())
write.vcf(out_vcf,file = opt$output)

cat(">>> Filtering complete <<<\n")
