# Read Q matrix from CellPhy *.bestModel output
read_cellphy_model <- function(path_to_cellphy_bestModel) {
  RateCats <- c("Zero", "AC", "AG", "AT", "CG", "CT", "GT")
  States <- c("AA", "CC", "GG", "TT", "AC",
              "AG", "AT", "CG", "CT", "GT")
  
  #Read bestmodel file
  txt <- read_file(file = as.character(path_to_cellphy_bestModel))
  
  #Digest rates matrix
  Rates <- txt %>%
    str_split(pattern = "\\+") %>%
    purrr::map(., 1) %>%
    unlist() %>%
    str_remove_all(pattern = "GT10") %>%
    str_remove_all(pattern = "[{}]") %>%
    str_split(pattern = "/") %>%
    unlist() %>%
    as.numeric()
  names(Rates) <- RateCats
  
  alpha <- Rates["AC"]
  beta <- Rates["AG"]
  gamma <- Rates["AT"]
  kappa <- Rates["CG"]
  lambda <- Rates["CT"]
  mu <- Rates["GT"]
  
  #Digest state frequencies matrix
  Freq <- txt %>%
    str_split(pattern = "\\+") %>%
    purrr::map(., 2) %>%
    unlist() %>%
    str_remove_all(pattern = "FU") %>%
    gsub("\\}.*", "", .) %>%
    str_remove_all(pattern = "\\{") %>%
    str_split(pattern = "/") %>%
    unlist() %>%
    as.numeric()
  names(Freq) <- States
  
  #Generate matrix columns c("AA","CC","GG","TT","AC","AG","AT","CG","CT","GT")
  AAvec <- c(NA,
             Rates["Zero"],
             Rates["Zero"],
             Rates["Zero"],
             alpha*Freq["AA"],
             beta*Freq["AA"],
             gamma*Freq["AA"],
             Rates["Zero"],
             Rates["Zero"],
             Rates["Zero"])
  AAvec <- unname(AAvec)
  
  CCvec <- c(Rates["Zero"],
             NA,
             Rates["Zero"],
             Rates["Zero"],
             alpha*Freq["CC"],
             Rates["Zero"],
             Rates["Zero"],
             kappa*Freq["CC"],
             lambda*Freq["CC"],
             Rates["Zero"])
  CCvec <- unname(CCvec)
  
  GGvec <- c(Rates["Zero"],
             Rates["Zero"],
             NA,
             Rates["Zero"],
             Rates["Zero"],
             beta*Freq["GG"],
             Rates["Zero"],
             kappa*Freq["GG"],
             Rates["Zero"],
             mu*Freq["GG"])
  GGvec <- unname(GGvec)
  
  TTvec <- c(Rates["Zero"],
             Rates["Zero"],
             Rates["Zero"],
             NA,
             Rates["Zero"],
             Rates["Zero"],
             gamma*Freq["TT"],
             Rates["Zero"],
             lambda*Freq["TT"],
             mu*Freq["TT"])
  TTvec <- unname(TTvec)
  
  ACvec <- c(alpha*Freq["AC"],
             alpha*Freq["AC"],
             Rates["Zero"],
             Rates["Zero"],
             NA,
             kappa*Freq["AC"],
             lambda*Freq["AC"],
             beta*Freq["AC"],
             gamma*Freq["AC"],
             Rates["Zero"])
  ACvec <- unname(ACvec)
  
  AGvec <- c(beta*Freq["AG"],
             Rates["Zero"],
             beta*Freq["AG"],
             Rates["Zero"],
             kappa*Freq["AG"],
             NA,
             mu*Freq["AG"],
             alpha*Freq["AG"],
             Rates["Zero"],
             gamma*Freq["AG"])
  AGvec <- unname(AGvec)
  
  ATvec <- c(gamma*Freq["AT"],
             Rates["Zero"],
             Rates["Zero"],
             gamma*Freq["AT"],
             lambda*Freq["AT"],
             mu*Freq["AT"],
             NA,
             Rates["Zero"],
             alpha*Freq["AT"],
             beta*Freq["AT"])
  ATvec <- unname(ATvec)
  
  CGvec <- c(Rates["Zero"],
             kappa*Freq["CG"],
             kappa*Freq["CG"],
             Rates["Zero"],
             beta*Freq["CG"],
             alpha*Freq["CG"],
             Rates["Zero"],
             NA,
             mu*Freq["CG"],
             lambda*Freq["CG"])
  CGvec <- unname(CGvec)
  
  CTvec <- c(Rates["Zero"],
             lambda*Freq["CT"],
             Rates["Zero"],
             lambda*Freq["CT"],
             gamma*Freq["CT"],
             Rates["Zero"],
             alpha*Freq["CT"],
             mu*Freq["CT"],
             NA,
             kappa*Freq["CT"])
  CTvec <- unname(CTvec)
  
  GTvec <- c(Rates["Zero"],
             Rates["Zero"],
             mu*Freq["GT"],
             mu*Freq["GT"],
             Rates["Zero"],
             gamma*Freq["GT"],
             beta*Freq["GT"],
             lambda*Freq["GT"],
             kappa*Freq["GT"],
             NA)
  GTvec <- unname(GTvec)
  
  Q <- matrix(c(AAvec, CCvec, GGvec, TTvec, ACvec,
                AGvec, ATvec, CGvec, CTvec, GTvec),
              nrow = 10, ncol = 10, byrow = F)
  Q[1,1] <- sum(Q[1,],na.rm = T) * -1
  Q[2,2] <- sum(Q[2,],na.rm = T) * -1
  Q[3,3] <- sum(Q[3,],na.rm = T) * -1
  Q[4,4] <- sum(Q[4,],na.rm = T) * -1
  Q[5,5] <- sum(Q[5,],na.rm = T) * -1
  Q[6,6] <- sum(Q[6,],na.rm = T) * -1
  Q[7,7] <- sum(Q[7,],na.rm = T) * -1
  Q[8,8] <- sum(Q[8,],na.rm = T) * -1
  Q[9,9] <- sum(Q[9,],na.rm = T) * -1
  Q[10,10] <- sum(Q[10,],na.rm = T) * -1
  
  colnames(Q) <- States
  rownames(Q) <- States
  
  return(Q)
}

# Compute conditional probabilities from Phred score vector (PL -> GP)
unphred<-function(phredscore){
  prob <- 10^(-phredscore/10)
  prob_norm <- prob/sum(prob)
  return(prob_norm)
}

# Read a vcf and generate a genotype probability matrix for each SNP site.
# This function assumes that the VCF contains diploid genotypes and
# has been loaded with vcfR::read.vcfR() to create a `vcfR-class` object.
compile_gt_states.snp <- function(vcf, mat_list = TRUE){
  print("Reading vcf data...")
  tidy_vcf <- vcfR2tidy(vcf,
                        gt_column_prepend = "",
                        allele.sep = ",",
                        format_fields = c("GT", "GQ", "PL"))
  df <- left_join(tidy_vcf$gt,
                  tidy_vcf$fix %>%
                    select(ChromKey, CHROM) %>%
                    unique(),
                  by = "ChromKey") %>%
    select(CHROM,POS,Indiv,GT,GQ,PL)
  df <- left_join(df,
                  tidy_vcf$fix %>%
                    select(CHROM,POS,REF,ALT,AC,AF),
                  by = c("CHROM", "POS"))
  
  print("Computing genotype probabilities...")
  df <- df %>%
    mutate(locus = paste(CHROM, POS, sep = "_"),
           gt_str = case_when(GT %in% c("0/0", "0|0") ~ paste0(REF, REF),
                              GT %in% c("1/0", "1|0", "0/1", "0|1") ~ paste0(REF, ALT),
                              GT %in% c("1/1","1|1") ~ paste0(ALT,ALT),
                              is.na(GT) ~ NA)) %>%
    rename(clone = Indiv) %>%
    select(locus, CHROM, POS, REF, ALT, AC, AF, clone, GT, gt_str, GQ, PL) %>%
    separate(PL, sep = ",", into = c("PL_ref", "PL_het", "PL_alt")) %>%
    mutate(GQ = as.integer(GQ),
           PL_ref = as.integer(PL_ref),
           PL_het = as.integer(PL_het),
           PL_alt = as.integer(PL_alt),
           AF = as.numeric(AF),
           AC = as.integer(AC)) %>%
    rowwise() %>%
    mutate(Pref = unphred(c(PL_ref, PL_het, PL_alt))[1],
           Phet = unphred(c(PL_ref, PL_het, PL_alt))[2],
           Palt = unphred(c(PL_ref, PL_het, PL_alt))[3])
  
  df <- df %>%
    arrange(CHROM, POS) %>%
    select(locus, REF, ALT, gt_str, clone, Pref, Phet, Palt) %>% 
    as.data.table()
  if (mat_list == FALSE){
    return(df) #Return a long data frame with genotype probabilities for each sample per variant
  } else {
    print("Assembling genotype state matrices for SNPs...")
    df.snp<-df %>%
      group_by(locus) %>%
      filter(REF %in% c("A", "C", "G", "T") & ALT %in% c("A", "C", "G", "T"))
    df.snp$gt_str <- str_split(df.snp$gt_str, pattern = "") %>% 
      pblapply(., function(x) {sort(x) %>% paste0(., collapse = "")}) %>%
      unlist()
    
    gt_list.snp<-df.snp %>% 
      group_by(locus) %>% 
      group_split() %>%
      pblapply(., function(x){
        ref <- head(x, n=1) %>% pull(REF)
        alt <- head(x, n=1) %>% pull(ALT)
        colnames(x)<-c("locus", "REF", "ALT", "gt_str", "clone",
                       paste0(c(ref, ref), collapse = ""),
                       paste0(sort(c(ref, alt)), collapse = ""),
                       paste0(c(alt, alt), collapse = ""))
        return(x %>% 
                 select(!c(REF, ALT, gt_str)) %>% 
                 as.data.table())
      }) %>%
      pblapply(.,function(x){
        cols_to_add<-setdiff(c("AA", "CC", "GG", "TT", "AC",
                               "AG", "AT", "CG", "CT", "GT"), colnames(x)[3:5])
        x[,cols_to_add] = 0
        
        x <- x %>% 
          select(locus,
                 clone,
                 c("AA", "CC", "GG", "TT", "AC",
                   "AG", "AT", "CG", "CT", "GT"))
        
        return(x)
      })
    names(gt_list.snp) <- df.snp %>% group_keys() %>% pull()
    
    return(gt_list.snp) 
  }
}

# Run stochastic character mapping for a single variant
runSCM_single<-function(locus_str, tree, gt_state_list, Qmat,
                        reduced = FALSE, reps = 100,
                        root_state = "equal", cores = 1) {
  if (length(locus_str) > 1){
    return(print("ERROR: >1 variant used for input"))
  }
  
  gt_matrix<-gt_state_list[[locus_str]] %>% #Subset input to genotype priors
    select(!locus) %>%
    column_to_rownames(var = "clone")
  
  if (reduced == TRUE) { #Optionally, remove all genotype states without evidence of presence
    gt_matrix <- gt_matrix %>% 
      select(where(~sum(.) != 0))
    Qmat <- Qmat[colnames(gt_matrix), colnames(gt_matrix)]
  }
  
  #Prep cores and print run info
  cl <- makeSOCKcluster(rep("localhost", cores))
  print(paste0("Performing stochastic character mapping on ", locus_str, "."))
  print(paste0("Cores: ", cores))
  print(paste0("Replicates: ", reps))
  print(paste0("Reduced Qmatrix: ", reduced))
  
  #Run simmap
  scm <- clusterApply(cl,
                      x=replicate(cores, as.matrix(gt_matrix), simplify = FALSE),
                      fun = make.simmap,
                      tree = tree,
                      Q = Qmat,
                      pi = root_state,
                      nsim = as.integer(round(reps/cores)))
  scm <- do.call("c", scm)
  
  #Reclass simmap if not properly classed by output above
  if(!("multiSimmap" %in% class(scm))) {
    class(scm) <- c("multiSimmap", class(scm))
  }
  
  stopCluster(cl)
  return(scm)
}

#Output a summary for SNP SCM, optionally plot
summarise_scm.snp<-function(multiSimmap, PPthreshold = 0.95, plot = FALSE,
                            legend = TRUE, title = NA, scale = TRUE){
  #load primary df from input
  scm_summary <- summary(multiSimmap) #longest step
  scm_summary <- scm_summary$ace %>% as.data.frame()
  tree <- multiSimmap[[1]]
  
  #annotate posteriors with consensus state and node numbers (renumbering tip names to node number)
  scm_summary$constate <- apply(scm_summary, 1, pull_consensus_state, threshold = PPthreshold)
  scm_summary$node <- rownames(scm_summary)
  scm_summary[(tree$Nnode+1):nrow(scm_summary),]$node<-1:(length(tree$tip.label))
  
  #pull vector of node consensus states
  genotype_states <- arrange(scm_summary, as.integer(node)) %>% 
    pull(constate)
  
  #format state posteriors df output 
  out<-scm_summary %>% 
    arrange(as.integer(node)) %>% 
    select(!c(node,constate))
  rownames(out) <- NULL
  cols_to_add<-setdiff(c("AA", "CC", "GG", "TT", "AC", 
                         "AG", "AT", "CG", "CT", "GT"),
                       colnames(out))
  out[,cols_to_add] = 0
  out <- out %>% 
    select("AA", "CC", "GG", "TT", "AC", 
           "AG", "AT", "CG", "CT", "GT")
  
  #pull consensus state posteriors (in same order as out)
  constate_posteriors<-apply(out, 1, max)
  
  #plotting
  if ( plot == T ){
    cols<-setNames(brewer.pal(n = 10, name = "Set3"),
                   c("AA", "CC", "GG", "TT", "AC", "AG", "AT", "CG", "CT", "GT"))
    plot(summary(multiSimmap),
         type = "phylogram",
         direction = "downwards",
         colors = cols,
         fsize = 1,
         ftype = "off",
         lwd = 1,
         cex=c(0.25, 0.25),
         cex=c(0.5, 0.5),
         mar=c(2, 0.1, 2, 0.1),
         outline=F)
    title(main = title)
    if ( legend == T ){
      legend("bottom",
             horiz = T,
             xpd = T,
             legend = c("AA", "CC", "GG", "TT", "AC", "AG", "AT", "CG", "CT", "GT"),
             xjust = 0.5,
             yjust = 0.5,
             pch = 22,
             pt.cex = 2,
             pt.bg = cols,
             inset = c(0,-0.02),
             bty = "n",
             cex = 1)
    }
    if ( scale == T ){
      add.scale.bar(length = 0.01)
    }
    if ( sum(detect_state_changes(tree$edge, genotype_states)) > 0) {
      edgelabels(text="",
                 edge = which(detect_state_changes(tree$edge,genotype_states) == 1),
                 col = "black",
                 frame = "none",
                 pch = 18,
                 cex = 1.5,
                 adj = c(0.5, 0.5))
    }
  }
  
  #output
  return(list("scm_summary" = out,
              "assigned" = detect_state_changes(tree$edge,genotype_states),
              "constate_PPs" = constate_posteriors,
              "QlogL" = multiSimmap[[1]]$logL[1]))
}

#Return the consensus state (>= threshold) for a single row of summary(multiSimmap)$ace
pull_consensus_state <- function(row, threshold) {
  if (max(row) >= threshold) {
    names(row)[which.max(row)]
  } else {
    NA
  }
} 

#Returns a vector (in same order as tree_edges, where tree_edges is <phylo>$edge) where state changes have occurred
detect_state_changes<-function(tree_edges, states){
  if ( is.na(states[tree_edges[1,][1]]) ){
    return(rep(0, length(states)))
  } # If root state is uncertain, return 0's for all branches & exits loop
  root_state<-states[1] #Pull high confidence root state
  
  apply(tree_edges,1,function(x){
    while ( is.na(states[x[1]]) ) {
      x[1] <- tree_edges[tree_edges[,2] == x[1]][1] # Pull the parent node to the current parent node
    } # If parent state uncertain, assign parent state to most recent high confidence state
    
    if ( is.na(states[x[2]]) ) {
      return(0)
    } # If daughter state is uncertain, no state change
    
    if ( states[x[1]] == states[x[2]] ){
      return(0)
    } # If parent state == daughter state, no state change
    
    else {
      return(1)
    }
  })
} 

