#!/usr/bin/env Rscript

# Test for evidence of clonal mixture using binomial mixture modeling and plot somatic VCF stats

# Load packages
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))

# Load functions (from Coorens et al. 2024 (Nature Protocols))
binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6, mode="Full"){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC) 
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}
dbinomtrunc = function(x, size, prob, minx=4) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}
estep = function(x,size,p.vector,prop.vector,ncomp, mode){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i])
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}
mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)   
}
em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}
apply_mix_model=function(NV,NR,plot=T,min_clonal_mut=35){
  peak_VAF=rep(0,ncol(NV))
  names(peak_VAF)=colnames(NV)
  autosomal=!grepl("X|Y",rownames(NV))

  out.df <- data.frame(sample = as.character(),
                       component = as.integer(),
                       prop = as.numeric(),
                       p = as.numeric())

  for(s in colnames(NV)){
    muts_include=NV[,s]>3&autosomal&!is.na(NV[,s])
    if(sum(muts_include,na.rm = TRUE)>5){
      NV_vec=NV[muts_include,s]
      NR_vec=NR[muts_include,s]
      res=binom_mix(NV_vec,NR_vec,mode="Truncated",nrange=1:3)
      for (i in 1:res$n){
        out.df <- rbind(out.df, data.frame(sample = s,
                                            component = i,
                                            prop = res$prop[i],
                                            p = res$p[i]))
      }
      
      if(plot){
        pdf(paste0(output_dir,s,"_binom_mix.pdf"))
        p=hist(NV_vec/NR_vec,breaks=20,xlim=c(0,1),col='gray',freq=F,xlab="Variant Allele Frequency",
               main=paste0(s,", (n=",length(NV_vec),")"))
        cols=c("red","blue","green","magenta","cyan")
        
        y_coord=max(p$density)-0.5
        y_intv=y_coord/5
        
        for (i in 1:res$n){
          depth=rpois(n=5000,lambda=median(NR_vec))
          sim_NV=unlist(lapply(depth,rbinom,n=1,prob=res$p[i]))
          sim_VAF=sim_NV/depth
          sim_VAF=sim_VAF[sim_NV>3]
          dens=density(sim_VAF)
          lines(x=dens$x,y=res$prop[i]*dens$y,lwd=2,lty='dashed',col=cols[i])
          y_coord=y_coord-y_intv/2
          text(y=y_coord,x=0.9,label=paste0("p1: ",round(res$p[i],digits=2)))
          segments(lwd=2,lty='dashed',col=cols[i],y0=y_coord+y_intv/4,x0=0.85,x1=0.95)
        }
        dev.off()
      }
#       peak_VAF[s]=max(res$p[(res$prop*length(res$Which_cluster))>min_clonal_mut])
    }
  }
  return(out.df)
}

# Read commandline arguments
cat("Reading commandline arguments...\n")
args <- commandArgs(trailingOnly = TRUE)
input_af <- args[1]
input_query_stats <- args[2]
input_sensitivity <- args[3]
input_psc <- args[4]
output_tsv <- args[5]
output_jpg <- args[6]

# Load & prep data
cat("Loading and preparing data...\n")
af <- read_tsv(input_af)
vaf.df <- read_tsv(input_query_stats)
vaf.df <- vaf.df %>% 
                mutate(AD = as.integer(round(DP * VAF)),
                       exp_AD = as.integer(round(DP * 0.5)))

sens <- read_tsv(input_sensitivity)
psc <- read_tsv(input_psc) %>%
        dplyr::select(3, 6, 10)
colnames(psc) <- c("sample", "nHets","varDP")
psc_sens.df <- left_join(psc, sens)

# BMM testing for sample mixture
nv <- vaf.df %>%
  dplyr::select(site,sample,AD) %>%
  pivot_wider(names_from = sample, values_from = AD) %>%
  as.data.frame() %>%
  column_to_rownames(var = "site")

nr <- vaf.df %>%
  dplyr::select(site,sample,DP) %>%
  pivot_wider(names_from = sample, values_from = DP) %>%
  as.data.frame() %>%
  column_to_rownames(var = "site")

bmm.df <- apply_mix_model(nv, nr, plot = FALSE)
bmm.df <- slice_max(bmm.df, prop, by = sample) %>%
                mutate(major_component = TRUE) %>%
                right_join(bmm.df) %>%
                mutate(major_component = ifelse(is.na(major_component),
                                                FALSE,
                                                major_component)) %>%
                arrange(sample, component)

write_tsv(bmm.df %>% filter(major_component == TRUE), output_tsv)

# Plot AFS
cat("Plotting allele frequency spectrum...\n")
g1 <- ggplot(af, aes(x = AF)) +
        geom_histogram(aes(y = after_stat(count) / 1000),
                        binwidth = 0.01, 
                        fill = "black", 
                        color = "transparent") +
        scale_x_continuous(limits = c(0, 1)) +
        theme_cowplot() +
        theme(axis.text = element_text(size = 6)) +
        labs(x = "AF", y = "Count (thousands)")

# Plot VAF distribution per sample
cat("Plotting VAF distribution per sample...\n")
g2 <- ggplot(vaf.df %>% filter(VAF > 0), 
             aes(x = sample, y = VAF)) +
        geom_hline(yintercept = seq(0, 1, 0.1),
                linetype = 1,
                color = "gray90") +
        geom_boxplot(outlier.alpha = 0, fill = "gray90") +
        geom_hline(yintercept = 0.5, linetype = 2, color = "#1d1d72") +
        scale_y_continuous(breaks = seq(0,1,0.1)) +
        theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, 
                                        hjust = 1, 
                                        vjust = 0.5, 
                                        size = 6)) +
        labs(x = "Sample", y = "VAF")

# Plot VAF x DP heat diagram
cat("Plotting VAF x DP heat diagram...\n")
g3 <- ggplot(vaf.df %>% filter(VAF > 0),
        aes(x = VAF, y = DP)) +
        geom_bin2d(aes(fill = after_stat(ndensity))) +
        scale_fill_viridis_c(option = "magma", guide = "none") +
        scale_y_continuous(breaks = seq(0, max(vaf.df %>% filter(VAF > 0) %>% pull(DP), by = 10))) +
        theme_cowplot() +
        theme(panel.background = element_rect(fill = "black"),
                axis.text = element_text(size = 6)) +
        labs(x = "VAF", y = "DP")

# Plot relationship between sensitivity and PSC nHets
cat("Plotting relationship between sensitivity and nHets...\n")
g4 <- ggplot(psc_sens.df, aes(x = nHets, y = sensitivity)) +
        geom_point(alpha = 0.75) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_x_continuous(limits = c(0, max(psc_sens.df$nHets))) +
        theme_cowplot() +
        theme(axis.text = element_text(size = 6)) +
        labs(x = "HET sites per sample (count)", y = "Sample sensitivity")

# Plot modelled BMM for each sample
## Create dataframe for plotting
sim_vaf.df <- split(bmm.df, bmm.df$sample) %>%
  lapply(function(x) {
    sample_name <- unique(x$sample)
    depth <- rpois(n = 5000, lambda = median(nr[, sample_name], na.rm = TRUE))
    lapply(1:nrow(x), function(y) {
        df <- x[y, ]

        sim_NV <- unlist(lapply(depth, rbinom, n = 1, prob = df$p))
        sim_VAF <- sim_NV / depth
        sim_VAF <- sim_VAF[sim_NV > 3]
        dens <- density(sim_VAF)

        out.df <- data.frame(VAF = dens$x,
                             dens_y = dens$y)
        out.df$sample <- sample_name
        out.df$component <- df$component
        out.df$prop <- df$prop
        out.df$p <- df$p
        out.df$density <- out.df$dens_y * out.df$prop
        return(out.df)
    }) %>%
        bind_rows() %>%
        return()
  }) %>% 
  bind_rows() %>% 
  tibble()

# Plot the modeled BMM for each sample
cat("Plotting VAF BMM results...\n")
g5 <- ggplot() +
        geom_histogram(data = vaf.df %>% filter(VAF > 0),
                        aes(x = VAF, y = after_stat(density)),
                        binwidth = 0.05, 
                        fill = "gray90", 
                        color = "black") +
        geom_ribbon(data = sim_vaf.df,
                        aes(x = VAF,
                            y = density,
                            ymin = density,
                            ymax = density,
                            group = as.factor(component),
                            color = ),
                        alpha = 0.5,
                        fill = "transparent",
                        color= "black",
                        linetype = 1,
                        linewidth = 0.5) +
        geom_rug(data = bmm.df %>% filter(major_component == TRUE),
                aes(x = p),
                linewidth = 2,
                outside = FALSE) +
        facet_wrap(~sample, scales = "free") +
        scale_x_continuous(limits = c(0, 1)) +
        theme_cowplot() +
        theme(strip.background = element_blank(),
              strip.text = element_text(size = 6, face = "bold"),
              axis.text = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "bottom") +
        labs(x = "VAF", y = "Density", color = "Major component")

# Merge plots and save
cat("Merging plots and saving...\n")
g <- ((g2 / g5 + plot_layout(heights = c(1, 3))) |
        (g3 / g4 / g1 + plot_layout(heights = c(2, 1, 1)))) +
        plot_layout(ncol = 2, widths = c(2, 1)) +
        plot_annotation(tag_levels = 'A') &
        theme(axis.title = element_text(size = 10))
ggsave(g,file = output_jpg, width = 15, height = 10, dpi = 300)
