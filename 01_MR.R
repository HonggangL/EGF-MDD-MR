#2026-09-11
#Honggang Lyu
#Rscript 01_MR.R

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ieugwasr))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(gwasvcf))
suppressMessages(library(RadialMR))
suppressMessages(library(openxlsx))
suppressMessages(library(R.utils))




#paths
datadir <- "data"
outdir <- "output"
dir.create(outdir, showWarnings = FALSE)

exposure <- "EGF"
exposure_file <- file.path(datadir, "mr_MetaBrain.csv")
outcome <- "MDD"
outcome_file <- file.path(datadir, "MDD.txt.gz")

result_dir <- outdir
sigstate <- 0

p_thresh <- 5e-6

#only needed if the LD proxy search block below is enabled
plink_path <- "plink"
ld_ref_path <- "path/to/1000G_EUR/g1000_eur"

cat("=================================================\n")
cat("                TwoSample MR analysis            \n")
cat("MR analysis config:\n")
cat(paste0("[1]exposure: ",exposure,"\n"))
cat(paste0("[2]exposurefile: ",exposure_file,"\n"))
cat(paste0("[3]outcome: ",outcome,"\n"))
cat(paste0("[4]outcome_file: ",outcome_file,"\n"))
cat(paste0("[5]output_dir: ",result_dir,"\n"))
cat(paste0("[6]p_thresh: ",p_thresh,"\n"))
cat(paste0("[7]plink_path: ",plink_path,"\n"))
cat(paste0("[8]ld_ref_path: ",ld_ref_path,"\n"))
cat("=================================================\n")


Result <- data.frame(
  Trait1 = rep(exposure,6),
  Trait2 = rep(outcome,6),
  p_thresh = rep(p_thresh,6)
)

# 1. Data Prepare (Clump)
exp_dat <- read_exposure_data(
  filename = exposure_file,
  clump = FALSE,
  sep = ",",
  eaf_col = "EAF",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  phenotype_col = "Symbol",
  units_col = "BrainRegion",
  #eaf_col = "EAF_1KG",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  gene_col = "Symbol",
  samplesize_col = "Samplesize",
  min_pval = 1e-200, 
  log_pval = FALSE, 
  chr_col = "CHR",
  pos_col = "BP"
)
# Outcome Data
out_dat <- read_outcome_data(
  filename = outcome_file, 
  snps = exp_dat$SNP, 
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  #eaf_col = "EAF_1KG",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N",
  min_pval = 1e-200, 
  log_pval = FALSE, 
  chr_col = "CHR",
  pos_col = "BP"
)

out_dat <- out_dat %>%
  mutate(outcome=!!outcome)

# # Identifying & printing exposure instruments missing from outcome GWAS
# missing_IVs <- exp_dat_clumped$SNP[!(exp_dat_clumped$SNP %in% out_dat$SNP)]
# missing_IVs <- exp_dat_clumped %>% filter(SNP %in% missing_IVs)
# print(paste0("Number of IVs missing from outcome GWAS: ", as.character(nrow(missing_IVs))))
# print("List of IVs missing from outcome GWAS:")
# for (i in 1:nrow(missing_IVs)) {
#   print(paste0(missing_IVs$SNP[i]))
# }
# 
# # Replacing missing instruments from outcome GWAS with proxies
# if(length(missing_IVs) == 0) {
#   print("All exposure IVs found in outcome GWAS.")
# } else {
#   print("Some exposure IVs missing from outcome GWAS.")
#   if (nrow(out_dat) - nrow(missing_IVs) < 20) {
#     print("Number of IVs is less than 20, searching for proxies.")
#     out_full <- fread(outcome_file)
#     for (i in 1:nrow(missing_IVs)) {
#       set_plink(plink_path)
#       proxies <-  get_ld_proxies( missing_IVs$SNP[i], ld_ref_path, searchspace = NULL,tag_kb = 5000,tag_nsnp = 5000,
#                                   tag_r2 = 0.8, threads = 10, out = tempfile() )
#       proxy_present = FALSE
#       if(length(proxies$SNP_B) == 0){
#         print(paste0("No proxy SNP available for ", missing_IVs$SNP[i]))
#       } else {
#         for (j in 1:nrow(proxies)) {
#           proxy_present <- proxies$SNP_B[j]  %in% out_full$SNP
#           if (proxy_present) {
#             proxy_SNP = proxies$SNP_B[j]
#             proxy_SNP_allele_1 = proxies$B1[j]
#             proxy_SNP_allele_2 = proxies$B2[j]
#             #original_SNP_allele_1 = proxies$A1[j]
#             #original_SNP_allele_2 = proxies$A2[j]
#             break
#           }
#         }
#       }
#       
#       if(proxy_present == TRUE) {
#         print(paste0("Proxy SNP found. ",  proxy_SNP," has replaced  ", missing_IVs$SNP[i]))
#         proxy_row <- out_dat[1, ]
#         proxy_row$SNP = missing_IVs$SNP[i]
#         proxy_row$beta.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "BETA"])
#         proxy_row$se.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "SE"])
#         
#         if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_1) {
#           proxy_row$effect_allele.outcome <- proxy_SNP_allele_1
#           proxy_row$other_allele.outcome <- proxy_SNP_allele_2
#         } else if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_2) {
#           proxy_row$effect_allele.outcome <- proxy_SNP_allele_2
#           proxy_row$other_allele.outcome <- proxy_SNP_allele_1
#         } else {
#           print(paste0("Mismatch in alleles for proxy SNP: ", proxy_SNP))
#           next  # 跳过后续的操作
#           #proxy_present == FALSE
#         }
#         
#         
#         proxy_row$pval.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "P"])
#         proxy_row$samplesize.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N"])
#         proxy_row$chr.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs$SNP[i], "chr.exposure"])
#         proxy_row$pos.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs$SNP[i], "pos.exposure"])
#         # if("Frq" %in% colnames(out_full)) proxy_row$eaf.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "Frq"])
#         if("Frq" %in% colnames(out_full)) {
#           proxy_row$eaf.outcome <- as.numeric(out_full[out_full$SNP == proxy_SNP, "Frq"])
#         } else {
#           proxy_row$eaf.outcome <- NA  # 或者你可以选择其他默认值来表示缺失的eaf.outcome
#         }
#         
#         out_dat <- rbind(out_dat, proxy_row)
#       }
#       
#       if(proxy_present == FALSE) {
#         #print("[",i,"]  ", paste0("  No proxy SNP available for ", missing_IVs$SNP[i], " in outcome GWAS."))
#         print(sprintf("[%d] No proxy SNP available for %s in outcome GWAS.", i, missing_IVs$SNP[i]))
#       }
#     }
#   }else {
#     print("Number of  IVs is 20 or more, skipping proxy search.")
#   } 
#   }
cat("Harmonise DATA!\n")

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = out_dat, 
  action = 2
)

#rm(exp_dat)

# Calculate the F-value
get_f<-function(dat,F_value=10){
  log<-is.na(dat$eaf.exposure)
  log<-unique(log)
  if(length(log)==1)
  {if(log==TRUE){
    print("Frq is not found, F-value can't be calculated!")
    return(dat)}
  }
  if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("BETA is not found, F-value can't be calculated!")
    return(dat)}
  if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("SE is not found, F-value can't be calculated!")
    return(dat)}
  if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("SampleSize is not found, F-value can't be calculated!")
    return(dat)}


  if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
    R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
    F<- (dat$samplesize.exposure-2)*R2/(1-R2)
    dat$R2<-R2
    dat$F<-F
    dat<-subset(dat,F>F_value)
    return(dat)
  }
}
dat1 <- get_f(dat)
# 
# cat(paste0("SNPs after f-stastic: ",as.character(nrow(dat1)),"\n"))

# Remove outliers
# cat("Remove outliers!\n")

# dat_radial <- format_radial(dat1$`beta.exposure`,dat1$`se.exposure`,dat1$`beta.outcome`,dat1$`se.outcome`,dat1$SNP)
# radial_ivw_obj <- ivw_radial(dat_radial,0.05,1,0.0001)  # Cochran Q
# radial_egger_obj <- egger_radial(dat_radial,0.05,1)   # Rucker's Q

# radial_ivw_dat <- radial_ivw_obj$data %>%
#   filter(Outliers=="Variant")
# radial_egger_dat <- radial_egger_obj$data %>%
#   filter(Outliers=="Variant")

# cat("SNPs qc is done!\n")
dat_qc <- dat

# dat_qc <- subset(dat1, SNP %in% radial_ivw_dat$SNP & SNP %in% radial_egger_dat$SNP)

# cat(paste0("SNPs after radial qc: ",as.character(nrow(dat_qc)),"\n"))

mr_modified <- function (dat, 
                         parameters = default_parameters(), 
                         method_list = subset(mr_method_list(), use_by_default)$obj) 
{
  library(TwoSampleMR)
  mr_raps_modified <- function (b_exp, b_out, se_exp, se_out,parameters) 
  {
    out <- try(suppressMessages(mr.raps::mr.raps(b_exp, b_out, se_exp, se_out,
                                                 over.dispersion = parameters$over.dispersion, 
                                                 loss.function = parameters$loss.function,
                                                 diagnosis = FALSE)),
               silent = T)
    
    # The estimated overdispersion parameter is very small. Consider using the simple model without overdispersion
    # When encountering such warning, change the over.dispersion as 'FASLE'
    
    if ('try-error' %in% class(out))
    {
      output = list(b = NA, se = NA, pval = NA, nsnp = NA)
    }
    else
    {
      output = list(b = out$beta.hat, se = out$beta.se, 
                    pval = pnorm(-abs(out$beta.hat/out$beta.se)) * 2, nsnp = length(b_exp))
    }
    return(output)
  }
  
  method_list_modified <- stringr::str_replace_all(method_list, "mr_raps","mr_raps_modified")
  
  mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"),function(x1)
  {
    x <- subset(x1, mr_keep)
    
    if (nrow(x) == 0) {
      message("No SNPs available for MR analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
      return(NULL)
    }
    else {
      message("Analysing '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
    }
    res <- lapply(method_list_modified, function(meth)
    {
      get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    }
    )
    
    methl <- mr_method_list()
    mr_tab <- data.frame(outcome = x$outcome[1], exposure = x$exposure[1], 
                         method = methl$name[match(method_list, methl$obj)], 
                         nsnp = sapply(res, function(x) x$nsnp), 
                         b = sapply(res, function(x) x$b), 
                         se = sapply(res, function(x) x$se), 
                         pval = sapply(res, function(x) x$pval))
    
    mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
    
    return(mr_tab)
  }
  )
  return(mr_tab)
}


res_tsmr<-mr_modified(dat_qc, method_list=c("mr_wald_ratio","mr_ivw","mr_raps","mr_egger_regression","mr_weighted_mode","mr_weighted_median"))
or_tsmr <- generate_odds_ratios(res_tsmr)

write.csv(or_tsmr, file.path(outdir, "mr.res.csv"), quote = F, row.names = F)


cat("Result table is built!\n")

result_prefix <- paste("tsmr",as.character(p_thresh),exposure,outcome,sep="_")

outprefix <- file.path(result_dir,result_prefix)

write.csv(or_tsmr,paste0(outprefix,".csv"),quote = F,row.names = F)

cat(paste0("MR analysis between ",exposure," and ",outcome, " is done!!\n"))
cat("=================================================\n")
cat(paste0("MR Ploting between ",exposure," and ",outcome, "\n"))


mr_scatter_plot2 <- function (mr_results, dat) 
{
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), 
                       function(d) {
                         d <- plyr::mutate(d)
                         if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
                           return(blank_plot("Insufficient number of SNPs"))
                         }
                         d <- subset(d, mr_keep)
                         index <- d$beta.exposure < 0
                         d$beta.exposure[index] <- d$beta.exposure[index] * -1
                         d$beta.outcome[index] <- d$beta.outcome[index] * -1
                         mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & 
                                           id.outcome == d$id.outcome[1])
                         mrres$a <- 0
                         if ("MR Egger" %in% mrres$method) {
                           temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
                           mrres$a[mrres$method == "MR Egger"] <- temp$b_i
                         }
                         if ("MR Egger (bootstrap)" %in% mrres$method) {
                           temp <- mr_egger_regression_bootstrap(d$beta.exposure, 
                                                                 d$beta.outcome, d$se.exposure, d$se.outcome, 
                                                                 default_parameters())
                           mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
                         }
                         ggplot2::ggplot(data = d, ggplot2::aes(x = beta.exposure, 
                                                                y = beta.outcome)) + ggplot2::geom_errorbar(ggplot2::aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome), 
                                                                                                            colour = "grey", width = 0) + ggplot2::geom_errorbarh(ggplot2::aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure), 
                                                                                                                                                                  colour = "grey", height = 0) + ggplot2::geom_point() + 
                           ggplot2::geom_abline(data = mrres, ggplot2::aes(intercept = a, slope = b, colour = method), show.legend = TRUE) + 
                           ggplot2::scale_colour_manual(values = c( "#33a02c","#1A2693","#B22580","#D1BC77","#6a3d9a","#D5564D",
                                                                    "#ffff99", "#b15928")) +
                           ggplot2::labs(colour = "MR Test", x = paste("SNP effect on", d$exposure[1]), y = paste("SNP effect on", d$outcome[1])) + 
                           ggplot2::theme(legend.position = "top",  legend.direction = "vertical") + 
                           ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))+
                           ggplot2::theme_bw()+
                           ggplot2::theme(legend.position  = "top")
                       })
  mrres
}




mr_leaveoneout_plot2<-function (leaveoneout_results) 
{
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + 
                                                ggplot2::geom_vline(xintercept = 0, linetype = "dotted") + 
                                                ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo,  xmax = up, linewidth = as.factor(tot), colour = as.factor(tot)),  height = 0) + 
                                                ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour = "grey") + 
                                                ggplot2::scale_colour_manual(values = c("black", "red4")) + 
                                                ggplot2::scale_linewidth_manual(values = c(0.3, 1)) + 
                                                ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 8),  axis.ticks.y = ggplot2::element_line(linewidth = 0),   axis.title.x = ggplot2::element_text(size = 8)) + 
                                                ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n'", 
                                                                                 d$exposure[1], "' on '", d$outcome[1], "'"))+
                                                ggplot2::theme_bw()+
                                                ggplot2::theme(legend.position = "none")
                                            })
  res
}

mr_forest_plot2<-function (singlesnp_results, exponentiate = FALSE) 
{
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), 
                     function(d) {
                       d <- plyr::mutate(d)
                       if (sum(!grepl("All", d$SNP)) < 2) {
                         return(blank_plot("Insufficient number of SNPs"))
                       }
                       levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
                       levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
                       am <- grep("All", d$SNP, value = TRUE)
                       d$up <- d$b + 1.96 * d$se
                       d$lo <- d$b - 1.96 * d$se
                       d$tot <- 0.01
                       d$tot[d$SNP %in% am] <- 1
                       d$SNP <- as.character(d$SNP)
                       nom <- d$SNP[!d$SNP %in% am]
                       nom <- nom[order(d$b)]
                       d <- rbind(d, d[nrow(d), ])
                       d$SNP[nrow(d) - 1] <- ""
                       d$b[nrow(d) - 1] <- NA
                       d$up[nrow(d) - 1] <- NA
                       d$lo[nrow(d) - 1] <- NA
                       d$SNP <- ordered(d$SNP, levels = c(am, "", nom))
                       xint <- 0
                       if (exponentiate) {
                         d$b <- exp(d$b)
                         d$up <- exp(d$up)
                         d$lo <- exp(d$lo)
                         xint <- 1
                       }
                       ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + 
                         ggplot2::geom_vline(xintercept = xint, linetype = "dotted") + 
                         ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = up, linewidth = as.factor(tot), colour = as.factor(tot)), height = 0) + 
                         ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                         ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in%  "")), colour = "grey") + 
                         ggplot2::scale_colour_manual(values = c("black",  "red4")) + 
                         ggplot2::scale_linewidth_manual(values = c(0.4,   1.0)) + 
                         ggplot2::theme(legend.position = "none",   axis.text.y = ggplot2::element_text(size = 8),   axis.ticks.y = ggplot2::element_line(linewidth = 0),     axis.title.x = ggplot2::element_text(size = 8)) + 
                         ggplot2::labs(y = "", x = paste0("MR effect size for\n'",    d$exposure[1], "' on '", d$outcome[1], "'"))+
                         ggplot2::theme_bw()+ggplot2::theme( legend.position =" none")
                     })
  res
}

mr_funnel_plot2 <- function(singlesnp_results)
{
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    am <- grep("All", d$SNP, value=TRUE)
    d$SNP <- gsub("All - ", "", d$SNP)
    am <- gsub("All - ", "", am)
    ggplot2::ggplot(subset(d, ! SNP %in% am), ggplot2::aes(y = 1/se, x=b)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(data=subset(d, SNP %in% am), ggplot2::aes(xintercept=b, colour = SNP)) +
      # ggplot2::scale_colour_brewer(type="qual") +
      ggplot2::scale_colour_manual(values = c("#33a02c","#1A2693","#B22580","#D1BC77","#6a3d9a","#D5564D",
                                              "#ffff99", "#b15928")) +
      ggplot2::labs(y=expression(1/SE[IV]), x=expression(beta[IV]), colour="MR Method") +
      ggplot2::theme(legend.position="top", legend.direction="vertical")+ggplot2::theme_bw()
  })
  res
}



plot_ivw = or_tsmr %>%  filter( nsnp >=2)

### scatter plot
for (i in unique(plot_ivw$exposure)){
  res = or_tsmr %>% filter(or_tsmr$exposure == i)
  dat = dat_qc %>% filter(dat_qc$exposure == i)
  p1 = mr_scatter_plot2(mr_results = res, dat = dat) #### 可更改 散点图样式
  pdf(file = paste0(outprefix,'.scatter.plot.pdf'), width = 8, height = 6,pointsize = 5);  # 导出 PDF 开始
  print( p1  );
  dev.off()
}

#### leave one out plot
for (i in unique(plot_ivw$exposure)){
  dat = dat_qc %>% filter(dat_qc$exposure == i)
  res_loo <- mr_leaveoneout(dat)
  p3 <- mr_leaveoneout_plot2(res_loo)
  pdf(file =paste0(outprefix,'.LOO.plot.pdf'), width = 8, height = 6,pointsize = 5);  # 导出 PDF 开始
  print( p3  );
  dev.off()
}

###   Funnel plot
for (i in unique(plot_ivw$exposure)){
  dat = dat_qc%>% filter(dat_qc$exposure == i)
  res_single <- mr_singlesnp(dat)
  p4 <- mr_funnel_plot2(res_single)
  pdf(file = paste0(outprefix,'.Funnel.plot.pdf'), width = 8, height = 6,pointsize = 5);  # 导出 PDF 开始
  print( p4  );
  dev.off()
}

### Forest plot
for (i in unique(plot_ivw$exposure)){
  dat = dat_qc %>% filter(dat_qc$exposure == i)
  res_single <- mr_singlesnp(dat)
  p2 <- mr_forest_plot2(res_single)
  pdf(file = paste0(outprefix,'.Forest.plot.pdf'), width = 8, height = 6,pointsize = 5);  # 导出 PDF 开始
  print( p2  );
  dev.off()
}

cat(paste0("MR Ploting between ",exposure," and ",outcome, " is done!!\n"))
