#2026-09-11
#Honggang Lyu
#Rscript 02_Forest.R

library(grid)
library(tidyverse)
library(forestploter)
library(ragg)
library(openxlsx)

#paths
datadir <- "data"
outdir <- "output"
dir.create(outdir, showWarnings = FALSE)

#BH-adjusted q-values were computed on the MR results of 01_MR.R as below
# res_adj <- dat %>%
#   group_by(exposure) %>%
#   mutate(pfdr = p.adjust(pval, method = "fdr")) %>%
#   mutate(pbh = p.adjust(pval, method = "BH")) %>%
#   ungroup()

dat = read.csv(file.path(datadir, "plot_EGF_mr.res.csv")) %>%
  rename(Methods = method, `P-value`=pval,Exposure = exposure , `q-value` =pfdr)

### 去除 NA
dat$Exposure <- ifelse(is.na(dat$Exposure), "  ", dat$Exposure)

### 增加空格
dat$Exposure <- ifelse(is.na(dat$Exposure),
                       dat$Exposures,
                       paste0("  ", dat$Exposure))


# NA to blank
dat$Methods <- ifelse(is.na(dat$Method), "", dat$Methods)
dat$nsnp <- ifelse(is.na(dat$nsnp), "", dat$nsnp)



dat$`P-value` <- ifelse(is.na(dat$`P-value`), " ",
                        format(dat$`P-value`, scientific = TRUE, digits = 2))
dat$`q-value` <- ifelse(is.na(dat$`q-value`), " ",
                        format(dat$`q-value`, scientific = TRUE, digits = 2))
# Create confidence interval column to display
dat$`OR (95% CI)` <- ifelse(is.na(dat$or), "",
                            sprintf("%.2f (%.2f to %.2f)",
                                    dat$or, dat$or_lci95, dat$or_uci95))


#dat$`OR.(95%.CI)` <- ifelse(is.na(dat$`OR.(95%.CI)`), " ", dat$`OR.(95%.CI)`)

dat[,15] <- paste(rep("  ", 12), collapse = " ")
dat[,16] <- paste(rep("  ", 1), collapse = " ")


#dat = dat[,-8]


or = as.numeric(dat$or)
or_lci95 = as.numeric(dat$or_lci95)
or_uci95 = as.numeric(dat$or_uci95)

# Define theme
tm <- forest_theme(base_size = 10,
                   ci_pch = 19,
                   ci_Theight = 0.2,
                   refline_col = "#AA3A3A",
                   arrow_type = "open",
                   ci_col = "#527192",
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # 自定义背景色、前景色。fontface:1常规，2粗体，3斜体，4粗斜体
                   core = list(bg_params = list(fill = c("#FFFFFF","#FFFFFF"), col=NA)) # "#f6f6f6"
)



p1 <- forest(dat[,c(1,16,2,3,6,14,15,12)],
             est = or,
             lower = or_lci95,
             upper = or_uci95,
             ci_column = 7,### plot 列
             ref_line = 1,
             arrow_lab = c("Protect Factor","Risk Factor"),
             xlim = c(0.85,1.35),
             ticks_at = c(0.85,1,1.15,1.35),
             theme = tm)
g = p1
g
#g <- add_border(g, part = "header", row=1,col = c(1,3:5,7), where = "bottom")
#g <- add_border(g, part = "header", row=1,col = c(1,3:5,7), where = "top")
g <- edit_plot(g, row = c(1,3,5,7,9), which = "background",
               gp = gpar(fill = "#f6f6f6"))


g <- edit_plot(g, row = c(1,3,5,7,9),col = 1, gp = gpar(col = "black", fontface = "bold"))

g




agg_tiff(filename = file.path(outdir, "mr_forest1.tiff"),
         width = 3400, height = 2800, units = "px", pointsize = 12,
         compression = "lzw",bg = "white", res = 400)
# Print plot

g
dev.off()
