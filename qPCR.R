# ================================== qPCR ================================== #


## Prep
## --------------------------------

setwd("/Users/mgeerts/Dropbox/SNPassay/")

sessionInfo()

suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(gginnards))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrepel))

for (package in c("openxlsx","ggpmisc","gginnards",
                  "ggplot2","ComplexHeatmap","RColorBrewer",
                  "reshape2","pheatmap", "ggpubr",
                  "ggrepel")) {
  cat(paste0(package, ':  ',packageVersion(package), '\n'))
}

th <- theme_minimal() + 
      theme(axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 14, face = 'bold'),
            strip.text = element_text(size = 12, face = 'bold'))

saveplot <- function(filetype, filename) {
  filename <- paste0("5_qPCR/", filename)
  if (filetype == "pdf") {
    filename <- paste0(filename, ".pdf")
    pdf(filename, width = 7, height = 5)
  } else {
    filename <- paste0(filename, ".png")
    png(filename, res = 3000, width = 7, height = 5, units = 'in')
  }
}


## A  Efficiency
## --------------------------------
## Efficiency is calculated only when all three replicates showed amplification
## 40 cycli


efficiency <- read.xlsx("1_Data/Suppl3_qPCR.xlsx", sheet="Efficiency")
efficiency$`DNA(fg)/reaction.` <- as.factor(efficiency$`DNA(fg)/reaction.`)
efficiency.na <- na.omit(efficiency)


## For a singleplex reaction, the efficiency of qPCR is calculated as follows:
##       Efficiency = 10^(-1/slope) - 1
## The slope is derived from a graph of Cycles to Threshold (Ct) values 
## plotted against the Log10 of the template amount. 
## A slope of -3.32 indicates an amplification efficiency of 100%.


eff.formula <- function(slope) {
  return(10^(-1/slope) - 1)
} 


## error bars = standard deviation
## gray lines = linear model


mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 7))
saveplot("png","A_efficiency")
ggplot(data=efficiency.na, aes(x=`DNA(fg)/reaction.`, y=avCq, group=qPCR, colour=qPCR)) + 
  geom_errorbar(aes(ymin=avCq-sdCq, ymax=avCq+sdCq), width=0.1, lty=2) + 
  geom_line(aes(group=qPCR)) +
  geom_smooth(method = lm, se = FALSE, formula = y ~ x, size=0.5, col="gray") + 
  stat_poly_eq(formula = y ~ x, 
               #geom="debug",
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,
               vjust=19,
               rr.digits = 2,
               coef.digits= 3,
               size = 3.5) + 
  geom_point(alpha=0.5, size=3) + 
  scale_colour_manual(values = mycolors) +
  xlab('DNA (fg) per qPCR reaction') + 
  ylab('Average Cq values') +
  facet_wrap(. ~ sample) + th
dev.off()


# max slope 
eff.formula(-3.29)*100 # 101.35

# min slope 
eff.formula(-3.47)*100 # 94.17


## B  Interference
## --------------------------------
## 40 cycli 

interference <- read.xlsx("1_Data/Suppl3_qPCR.xlsx", sheet="Interference")
# interference.na <- na.omit(interference)
interference[is.na(interference)] <- 0

## if cut-off to 35 cycli: 
# interference[,4:7][interference[,4:7]>35] <- 0

saveplot("png","B_interference")
ggplot(data=interference, aes(x=qPCR, y=avCq, colour=qPCR)) + 
  geom_point(size=5, alpha=0.5, show.legend = F) + 
  scale_colour_manual(values = mycolors) + th + 
  geom_label_repel(aes(label=ifelse(avCq>1,sample,'')), size=2) + 
  theme(legend.position = "none")
dev.off()


## C  Analytical sensitivity
## --------------------------------
## Analytical sensitivity is calculated only when all three replicates showed amplification
## Only MBA

anal.sens <- read.xlsx("1_Data/Suppl3_qPCR.xlsx", sheet="Analytical sensitivity", sep.names = " ")
anal.sens$`DNA(fg)/reaction` <- round(anal.sens$`DNA(fg)/reaction`,0)
anal.sens$`DNA(fg)/reaction` <- as.factor(anal.sens$`DNA(fg)/reaction`)
anal.sens$`estimated CN` <- as.factor(anal.sens$`estimated CN`)
anal.sens.MBA <- anal.sens[which(anal.sens$sample=="MBA in heparine blood"),]
anal.sens.MBA.na <- na.omit(anal.sens.MBA)

## if cut-off to 35 cycli: 
# anal.sens[,5:9][anal.sens[,5:9]>35] <- NA

mycolors = brewer.pal(name="Oranges", n = 9)[c(4,6,9)]
saveplot("png", "C_analytical.sensitivity")
ggplot(data=anal.sens.MBA.na[-12,], aes(x=`DNA(fg)/reaction`, y=avCq, colour=`estimated CN`, shape=qPCR)) + 
  geom_errorbar(aes(ymin=avCq-sdCq, ymax=avCq+sdCq), width=0.1, lty=2) + 
  geom_line(aes(group=qPCR)) +
  geom_point(size=5, alpha=0.5, aes(shape=factor(qPCR))) + 
  scale_colour_manual(values = mycolors) + th +
  xlab('DNA (fg) per qPCR reaction') + 
  ylab('Average Cq values') 
dev.off()


## D  Specificity
## --------------------------------


## -- some stats


spec <- read.xlsx("1_Data/Suppl3_qPCR.xlsx", sheet="Specificity", sep.names = " ")
table(spec$subspecies)
table(spec$species)
table(spec$subgenus)
nrow(spec)
length(unique(spec$sample))
table(spec[!duplicated(spec$sample),"subspecies"])

table(spec[which(spec$`extraction method`=="maxwell"),]$species)
table(spec[which(spec$`extraction method`=="phenol chloroform"),]$species)
table(spec$`extraction method`)


## -- all raw Cq values


spec[is.na(spec)] <- 0
spec[which(spec$subspecies!="Tbg1"),c('sample','g-qPCR1-Mini1','g-qPCR2-Mini2','g-qPCR3-Mini3','Tbg-qPCR-TgsGP')]

pcrs <- 'g-qPCR1-Mini1|g-qPCR2-Mini2|g-qPCR3-Mini3|Tbg-qPCR-TgsGP|Tbg-qPCR-18SB|Tbg-qPCR-GPI-PLC'
spec.f <- as.matrix(spec[,grep(pcrs, colnames(spec))])
rownames(spec.f) <- paste(spec$species, spec$sample)


col_fun = circlize::colorRamp2(c(0,15,30,30.1), c("white", "yellow", "red","blue"))

png("5_qPCR/D_specificity.rawCq.png", res = 3000, width = 10, height = 15, units = 'in')
Heatmap(spec.f, 
        col = col_fun,
        show_column_names = T, 
        column_names_rot = 0,
        show_row_names = T,
        column_names_gp = grid::gpar(fontsize = 7), 
        column_names_centered = T,
        row_names_gp = grid::gpar(fontsize = 5),
        show_column_dend = F, 
        show_row_dend = T,
        heatmap_legend_param = list(title = "Cq"),
        row_split = spec$extraction.method,
        row_names_side = "left",
        column_split = c("g-qPCR1","g-qPCR2","g-qPCR3","Tbg-qPCR","Tbg-qPCR","Tbg-qPCR"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          #if(spec.f[i, j] > 0)
          grid.text(sprintf("%.1f", spec.f[i, j]), x, y, gp = gpar(fontsize = 7))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          } 
        ) + 
  Heatmap(as.factor(spec$subspecies), name = " ", width = grid::unit(5, "mm"), 
          col=4:5) + 
  Heatmap(as.factor(spec$subgenus), name = " ", width = grid::unit(5, "mm"), 
          col=2:3)
dev.off()


## -- Normalized Cq values (only on Tbg1 strains)


spec.tbg <- spec[which(spec$species=="Tbg1"),]
pc <- which(spec.tbg$`extraction method`=="phenol chloroform")
maxwell <- which(spec.tbg$`extraction method`=="maxwell")
spec.tbg.m <- melt(spec.tbg[,c(1,2,grep("Δ", colnames(spec.tbg)))])

spec.tbg.f <- data.matrix(spec.tbg[,grep("Δ", colnames(spec.tbg))])
rownames(spec.tbg.f) <- paste(spec.tbg$species, spec.tbg$sample)
apply(spec.tbg.f, 2, median)

saveplot("png","D_specificity.normCq")
Heatmap(spec.tbg.f, 
        show_column_names = T, 
        column_names_rot = 0,
        show_row_names = T,
        column_names_gp = grid::gpar(fontsize = 7), 
        column_names_centered = T,
        row_names_gp = grid::gpar(fontsize = 5),
        show_column_dend = F, 
        show_row_dend = T,
        row_split = spec.tbg$extraction.method,
        heatmap_legend_param = list(title = "Cq"),
        row_names_side = "left",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", spec.tbg.f[i, j]), x, y, gp = gpar(fontsize = 7))
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
        } 
) 
dev.off()


## -- Normalized Cq value distributions considering extraction method


saveplot("png","D_specificity.normCq.distribution")
ggplot(spec.tbg.m, aes(x=variable, y=value, col=`extraction method`, fill=`extraction method`)) + 
  geom_boxplot(alpha=0.5) + 
  ylab('ΔCq values') + xlab('qPCR assay') +
  geom_hline(yintercept = 0, linetype="dashed", col="blue", size=0.9) + 
  scale_y_continuous(breaks=seq(-5,10,1)) + 
  facet_grid(. ~ variable, scales = "free") + 
  scale_fill_manual(values = brewer.pal("Dark2", n=4)) + 
  scale_colour_manual(values = brewer.pal("Dark2", n=4)) +
  th + theme(strip.text.x = element_blank()) 
dev.off()


## -- Compare phenol vs maxwell with 4 strains


dupl <- data.frame(table(spec.tbg$sample))
spec.tbg.dupl <- spec.tbg[which(spec.tbg$sample %in% dupl[dupl$Freq >=2, ]$Var1),]
spec.tbg.dupl <- spec.tbg.dupl[-grep('163AT', spec.tbg.dupl$sample),]
pc <- which(spec.tbg.dupl$`extraction method`=="phenol chloroform")
maxwell <- which(spec.tbg.dupl$`extraction method`=="maxwell")

2^abs(median(spec.tbg.dupl[pc,"ΔTbg-qPCR"]) - median(spec.tbg.dupl[pc,"Δg-qPCR2"]))

2^abs(median(spec.tbg.dupl[maxwell,"ΔTbg-qPCR"]) - median(spec.tbg.dupl[maxwell,"Δg-qPCR2"]))

spec.tbg.dupl.m <- melt(spec.tbg.dupl[,c(1:2, grep(pcrs, colnames(spec.tbg.dupl)))])


ggplot(spec.tbg.dupl.m, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=`extraction method`), alpha=0.8) + 
  geom_point(position=position_dodge(width=0.75), aes(group=`extraction method`)) +
  geom_text(position=position_dodge(width=0.75),
            aes(label=sample, group=`extraction method`),
            hjust = 1.1, size=3) + th

ggarrange(
  ggplot(spec.tbg.dupl[pc,], aes(x=`Δg-qPCR1`, y=spec.tbg.dupl[maxwell,"Δg-qPCR1"])) + 
    geom_point(size=5, alpha=0.5) + ylim(-6,3) +
    geom_text(aes(label=sample, vjust=-1.5)) + 
    xlab("phenol chloroform") + ylab("maxwell") + 
    geom_smooth(method = "lm", se = FALSE),
  ggplot(spec.tbg.dupl[pc,], aes(x=`Δg-qPCR2`, y=spec.tbg.dupl[maxwell,"Δg-qPCR2"])) + 
    geom_point(size=5, alpha=0.5) + ylim(-6,3) +
    geom_text(aes(label=sample, vjust=-1.5)) + 
    xlab("phenol chloroform") + ylab("maxwell") + 
    geom_smooth(method = "lm", se = FALSE),
  ggplot(spec.tbg.dupl[pc,], aes(x=`Δg-qPCR3`, y=spec.tbg.dupl[maxwell,"Δg-qPCR3"])) + 
    geom_point(size=5, alpha=0.5) + ylim(-6,3) +
    geom_text(aes(label=sample, vjust=-1.5)) + 
    xlab("phenol chloroform") + ylab("maxwell") + 
    geom_smooth(method = "lm", se = FALSE),
  ncol=3, nrow=1
)


## -- Compare Normalized Cq values of g-qPCRs with qTgsGP


spec.tbg <- na.omit(spec.tbg)
pc <- which(spec.tbg$`extraction method`=="phenol chloroform")
maxwell <- which(spec.tbg$`extraction method`=="maxwell")

sens <- function(samples, pcr) {
  if (samples == "all") {
    Cq <- median(spec.tbg[,"ΔTbg-qPCR"]) - median(spec.tbg[,pcr])
    Cq2 <- 2^abs(median(spec.tbg[,"ΔTbg-qPCR"]) - median(spec.tbg[,pcr]))
  } else if (samples == "pc") {
    Cq <- median(spec.tbg[pc,"ΔTbg-qPCR"]) - median(spec.tbg[pc,pcr])
    Cq2 <- 2^abs(median(spec.tbg[pc,"ΔTbg-qPCR"]) - median(spec.tbg[pc,pcr]))
  } else if (samples == "maxwell") {
    Cq <- median(spec.tbg[maxwell,"ΔTbg-qPCR"]) - median(spec.tbg[maxwell,pcr])
    Cq2 <- 2^abs(median(spec.tbg[maxwell,"ΔTbg-qPCR"]) - median(spec.tbg[maxwell,pcr]))
  } else {
    print('error')
  }
  df <- cbind(Cq, Cq2)
  colnames(df) <- c("Cq", "sensitivity")
  return(df)
}


# all samples
apply(spec.tbg[,-c(1:5)], 2, median)
nrow(spec.tbg)
rbind(sens("all","Δg-qPCR1"),
      sens("all","Δg-qPCR2"),
      sens("all", "Δg-qPCR3"))

# phenol chloroform extractions
apply(na.omit(spec.tbg[pc,-c(1:5)]), 2, median)
nrow(spec.tbg[pc,])
rbind(sens("pc","Δg-qPCR1"),
      sens("pc","Δg-qPCR2"),
      sens("pc", "Δg-qPCR3"))

# maxwell extractions
apply(spec.tbg[maxwell,-c(1:5)], 2, median)
nrow(spec.tbg[maxwell,])
rbind(sens("maxwell","Δg-qPCR1"),
      sens("maxwell","Δg-qPCR2"),
      sens("maxwell", "Δg-qPCR3"))


## Correlation Cq - CN
## --------------------------------

correlation <- read.xlsx("1_Data/Suppl3_qPCR.xlsx", sheet="correlation")

ggplot(correlation, aes(x=CN.MSC1, y=`Δg-qPCR1`)) + 
  geom_point(size=4,alpha=0.5, show.legend = T) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab('estimated copy number') + ylab('Cq value')
ggplot(correlation, aes(x=CN.MSC2, y=`Δg-qPCR2`)) + 
  geom_point(size=4,alpha=0.5, show.legend = T) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab('estimated copy number') + ylab('Cq value')
ggplot(correlation, aes(x=CN.MSC2, y=`Δg-qPCR3`)) + 
  geom_point(size=4,alpha=0.5, show.legend = T) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab('estimated copy number') + ylab('Cq value')
ggplot(correlation, aes(x=CN.TgsGP, y=`ΔTbg-qPCR`)) + 
  geom_point(size=4,alpha=0.5, show.legend = T) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab('estimated copy number') + ylab('Cq value')












