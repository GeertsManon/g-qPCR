# ================================== kDNA ================================== #


utils::Sweave("vignette.Rnw")
pdflatex("vignette.tex")


## Prep
## --------------------------------

sessionInfo()

suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(rKOMICS))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(adegenet))

for (package in c("openxlsx","rKOMICS","ggplot2",
                  "circlize","ComplexHeatmap","RColorBrewer",
                  "reshape2","adegenet")) {
  cat(paste0(package, ':  ',packageVersion(package), '\n'))
}

th <- theme(axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 6),
            axis.title = element_text(size = 14, face = 'bold'),
            strip.text = element_text(size = 12, face = 'bold.italic'))

saveplot <- function(filetype, filename) {
  filename <- paste0("3_kDNA/1_Minicircles/5_Figures/", filename)
  if (filetype == "pdf") {
    filename <- paste0(filename, ".pdf")
    pdf(filename, width = 7, height = 5)
  } else {
    filename <- paste0(filename, ".png")
    png(filename, res = 3000, width = 7, height = 5, units = 'in')
  }
}


## Isolate info
## --------------------------------

setwd("~/Dropbox/SNPassay/")
samplesheet <- read.xlsx("1_Data/Supplementary Table 1.xlsx", sheet=1)
omicsdata <- read.xlsx("1_Data/Supplementary Table 2.xlsx", sheet=1)[1:38,]
samplesheet$analysis.code == omicsdata$strain
meta <- merge(samplesheet,omicsdata, by.x="analysis.code", by.y = "strain")
# excl <- which(meta$alternate.allele.frequency=="NOK")
# meta <- meta[-excl,]
nrow(meta)
meta$typed.as <- gsub('T.b. gambiense II', 'Tbg2', meta$typed.as)
meta$typed.as <- gsub('T.b. brucei', 'Tbb', meta$typed.as)
meta$typed.as <- gsub('T.b. rhodesiense', 'Tbr', meta$typed.as)
meta$typed.as <- gsub('T.b. gambiense', 'Tbg1', meta$typed.as)
meta$typed.as <- gsub('uncertain', 'unc', meta$typed.as)
table(meta$typed.as)
subspecies <- as.factor(meta$typed.as)

Tbg1 <- which(subspecies=="Tbg1")
Tbg2 <- which(subspecies=="Tbg2")
Tbb <- which(subspecies=="Tbb")
Tbr <- which(subspecies=="Tbr")
unc <- which(subspecies=="unc")
all <- c(Tbg1, Tbg2, Tbb, Tbr, unc)


## A  Inspect MC length
## --------------------------------

MSClength <- msc.length("3_kDNA/1_Minicircles/1_MinicirclesFasta/all.minicircles.fasta", 
                        samples$analysis.code, asubspecies)



saveplot("png", "A_MSClength.beforefiltering")

MSClength$plot + 
  geom_vline(xintercept=800, linetype="dashed", color = "red") + 
  geom_vline(xintercept=1200, linetype="dashed", color = "red") + 
  xlim(0,3000) + th + 
  theme(axis.text.x = element_text(size=10))

dev.off()

length(MSClength$length)
length(which(MSClength$length>800 & MSClength$length<1200))/length(MSClength$length)*100


## B  Filter MC
## --------------------------------

fastas <- list.files("3_kDNA/1_Minicircles/1_MinicirclesFasta", pattern = "minicircles.fasta", full.names = T)
fastas <- fastas[-grep("all", fastas)]
# fastas <- fastas[-excl]
preproc <- preprocess(fastas, 
                      groups = subspecies,
                      circ = T,
                      min = 800, max = 1200, 
                      writeDNA = F)

mean(preproc$N_MC$afterfiltering/preproc$N_MC$beforefiltering*100)

saveplot("png", "B_Preproc")
preproc$plot + th
dev.off()

preproc$N_MC$beforefiltering
preproc$summary
results <- preproc$N_MC


## C  Inspect filtered MC length
## --------------------------------

MSClength.f <- msc.length("3_kDNA/1_Minicircles/2_MinicirclesCircFasta/all.minicircles.circ.fasta",
                          meta$analysis.code, subspecies)
length(MSClength.f$length)

saveplot("png", "C_MSClength.afterfiltering")
MSClength.f$plot + th + 
  theme(axis.text.x = element_text(size=10))
dev.off()

########################## --> concatenate all *.minicircles.circ.fasta files
########################## --> perform VSEARCH


## D  UC info
## --------------------------------

ucs <- list.files("3_kDNA/1_Minicircles/3_UC/", pattern = "uc", full.names = T)
ucs <- ucs[c(2:9,1)]
ucs.info <- msc.uc(ucs)

ucs.info$MSCs
ucs.info$`perfect alignments`

saveplot("png", "D_UCinfo")
ucs.info$plots$`MSCs and perfect aligments` + th
dev.off()


## E  MSC matrix
## --------------------------------

matrices <- msc.matrix(ucs, meta$analysis.code, subspecies)
# matrices <- lapply(matrices, function (x) x[,all])
# samples <- meta$analysis.code[all]
# subspecies <- subspecies[all]

id <- "id98"

saveplot("png", "E_MSCheatmap.id98")
msc.heatmap(matrices[[id]], groups = subspecies, samples = meta$analysis.code)
dev.off()


## F  Richness
## --------------------------------

richness <- msc.richness(matrices, meta$analysis.code, subspecies)

saveplot("png", "F_MSCrichness")
richness$plot + th + ylab("Number of MSCs for \ndifferent percent identities")
dev.off()

results <- cbind(results,richness$table[,id])

rich.tbgI <- apply(richness$table[Tbg1,-c(1,2)], 2, mean)
mean(rich.tbgI)
apply(richness$table[Tbg1,-c(1,2)], 2, median)
rich.tbgII <- apply(richness$table[Tbg2,-c(1,2)], 2, mean)
rich.tbb <- apply(richness$table[Tbb,-c(1,2)], 2, mean)
rich.tbr <- apply(richness$table[Tbr,-c(1,2)], 2, mean)
rich.unc <- apply(richness$table[unc,-c(1,2)], 2, mean)
rich.nTbg1 <- apply(richness$table[-Tbg1,-c(1,2)], 2, mean)

mean(rich.nTbg1/rich.tbgI)

mean(apply(rbind(rich.tbgII/rich.tbgI,
            rich.tbb/rich.tbgI,
            rich.tbr/rich.tbgI),
      2, mean))

mean(rich.tbgII/rich.tbgI)
mean(rich.tbgII/rich.tbb)
mean(rich.tbgII/rich.tbr)

Li1.5 <- grep("LiTat-1-5", richness$table$samples)
rich.Li1.5 <- richness$table[grep("LiTat-1-5", richness$table$samples),-c(1,2)]
rich.all <- apply(richness$table[,-c(1,2)], 2, mean)

mean(as.numeric(rich.tbgI/rich.Li1.5))
mean(as.numeric(rich.all/rich.Li1.5))

mean(richness$table[Tbg1,"id100"] - richness$table[Tbg1,"id70"])
mean(richness$table[Tbg2,"id100"] - richness$table[Tbg2,"id70"])
mean(richness$table[Tbb,"id100"] - richness$table[Tbb,"id70"])
mean(richness$table[Tbr,"id100"] - richness$table[Tbr,"id70"])
mean(richness$table[unc,"id100"] - richness$table[unc,"id70"])


## G  Similarity 
## --------------------------------

sim <- msc.similarity(matrices, meta$analysis.code, subspecies)

matrices.sub <- lapply(matrices, "[", TRUE, -c(unc))
matrices.sub <- lapply(matrices.sub, function (x) x[-which(rowSums(x)==0),])
sim.sub <- msc.similarity(matrices.sub, meta$analysis.code[-unc], subspecies[-unc])

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 7))

saveplot("png","G_MSCsimilarity.4groups")
sim.sub$relfreq.plot +
  xlab("MPI") +
  th + ggtitle("") + 
  theme(axis.title = element_text(face = "bold"), 
        plot.caption = element_text(face = "italic", color = "darkgray"),
        legend.title = element_blank()) +
  scale_fill_manual(values = mycolors) 
dev.off()


subspecies2 <- subspecies
levels(subspecies2)[levels(subspecies2)!="Tbg1"] <- "non-Tbg1"
sim <- msc.similarity(matrices, meta$analysis.code, subspecies2)

saveplot("png","G_MSCsimilarity.2groups")
sim$relfreq.plot +
  xlab("Minimum percent identity") + ylab("Proportion of unique \nand shared MSCs") + 
  th + ggtitle("") + 
  theme(axis.title = element_text(face = "bold"), 
        plot.caption = element_text(face = "italic", color = "darkgray"),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(values = mycolors) 
dev.off()


## H  PCA
## --------------------------------

gl <- new('genlight', as.list(as.data.frame(matrices[[id]])))
gl@loc.names <- rownames(matrices[[id]])
genogl.pca <- glPca(gl)

ggplot(as.data.frame(genogl.pca$scores), aes(x=PC1, y=PC2, fill=meta$typed.as, color=meta$typed.as)) + 
  geom_point(size=10, alpha=0.5) + 
  theme(legend.text = element_text(size = 12)) + th

saveplot("png", "H_MSCpca")
msc.pca(clustmatrix = matrices[[id]], meta$analysis.code, subspecies)$plot
dev.off()

saveplot("png", "H_MSCpca.Tbg1")
msc.pca(clustmatrix = matrices[[id]], meta$analysis.code[Tbg1], as.factor(meta$typed.as)[Tbg1])$plot
dev.off()

msc.pca(clustmatrix = matrices[["id70"]], meta$analysis.code, as.factor(meta$typed.as))
msc.pca(clustmatrix = matrices[["id100"]], meta$analysis.code, as.factor(meta$typed.as))


## I  Tbg-spec MSC
## --------------------------------

lapply(matrices, function(x) msc.subset(x, subset = which(subspecies=="Tbg1"))$sum)

subs.id98 <- msc.subset(matrices[[id]], subset = which(subspecies=="Tbg1"))
subs.id98$sum

saveplot("png", "I_MSCsubset")
msc.heatmap(subs.id98$matrix, meta$analysis.code, subspecies)
dev.off()

clusts <- names(which(subs.id98$freq==length(which(subspecies=="Tbg1"))))
length(clusts); clusts
names(clusts) <- c("MSC1","MSC3","MSC2")

msc.heatmap(subs.id98$matrix[clusts,], meta$analysis.code, subspecies)

seq <- msc.seqs(fastafile = "3_kDNA/1_Minicircles/2_MinicirclesCircFasta/all.minicircles.circ.fasta",
                ucfile = ucs[grep("98", ucs)], 
                clustnumbers = clusts, writeDNA = F)

TgsMSC <- seq$summary$`cluster strain`

uc.id70 <- read.uc(ucs[grep("70", ucs)])
uc.id98 <- read.uc(ucs[grep("98", ucs)])

TgsMSC.id70 <- list()
for (i in 1:length(TgsMSC)) {
  TgsMSC.id70[[i]] <- c(TgsMSC[i],
                        uc.id70$hits[grep(TgsMSC[i], uc.id70$hits$V10),]$V9,
                        uc.id70$hits[grep(TgsMSC[i], uc.id70$hits$V9),]$V10,
                        uc.id98$hits[grep(TgsMSC[i], uc.id98$hits$V10),]$V9,
                        uc.id98$hits[grep(TgsMSC[i], uc.id98$hits$V9),]$V10)
}
names(TgsMSC.id70) <- TgsMSC
TgsMSC.id70


# seq <- msc.seqs(fastafile = "all.minicircles.circ.BGI.fasta",
#                 ucfile = ucs[grep("70", ucs)], 
#                 clustnumbers = "C384", writeDNA = TRUE)
# View(as.data.frame(seq))


## J  Depth
## --------------------------------

depth <- msc.depth(depthstats = list.files("3_kDNA/1_Minicircles/4_DepthStats/", 
                                           pattern="stats.txt", 
                                           full.names = T), 
                   groups = subspecies,
                   HCN = meta$haploid.copy.number)


saveplot("png", "J_MSCdepth")
depth$CN + th
dev.off()

totalMC <- averageCN <- medianCN <- maxCN <- minCN <- vector()
for (i in 1:nrow(meta)) {
  totalMC[i] <- sum(depth$all[depth$all$sample==meta$analysis.code[i],"CN"])
  averageCN[i] <- mean(depth$all[depth$all$sample==meta$analysis.code[i],"CN"])
  medianCN[i] <- median(depth$all[depth$all$sample==meta$analysis.code[i],"CN"])
  maxCN[i] <- max(depth$all[depth$all$sample==meta$analysis.code[i],"CN"])
  minCN[i] <- min(depth$all[depth$all$sample==meta$analysis.code[i],"CN"])
}

results <- cbind(results, totalMC, averageCN, medianCN, maxCN, minCN)
Li1.5 <- grep("LiTat-1-5", results$samples)
results[Li1.5,]
mean(results[, "maxCN"])
mean(results[which(results$groups=="T.b. gambiense"),"averageCN"][-Li1.5])
mean(results[,"averageCN"])
results[which(results$groups=="T.b. gambiense II"),]

write.xlsx(results, "3_kDNA/1_Minicircles/results.rkomics.xlsx", sheetName="results")

depth$all[depth$all$sample=="ABBA","CN"]
mean(depth$all[depth$all$sample=="57AT","CN"])


TgsMSC.id98 <- list()
for (i in 1:length(TgsMSC)) {
  TgsMSC.id98[[i]] <- c(TgsMSC[i],
                        uc.id98$hits[grep(TgsMSC[i], uc.id98$hits$V10),]$V9,
                        uc.id98$hits[grep(TgsMSC[i], uc.id98$hits$V9),]$V10)
}
names(TgsMSC.id98) <- TgsMSC
TgsMSC.id98

CN.clusts <- list()
for (j in 1:length(TgsMSC.id98)) {
  CN.clust <- vector()
  for (i in 1:length(TgsMSC.id98[[j]])) {
    d <- grep(TgsMSC.id98[[j]][[i]], rownames(depth$all))
    if (length(d) == 0) {
      CN.clust[i] <- NA
    } else {
      CN.clust[i] <- depth$all[d,]$CN
    }
  }
  names(CN.clust) <- TgsMSC.id98[[j]]
  CN.clusts[[j]] <- CN.clust
}
names(CN.clusts) <- names(TgsMSC.id98)

do.call(rbind, CN.clusts)
df <- melt(CN.clusts)

df[,3] <- c(names(CN.clusts[[1]]), names(CN.clusts[[2]]), 
            names(CN.clusts[[3]]))

write.xlsx(df, "3_kDNA/1_Minicircles/4_DepthStats/DepthPerMinicircle.xlsx", sheetName="CN")


## K  shared MSC
## --------------------------------

display.brewer.pal(n = 4, name = 'Paired')
rainbow(5)

sharedMSC <- function(matrix) {
  matrix.t <- t(matrix)
  clusts <- apply(matrix.t, 1, function(x) sum(x>=1))
  shared.abs <- list()
  for (k in 1:nrow(matrix.t)) {
    shared <- vector()
    clusts.row.k <- names(which(matrix.t[k,]>=1))
    #clusts.row.k <- colnames(t.id97[k,which(t.id97[k,]>=1)])
    for (i in 1:nrow(matrix.t)) {
      clusts.row.i <- names(which(matrix.t[i,]>=1))
      #clusts.row.i <- colnames(t.id97[i,which(t.id97[i,]>=1)])
      shared[i] <- length(intersect(clusts.row.k,clusts.row.i))
    }
    names(shared) <- rownames(matrix.t)
    shared.abs[[k]] <- shared
    names(shared.abs)[k] <- rownames(matrix.t)[k]
  }
  shared.abs <- do.call(rbind, shared.abs)
  return(shared.abs)
}


heatmap(sharedMSC(matrices[[id]]))
heatmap(sharedMSC(subs.id98$matrix))

ncols <- length(levels(as.factor(meta$typed.as)))+1
col_fun = circlize::colorRamp2(c(0,max(sharedMSC(matrices[[id]]))/2,max(sharedMSC(matrices[[id]]))), c("white", "red", "black"))
captions <- vector()
for (i in 1:length(levels(as.factor(meta$typed.as)))) {
  captions[i] <- paste0(i, ": ", levels(as.factor(meta$typed.as))[i])
}

saveplot("png","K_MSCshared.allMSC")
draw(
  Heatmap(sharedMSC(matrices[[id]]), col=col_fun, show_column_names = T, show_row_names = T,
          column_names_gp = grid::gpar(fontsize = 5),
          row_names_gp = grid::gpar(fontsize = 5),
          show_column_dend = T, show_row_dend = T,
          heatmap_legend_param = list(title = "N shared MSC"),
          top_annotation = HeatmapAnnotation(subspecies = as.numeric(as.factor(meta$typed.as)), 
                                             col = list(subspecies = colorRamp2(c(1,2,3,4,5), 2:6)))) + 
    Heatmap(as.factor(meta$typed.as), name = " ", width = grid::unit(5, "mm"), 
            col=2:6),
  show_annotation_legend = FALSE) 
dev.off()

saveplot("png","K_MSCshared.Tbg1MSC")
draw(
  Heatmap(sharedMSC(subs.id98$matrix), col=col_fun, show_column_names = T, show_row_names = T,
          column_names_gp = grid::gpar(fontsize = 5),
          row_names_gp = grid::gpar(fontsize = 5),
          show_column_dend = T, show_row_dend = T,
          heatmap_legend_param = list(title = "N shared MSC"),
          top_annotation = HeatmapAnnotation(subspecies = as.numeric(as.factor(meta$typed.as)), 
                                             col = list(subspecies = colorRamp2(c(1,2,3,4,5), 2:6)))) + 
    Heatmap(as.factor(meta$typed.as), name = " ", width = grid::unit(5, "mm"), 
            col=2:6),
  show_annotation_legend = FALSE) 
dev.off()


## J High copy minicircles
## --------------------------------

nrow(depth$all)
nrow(depth$all[which(depth$all[,"CN"]>5),])
selection <- rownames(depth$all[which(depth$all[,"CN"]>5),])

dna <- ape::read.dna("3_kDNA/1_Minicircles/1_MinicirclesFasta/all.minicircles.fasta", "fasta")
dna2 <- dna[selection]
ape::write.FASTA(dna2, 
               file = "3_kDNA/1_Minicircles/1_MinicirclesFasta/high.copy.numbers.fasta")

## perform VSEARCH

ucs <- list.files("3_kDNA/1_Minicircles/1_MinicirclesFasta", pattern = "uc", full.names = T)
ucs <- ucs[c(2:7,1)]
ucs.info <- msc.uc(ucs)

matrices <- msc.matrix(ucs, meta$analysis.code, as.factor(meta$typed.as))
heatmaps <- lapply(matrices, function(x) msc.heatmap(x, 
                                                     groups = subspecies, samples = meta$analysis.code ))


lapply(matrices, function(x) msc.subset(x, subset = which(subspecies=="T.b. gambiense"))$sum)

subs.id98 <- msc.subset(matrices[[id]], subset = which(subspecies=="T.b. gambiense"))
subs.id98$sum

msc.heatmap(subs.id98$matrix, meta$analysis.code, subspecies)
max(subs.id98$freq)
clusts <- names(which(subs.id98$freq>=10))
length(clusts)

msc.heatmap(subs.id98$matrix[clusts,], meta$analysis.code, subspecies)

seq <- msc.seqs(fastafile = "3_kDNA/1_Minicircles/1_MinicirclesFasta/high.copy.numbers.fasta",
                ucfile = ucs[grep("98", ucs)], 
                clustnumbers = "C23", writeDNA = T)



TgsMSC <- seq$summary$`cluster strain`

uc.id70 <- read.uc(ucs[grep("70", ucs)])
uc.id98 <- read.uc(ucs[grep("98", ucs)])

TgsMSC.id70 <- list()
for (i in 1:length(TgsMSC)) {
  TgsMSC.id70[[i]] <- c(TgsMSC[i],
                        uc.id70$hits[grep(TgsMSC[i], uc.id70$hits$V10),]$V9,
                        uc.id70$hits[grep(TgsMSC[i], uc.id70$hits$V9),]$V10,
                        uc.id98$hits[grep(TgsMSC[i], uc.id98$hits$V10),]$V9,
                        uc.id98$hits[grep(TgsMSC[i], uc.id98$hits$V9),]$V10)
}
names(TgsMSC.id70) <- TgsMSC
TgsMSC.id70



