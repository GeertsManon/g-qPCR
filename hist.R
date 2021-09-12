#!/usr/bin/env Rscript 

###############################################################
#
#
# VIS variant filtration
# ======================
#
#	$ module load R
# if error message LC_CTYPE failed, using "C" do:
#	
#	$ unset LC_CTYPE
#
#	$ ./hist.R --f QUAL.values --o QUAL.hist.pdf > QUAL.summary
#
#
###############################################################


# -------------- Options


suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option ( c("--f"), help = "filter: QUAL, DP, FS, MQ, DP, MQRankSum or ReadPosRanksum"),

  make_option ( c("--o"), help = "output file name: <filter>.hist.pdf")
  
)


parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list,
  description = "\nVIS variant filtration"
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options



# -------------- Run


base.h <- function(values) {
  abline(v=mean(values), col="blue", lty=1, lwd=3)
  abline(v=median(values), col="blue", lty=2, lwd=2)
  lines(density(values), lwd = 1, col=2)
}

pdf(opt$o, width=12, height=7)
par(mfrow=c(2,1))

values <- read.csv(opt$f)[,1]
hist(values, prob=T, main="", breaks = 100)
base.h(values)

if ( opt$f == "QUAL.values") {
  hist(values, prob=T, main="", breaks = 1000, xlim=c(0,500000))
  base.h(values)
  abline(v=500, col="green", lty=1, lwd=2)
}

if ( opt$f == "DP.values") {
  hist(values, prob=T, main="", breaks = 1000, xlim=c(0,25000))
  base.h(values)
  abline(v=5, col="green", lty=1, lwd=2)
}

if ( opt$f == "FS.values") {
  hist(values, prob=T, main="", breaks = 1000, xlim=c(0,200))
  base.h(values)
  abline(v=60, col="green", lty=1, lwd=2)
}

if ( opt$f == "QD.values") {
  abline(v=2, col="green", lty=1, lwd=2)
}

if ( opt$f == "MQ.values") {
  abline(v=40, col="green", lty=1, lwd=2)
}

if ( opt$f == "MQRankSum.values") {
  abline(v=-12.5, col="green", lty=1, lwd=2)
}

if ( opt$f == "ReadPosRanksum.values") {
  abline(v=-8, col="green", lty=1, lwd=2)
}


invisible(dev.off())

# -------------- Output


cat("mean:", mean(values), "\n")
cat("median:", median(values), "\n")
cat("max:", max(values), "\n")
cat("min:", min(values), "\n")
