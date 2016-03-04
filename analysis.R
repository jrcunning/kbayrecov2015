# =================================================================================================
# DATA PREPARATION
# =================================================================================================
# • Load libraries --------------------------------------------------------------------------------
library(lme4); library(MASS); library(reshape2); library(lattice); library(lmerTest); library(lsmeans)
library(pbkrtest); library(scales); library(merTools); library(devtools); library(pBrackets)
## SPIDA package available at http://r-forge.r-project.org/projects/spida/
#system(paste("svn checkout svn://svn.r-forge.r-project.org/svnroot/spida/"))
#devtools::install("spida/pkg")
library(spida)

addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}

# -------------------------------------------------------------------------------------------------
# • Load data -------------------------------------------------------------------------------------
# Use steponeR to import data and calculate S/H ratios
source_url("https://raw.githubusercontent.com/jrcunning/steponeR/master/steponeR.R")
# Get list of plate files to read in
Mcap.plates <- list.files(path=".", pattern="txt$", full.names=T)
Mcap <- steponeR(files=Mcap.plates, delim="\t", target.ratios=c("C.Mc", "D.Mc"),
                 fluor.norm=list(C=2.26827, D=0, Mc=0.84815),
                 copy.number=list(C=33, D=3, Mc=1),
                 ploidy=list(C=1, D=1, Mc=2), 
                 extract=list(C=0.813, D=0.813, Mc=0.982))
Mcap <- Mcap$result

# If C or D only detected in one technical replicate, set its ratio to NA (becomes zero)
Mcap$D.Mc[which(Mcap$D.reps==1)] <- NA
Mcap$C.Mc[which(Mcap$C.reps==1)] <- NA

# Remove positive controls
Mcap <- Mcap[grep("+", Mcap$Sample.Name, fixed=T, invert=T), ]

# Parse sample names and dates
sample.names <- rbind.fill(lapply(strsplit(as.character(Mcap$Sample.Name), split="_|-"), 
                                  function(X) data.frame(t(X))))
colnames(sample.names) <- c("sample", "date")
Mcap <- cbind(sample.names, Mcap[, -1])
Mcap$date <- as.Date(Mcap$date, format="%m.%e.%y")
Mcap$fdate <- factor(Mcap$date)

# Check for date errors
levels(Mcap$fdate)
Mcap[which(Mcap$date=="2011-12-17"), ]
Mcap[which(Mcap$date=="2015-01-20"| Mcap$date=="2015-02-11"),]

# Check Mcap CT distribution
hist(Mcap$Mc.CT.mean, breaks=seq(0,40,1))

bad <- Mcap[which(Mcap$Mc.CT.mean >= 30), ]
nrow(bad)

bad <- with(bad, bad[order(Mc.CT.mean, decreasing=T), ])
bads <- bad[, c("sample", "date", "File.Name", "C.CT)]"

table(bad$date)
table(Mcap$date)

nrow(Mcap)


# Calculate total S/H ratio and D/C ratio and propD
colnames(Mcap)[which(colnames(Mcap) %in% c("C.Mcap", "D.Mcap"))] <- c("C.SH", "D.SH")  # Rename cols
Mcap$C.SH[is.na(Mcap$C.SH)] <- 0
Mcap$D.SH[is.na(Mcap$D.SH)] <- 0
Mcap$tot.SH <- Mcap$C.SH + Mcap$D.SH  # Add C and D to get total SH
Mcap$logDC <- log(Mcap$D.SH / Mcap$C.SH)  # Calculate logDC ratio
Mcap$propD <- Mcap$D.SH / (Mcap$D.SH + Mcap$C.SH)

# Identify symbiont clades present (C=C only, CD=C > D, DC=D > C, D=D only)
Mcap$syms <- factor(ifelse(Mcap$C.SH > Mcap$D.SH, ifelse(Mcap$D.SH!=0, "CD", "C"), 
                           ifelse(Mcap$D.SH > Mcap$C.SH, ifelse(Mcap$C.SH!=0, "DC", "D"), NA)),
                    levels=c("C", "CD", "DC", "D"))
# Identify dominant symbiont clade
Mcap$dom <- factor(substr(as.character(Mcap$syms), 1, 1))
# Set zeros to NA to facilitate log transformation
Mcap$C.SH[which(Mcap$C.SH==0)] <- NA
Mcap$D.SH[which(Mcap$D.SH==0)] <- NA
Mcap$tot.SH[which(Mcap$tot.SH==0)] <- NA

#table(is.na(Mcap$tot.SH))  # 5 SAMPLES HAVE NO DATA

# Assign visual ID and reef location metadata
Mcap$vis <- factor(ifelse(as.numeric(as.character(Mcap$sample)) %% 2 == 0, "not bleached", "bleached"))
Mcap$reef <- factor(ifelse(as.numeric(as.character(Mcap$sample)) <= 50, "HIMB",
                           ifelse(as.numeric(as.character(Mcap$sample)) >= 101, "44", "25")))

# Replace date values, create date factor, POSIX date, and numeric days
Mcap$fdate <- revalue(Mcap$date, c("10.24"="20141024", "11.04"="20141104", "11.24"="20141124",
                                   "12.16"="20141216", "01.14"="20150114", "05.06"="20150506"))
Mcap$date <- as.Date(Mcap$fdate, format="%Y%m%d")
Mcap$time <- as.POSIXct(Mcap$date)
Mcap$days <- as.numeric(Mcap$date) - min(as.numeric(Mcap$date))

# Replace "sample" column name with "colony"
colnames(Mcap)[which(colnames(Mcap)=="sample")] <- "colony"