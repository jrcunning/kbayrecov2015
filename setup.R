# =================================================================================================
# DATA PREPARATION
# =================================================================================================
# • Load libraries --------------------------------------------------------------------------------
library(lme4); library(MASS); library(reshape2); library(lattice); library(lmerTest); library(lsmeans)
library(scales); library(merTools); library(devtools); library(pBrackets); library(lattice)
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
Mcap.plates <- list.files(path="b_16-22_t_05/", pattern="txt$", full.names=T)
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
# Remove NECs
Mcap <- Mcap[grep("NEC", Mcap$Sample.Name, fixed=T, invert=T), ]

# Parse sample names and dates
Mcap$plate <- unlist(lapply(strsplit(Mcap$File.Name, split="_"), "[", 4))
sample.names <- rbind.fill(lapply(strsplit(as.character(Mcap$Sample.Name), split="_|-"), 
                                  function(X) data.frame(t(X))))
colnames(sample.names) <- c("sample", "date")
Mcap <- cbind(sample.names, Mcap[, -1])
Mcap$date <- as.Date(Mcap$date, format="%m.%e.%y")
Mcap$fdate <- factor(Mcap$date)
Mcap$days <- as.numeric(Mcap$date) - min(as.numeric(Mcap$date))

# Check for date errors
levels(Mcap$fdate)
Mcap[Mcap$date=="2015-12-04",]
# Check for missing dates
Mcap[is.na(Mcap$date), ]

# Replace "sample" column name with "colony"
colnames(Mcap)[which(colnames(Mcap)=="sample")] <- "colony"

# Make new column to indicate if sample failed (host assay did not amplify in both technical replicates)
Mcap$fail <- ifelse(Mcap$Mc.reps < 2, TRUE, FALSE)
table(Mcap$fail)
fails <- Mcap[Mcap$fail==TRUE, ]
#write.table(fails, file="fails.tsv", sep="\t", row.names=F)
# Remove failed samples
Mcap <- Mcap[which(Mcap$fail==FALSE), ]

# Calculate total S/H ratio and D/C ratio and propD
colnames(Mcap)[which(colnames(Mcap) %in% c("C.Mc", "D.Mc"))] <- c("C.SH", "D.SH")  # Rename cols
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

# Assign visual ID and reef location metadata
Mcap$vis <- factor(ifelse(as.numeric(as.character(Mcap$colony)) %% 2 == 0, "not bleached", "bleached"))
Mcap$reef <- cut(as.numeric(as.character(Mcap$colony)), 
                 breaks=c(1,51,101,201,251), labels=c("HIMB", "25", "44", "42"))


# # Identify overall dominant symbiont clade across time points based on mean proportion clade D
meanpropD <- aggregate(Mcap$propD, by=list(colony=Mcap$colony), FUN=mean, na.rm=T)
meanpropD$tdom <- factor(ifelse(meanpropD$x > 0.5, "D", "C"))
rownames(meanpropD) <- meanpropD$colony
Mcap$tdom <- meanpropD[as.character(Mcap$colony), "tdom"]

# Identify samples in which no symbionts were detected
Mcap[which(Mcap$tot.SH==0), ]
Mcap[which(Mcap$tot.SH==0), "tot.SH"] <- 2e-5

# Manually remove inappropriate data points
Mcap <- Mcap[!(Mcap$colony=="77" & Mcap$date=="2015-12-17"),] #Photo reveals variation within colony, this sample probably came from a bleached tip while other parts of the colony appeared to be recovering, not representative of full colony recovery
Mcap <- Mcap[!(Mcap$colony=="130" & Mcap$date=="2015-09-14"),] #Photo reveals different bleached colony very close to actual colony that is assumed to be mistaken for colony 130 on this timepoint. Bleached colony died by next timepoint.
Mcap <- Mcap[!(Mcap$colony=="58" & Mcap$date=="2015-11-04"),]
Mcap <- Mcap[!(Mcap$colony=="31" & Mcap$date=="2015-10-01"),]
Mcap <- Mcap[!(Mcap$colony=="31" & Mcap$date>"2015-12-03"), ] # This colony is dead in pics, likely continued sampling neighboring colony
Mcap <- Mcap[!(Mcap$colony=="123" & Mcap$date=="2016-01-20"),]

# Filter duplicates -----------------------
filter.dups <- function(data) {
  keep <- data.frame()  # Create empty data frame to receive runs to keep
  # Identify duplicated samples (i.e., samples run multiple times)
  dups <- unique(rbind(data[duplicated(interaction(data$colony, data$date)), ], 
                       data[rev(duplicated(rev(interaction(data$colony, data$date)))), ]))
  dups <- droplevels(dups)
  # Create list of data frames for each duplicated sample
  dups.list <- split(dups, f=interaction(dups$colony, dups$date))
  dups.list <- dups.list[!unlist(lapply(dups.list, empty))]
  # Loop through each duplicated sample and select which run to keep
  for (sample in dups.list) {
    # If any of the runs are 12-17 time point, delete the first run
    if ("2015-12-17" %in% sample$date) {  
      sample <- sample[which.max(sample$plate), ]
    }
    # If multiple runs of the original or re-extracted sample still exist, keep the run with
    #   the lowest average standard deviation of technical replicates of each target
    sample <- sample[which.min(rowMeans(sample[, c("Mc.CT.sd", "C.CT.sd", "D.CT.sd")], na.rm=T)), ]
    # Collect runs to keep
    keep <- merge(keep, sample, all=T)
  }
  # Filter out those runs not kept from data
  nondups <- data[setdiff(rownames(data), rownames(dups)), ]
  result <- merge(nondups, keep, all=T)
  return(result)
}

Mcap.f <- filter.dups(Mcap)  # filter duplicate sample runs

# Filter out samples with high Mc.CT values
boxplot(Mcap.f$Mc.CT.mean)
thresh <- boxplot.stats(Mcap.f$Mc.CT.mean)$stats[5]
Mcap.ff <- Mcap.f[which(Mcap.f$Mc.CT.mean <= thresh), ]
#Mcap.ff <- Mcap.ff[which(Mcap.ff$tot.SH < 1), ]
range(Mcap.ff$tot.SH)

# IMPORT YEAR1 DATA/ read in coral condition data and merge with Mcap.ff.all
Mcap2014.f <- read.csv("Mcapf_year1.csv")
Mcap2014.ff <- read.csv("Mcapff_year1.csv")

Mcap2014.ff <- Mcap2014.ff[, c("colony","date","C.SH","D.SH","tot.SH","propD","syms","dom","tdom","vis","reef")]
Mcap2014.ff$colony <- as.factor(as.character(Mcap2014.ff$colony))
Mcap2014.ff$date <- as.Date(Mcap2014.ff$date)
Mcap2015.ff <- Mcap.ff[, c("colony","date","C.SH","D.SH","tot.SH","propD","syms","dom","tdom","vis","reef")]

Mcap.ff.all <- rbind(Mcap2014.ff, Mcap2015.ff)
Mcap.ff.all$syms <- factor(Mcap.ff.all$syms, levels=c("C", "CD", "DC", "D"))
Mcap.ff.all$days <- as.numeric(Mcap.ff.all$date - as.Date("2014-10-24"))
levels(Mcap.ff.all$syms)

condition <- read.csv("coralcondition.csv", header=TRUE)
condition$date <- as.Date(condition$date, format="%m/%e/%y")
condition$date
condition$colony <- as.factor(condition$colony)
condition$score <- as.integer(as.character(condition$score))
Mcap.ff.all <- merge(condition, Mcap.ff.all, all.y=T)

depth <- read.csv("depth.csv", header=TRUE)
depth$colony <- as.factor(depth$colony)
Mcap.ff.all <- merge(depth, Mcap.ff.all, all.y=T)

# FUNCTION TO TO DETERMINE SHUFFLERS USING VARIOUS FITTING TECHNIQUES
plotpropD <- function(col, method=c("locfit", "locfitr", "loess", "loess.s", "gam")) {
  method <- factor(method, levels=c("locfit", "locfitr", "loess", "loess.s", "gam"))
  df <- Mcap.ff.all[Mcap.ff.all$colony==col, ]
  plot(propD ~ days, df, ylim=c(0,1), main=col)
  dayrange <- as.numeric(seq(min(df$days), max(df$days), length.out=100))
  shuffle <- as.data.frame(setNames(replicate(5, NA, simplify=F), levels(method)))
  success <- as.data.frame(setNames(replicate(5, NA, simplify=F), levels(method)))
  rownames(shuffle) <- col
  colors <- brewer.pal(5, "Dark2")
  # Fit Locfit
  if ("locfit" %in% method) {
    tryCatch({
      locf <- locfit(propD ~ lp(days, nn=1), data=df, family="betar", lfproc=locfit.raw)
      locfpred <- predict(locf, newdata=data.frame(days=dayrange))
      lines(locfpred ~ dayrange, lty=1, col=colors[which(levels(method)=="locfit")])
      CtoD <- dayrange[which(diff(sign(locfpred-0.5))>0)]
      DtoC <- dayrange[which(diff(sign(locfpred-0.5))<0)]
      points(c(CtoD, DtoC), rep(0.5, length(c(CtoD, DtoC))), pch="*", cex=2, col=colors[which(levels(method)=="locfit")])
      shuff <- ifelse(length(c(CtoD, DtoC))!=0, "shuff", "noshuff")
      shuffle$locfit <- shuff
      success$locfit <- "fitted"
    },
    error=function(e) {
      success$locfit <<- "failed"
      shuffle$locfit <<- "noshuff"
    },
    warning=function(w) {
      success$locfit <<- "unreliable"
      shuffle$locfit <<- "noshuff"
    })
  } 
  if ("locfitr" %in% method) {
    tryCatch({
      df$propD[df$propD==0] <- 0.0000001
      df$propD[df$propD==1] <- 0.9999999
      locf2 <- locfit(propD ~ lp(days), data=df, family="betar", lfproc=locfit.robust)
      locf2pred <- predict(locf2, newdata=data.frame(days=dayrange))
      lines(locf2pred ~ dayrange, lty=2, col=colors[which(levels(method)=="locfitr")])
      CtoD <- dayrange[which(diff(sign(locf2pred-0.5))>0)]
      DtoC <- dayrange[which(diff(sign(locf2pred-0.5))<0)]
      points(c(CtoD, DtoC), rep(0.5, length(c(CtoD, DtoC))), pch="*", cex=2, col=colors[which(levels(method)=="locfitr")])
      shuff <- ifelse(length(c(CtoD, DtoC))!=0, "shuff", "noshuff")
      shuffle$locfitr <- shuff
      success$locfitr <- "fitted"
    },
    error=function(e) {
      success$locfitr <<- "failed"
      shuffle$locfitr <<- "noshuff"
    },
    warning=function(w) {
      success$locfitr <<- "unreliable"
      shuffle$locfitr <<- "noshuff"
    })
  } 
  if ("loess" %in% method) {
    tryCatch({
      loess <- loess(propD ~ days, data=df, family="symmetric")
      loesspred <- predict(loess, newdata=data.frame(days=dayrange))
      lines(loesspred ~ dayrange, lty=3, col=colors[which(levels(method)=="loess")])
      CtoD <- dayrange[which(diff(sign(loesspred-0.5))>0)]
      DtoC <- dayrange[which(diff(sign(loesspred-0.5))<0)]
      points(c(CtoD, DtoC), rep(0.5, length(c(CtoD, DtoC))), pch="*", cex=2, col=colors[which(levels(method)=="loess")])
      shuff <- ifelse(length(c(CtoD, DtoC))!=0, "shuff", "noshuff")
      shuffle$loess <- shuff
      success$loess <- "fitted"
    }, 
    error=function(e) {
      success$loess <<- "failed"
      shuffle$loess <<- "noshuff"
    },
    warning=function(w) {
      success$loess <<- "unreliable"
      shuffle$loess <<- "noshuff"
    })
  } 
  if ("loess.s" %in% method) {
    tryCatch({
      loesssmooth <- loess.smooth(df$days, df$propD, family="symmetric", evaluation = 100)
      lines(loesssmooth, lty=4, col=colors[which(levels(method)=="loess.s")])
      CtoD <- dayrange[which(diff(sign(loesssmooth$y-0.5))>0)]
      DtoC <- dayrange[which(diff(sign(loesssmooth$y-0.5))<0)]
      points(c(CtoD, DtoC), rep(0.5, length(c(CtoD, DtoC))), pch="*", cex=2, col=colors[which(levels(method)=="loess.s")])
      shuff <- ifelse(length(c(CtoD, DtoC))>0, "shuff", "noshuff")
      shuffle$loess.s <- shuff
      success$loess.s <- "fitted"
    }, 
    error=function(e) {
      success$loess.s <<- "failed"
      shuffle$loess.s <<- "noshuff"
    },
    warning=function(w) {
      success$loess.s <<- "unreliable"
      shuffle$loess.s <<- "noshuff"
    })
  } 
  if ("gam" %in% method) {
    tryCatch({
      gam <- gam(propD ~ s(days), family="gaussian", data=df)
      gampred <- predict(gam, newdata=data.frame(days=dayrange), type="response")
      lines(gampred ~ dayrange, lty=5, col=colors[which(levels(method)=="gam")])
      CtoD <- dayrange[which(diff(sign(gampred-0.5))>0)]
      DtoC <- dayrange[which(diff(sign(gampred-0.5))<0)]
      points(c(CtoD, DtoC), rep(0.5, length(c(CtoD, DtoC))), pch="*", cex=2, col=colors[which(levels(method)=="gam")])
      shuff <- ifelse(length(c(CtoD, DtoC))>0, "shuff", "noshuff")
      shuffle$gam <- shuff
      success$gam <- "fitted"
    }, 
    error=function(e) {
      success$gam <<- "failed"
      shuffle$gam <<- "noshuff"
    },
    warning=function(w) {
      success$gam <<- "unreliable"
      shuffle$gam <<- "noshuff"
    })
  }
  success <- success[which(!is.na(success))]
  legend(x=par("usr")[1], y=1.4, legend=paste(colnames(success), "(", success, ")", sep=""), 
         lty=c(1,2,3,4,5)[method], col=colors[method], xpd=NA, bty="n", cex=0.8)
  return(shuffle)
}