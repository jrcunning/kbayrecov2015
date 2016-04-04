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
Mcap.ff <- Mcap.ff[which(Mcap.ff$tot.SH < 1), ]
range(Mcap.ff$tot.SH)

# IMPORT YEAR1 DATA/ read in coral condition data and merge with Mcap.ff.all
Mcap2014.f <- read.csv("Mcapf_year1.csv")
Mcap2014.ff <- read.csv("Mcapff_year1.csv")

Mcap2014.ff <- Mcap2014.ff[, c("colony","date","C.SH","D.SH","tot.SH","propD","syms","dom","vis","reef")]
Mcap2014.ff$colony <- as.factor(as.character(Mcap2014.ff$colony))
Mcap2014.ff$date <- as.Date(Mcap2014.ff$date)
Mcap2015.ff <- Mcap.ff[, c("colony","date","C.SH","D.SH","tot.SH","propD","syms","dom","vis","reef")]

Mcap.ff.all <- rbind(Mcap2014.ff, Mcap2015.ff)
Mcap.ff.all$syms <- factor(Mcap.ff.all$syms, levels=c("C", "CD", "DC", "D"))
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


# Clustering colonies ----------
Mcap2015.ff <- Mcap.ff.all[Mcap.ff.all$date > "2015-05-07", c("colony","date","C.SH","D.SH","tot.SH","propD","syms","dom","vis","reef","score")]
Mcap2015.ff$logtot <- log(Mcap2015.ff$tot.SH)
Mcap2015.ff$asinsqrtpropD <- asin(sqrt(Mcap2015.ff$propD))
Mcap2015.ff <- Mcap2015.ff[!Mcap2015.ff$colony %in% c(127,233,27),]

Mcapm <- melt(Mcap2015.ff, id.vars=c("colony","vis","reef","date"), measure.vars=c("logtot", "asinsqrtpropD", "score"))
Mcapd <- dcast(Mcapm, colony ~ date + variable)
rownames(Mcapd) <- Mcapd$colony
Mcapd <- Mcapd[,-1]
head(Mcapd)

library(fpc)
pamk.best <- pamk(Mcapd, scaling=F, criterion="asw")
pamk.best$nc

Mcapdclust <- pam(Mcapd, k=3, metric="manhattan")

#plot(Mcapdclust)

Mcap2015.ff$cluster <- as.factor(Mcapdclust$clustering[as.character(Mcap2015.ff$colony)])
#par(mfrow=c(1,3), mar=c(2,2,2,1))
for (i in 1:nlevels(Mcap2015.ff$cluster)) {
  df <- Mcap2015.ff[Mcap2015.ff$cluster==i, ]
  df <- df[order(df$date), ]
  plot(NA, xlim=range(df$date), ylim=c(-11,0))
  for (c in levels(df$colony)) {
    lines(logtot ~ date, data=df[df$colony==c, ], type="o",
          pch=21, bg=c("blue","lightblue","pink","red")[syms])
  }
}

# Prepare data for clustering
# CHANGE 5-07 to 5-04 to INCLUDE MAY 2015 TIME POINT
Mcap2015.ff <- Mcap.ff.all[Mcap.ff.all$date > "2015-05-07", c("colony","date","C.SH","D.SH","tot.SH","propD","syms","dom","vis","reef","score", "mortality")]
Mcap2015.ff$logtot <- log(Mcap2015.ff$tot.SH)
Mcap2015.ff$C.SH[Mcap2015.ff$C.SH==0] <- 2e-5
Mcap2015.ff$logC <- log(Mcap2015.ff$C.SH)
Mcap2015.ff$D.SH[Mcap2015.ff$D.SH==0] <- 2e-5
Mcap2015.ff$logD <- log(Mcap2015.ff$D.SH)
Mcap2015.ff$asinsqrtpropD <- asin(sqrt(Mcap2015.ff$propD))
# Exclude certain colonies with not enough data
Mcap2015.ff <- Mcap2015.ff[!Mcap2015.ff$colony %in% c(127,233,27),]
# Exclude reef 42 and dates after 2015-12-17 (focus on bleaching and recovery period)
Mcap2015.ff <- Mcap2015.ff[Mcap2015.ff$reef!=42 & Mcap2015.ff$date < "2015-12-18", ]
# Reformat data for clustering and choose variables to use
Mcapm <- melt(Mcap2015.ff, id.vars=c("colony","vis","reef","date"), measure.vars=c("logtot", "asinsqrtpropD", "score"))
Mcapd <- dcast(Mcapm, colony ~ date + variable)
rownames(Mcapd) <- Mcapd$colony
Mcapd <- Mcapd[,-1]
#test to see which colonies dont overlap
#vegdist(Mcapd, na.rm=T)

# Determine number of clusters, cluster, and assign cluster to each colony
pamk(Mcapd, scaling=T)$nc
Mcapdclust <- pamk(Mcapd, scaling=T, k=5)
plot(Mcapdclust$pamobject)
Mcap2015.ff$cluster <- as.factor(Mcapdclust$pamobject$clustering[as.character(Mcap2015.ff$colony)])

# Plot individual colonies in each cluster
for (i in 1:nlevels(Mcap2015.ff$cluster)) {
  df <- Mcap2015.ff[Mcap2015.ff$cluster==i, ]
  df <- df[order(df$date), ]
  plot(NA, xlim=range(df$date), ylim=c(-11,1))
  for (c in levels(df$colony)) {
    lines(logtot ~ date, data=df[df$colony==c, ], type="o",
          pch=21, bg=c("blue","lightblue","pink","red")[syms])
  }
}

# Calculate average SH ratio and propD for each cluster over time
clustavg <- with(Mcap2015.ff, 
                 aggregate(cbind(logtot, asinsqrtpropD, score), by=list(cluster, date), FUN=mean, na.rm=T))
colnames(clustavg) <- c("cluster", "date", "logtot", "asinsqrtpropD", "score")
# Plot cluster mean SH and propD over time
xyplot(logtot ~ date, groups=~cluster, data=clustavg, type="o")
xyplot(asinsqrtpropD ~ date, groups=~cluster, data=clustavg, type="o")
xyplot(score ~ date, groups=~cluster, data=clustavg, type="o")

# Plot cluster mean SH over time and color points by mean propD
rbPal <- colorRampPalette(c('blue','red'))
clustavg$color <- rbPal(10)[as.numeric(cut(clustavg$asinsqrtpropD, breaks = 10))]
with(clustavg, {
  plot(NA, xlim=range(date), ylim=range(logtot))
  for (c in rev(1:nlevels(cluster))) {
    lines(logtot ~ date, clustavg[cluster==c, ], type="o", pch=c(21,22,23,24,25)[c], 
          cex=1.5, lty=c(1,2,3,4,5)[c], lwd=1.5, bg=color)
  }
})


# GROUP MANUALLY BASED ON SET RULES --------
# Rule 1: Colony "bleached" if 1 or more samples had total log(S/H) less than -5.1
#nlow <- aggregate(log(Mcap.ff$tot.SH), by=list(colony=Mcap.ff$colony), FUN=function(x) table(x < -6)[2])
#Mcap.ff$bleach <- ifelse(Mcap.ff$colony %in% nlow[nlow$x >= 1, "colony"], "bleach", "notbleached")
# Alt. Rule 1: Colony "bleached" if was rated as a 1 ever
Mcap.ff$bleach <- ifelse(Mcap.ff$colony %in% Mcap.ff.all[Mcap.ff.all$score==1, "colony"], "bleach", "notbleached")
# Rule 2: Colony ends with either C or D based on dominant clade on 2016-02-11
#Mcap.ff$endwith <- ifelse(Mcap.ff$colony %in% Mcap.ff[Mcap.ff$date=="2016-02-11" & Mcap.ff$dom=="C", "colony"], "C", "D")
# Alt Rule 2: Colony ends with either C or D based on dominant average of 2016-01-20 and 2016-02-11
#avgendpropD <- with(Mcap.ff[Mcap.ff$date=="2016-01-20" | Mcap.ff$date=="2016-02-11", ], 
#                    aggregate(propD, by=list(colony=colony), FUN=mean))
#Mcap.ff$endwith <- ifelse(Mcap.ff$colony %in% avgendpropD[avgendpropD$x >= 0.5, "colony"], "D", "C")
# Rule 3: If colony did not bleach and ended with D, did it start with C on 2015-08-11?
# Mcap.ff$startwith <- ifelse(Mcap.ff$colony %in% Mcap.ff[Mcap.ff$date=="2015-08-11" & Mcap.ff$dom=="C", "colony"], "C", 
#                             ifelse(Mcap.ff$colony %in% Mcap.ff[Mcap.ff$date=="2015-08-11" & Mcap.ff$dom=="D", "colony"], "D", NA))
# Mcap.ff$startdiff <- ifelse(Mcap.ff$bleach=="notbleached" & Mcap.ff$endwith=="D" & Mcap.ff$startwith=="C", "C", "NA")
# Rule 4: Colony shuffles if has significant (p<0.2) change in propD over time from regression analysis
#xyplot(asin(sqrt(propD)) ~ date | colony, type="r", data=Mcap.ff)
modcoefs <- summary(lm(asin(sqrt(propD)) ~ date * colony, data=Mcap.ff))$coef
modcoefs <- modcoefs[grep(pattern=":", x=rownames(modcoefs)), ]
sigs <- modcoefs[modcoefs[,"Pr(>|t|)"] < 0.4, ]
sigcols <- gsub("[^0-9]", "", rownames(sigs))
Mcap.ff$shuff <- ifelse(Mcap.ff$colony %in% sigcols, "shuff", "noshuff")
# Rule 5: Nobleach & noshuff colonies get separated into C or D
Mcap.ff$nbnsdom <- ifelse(Mcap.ff$bleach=="notbleached" & Mcap.ff$shuff=="noshuff" & Mcap.ff$tdom=="C", "C", "NA")
# Identify and count colonies in each group
#Mcap.ff$group <- interaction(Mcap.ff$bleach, Mcap.ff$endwith, Mcap.ff$startdiff)
Mcap.ff$group <- interaction(Mcap.ff$bleach, Mcap.ff$shuff, Mcap.ff$nbnsdom)
aggregate(Mcap.ff$colony, by=list(Mcap.ff$group), FUN=function(x) unique(as.character(x)))
aggregate(Mcap.ff$colony, by=list(Mcap.ff$group), FUN=function(x) length(unique(as.character(x))))
# Plot individual colonies in each group
#xyplot(log(tot.SH) ~ date | group, groups=~colony, data=Mcap.ff, type="o")
#xyplot(propD ~ date | group, groups=~colony, data=Mcap.ff, type="o")
# Aggregate data for each group on each date (mean and SD)
mcdf <- with(Mcap.ff, {
  data.frame(aggregate(cbind(logtot=log(tot.SH), propD), 
                       by=list(bleach=bleach, nbnsdom=nbnsdom, date=date, shuff=shuff), 
                       FUN=mean, na.rm=T),
             logtotSD=aggregate(log(tot.SH), 
                                by=list(bleach=bleach, nbnsdom=nbnsdom, date=date, shuff=shuff),
                                FUN=sd, na.rm=T)$x,
             logtotSE=aggregate(log(tot.SH), 
                                by=list(bleach=bleach, const=const, date=date, shuff=shuff),
                                FUN=function(x) sd(x, na.rm=T)/sqrt(length(x)))$x)
})
mcdf$group <- interaction(mcdf$bleach, mcdf$shuff, mcdf$nbnsdom)
mcdf <- droplevels(mcdf)
#xyplot(logtot ~ date, groups=~group, data=mcdf, type="o")
#xyplot(propD ~ date, groups=~group, data=mcdf, type="o")
rbPal <- colorRampPalette(c('blue','red'))
mcdf$color <- rbPal(10)[as.numeric(cut(mcdf$propD, breaks = 10))]
par(mfrow=c(2,1), mgp=c(2,0.6,0.1))
with(mcdf, {
  par(mar=c(2,3,1,1.5))
  plot(NA, xlim=range(date), ylim=range(logtot), xlab="", ylab="ln S/H", xaxt="n", bty="n", cex.axis=0.8)
  text(par("usr")[1],par("usr")[4], expression(bold(A.)), xpd=NA, adj=2.8, cex=1.2)
  for (c in c(1,2,3)) {
    df=mcdf[group==levels(group)[c], ]
    arrows(df$date, df$logtot + df$logtotSE, df$date, df$logtot - df$logtotSE, code=3, angle=90, length=0.025, xpd=NA)
    lines(logtot ~ date, df, type="o", pch=c(21,22,23)[c], 
          cex=2, lwd=1, bg=color, lty=c(1,2,3)[c], xpd=NA)
  }
  legend("bottomright", pch=c(21,22,23), lty=c(1,2,3), bty="n", inset=0.075, cex=0.8, xpd=NA,
         legend=c("Group 1 (n=12)", "Group 2 (n=22)", "Group 3 (n=22)"), )
  par(mar=c(3,3,0,1.5))
  plot(NA, xlim=range(date), ylim=range(logtot), xlab="", ylab="ln S/H", xaxt="n", bty="n", cex.axis=0.8)
  axis(side=1, at=as.Date(c("2015-08-01", "2015-09-01", "2015-10-01", "2015-11-01", 
                                       "2015-12-01", "2016-01-01", "2016-02-01")),
       labels=c("Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb"), cex.axis=0.8, line=1)
  text(par("usr")[1],par("usr")[4], expression(bold(B.)), xpd=NA, adj=2.8, cex=1.2)
  for (c in c(4,5)) {
    df=mcdf[group==levels(group)[c], ]
    arrows(df$date, df$logtot + df$logtotSE, df$date, df$logtot - df$logtotSE, code=3, angle=90, length=0.025, xpd=NA)
    lines(logtot ~ date, df, type="o", pch=c(21,22,23,24,25)[c], 
          cex=2, lwd=1, bg=color, lty=c(1,2,3,4,5)[c], xpd=NA)
  }
  legend("bottomright", pch=c(24,25), lty=c(4,5), bty="n", inset=0.075, cex=0.8, xpd=NA,
         legend=c("Group 4 (n=10)", "Group 5 (n=2)"))
  # Plot color bar
  x <- quantile(par("usr")[1:2], probs=seq(0.65, 1, length.out=11))
  y <- rep(quantile(par("usr")[3:4], 1) * -0.6, 2) - c(0.4, 1)
  #rect(x[1], y[1], x[7], y[2], xpd=T)
  for (i in 1:10) {
    rect(x[i], y[1], x[i + 1], y[2], xpd=NA,
         border=NA, col=c(rbPal(10))[i])
  }
  text(x[1]-8,y[1]-0.25,"C", xpd=NA, cex=0.8)
  text(x[10]+14,y[1]-0.25,"D", xpd=NA, cex=0.8)
  text(x[5]+8,y[1]-0.4, "mean prop. D", xpd=NA, pos=1, cex=0.7)
})



# Which colonies "bleached" first year and did "notbleach" second year (according to -6 threshold)
Mcap.ff[Mcap.ff$vis=="bleached" & Mcap.ff$bleach=="notbleached", c("colony", "date", "tot.SH", "propD")]
plotcolonyXL("109")




plotpropD <- function(colony) {
  plot(propD ~ days, data=Mcap.ff[Mcap.ff$colony==colony, ])
  abline(lm(propD ~ days, data=Mcap.ff[Mcap.ff$colony==colony, ]))
}

xyplot(propD ~ days | reef, groups=~colony, type="o",
       data=Mcap.ff[Mcap.ff$colony %in% c(125, 123, 11, 54, 71, 31, 201, 119, 215, 40, 69), ])
head(Mcap.ff)



shufflers <- c(11,31,40,54,71,119,125,223,227)
nonshufflers <- as.numeric(as.character(levels(Mcap2015.ff$colony)[!levels(Mcap2015.ff$colony) %in% shufflers]))
par(mfrow=c(1,2))
for (g in c("shufflers", "nonshufflers")) {
  g <- get(g)
  df <- Mcap2015.ff[as.character(Mcap2015.ff$colony) %in% as.character(g), ]
  df <- df[order(df$date), ]
  plot(NA, xlim=range(df$date), ylim=c(-11,0))
  for (c in levels(df$colony)) {
    lines(logtot ~ date, data=df[df$colony==c, ], type="o",
          pch=21, bg=c("blue","lightblue","pink","red")[syms])
  }
}

Mcap.ff.all$shuff <- ifelse(Mcap.ff.all$colony %in% shufflers, 1, 0)
anova(lm(shuff ~ vis, data=Mcap.ff))
plot(factor(shuff) ~ interaction(factor(vis), factor(reef)), data=Mcap.ff)
minsh <- aggregate(data.frame(min=log(Mcap.ff$tot.SH)), by=list(colony=Mcap.ff$colony), FUN=min)
Cs <- unique(Mcap2014.f[Mcap2014.f$tdom=="C", "colony"])
minsh$shuff <- ifelse(minsh$colony %in% shufflers, 1, 0)
minsh <- minsh[minsh$colony %in% Cs, ]
plot(shuff ~ min, data=minsh)
lr <- glm(shuff~min,data=minsh,family=binomial)
lr
predict(lr,data.frame(x=-6),type="response")


# • Analysis: Symbiodinium community structure --------------------------
# Proportion of samples with C only, D only, and C+D mixtures
symtab <- table(Mcap.f$syms)
symtab
samples <- c(symtab[1], symtab[2] + symtab[3], symtab[4])
prop.table(samples)
# Proportion clade D in samples with C+D mixtures
propD <- Mcap.f$propD[which(Mcap.f$propD > 0 & Mcap.f$propD < 1)]
hist(propD)
range(propD)
# Percent of samples with >10% non-dominant symbiont (between 10% and 90% clade D)
sum(prop.table(hist(propD, plot=F)$counts)[2:9])
# Proportion of colonies with C only, D only, and C+D mixtures, aggregated over time
colonies <- aggregate(Mcap.f$syms, by=list(colony=Mcap.f$colony), FUN=paste, collapse="")
colonies$C[grep("C", colonies$x)] <- "C"
colonies$D[grep("D", colonies$x)] <- "D"
colonies$present <- ifelse(is.na(colonies$C), ifelse(is.na(colonies$D), "none", "D only"), ifelse(is.na(colonies$D), "C only", "C+D"))
prop.table(table(colonies$present))
# Summarize clade composition of samples and colonies
clades <- data.frame(Colonies=matrix(prop.table(table(colonies$present))), Samples=prop.table(samples))
clades
# Proportion of colonies with overall C or D dominance (most abundant over time)
propdom <- prop.table(table(Mcap.f[which(Mcap.f$fdate=="2015-08-11"), "tdom"]))
propdom
with(Mcap.f[which(Mcap.f$fdate=="2015-10-01"), ], table(interaction(tdom, reef)))

# quick plot of S/H over time for all colonies --------
xyplot(log(tot.SH) ~ days | vis + reef, groups= ~ colony, data=Mcap.ff[order(Mcap.ff$days), ],
       type="o", ylim=c(-11,1))

xyplot(log(tot.SH) ~ date | vis + reef, groups= ~ colony, data=Mcap.ff.all[order(Mcap.ff.all$date), ],
       type="o", ylim=c(-11,1))

# Analysis of mean SH over time at each reef ----------
# filter out 12.17 because data are bad
# Mcap <- Mcap[which(Mcap$date!="2015-12-17" & !(Mcap$date=="2015-10-21" & Mcap$Mc.CT.mean>30)), ]
# Mcap <- Mcap[which(Mcap$Mc.CT.mean < 30), ]

happy <- aggregate(log(Mcap.ff$tot.SH), by=list(Mcap.ff$date, Mcap.ff$vis), FUN=mean, na.rm=T)
colnames(happy) <- c("date", "vis", "logSH")
yay <- (subset(happy, vis=="bleached"))
cool <- (subset(happy, vis=="not bleached"))
plot(yay$date, yay$logSH, type="b", xlab="Date", ylab="Value of logSH", main="Average Mcap SH Values", xaxt="n")
lines(cool$date, cool$logSH, type="b")
dates <- as.Date(c("2015-08-11", "2015-09-14", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04", "2016-01-16", "2016-02-11"))
axis(side=1, at=dates, labels=as.character(dates))

rainbow <- aggregate(log(Mcap.ff$tot.SH), by=list(Mcap.ff$date, Mcap.ff$reef, Mcap.ff$vis), FUN=mean, na.rm=T)
colnames(rainbow) <- c("date", "reef", "vis", "logSH")
thing1 <- subset(rainbow, (vis=="bleached" & reef=="HIMB"))
thing2 <- subset(rainbow, (vis=="bleached" & reef=="25"))
thing3 <- subset(rainbow, (vis=="bleached" & reef=="44"))
thing4 <- subset(rainbow, (vis=="bleached" & reef=="42"))
thing5 <- subset(rainbow, (vis=="not bleached" & reef=="HIMB"))
thing6 <- subset(rainbow, (vis=="not bleached" & reef=="25"))
thing7 <- subset(rainbow, (vis=="not bleached" & reef=="44"))
thing8 <- subset(rainbow, (vis=="not bleached" & reef=="42"))

plot(thing1$date, thing1$logSH, type="b", xlab="Date", ylab="Value of logSH", main="Average Mcap SH Values", xaxt="n", col="darkorange", pch=5, ylim=c(-10,0))
lines(thing2$date, thing2$logSH, type="b", col="magenta", pch=0)
lines(thing3$date, thing3$logSH, type="b", col="turquoise", pch=1)
lines(thing4$date, thing4$logSH, type="b", col="slateblue", pch=2)
lines(thing5$date, thing5$logSH, type="b", col="darkorange", pch=5, lty=2)
lines(thing6$date, thing6$logSH, type="b", col="magenta", pch=0, lty=2)
lines(thing7$date, thing7$logSH, type="b", col="turquoise", pch=1, lty=2)
lines(thing8$date, thing8$logSH, type="b", col="slateblue", pch=2, lty=2)
dates <- as.Date(c("2015-08-11", "2015-09-14", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04", "2016-01-16", "2016-02-11"))
axis(side=1, at=dates, labels=as.character(dates))
legend("bottomright", legend=c("Reef HIMB", "Reef 25", "Reef 44", "Reef 42"), lty=c(1,1,1,1), col=c("darkorange","magenta","turquoise","slateblue"),inset=.05)
legend("bottomright", legend=c("Bleached", "Not Bleached"), lty=c(1,2), col=c("black","black"), inset=c(.05,.25))


# Plot abundance trajectory for a single reef ------
reefHIMB <- Mcap.ff[Mcap.ff$reef=="HIMB", ] 
reef25 <- Mcap.ff[Mcap.ff$reef=="25", ]
reef42 <- Mcap.ff[Mcap.ff$reef=="42", ]
reef44 <- Mcap.ff[Mcap.ff$reef=="44", ]

plot(reefHIMB$date, log(reefHIMB$tot.SH), 
     pch=21, type="b", bg=c("pink","purple")[reefHIMB$vis], 
     lines(reefHIMB$colony)
     )

plot(reef25$date, log(reef25$tot.SH), 
     pch=21, type="b", bg=c("pink","purple")[reef25$vis], 
     lines(reef25$colony)
)
plot(reef42$date, log(reef42$tot.SH), 
      pch=21, type="b", bg=c("pink","purple")[reef42$vis], 
      lines(reef42$colony)
)
plot(reef44$date, log(reef44$tot.SH), 
     pch=21, type="b", bg=c("pink","purple")[reef44$vis], 
     lines(reef44$colony)
)

xyplot(log(tot.SH) ~ days | vis + reef, groups= ~ colony, data=Mcap.ff[order(Mcap.ff$days), ],
       type="o", ylim=c(-11,1))

# Plot abundance trajectory for a single colony -----
plotcolony <- function(colony) {
  df <- Mcap.ff[Mcap.ff$colony==colony, ]
  df <- df[order(df$date), ]
  plot(df$date, log(df$tot.SH), type="b", pch=21, cex=2, bg=c("blue","lightblue","pink","red")[df$syms], ylim=c(-9,1), xaxt="n")
  dates <- as.Date(c("2015-08-11", "2015-09-14", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04", "2016-01-16", "2016-02-11"))
  axis(side=1, at=dates, labels=as.character(dates))
  abline(h=-1, lty=2)
}

plotcolonyXL <- function(colony) {
  df <- Mcap.ff.all[Mcap.ff.all$colony==colony, ]
  df <- df[order(df$date), ]
  par(mar = c(5,5,2,5))
  plot(df$date, log(df$tot.SH), type="b", pch=21, cex=2, bg=c("blue","lightblue","pink","red")[df$syms], ylim=c(-11,1), xaxt="n", xlab="Date", ylab="log SH")
  dates <- as.Date(c("2014-10-24", "2014-11-04", "2014-11-24", "2014-12-16", "2015-01-14", "2015-05-06", "2015-08-11", "2015-09-14", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04", "2016-01-16", "2016-02-11"))
  axis(side=1, at=dates, labels=as.character(dates))
  abline(h=-1, lty=2)
  par(new = T)
  plot(df$date, df$score, type="b", pch=16, axes=F, xlab=NA, ylab=NA, ylim=c(1,3))
  axis(side=4, at = c(1,2,3), labels= c(1,2,3))
  mtext(side = 4, line = 3, 'Visual Score')
}

plotcolonyXXL <- function(colony) {
  df <- Mcap.ff.all[Mcap.ff.all$colony==colony, ]
  df <- df[order(df$date), ]
  par(mar = c(5,5,2,5))
  plot(df$date, log(df$tot.SH), type="b", pch=21, cex=2, bg=c("blue","lightblue","pink","red")[df$syms], ylim=c(-11,1), xaxt="n", xlab="Date", ylab="log SH")
  dates <- as.Date(c("2014-10-24", "2014-11-04", "2014-11-24", "2014-12-16", "2015-01-13", "2015-02-10", "2015-03-10", "2015-05-06", "2015-06-05", "2015-07-14", "2015-08-11", "2015-09-14", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04", "2016-01-16", "2016-02-11"))
  axis(side=1, at=dates, labels=as.character(dates))
  abline(h=-1, lty=2)
  par(new = T)
  plot(df$date, df$depth, type="l", pch=16, axes=F, xlab=NA, ylab=NA, ylim=c(1,10))
  axis(side=4, at = c(1,2,3,4,5,6,7,8,9,10), labels= c(1,2,3,4,5,6,7,8,9,10))
  mtext(side = 4, line = 3, 'Depth')
}

plotcolonyXXL("11") #2
plotcolonyXXL("31") #2
plotcolonyXXL("40") #5
plotcolonyXXL("54") #8
plotcolonyXXL("71") #5
plotcolonyXXL("78") #6
plotcolonyXXL("80") #7
plotcolonyXXL("119") #5
plotcolonyXXL("125") #7
plotcolonyXXL("130") #6
plotcolonyXXL("215") #2
plotcolonyXXL("223") #2
plotcolonyXXL("227") #4


#REEF HIMB
plotcolonyXL("3") #C, bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("4") #D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("8") #D>C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("11") #C, severe bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("19") #C, bleaching, no shuffling, NO: SH recovers, score does not (all 1)
plotcolonyXL("20") #D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("25") #C, no bleaching, no shuffling, NO: SH seems to not bleach, score is ones
plotcolonyXL("26") #D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("31") #C, severe bleaching, odd shuffling to D, YES score corresponds to SH *SHUFFLING
plotcolonyXL("32") #D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("40") #C>D, no bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("43") #C, severe bleaching, no shuffling, NO: SH recovers, score does not
plotcolonyXL("44") #C, no bleaching, no shuffling, YES score corresponds to SH

#REEF 25
plotcolonyXL("51") #C, bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("52") #D>C, no bleaching, for 9/14 shuffled to C>D then returned to D>C, YES score corresponds to SH
plotcolonyXL("53") #C, bleaching, no shuffling, YES score corresponds to SH (although 12.04.15 is wierd)
plotcolonyXL("54") #C, no bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("57") #C, severe bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("58") #D>C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("65") #C, baby bleaching, no shuffling, NO: SH does not appear to change, vis shows bleaching
plotcolonyXL("66") #C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("69") #C, severe bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("70") #D>C, no bleaching, baby shuffling to D on 2/11, YES score corresponds to SH *(baby)SHUFFLING
plotcolonyXL("71") #C, bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("72") #C>D, no bleaching, random shuffling, YES score corresponds to SH *SHUFFLING
plotcolonyXL("77") #C, severe bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("78") #D>C, no bleaching, no shufflilng, YES score corresponds to SH
plotcolonyXL("79") #C, bleaching, random shuffling, YES score corresponds to SH *SHUFFLING
plotcolonyXL("80") #D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("83") #C, bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("84") #C, no bleaching, no shuffling, YES score corresponds to SH

#REEF 44
plotcolonyXL("109") #C, bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("110") #C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("112") #D>C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("119") #C>D, NO bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("120") #D>C, no bleaching, no shuffling, YES score correponds to SH
plotcolonyXL("121") #C, NO bleaching, no shuffling, NO score shows bleaching, symbionts don't
plotcolonyXL("122") #C>D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("123") #C, bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("124") #C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("125") #C>D, bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("126") #D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("127") # only one data point, C
plotcolonyXL("130") #D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("131") #C, bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("132") #C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("137") #C, no bleaching, no shuffling, NO score shows bleaching, SH does not
plotcolonyXL("138") #D>C, no bleaching, shuffling to D, YES score corresponds to SH *SHUFFLING

#REEF 42
plotcolonyXL("201") #C, bleaching, shuffling, YES score corresponds to SH *SHUFFLING
plotcolonyXL("202") #C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("207")
plotcolonyXL("211") #C, NO bleaching, no shuffling, NO score shows bleaching, but SH doesn't 
plotcolonyXL("212") #D>C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("215") #C, severe bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("216") #D, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("223") #C, bleaching, shuffling to D>C, YES score corresponds to SH *SHUFFLING
plotcolonyXL("224") #C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("227") #D, bleaching, random shuffling, YES score corresponds to SH *SHUFFLING
plotcolonyXL("228") #C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("233")
plotcolonyXL("234") #D, no bleaching, random shuffling but ended at D, YES score corresponds to SH *SHUFFLING
plotcolonyXL("237") #C, no bleaching, no shuffling, NO score shows bleaching, SH doesn't
plotcolonyXL("238") #D>C, no bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("239") #C, bleaching, no shuffling, YES score corresponds to SH
plotcolonyXL("240") # D, no bleaching, no shuffling, YES score corresponds to SH

#SHUFFLING ---------

# Create matrix for image function
clades <- melt(Mcap.f, id.vars=c("colony", "date", "vis", "reef", "tdom"), measure.vars="syms",
               factorsAsStrings=FALSE)
head(clades)
clades$value <- as.numeric(factor(clades$value))
clades <- unique(clades)
clades <- dcast(clades, vis + colony + reef + tdom ~ date, drop=T)
head(clades)
clades[is.na(clades)] <- -1  # Recode missing values as -1

#clades[which(clades$colony=="129"), 8:10] <- -2  # Recode mortality as -2
clades.m0 <- clades[with(clades, order(clades[, 12], clades[, 5], clades[, 7], 
                                       clades[, 8], clades[, 9], clades[, 10])), ]
clades.m <- as.matrix(clades.m0[,5:12])
rownames(clades.m) <- as.character(clades.m0$colony)
# Plot figure
library(RColorBrewer)
par(mfrow=c(1,1), mar=c(3,5,2,2), bg="white")
image(x=seq(1, ncol(clades.m)), y=seq(1, nrow(clades.m)), z=t(clades.m), 
      xaxt="n", yaxt="n", xlab="", ylab="",
      breaks=c(-2,-1.1,0,1,2,3,4,5),
      col=c("black", "white", rev(brewer.pal(11, "RdYlBu")[c(2,1,3,9,11)])))
# Plot date axis
#axis(side=3, at=seq(1:8), labels=FALSE, cex.axis=0.75, par("tck"=-0.025), xpd=T)
text(0.5:7.5, par("usr")[4], xpd=T, cex=0.6, pos=4,
     labels=levels(Mcap$fdate), srt=45, adj=-0.1)
# # Plot Bleached vs. Not Bleached rectangles
# rect(par("usr")[1] - 2.25, par("usr")[3], par("usr")[1] - 1.25, par("usr")[4], xpd=T)
# text(par("usr")[1] - 1.75, quantile(par("usr")[3:4])[c(2, 4)], labels=c("Not Bleached", "Bleached"),
#      srt=90, xpd=2)
# Plot colony numbers
text(0, 1:nrow(clades.m), labels=rownames(clades.m), xpd=T, cex=0.5)

# get shufflers
Mcap[which(Mcap$colony %in% c(71,54,40,119)), ]

head(clades.m)
# Plot Row Side Colors
reefcols <- c("#bebada", "#8dd3c7", "#d9d9d9")
for (i in 1:nrow(clades.m0)) {
  reef <- clades.m0$reef[i]
  rect(par("usr")[1] - 1.25, par("usr")[3] + 1 * (i - 1), 
       par("usr")[1] - 0.25, par("usr")[3] + 1 * (i - 1) + 1, col=reefcols[as.numeric(reef)],
       xpd=T, border=NA)
}
rect(par("usr")[1] - 1.25, par("usr")[3], par("usr")[1] - 0.25, par("usr")[4], xpd=T)
lines(x=c(par("usr")[1] - 2.25, par("usr")[1] - 0.25), y=rep(quantile(par("usr")[3:4], 0.5), 2), xpd=T)
# Plot tdom side boxes
breaks <- c(0, which(diff(as.numeric(clades.m0$tdom))!=0), length(clades.m0$tdom))
doms <- clades.m0$tdom[breaks]
for (i in 2:length(breaks)) {
  rect(par("usr")[2] + 0.25, par("usr")[3] + breaks[i-1], par("usr")[2] + 0.75, par("usr")[3] + breaks[i], xpd=T)
}
for (i in 1:(length(breaks)-1)) {
  text(par("usr")[2] + 0.5, (breaks[i] + breaks[i+1]) / 2, paste(doms[i], "dominant"), xpd=T, srt=90,
       cex=0.75)
}
# Plot Row Side Color Key
for (i in 1:3) {
  rect(par("usr")[1] - 1.25, quantile(par("usr")[3:4], 0) * -1.05 - ((i - 1) * 1),
       par("usr")[1] - 0.25, quantile(par("usr")[3:4], 0) * -1.05 - ((i - 1) * 1) - 1, xpd=T,
       border=NA, col=reefcols[i])
}
rect(par("usr")[1] - 1.25, quantile(par("usr")[3:4], 0) * -1.05, 
     par("usr")[1] - 0.25, quantile(par("usr")[3:4], 0) * -1.05 - 3, xpd=T)
axis(side=2, xpd=T, pos=par("usr")[1] - 1.25, lwd=0, lwd.ticks=0,
     at=c(-1, -2, -3), labels=c("Rf 25", "Rf 44", "HIMB"), las=2, cex.axis=0.6, mgp=c(0,0.4,0))
# Plot Heatmap Key
x <- quantile(par("usr")[1:2], probs=seq(0, 1, length.out=7))
y <- rep(quantile(par("usr")[3:4], 0) * -1.05, 2) - c(0, 1)
rect(x[1], y[1], x[7], y[2], xpd=T)
for (i in 1:6) {
  rect(x[i], y[1], x[i + 1], y[2], xpd=T,
       border=NA, col=c(brewer.pal(11, "RdYlBu")[c(1,3,9,11)], "white", "black")[i])
}


text(xpd=T, y=y[1] - 0.75, pos=1, cex=0.6,
     x=seq(par("usr")[1], par("usr")[2], length=7)[-7] + 0.5, 
     labels=c("D only", "D > C", "C > D", "C only", "no data", "dead"))
text(xpd=T, y=quantile(par("usr")[3:4], 0) * -1.05 - 2.5, pos=1, cex=0.9,
     x=quantile(par("usr")[1:2], 0.5),
     labels=expression(italic(Symbiodinium)~clades))





# MODEL TRAJECTORIES ------
# Model trajectories of symbiont populations over time using mixed model
#   Build piecewise polynomial model with knot at 82 days (January time point)
#   From October to January, fit a quadratic polynomial (1st element of degree=2)
#   From January to May, fit a linear model (2nd element of degree=1)
#   Function is continuous at time=82 days (smooth=0)
Mcap.ff$days <- as.numeric(Mcap.ff$date - as.Date("2015-08-11"))
#offset <- 0  # optional to center "days" axis at any point
sp <- function(x) gsp(x, knots=c(51,85,128), degree=c(2,3,2,2))
#sp <- function(x) cs(x)
#sp <- function(x) bs(x, knots=c(71,128))
# Build full model with fixed effects of vis, tdom, reef, and time, random effect of colony
Mcapdf <- Mcap.ff[Mcap.ff$reef!="42",]
#mod.all.full <- lmerTest::lmer(log(Mcapdf$tot.SH) ~ sp(Mcapdf$days) * Mcapdf$vis * Mcapdf$reef + (sp(Mcapdf$days) | Mcapdf$colony))
#mod.all.full <- lmerTest::lmer(log(tot.SH) ~ poly(days, 3) * vis * reef + (1 | colony), data=Mcap.ff[Mcap.ff$reef!="42",])
#plot(Effect(c("days", "vis", "reef"), mod.all.full, xlevels=list(days=unique(Mcap.ff$days))), 
#     multiline=T, z.var="reef", ci.style="bars")
# Test significance of fixed effects by backwards selection
#modselect <- step(mod.all.full, lsmeans.calc=F, difflsmeans.calc=F, alpha.fixed=0.05)
#modselect
#summary(mod.all.full)
# Rebuild model omitting non-significant fixed effects
#mod.all <- mod.all.full
# Identify outliers with standardized residuals > 2.5
#out <- abs(residuals(mod.all)) > sd(residuals(mod.all)) * 2.5
#Mcap.ff[out, ]  # outlying data points
# Refit model without outliers
#Mcapdf <- Mcapdf[!out, ]
mod.all <- lmerTest::lmer(log(tot.SH) ~ sp(days) * vis * reef + (1 | colony), data=Mcapdf)
#mod.all <- mod.all.full
# Print and save ANOVA table for model
anovatab <- anova(mod.all)
#write.csv(round(anovatab, digits=3), file="output/Table1.csv")
# pseudo-r2 value-- squared correlation between fitted and observed values
summary(lm(model.response(model.frame(mod.all)) ~ fitted(mod.all)))$r.squared

# Plotting function
plotreefs <- function(mod, n) {
  dat <- get(as.character(summary(mod)$call$data))
  dat <- droplevels(dat)
  levs <- expand.grid(reef=levels(dat$reef), vis=levels(dat$vis), days=as.numeric(levels(as.factor(dat$days))))
  datlevs <- list(interaction(dat$reef, dat$vis, dat$days))
  datsumm <- data.frame(levs,
                        mean=with(dat, aggregate(log(tot.SH), by=datlevs, FUN=mean)$x),
                        sd=with(dat, aggregate(log(tot.SH), by=datlevs, FUN=sd)$x),
                        se=with(dat, aggregate(log(tot.SH), by=datlevs, FUN=function(x) sd(x)/sqrt(length(x)))$x),
                        conf95=with(dat, aggregate(log(tot.SH), by=datlevs, FUN=function(x) sd(x)/sqrt(length(x)) * qt(0.975, length(x)-1))$x))
  datlist <- split(datsumm, f=datsumm$reef)
  datlist <- lapply(datlist, function(x) rev(split(x, f=x$vis)))
  pred <- expand.grid(days=seq_len(max(dat$days)), reef=levels(dat$reef), vis=levels(dat$vis))
  bootfit <- bootMer(mod, FUN=function(x) predict(x, pred, re.form=NA), nsim=n)
  # Extract 90% confidence interval on predicted values
  pred$fit <- predict(mod, pred, re.form=NA)
  pred$lci <- apply(bootfit$t, 2, quantile, 0.05)
  pred$uci <- apply(bootfit$t, 2, quantile, 0.95)
  predlist <- split(pred, f=pred$reef)
  predlist <- lapply(predlist, function(x) rev(split(x, f=x$vis)))
  par(mgp=c(1.75,0.4,0), oma=c(0,0,0,0))
  par(mar=c(0,3,0.3,1))
  layout(mat=matrix(seq_len(nlevels(dat$reef)+1)))
  for (reef in levels(dat$reef)) {
    with(datlist[[reef]], {
      # Create plot frame for each reef
      plot(NA, xlim=range(dat$days), ylim=c(-9,-1), xaxt="n", bty="n", tck=-0.03, ylab="ln S/H")
      title(paste("Reef", reef), line=-0.9, adj=0, outer=F)
      # Plot model fit line and shaded CI for bleached and/or not bleached corals
      with(predlist[[reef]], {
        lapply(predlist[[reef]], function(vis) {
          addpoly(vis$days, vis$lci, vis$uci, col=alpha(reefcols[[reef]], 0.4), xpd=NA)
          lines(vis$days, vis$fit, lty=vislty[[vis$vis[1]]])
        })
      })
      # Plot raw data +/- standard error
      lapply(datlist[[reef]], function(vis) {
        arrows(vis$days, vis$mean + vis$se, vis$days, vis$mean - vis$se, code=3, angle=90, length=0.03, xpd=NA)
        points(vis$days, vis$mean, pch=vispch[[vis$vis[1]]], bg=visbg[[vis$vis[1]]])
      })
    })
    #rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
  }
  axis(side=1, at=as.numeric(as.Date(c("2015-08-01", "2015-09-01", "2015-10-01", "2015-11-01", 
                                       "2015-12-01", "2016-01-01", "2016-02-01")) - as.Date("2015-08-11")),
       labels=c("Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb"))
  return(list(predlist=predlist, datlist=datlist))
}
# Plot
reefcols <- list(`25`="#bebada", `44`="#8dd3c7", HIMB="#d9d9d9", `42`="green")
vislty <- list("bleached"=2, "not bleached"=1)
vispch <- list("bleached"=24, "not bleached"=21)
visbg <- list("bleached"="white", "not bleached"="black")
modelplot <- plotreefs(mod.all, 99)

Mcap42 <- Mcap.ff[Mcap.ff$reef==42,]
sp2 <- function(x) gsp(x, knots=c(85,128), degree=c(2,2,2))
mod.42 <- lmerTest::lmer(log(tot.SH) ~ sp2(days) * vis + (1 | colony), data=Mcap42)
plot42 <- function(mod, n) {
  dat <- get(as.character(summary(mod)$call$data))
  dat <- droplevels(dat)
  levs <- expand.grid(vis=levels(dat$vis), days=as.numeric(levels(as.factor(dat$days))))
  datlevs <- list(interaction(dat$vis, dat$days))
  datsumm <- data.frame(levs,
                        mean=with(dat, aggregate(log(tot.SH), by=datlevs, FUN=mean)$x),
                        sd=with(dat, aggregate(log(tot.SH), by=datlevs, FUN=sd)$x),
                        se=with(dat, aggregate(log(tot.SH), by=datlevs, FUN=function(x) sd(x)/sqrt(length(x)))$x),
                        conf95=with(dat, aggregate(log(tot.SH), by=datlevs, FUN=function(x) sd(x)/sqrt(length(x)) * qt(0.975, length(x)-1))$x))
  datlist <- split(datsumm, f=datsumm$vis)
  pred <- expand.grid(days=seq(min(dat$days), max(dat$days)), vis=levels(dat$vis))
  bootfit <- bootMer(mod, FUN=function(x) predict(x, pred, re.form=NA), nsim=n)
  # Extract 90% confidence interval on predicted values
  pred$fit <- predict(mod, pred, re.form=NA)
  pred$lci <- apply(bootfit$t, 2, quantile, 0.05)
  pred$uci <- apply(bootfit$t, 2, quantile, 0.95)
  predlist <- split(pred, f=pred$vis)
 
      plot(NA, xlim=c(0, max(dat$days)), ylim=c(-9,-1), xaxt="n", bty="n", tck=-0.03, ylab="ln S/H")
      title("Reef 42", line=-0.9, adj=0, outer=F)
      # Plot model fit line and shaded CI for bleached and/or not bleached corals
      with(predlist, {
        lapply(predlist, function(vis) {
          addpoly(vis$days, vis$lci, vis$uci, col=alpha("green", 0.4), xpd=NA)
          lines(vis$days, vis$fit, lty=vislty[[vis$vis[1]]])
        })
      })
      # Plot raw data +/- standard error
      lapply(datlist, function(vis) {
        arrows(vis$days, vis$mean + vis$se, vis$days, vis$mean - vis$se, code=3, angle=90, length=0.03, xpd=NA)
        points(vis$days, vis$mean, pch=vispch[[vis$vis[1]]], bg=visbg[[vis$vis[1]]])
      })
    #rect(xleft = 0, ybottom = -6, xright = 82, ytop = -1, lty = 3, border="black")
  return(list(predlist=predlist, datlist=datlist))
}
plot42(mod.42, 99)


# MORTALITY
plot(Mcap.ff.all$mortality ~ Mcap.ff.all$colony)
hist(Mcap.ff.all$mortality)
