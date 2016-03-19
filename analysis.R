# =================================================================================================
# DATA PREPARATION
# =================================================================================================
# • Load libraries --------------------------------------------------------------------------------
library(lme4); library(MASS); library(reshape2); library(lattice); library(lmerTest); library(lsmeans)
library(pbkrtest); library(scales); library(merTools); library(devtools); library(pBrackets); library(lattice)
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

# Parse sample names and dates
sample.names <- rbind.fill(lapply(strsplit(as.character(Mcap$Sample.Name), split="_|-"), 
                                  function(X) data.frame(t(X))))
colnames(sample.names) <- c("sample", "date")
Mcap <- cbind(sample.names, Mcap[, -1])
Mcap$date <- as.Date(Mcap$date, format="%m.%e.%y")
Mcap$fdate <- factor(Mcap$date)
Mcap$days <- as.numeric(Mcap$date) - min(as.numeric(Mcap$date))

# Check for date errors
levels(Mcap$fdate)
Mcap[which(Mcap$date=="2011-12-17"), ]
Mcap[which(Mcap$date=="2015-01-20"| Mcap$date=="2015-02-11"),]

# Check Mcap CT distribution
hist(Mcap$Mc.CT.mean, breaks=seq(0,40,1))
hist(Mcap[which(Mcap$date==as.Date("2015-12-17")), "Mc.CT.mean"], breaks=seq(0,40,1))
hist(Mcap[which(Mcap$date!=as.Date("2015-12-17")), "Mc.CT.mean"], breaks=seq(0,40,1))

goodMcap <- Mcap[which(Mcap$Mc.CT.mean < 30 & Mcap$date!=as.Date("2015-12-17")), ]

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
# Set zeros to NA to facilitate log transformation
Mcap$C.SH[which(Mcap$C.SH==0)] <- NA
Mcap$D.SH[which(Mcap$D.SH==0)] <- NA
Mcap$tot.SH[which(Mcap$tot.SH==0)] <- NA

#table(is.na(Mcap$tot.SH))  # 5 SAMPLES HAVE NO DATA

# Assign visual ID and reef location metadata
Mcap$vis <- factor(ifelse(as.numeric(as.character(Mcap$sample)) %% 2 == 0, "not bleached", "bleached"))
Mcap$reef <- cut(as.numeric(as.character(Mcap$sample)), 
                 breaks=c(1,51,101,201,251), labels=c("HIMB", "25", "44", "42"))

# Replace "sample" column name with "colony"
colnames(Mcap)[which(colnames(Mcap)=="sample")] <- "colony"

# # Identify overall dominant symbiont clade across time points based on mean proportion clade D
meanpropD <- aggregate(Mcap$propD, by=list(colony=Mcap$colony), FUN=mean, na.rm=T)
meanpropD$tdom <- factor(ifelse(meanpropD$x > 0.5, "D", "C"))
rownames(meanpropD) <- meanpropD$colony
Mcap$tdom <- meanpropD[as.character(Mcap$colony), "tdom"]
head(Mcap)

# • Analysis: Symbiodinium community structure --------------------------
# Proportion of samples with C only, D only, and C+D mixtures
symtab <- table(Mcap$syms)
symtab
samples <- c(symtab[1], symtab[2] + symtab[3], symtab[4])
prop.table(samples)
# Proportion clade D in samples with C+D mixtures
propD <- Mcap$propD[which(Mcap$propD > 0 & Mcap$propD < 1)]
hist(propD)
range(propD)
# Percent of samples with >10% non-dominant symbiont (between 10% and 90% clade D)
sum(prop.table(hist(propD, plot=F)$counts)[2:9])
# Proportion of colonies with C only, D only, and C+D mixtures, aggregated over time
colonies <- aggregate(Mcap$syms, by=list(colony=Mcap$colony), FUN=paste, collapse="")
colonies$C[grep("C", colonies$x)] <- "C"
colonies$D[grep("D", colonies$x)] <- "D"
colonies$present <- ifelse(is.na(colonies$C), ifelse(is.na(colonies$D), "none", "D only"), ifelse(is.na(colonies$D), "C only", "C+D"))
prop.table(table(colonies$present))
# Summarize clade composition of samples and colonies
clades <- data.frame(Colonies=matrix(prop.table(table(colonies$present))), Samples=prop.table(samples))
clades
# Proportion of colonies with overall C or D dominance (most abundant over time)
propdom <- prop.table(table(Mcap[which(Mcap$fdate=="2015-08-11"), "tdom"]))
propdom
with(Mcap[which(Mcap$fdate=="2015-10-01"), ], table(interaction(tdom, reef)))

# # filter data
# Mcap <- Mcap[which(Mcap$fdate!="2015-12-17"), ]
# boxplot.stats(Mcap$Mc.CT.mean)$stats
# Mcap <- Mcap[which(Mcap$Mc.CT.mean < 30), ]
# Mcap <- Mcap[which(Mcap$Mc.CT.mean > 20), ]

# quick plot of S/H over time
xyplot(log(tot.SH) ~ days | vis + reef, groups= ~ colony, data=Mcap[order(Mcap$days), ],
       type="o", ylim=c(-10,0))

Mcap[which(Mcap$days==71 & Mcap$tdom=="D"), ]


# what is average symbiont to host cell ratio (log-transformed) in bleached and non-bleached colonies on august 11?
aug11 <- subset(Mcap, date=="2015-08-11")
aug11b <- subset(aug11, vis=="bleached")
mean(log(aug11b$tot.SH))
exp(mean(log(aug11b$tot.SH)))  # 0.0444

aug11nb <- subset(aug11, vis=="not bleached")
mean(log(aug11nb$tot.SH))
exp(mean(log(aug11nb$tot.SH))) # 0.06499

# what is average symbiont to host cell ratio (log-transformed) in bleached and non-bleached colonies on september 14th?
sept14 <- subset(Mcap, date=="2015-09-14")
sept14b <- subset(sept14, (vis=="bleached") & (tot.SH!="NA"))
mean(log(sept14b$tot.SH))
exp(mean(log(sept14b$tot.SH))) # 0.01447

sept14nb <- subset(sept14, vis=="not bleached")
mean(log(sept14nb$tot.SH))
exp(mean(log(sept14nb$tot.SH))) # 0.09232

# what is average symbiont to host cell ratio (log-transformed) in bleached and non-bleached colonies on october 1st?
oct01 <- subset(Mcap, date=="2015-10-01")
oct01b <- subset(oct01, (vis=="bleached") & (tot.SH!="NA"))
mean(log(oct01b$tot.SH))
exp(mean(log(oct01b$tot.SH))) # 0.010331

oct01nb <- subset(oct01, (vis=="not bleached") & (tot.SH!="NA"))
mean(log(oct01nb$tot.SH))
exp(mean(log(oct01nb$tot.SH))) # 0.07662

# what is average symbiont to host cell ratio (log-transformed) in bleached and non-bleached colonies on october 21st?
oct21 <- subset(Mcap, date=="2015-10-21")
oct21b <- subset(oct21, (vis=="bleached") & (tot.SH!="NA"))
mean(log(oct21b$tot.SH))
exp(mean(log(oct21b$tot.SH))) #0.03430

oct21nb <- subset(oct21, (vis=="not bleached") & (tot.SH!="NA"))
mean(log(oct21nb$tot.SH))
exp(mean(log(oct21nb$tot.SH))) #0.2295

# what is average symbiont to host cell ratio (log-transformed) in bleached and non-bleached colonies on november 4th?
nov04 <- subset(Mcap, date=="2015-11-04")
nov04b <- subset(nov04, vis=="bleached" & tot.SH!="NA")
mean(log(nov04b$tot.SH))
exp(mean(log(nov04b$tot.SH))) #0.03812

nov04nb <- subset(nov04, vis=="not bleached" & tot.SH!="NA")
mean(log(nov04nb$tot.SH))
exp(mean(log(nov04nb$tot.SH))) #0.1835

# what is average symbiont to host cell ratio (log-transformed) in bleached and non-bleached colonies on december 4th?
dec04 <- subset(Mcap, date=="2015-12-04")
dec04b <- subset(dec04, (vis=="bleached") & (tot.SH!="NA"))
mean(log(dec04b$tot.SH))
exp(mean(log(dec04b$tot.SH))) #0.05548

dec04nb <- subset(dec04, (vis=="not bleached") & (tot.SH!="NA"))
mean(log(dec04nb$tot.SH))
exp(mean(log(dec04nb$tot.SH))) #0.06052

# what is average symbiont to host cell ratio (log-transformed) in bleached and non-bleached colonies on january 20th?
jan20 <- subset(Mcap, date=="2016-01-20")
jan20b <- subset(jan20, (vis=="bleached") & (tot.SH!="NA"))
mean(log(jan20b$tot.SH))
exp(mean(log(jan20b$tot.SH))) #0.060698

jan20nb <- subset(jan20, (vis=="not bleached") & (tot.SH!="NA"))
mean(log(jan20nb$tot.SH))
exp(mean(log(jan20nb$tot.SH))) #0.09441

# what is average symbiont to host cell ratio (log-transformed) in bleached and non-bleached colonies on february 11th?
feb11 <- subset(Mcap, date=="2016-02-11")
feb11b <- subset(feb11, (vis=="bleached") & (tot.SH!="NA"))
mean(log(feb11b$tot.SH))
exp(mean(log(feb11b$tot.SH))) #0.1199

feb11nb <- subset(feb11, (vis=="not bleached") & (tot.SH!="NA"))
mean(log(feb11nb$tot.SH))
exp(mean(log(feb11nb$tot.SH))) #0.1110

#---------- aggregate party time
# filter out 12.17 because data are bad
Mcap <- Mcap[which(Mcap$date!="2015-12-17" & !(Mcap$date=="2015-10-21" & Mcap$Mc.CT.mean>30)), ]
Mcap <- Mcap[which(Mcap$Mc.CT.mean < 30), ]

happy <- aggregate(log(Mcap$tot.SH), by=list(Mcap$date, Mcap$vis), FUN=mean, na.rm=T)
colnames(happy) <- c("date", "vis", "logSH")
yay <- (subset(happy, vis=="bleached"))
cool <- (subset(happy, vis=="not bleached"))
plot(yay$date, yay$logSH, type="b", xlab="Date", ylab="Value of logSH", main="Average Mcap SH Values", xaxt="n")
lines(cool$date, cool$logSH, type="b")
dates <- as.Date(c("2015-08-11", "2015-09-14", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04", "2016-01-16", "2016-02-11"))
axis(side=1, at=dates, labels=as.character(dates))

rainbow <- aggregate(log(Mcap$tot.SH), by=list(Mcap$date, Mcap$reef, Mcap$vis), FUN=mean, na.rm=T)
colnames(rainbow) <- c("date", "reef", "vis", "logSH")
thing1 <- subset(rainbow, (vis=="bleached" & reef=="HIMB"))
thing2 <- subset(rainbow, (vis=="bleached" & reef=="25"))
thing3 <- subset(rainbow, (vis=="bleached" & reef=="44"))
thing4 <- subset(rainbow, (vis=="bleached" & reef=="42"))
thing5 <- subset(rainbow, (vis=="not bleached" & reef=="HIMB"))
thing6 <- subset(rainbow, (vis=="not bleached" & reef=="25"))
thing7 <- subset(rainbow, (vis=="not bleached" & reef=="44"))
thing8 <- subset(rainbow, (vis=="not bleached" & reef=="42"))

plot(thing1$date, thing1$logSH, type="b", xlab="Date", ylab="Value of logSH", main="Average Mcap SH Values", xaxt="n", col="darkorange", pch=5, ylim=c(-8,-1))
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


# Plot abundance trajectory for a single colony
plotcolony <- function(colony) {
  df <- Mcap[Mcap$colony==colony, ]
  df <- df[order(df$date), ]
  plot(df$date, log(df$tot.SH), type="b", pch=21, cex=2, bg=c("blue","lightblue","pink","red")[df$syms], ylim=c(-9,1), xaxt="n")
  dates <- as.Date(c("2015-08-11", "2015-09-14", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04", "2016-01-16", "2016-02-11"))
  axis(side=1, at=dates, labels=as.character(dates))
  abline(h=-1, lty=2)
}

plotcolony(colony="132") # C, no bleaching, no shuffling
plotcolony("77") # C, bleaching, no shuffling
plotcolony("40") # C>D, no bleaching, shuffling to D>C
plotcolony("8") # D, no bleaching, no shuffling
plotcolony("71") # C, bleaching, shuffling to D>C

plotcolony("3") # C, bleaching, no shuffling, 2 to majority 1 to 2
plotcolony("4") # D, no bleaching, no shuffling, all 3
plotcolony("8") # D>C, no bleaching, no shuffling, 2 to majority 3
plotcolony("11") # C, bleaching, shuffling to D>C, 2 to majority 1 to 2
plotcolony("19") # C, bleaching, no shuffling, 2 to majority 1
plotcolony("20") # D, no bleaching, no shuffling, all 3
plotcolony("25") # C, no bleaching, no shuffling, 2 to majority 1, *noted as bleached in pictures
plotcolony("26") # D, no bleaching, no shuffling, all 3
plotcolony("31") # C, bleaching, shuffling to D, 1 until 12/17 is 2, *not a lot of data
plotcolony("32") # D, no bleaching, no shuffling, all 3
plotcolony("40") # C>D, no bleaching, shuffling to D>C, 2 to majority 3 *no bleaching but shuffling
plotcolony("43") # C, bleaching, no shuffling, 2 to majority 1
plotcolony("44") # C, no bleaching, no shuffling, all 3

plotcolony("51") # C, bleaching, no shuffling, 2 to majority 1
plotcolony("52") # D>C, no bleaching, for 9/14 shuffled to C>D then returned to D>C, all 3
plotcolony("53") # C, baby bleaching, 2 to majority 1 to 2
plotcolony("54") # C, no bleaching, shuffling to D>C, all 3 *no bleaching but shuffling
plotcolony("57") # C, bleaching, no shuffling, 3 to majority 1 to 2
plotcolony("58") # D>C, no bleaching, no shuffling, half 2 to half 3 *unrealistic data point
plotcolony("65") # C, baby bleaching, no shuffling, 2 to 1 to majority 2 to 3
plotcolony("66") # C, no bleaching, no shuffling, all 3
plotcolony("69") # C, bleaching, shuffling to D>C, 2 to 1 to 2
plotcolony("70") # D>C, no bleaching, baby shuffling to D on 2/11, all 3
plotcolony("71") # C, bleaching, shuffling to D>C, majority 1 to 2
plotcolony("72") # C>D, no bleaching, random shuffling interspersed to C, all 3
plotcolony("77") # C, bleaching, no shuffling, 2 to majority 1 to 2
plotcolony("78") # D>C, no bleaching, no shufflilng, all 3
plotcolony("79") # C, bleaching, interspersed shuffling to C>D ending in all C, 2 to 1 to 2
plotcolony("80") # D, no bleaching, no shuffling, all 3
plotcolony("83") # C, bleaching, no shuffling, 2 to majority 1 to 2
plotcolony("84") # C, no bleaching, no shuffling, all 3

plotcolony("109") # C, bleaching, no shuffling, 2 to 1 to 2
plotcolony("110") # C, no bleaching (tiny dip maybe?), no shuffling, all 3
plotcolony("112") # D>C, no bleaching, no shuffling, all 3
plotcolony("119") # C>D, bleaching, possible shuffling to D>C, majority 2 to 3, *odd data
plotcolony("120") # D>C, no bleaching, no shuffling, all 3
plotcolony("121") # C, baby bleaching, no shuffling, majority 1 to 2 to 3
plotcolony("122") # C>D, no bleaching, no shuffling, 3 to one 2 to 3
plotcolony("123") # C, bleaching, no shuffling, majority 1 to 3 *odd data point
plotcolony("124") # C, no bleaching, no shuffling, all 3
plotcolony("125") # C>D, bleaching, shuffling to D>C, 2 to 1 to 2
plotcolony("126") # D, no bleaching but very up and down, no shuffling, all 3
plotcolony("127") # only one data point, C
plotcolony("130") # D, no bleaching, no shuffling, 2 to 3 *unrealistic data point
plotcolony("131") # C, bleaching, no shuffling, 2 to 1 to 2
plotcolony("132") # C, no bleaching, no shuffling, all 3
plotcolony("137") # C, bleaching, no shuffling, majority 1 to 2
plotcolony("138") # D>C, no bleaching, shuffling to D, all 3 *shuffling without bleaching

plotcolony("201") # C, bleaching, shuffling back and forth to D>C, 1 to 2 *odd data
plotcolony("202") # C, no bleaching, no shuffling, all 3
plotcolony("211") # C, bleaching, no shuffling, 1 to 2
plotcolony("212") # D>C, no bleaching, maybe shuffling on 1/20 to C>D, all 3
plotcolony("215") # NO DATA
plotcolony("216") # D, no bleaching, no shuffling, all 3
plotcolony("223") # C, bleaching, shuffling to D>C, 1 to 2
plotcolony("224") # C, no bleaching, no shuffling, all 3
plotcolony("227") # D, severe bleaching, various inconclusive shuffling, 1 to 2
plotcolony("228") # C, no bleaching, no shuffling, all 3
plotcolony("233")
plotcolony("234") # D, no bleaching, inconclusive shuffling but ended at D, all 3
plotcolony("237") # C, no shuffling, 1 to 2, *odd data
plotcolony("238") # D>C, no bleaching, no shuffling, all 3
plotcolony("239") # C, bleaching, no shuffling, 1 to 2
plotcolony("240") # D, no bleaching, no shuffling, all 3

# identify duplicates
howmany <- table(interaction(Mcap$colony, Mcap$date))
howmany

dups <- howmany[which(howmany==2)]
dups
# WITH 12/17: IDENTIFY BEST BY LATER PLATE NUMBER




#SHUFFLING
# Create matrix for image function
Mcap.f <- Mcap
clades <- melt(Mcap.f, id.vars=c("colony", "date", "vis", "reef", "tdom"), measure.vars="syms",
               factorsAsStrings=FALSE)
head(clades)
clades$value <- as.numeric(factor(clades$value))
clades <- unique(clades)
clades <- dcast(clades, vis + colony + reef + tdom ~ date, drop=T)
head(clades)
clades[is.na(clades)] <- -1  # Recode missing values as -1

#clades[which(clades$colony=="129"), 8:10] <- -2  # Recode mortality as -2
clades.m0 <- clades[with(clades, order(clades[, 5], clades[, 6], clades[, 7], 
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
