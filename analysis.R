# LOAD, MERGE, AND QC DATA
source("setup.R")

# # FILTER OUT COLONIES WITH <= 4 OBSERVATIONS -----
nobs <- aggregate(data.frame(obs=Mcap.ff.all$colony), by=list(colony=Mcap.ff.all$colony), FUN=length)
Mcap.ff.all <- droplevels(Mcap.ff.all[Mcap.ff.all$colony %in% nobs[nobs$obs > 4, "colony"], ])
Mcap.ff.all <- droplevels(Mcap.ff.all[Mcap.ff.all$reef!="42", ])

# WAS VISUAL BLEACHING SCORE THE SAME IN BOTH YEARS?
Mcap.ff.all$year <- ifelse(Mcap.ff.all$date < as.Date("2015-05-07"), "y1", "y2")
bscore <- aggregate(data.frame(minscore=Mcap.ff.all$score), 
                    by=list(colony=Mcap.ff.all$colony, year=Mcap.ff.all$year), FUN=min, na.rm=T)
dcast(bscore, colony ~ year, value.var="minscore")

# IDENTIFY COLONIES THAT SHUFFLED SYMBIONTS -----
res <- ldply(levels(Mcap.ff.all$colony), plotpropD)
rownames(res) <- unlist(levels(Mcap.ff.all$colony))
apply(res, 2, table)  # Count number of shufflers determined by each fitting technique
Mcap.ff.all$shuff <- res[Mcap.ff.all$colony, "loess"]

# GROUP COLONIES BY SHUFFLING, BLEACHING, AND DOMINANT CLADE -----
Mcap.ff.all$bleach1 <- ifelse(Mcap.ff.all$colony %in% Mcap.ff.all[Mcap.ff.all$score==1 & Mcap.ff.all$year=="y1", "colony"], "bleach", "notbleached")
Mcap.ff.all$bleach2 <- ifelse(Mcap.ff.all$colony %in% Mcap.ff.all[Mcap.ff.all$score==1 & Mcap.ff.all$year=="y2", "colony"], "bleach", "notbleached")
Mcap.ff.all$group <- as.character(droplevels(interaction(Mcap.ff.all$bleach2, Mcap.ff.all$shuff)))
Mcap.ff.all[Mcap.ff.all$group=="notbleached.noshuff", "group"] <- ifelse(Mcap.ff.all[Mcap.ff.all$group=="notbleached.noshuff", "tdom"]=="C", "notbleached.noshuff.C", "notbleached.noshuff.D")
Mcap.ff.all$group <- factor(Mcap.ff.all$group)
#Mcap.ff.all[Mcap.ff.all$colony=="207", "group"] <- "bleach.shuff"  # assume this colony was C-dominated prior to bleaching
# ADDITIONAL GROUPING FACTOR BY SURVIVAL OF SECOND BLEACHING EVENT
numberofsamplesafter12042015 <- aggregate(data.frame(n=Mcap.ff.all$date>="2015-12-04"), by=list(colony=Mcap.ff.all$colony), FUN=function(x) table(x)[2])
Mcap.ff.all$survival <- ifelse(Mcap.ff.all$colony %in% numberofsamplesafter12042015[numberofsamplesafter12042015$n>=2, "colony"], TRUE, FALSE)
Mcap.ff.all$group2 <- interaction(Mcap.ff.all$survival, Mcap.ff.all$group)

# IDENTIFY AND COUNT COLONIES IN EACH GROUP
cols <- aggregate(Mcap.ff.all$colony, by=list(Mcap.ff.all$group), FUN=function(x) unique(as.character(x)))
ncols <- aggregate(Mcap.ff.all$colony, by=list(Mcap.ff.all$group), FUN=function(x) length(unique(as.character(x))))
cols <- aggregate(Mcap.ff.all$colony, by=list(Mcap.ff.all$group2), FUN=function(x) unique(as.character(x)))
ncols <- aggregate(Mcap.ff.all$colony, by=list(Mcap.ff.all$group2), FUN=function(x) length(unique(as.character(x))))

#Plot Bleaching vs non bleaching dominant symbiont clade
eight11 <- Mcap.ff.all[Mcap.ff.all$date=="2015-08-11",]
st <- table(eight11$dom, eight11$bleach2)
bars <- barplot(st, col=c("blue","red"), width=1, xlim=c(0,6), ylab="Number of Colonies", names.arg=c("Bleached", "Not Bleached"))
text(0.7, 22, labels="n=21", xpd=NA)
text(1.9, 26, labels="n=25", xpd=NA)
legend("topleft",legend=c("D","C"), bty="n", pt.cex=2, pch=22, pt.bg=c("red","blue"))

#Pie chart function for proportion D at a specific timepoint
h <- Mcap.ff.all[(Mcap.ff.all$colony=="71" & Mcap.ff.all$date=="2015-10-21"),]
htable <- c(h$propD, 1-h$propD)
pie(htable)

pieintheface <- function(x,y) {
  h <- Mcap.ff.all[(Mcap.ff.all$colony==x & Mcap.ff.all$date==y),]
  htable <- c(h$propD, 1-h$propD)
  lbls <- c("Clade D","Clade C")
  pct <- round(htable/sum(htable)*100)
  lbls <- paste(lbls,pct)
  lbls <- paste(lbls,"%",sep="")
  pie(htable, col=c("red","blue"), labels=lbls, main=y)
}

pieintheface("71", "2016-02-11")

plotcolony <- function(colony) {
  df <- Mcap.ff.all[Mcap.ff.all$colony==colony, ]
  df <- df[order(df$date), ]
  par(mar=c(5,3,1,1))
  plot(df$date, log(df$tot.SH), type="b", pch=21, cex=2, bg=c("blue","lightblue","pink","red")[df$syms], ylim=c(-9,1), xlab="", ylab="Log SH", xaxt="n")
  dates <- as.Date(c("2014-10-24","2014-11-04","2014-11-24","2014-12-16","2015-01-14","2015-05-06","2015-08-11", "2015-09-14", "2015-10-01", "2015-10-21", "2015-11-04", "2015-12-04","2015-12-17", "2016-01-20", "2016-02-11","2016-03-31"))
  axis(side=1, at=dates, labels=FALSE)
  text(x=dates, y=par("usr")[3]-.2, srt=45, labels=as.character(dates), xpd=NA, pos=2)
  legend("topleft", legend=c("C","C>D","D>C","D"), pch=21, pt.cex=2, pt.bg=c("blue","lightblue","pink","red"))
}

plotcolony(71)

# PLOT SYMBIONT ABUNDANCE AND COMPOSITION FOR INDIVIDUAL COLONIES -----
# XYPLOT ALL COLONIES IN EACH GROUP
xyplot(log(tot.SH) ~ date | group, groups=~colony, ylim=c(-11,1), data=Mcap.ff.all, type="o", cex=0.25)
xyplot(propD ~ date | group, groups=~colony, ylim=c(-0.1,1.1), data=Mcap.ff.all, type="o", cex=0.25)

# XYPLOT INDIVIDUAL COLONIES BY GROUP, RAW DATA
for (g in levels(Mcap.ff.all$group)) {
  df <- subset(Mcap.ff.all, group==g)
  print(doubleYScale(
    # Plot total S/H with GAM fit
    xyplot(log(tot.SH) ~ date | colony, ylim=c(-11,1), data=df, type="o", cex=0.25, main=g),
    # Plot propD with locfit
    xyplot(propD ~ date | colony, ylim=c(-0.1, 1.1), data=df, type="o", cex=0.25)
  ))
}

# XYPLOT INDIVIDUAL COLONIES BY GROUP, FITTED RESPONSES
for (g in levels(Mcap.ff.all$group)) {
  df <- subset(Mcap.ff.all, group==g)
  print(doubleYScale(
    # Plot total S/H with GAM fit
    xyplot(log(tot.SH) ~ days | colony, ylim=c(-11,1), data=df, main=g, panel = function(x, y, ...) {
      panel.xyplot(x, y, cex=0.5, ...)
      dayrange <- seq(min(x), max(x), 1)
      tryCatch({
        m <- gam(y ~ s(x), family="gaussian")
        p <- predict(m, newdata=data.frame(x=dayrange))
        panel.lines(p ~ dayrange)
      },
      error=function(e) {
        m <- gam(y ~ s(x, k=3), family="gaussian")
        p <- predict(m, newdata=data.frame(x=dayrange))
        panel.lines(p ~ dayrange)
      },
      warning=function(w) print(w))
    }),
    # Plot propD with locfit
    xyplot(propD ~ days | colony, ylim=c(-0.1, 1.1), data=df, panel = function(x, y, ...) {
      panel.xyplot(x, y, cex=0.25, ...)
      dayrange <- seq(min(x), max(x), 1)
      tryCatch({
        m <- locfit(y ~ lp(x, nn=1), family="betar", lfproc=locfit.raw)
        p <- predict(m, newdata=data.frame(x=dayrange))
        panel.lines(p ~ dayrange)
        CtoD <- dayrange[which(diff(sign(p-0.5))>0)]
        DtoC <- dayrange[which(diff(sign(p-0.5))<0)]
        panel.xyplot(c(CtoD, DtoC), rep(0.5, length(c(CtoD, DtoC))), pch="*", cex=2, col="red")
      },
      error=function(e) print(e),
      warning=function(w) print(w))
    })
  ))
}

# MODEL SYMBIONT ABUNDANCE AND COMPOSITION FOR EACH GROUP -----
# Exclude groups that didn't quite make it
df <- Mcap.ff.all[as.numeric(Mcap.ff.all$group2) %in% c(2,4,6,8,10), ]
df <- droplevels(df)
# FIT PROPD GAMM BY GROUP
xyplot(propD ~ days | group2, data=df)
propDmod <- gamm4(propD ~ group2 + s(days, by=group2), random=~(1|colony), data=df)
# FIT TOTSH GAMM BY GROUP
xyplot(log(tot.SH) ~ days | group, data=df)
totSHmod <- gamm4(log(tot.SH) ~ group2 + s(days, by=group2), random=~(1|colony), data=df)
# GET FITTED VALUES FOR EACH GROUP
newdata <- expand.grid(days=seq(0,524,1), group2=levels(df$group2))
newdata$tot.SH <- predict(totSHmod$gam, newdata)
newdata$propD <- predict(propDmod$gam, newdata)
newdata$predse <- predict(totSHmod$gam, newdata, se.fit=T)$se.fit
# PLOT FITTED VALUES FOR EACH GROUP
xyplot(tot.SH ~ days, groups=~group2, newdata)
xyplot(propD ~ days, groups=~group2, newdata, ylim=c(0,1))
doubleYScale(xyplot(tot.SH ~ days | group2, newdata, type="l"),
             xyplot(propD ~ days | group2, newdata, type="l", ylim=c(-0.1,1.1)))

# PLOT FITTED RESPONSES FOR EACH GROUP, MULTIPANEL SHUFFLERS vs. NONSHUFFLERS -----
rbPal <- colorRampPalette(c('dodgerblue','red'))
newdata$color <- rbPal(100)[as.numeric(cut(newdata$propD, breaks = 100))]
par(mfrow=c(2,1), mar=c(1,3,1,2), mgp=c(1.5,0.4,0), tcl=-0.25)
plot(NA, ylim=c(-7,0), xlim=range(newdata$days), xaxs="i", xaxt="n", yaxt="n", ylab="")
axis(side=2, at=seq(-7,-1,1), cex.axis=0.75)
dateticks <- seq.Date(as.Date("2014-11-01"), as.Date("2016-02-01"), by="month")
axis(side=1, at=as.numeric(dateticks-as.Date("2014-10-24")), labels=NA)
for (group2 in levels(newdata$group2)[c(1,3,4)]) {
  df <- newdata[newdata$group2==group2, ]
  addpoly(df$days, df$tot.SH - 1.96*df$predse, df$tot.SH + 1.96*df$predse, col=alpha("gray", 0.7))
}
points(tot.SH ~ days, newdata[as.numeric(newdata$group2) %in% c(1,3,4), ], pch=21, col=color, bg=color)
text(par("usr")[1], quantile(par("usr")[3:4], 0.9), pos=4,
     expression(bold("A. Non-shuffling colonies")))
gradient.rect(quantile(par("usr")[1:2], 0.1), quantile(par("usr")[3:4], 0.1),
              quantile(par("usr")[1:2], 0.35), quantile(par("usr")[3:4], 0.175),
              col=rbPal(100), border=NA)
text(quantile(par("usr")[1:2], c(0.1, 0.35)), rep(quantile(par("usr")[3:4], 0.1375), 2), pos=c(2,4), labels=c("C", "D"), cex=0.75)
par(mar=c(2,3,0,2))
plot(NA, ylim=c(-7,0), xlim=range(newdata$days), xaxs="i", xlab="", ylab="", xaxt="n", yaxt="n", xpd=NA)
axis(side=2, at=seq(-7,-1,1), cex.axis=0.75)
mtext(side=2, text="Symbiont abundance (ln S/H)", line=-1.5, outer=T)
dateticks <- seq.Date(as.Date("2014-11-01"), as.Date("2016-02-01"), by="month")
axis(side=1, at=as.numeric(dateticks-as.Date("2014-10-24")), labels=format(dateticks, "%b"), cex.axis=0.75)
for (group2 in levels(df$group2)[c(2,5)]) {
  df <- newdata[newdata$group2==group2, ]
  addpoly(df$days, df$tot.SH - 1.96*df$predse, df$tot.SH + 1.96*df$predse, col=alpha("gray", 0.7))
}
points(tot.SH ~ days, newdata[as.numeric(newdata$group2) %in% c(2,5), ], pch=21, col=color, bg=color)
gradient.rect(quantile(par("usr")[1:2], 0.1), quantile(par("usr")[3:4], 0.1),
              quantile(par("usr")[1:2], 0.35), quantile(par("usr")[3:4], 0.175),
              col=rbPal(100), border=NA)
text(quantile(par("usr")[1:2], c(0.1, 0.35)), rep(quantile(par("usr")[3:4], 0.1375), 2), pos=c(2,4), labels=c("C", "D"), cex=0.75)
text(par("usr")[1], quantile(par("usr")[3:4], 0.9), pos=4,
     expression(bold("B. Shuffling colonies")))

# DOES SHUFFLING DEPEND ON REEF?
MC <- unique(Mcap.ff.all[, c("colony", "reef","shuff","bleach", "group", "depth")])
MC$shuff <- factor(MC$shuff)

MCb <- droplevels(MC[MC$bleach=="bleach", ])
plot(MCb$group ~ MCb$reef)
chisq.test(MCb$reef, MCb$group)

MCnb <- droplevels(MC[MC$bleach=="notbleached", ])
plot(MCnb$group ~ MCnb$reef)
chisq.test(MCnb$reef, MCnb$group)

# DOES SHUFFLING DEPEND ON DEPTH?
plot(shuff ~ depth, MCb)
chisq.test(MCb$depth, MCb$shuff)
plot(as.numeric(shuff) ~ depth, MCnb)


# DOES SHUFFLING RELATE TO BLEACHING SEVERITY?
bsev <- aggregate(data.frame(minscore=Mcap.ff.all$tot.SH), 
                    by=list(colony=Mcap.ff.all$colony), FUN=min, na.rm=T)
bsev <- merge(MCb, bsev)

plot(bsev$shuff ~ log(bsev$minscore))
plot(as.numeric(bsev$shuff) ~ log(bsev$minscore))

mod <- glm(shuff ~ depth * log(minscore), family="binomial", data=bsev)
plot(mod)
anova(mod, test="Chisq")
plot(effect("depth:log(minscore)", mod), x.var="minscore")
plot(effect("depth", mod))






# HEATMAP OF SHUFFLING ---------

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
# MODEL TRAJECTORIES USING SPIDA ------
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


# MODEL TRAJECTORIES FOR BOTH YEARS
Mcap.ff.all$days <- as.numeric(Mcap.ff.all$date - as.Date("2014-10-24"))
Mcap.ff.all$days <- scale(Mcap.ff.all$days)
points <- unique(Mcap.ff.all$days)
knots <- points[c(5,7,9,11,13)]
sp <- function(x) gsp(x, knots=knots, degree=c(2,1,2,3,2,1), smooth=c(0,1,1,1,1))
Mcapdf <- Mcap.ff.all[Mcap.ff.all$reef!="42",]
Mcapdf$reef <- factor(Mcapdf$reef, levels=c("44","25","HIMB"))
Mcapdf$batch <- Mcapdf$days < 195
mod.all <- lmerTest::lmer(log(tot.SH) ~ sp(days) * vis * reef + (1 | colony), data=Mcapdf)
#anova(mod.all)
summary(lm(model.response(model.frame(mod.all)) ~ fitted(mod.all)))$r.squared

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
  pred <- expand.grid(days=seq(min(dat$days), max(dat$days), length.out=475), reef=levels(dat$reef), vis=levels(dat$vis))
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
      plot(NA, xlim=range(dat$days), ylim=c(-9,0), xaxt="n", bty="n", tck=-0.03, ylab="ln S/H")
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
  axdates <- seq.Date(as.Date("2014-11-01"), as.Date("2016-02-01"), by="month")
  unscaled <- Mcap.ff.all$days * attr(Mcap.ff.all$days, 'scaled:scale') + attr(Mcap.ff.all$days, 'scaled:center')
  axdays <- as.numeric(axdates - as.Date("2014-10-24"))
  tdays <- scales::rescale(axdays, to=range(dat$days), from=range(unscaled))

  axis(side=1, at=tdays, labels=format(axdates, format="%b\n%y"), padj=1)
  return(list(predlist=predlist, datlist=datlist))
}
pl <- plotreefs(mod.all, 99)

