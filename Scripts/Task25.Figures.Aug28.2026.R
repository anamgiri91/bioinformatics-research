# ================================================================
# Task 25 - Figures for results.md
# Date: Aug 28, 2026
#
# Five figures, all built from Task 20-23 output CSVs plus the chr22
# beta matrix. Base graphics only, so this runs anywhere R does.
# ================================================================

options(stringsAsFactors = FALSE)
pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d); stop("no dir") }
repo <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                 path.expand("~/Desktop/bioinformatics-research"))
src  <- pick_dir("/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
                 path.expand("~/Downloads/tcga_brca_380355cg"))
R <- file.path(repo, "Results")
ST <- c("L","LM","M","HM","H","R")
PAL <- c(L="#4E79A7", LM="#76B7B2", M="#B0B0B0", HM="#F28E2B", H="#E15759", R="#59A14F")

png_open <- function(f, w=1800, h=1100) {
  png(file.path(R, f), width=w, height=h, res=190)
  par(mar=c4 <- c(4.5,4.8,3.2,1.4), mgp=c(2.9,0.8,0), las=1,
      cex.axis=0.85, cex.lab=0.95, cex.main=1.05, font.main=1)
}
done <- function(f) { dev.off(); cat("  written:", file.path(R,f), "\n") }

# ---------------- Fig 1: flag rate by state and method -----------
f <- "Fig1_FlagRateByState.png"; png_open(f, 2100, 1000)
par(mfrow=c(1,2))
for (ti in c("Normal","Tumor")) {
  d <- read.csv(file.path(R, sprintf("Task23_StateFlagTable_extref_%s.csv", ti)))
  d$rate <- 100*(d$count.neg1+d$count.pos1)/d$n.evaluated
  d <- d[match(ST, d$methy.state),]
  bp <- barplot(rbind(100*d$count.neg1/d$n.evaluated, 100*d$count.pos1/d$n.evaluated),
    beside=FALSE, names.arg=ST, col=c("#8FAADC","#E8A0A0"), border=NA,
    ylab="% of cells flagged", main=paste0(ti, " - external reference, all 380,355 CpGs"),
    ylim=c(0, max(d$rate)*1.25))
  text(bp, d$rate+max(d$rate)*0.05, sprintf("%.1f%%", d$rate), cex=0.8)
  legend("topleft", c("hypo (-1)","hyper (+1)"), fill=c("#8FAADC","#E8A0A0"),
         border=NA, bty="n", cex=0.85)
}
done(f)

# ---------------- Fig 2: magnitude of flagged differences --------
f <- "Fig2_FlagMagnitudeByState.png"; png_open(f, 2100, 1000)
par(mfrow=c(1,2))
m <- read.csv(file.path(R,"Task21_Chr22_ExtFlagMagnitude.csv"))
for (ti in c("Normal","Tumor")) {
  d <- m[m$tissue==ti,]; d <- d[match(ST, d$methy.state),]
  bp <- barplot(d$median.abs.dev, names.arg=ST, col=PAL[ST], border=NA,
    ylab=expression("median |"*Delta*beta*"| of flagged cells"),
    main=paste0(ti, " - chr22: how big is a flagged difference?"),
    ylim=c(0, max(m$median.abs.dev)*1.3))
  abline(h=0.10, lty=2, col="grey35"); abline(h=0.15, lty=3, col="grey55")
  text(bp, d$median.abs.dev+0.012, sprintf("%.3f", d$median.abs.dev), cex=0.8)
  text(par("usr")[2], 0.104, "0.10", adj=c(1.1,0), cex=0.72, col="grey35")
  text(par("usr")[2], 0.154, "0.15 (epimutacions floor)", adj=c(1.05,0), cex=0.72, col="grey55")
}
done(f)

# ---------------- Fig 3: threshold geometry ----------------------
f <- "Fig3_ThresholdGeometry.png"; png_open(f, 2100, 1000)
par(mfrow=c(1,2))
for (ti in c("Normal","Tumor")) {
  g <- read.csv(file.path(R, sprintf("Task22_Chr22_ThresholdGeometry_%s.csv", ti)))
  g <- g[match(ST, g$methy.state),]
  plot(NA, xlim=c(0,1), ylim=c(0.5, length(ST)+0.5), yaxt="n",
       xlab=expression(beta), ylab="",
       main=paste0(ti, " - where the outlier zones sit"))
  axis(2, at=seq_along(ST), labels=ST)
  for (i in seq_along(ST)) {
    y <- i
    rect(0, y-0.30, g$N.threshold[i], y+0.30, col="#8FAADC", border=NA)
    rect(g$N.threshold[i], y-0.30, g$P.threshold[i], y+0.30, col="grey92", border=NA)
    rect(g$P.threshold[i], y-0.30, 1, y+0.30, col="#E8A0A0", border=NA)
    text(g$P.threshold[i], y, sprintf(" %.3f", 1-g$P.threshold[i]),
         adj=c(0,0.5), cex=0.7, col="grey20")
  }
  legend("bottomright", c("hypo zone","normal","hyper zone"),
         fill=c("#8FAADC","grey92","#E8A0A0"), border=NA, bty="n", cex=0.78)
  mtext("labels = width of the hyper zone", side=3, line=0.1, cex=0.68, col="grey40")
}
done(f)

# ---------------- Fig 4: stability under noise -------------------
f <- "Fig4_FlagStability.png"; png_open(f, 2100, 1000)
par(mfrow=c(1,2))
for (ti in c("Normal","Tumor")) {
  s <- read.csv(file.path(R, sprintf("Task22_Chr22_FlagStability_%s.csv", ti)))
  cn <- grep("^pct.unstable", names(s), value=TRUE)
  sts <- sub("pct.unstable.", "", cn)
  plot(NA, xlim=range(s$noise.sd), ylim=c(0, 50), log="x",
       xlab="added measurement noise (sd, beta units)",
       ylab="% of baseline flags that change call",
       main=paste0(ti, " - are the flags reproducible?"))
  abline(v=c(0.01,0.03), col="grey85", lty=3)
  for (k in seq_along(cn))
    lines(s$noise.sd, s[[cn[k]]], type="b", pch=19, cex=0.7,
          col=PAL[sts[k]], lwd=2)
  legend("topleft", sts, col=PAL[sts], lwd=2, pch=19, bty="n", cex=0.8, ncol=2)
  mtext("grey band = 450k technical replicate sd", side=3, line=0.1,
        cex=0.68, col="grey40")
}
done(f)

# ---------------- Fig 5: the N37 cluster -------------------------
# NOTE: the probes in this window are NOT evenly spaced - there is a
# 6.4 kb stretch with no 450k coverage in the middle. Joining points
# across it draws a peak that was never measured, so the series is
# broken wherever consecutive probes are more than 1 kb apart.
f <- "Fig5_N37_LocalCluster.png"; png_open(f, 1900, 1050)
suppressWarnings({
  rd <- function(p) if (requireNamespace("data.table", quietly=TRUE))
    as.data.frame(data.table::fread(p, sep="\t", header=TRUE, showProgress=FALSE))
    else read.table(p, header=TRUE, sep="\t")
  d <- rd(file.path(src,"Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt"))
})
d <- d[d$Chromosome=="chr22",]; d <- d[order(d$Start),]
w <- d[d$Start >= 16600500 & d$Start <= 16611000,]
smp <- paste0("N",1:53)
B <- as.matrix(w[,smp]); storage.mode(B) <- "numeric"
pos <- w$Start

# break every series wherever the probe gap exceeds 1 kb
brk <- which(diff(pos) > 1000)
ins <- function(v) { out <- v; for (b in rev(brk)) out <- append(out, NA, after=b); out }
posb <- ins(pos)

plot(NA, xlim=range(pos), ylim=c(0, max(B, na.rm=TRUE)*1.05),
     xlab="chr22 position (bp)", ylab=expression(beta),
     main="N37: one focal hypermethylation event, counted ten times")
for (b in brk)
  rect(pos[b], -1, pos[b+1], 2, col="grey96", border=NA)
for (j in seq_len(ncol(B)))
  lines(posb, ins(B[,j]), col="#C8C8C880", lwd=1)
lines(posb, ins(apply(B,1,median)), col="grey25", lwd=2.4, lty=2)
lines(posb, ins(B[,"N37"]), col="#D62728", lwd=3)
points(pos, B[,"N37"], col="#D62728", pch=19, cex=0.9)
rug(pos, col="grey30", lwd=1.4)
for (b in brk)
  text(mean(pos[c(b,b+1)]), max(B,na.rm=TRUE)*0.5,
       sprintf("no probes\n%s bp", format(pos[b+1]-pos[b], big.mark=",")),
       col="grey55", cex=0.72)
legend("topright", c("N37","cohort median","other 52 samples"),
       col=c("#D62728","grey25","#C8C8C8"), lwd=c(3,2.4,1), lty=c(1,2,1),
       bty="n", cex=0.85, bg="white")
mtext(sprintf("%d CpGs across %s bp; ticks mark actual probe positions",
      nrow(w), format(diff(range(pos)), big.mark=",")),
      side=3, line=0.1, cex=0.7, col="grey40")
done(f)

cat("\nTASK 25 COMPLETE - 5 figures\n")
