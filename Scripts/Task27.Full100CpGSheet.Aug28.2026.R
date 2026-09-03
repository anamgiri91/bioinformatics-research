# ================================================================
# Task 27 - The complete 100-CpG sheet, site by site
# Date: Aug 28, 2026
#
# Tasks 18-24 report the 100-CpG window aggregated BY STATE. That
# hides which individual sites carry the flags, and the supervisor
# asked for the whole sheet in the report rather than the summary.
#
# This emits one row per CpG - all 100 - carrying, for both tissues:
#   position, gap to the previous probe, methylation state,
#   cohort median / quartiles of beta,
#   and the count of -1 and +1 flags from each of bio, self and ext.
#
# PRIVACY: counts and quartiles only. No per-sample value appears, so
# no row resolves to an individual. Cohort MIN and MAX are deliberately
# excluded - at n = 53 those ARE single individuals' beta values, which
# is the same reasoning that makes the referenceMeth() tables private
# (CLASSIFICATION.md). Median/Q1/Q3 collapse across samples and stay.
#
# Two figures are produced from the same table:
#   Fig6 - the window as a map: methylation landscape, state track,
#          and where each method's flags actually fall
#   Fig7 - flags per site per method, which makes the n = 53
#          degeneracy visible site by site rather than in aggregate
# ================================================================

options(stringsAsFactors = FALSE)

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d); stop("no dir") }
src_dir <- pick_dir(
  "/home/s_s355/research3.data/TCGA/filter.BRCA.UCEC.Aug17.2022/sorted.TCGA.380355cg.files",
  path.expand("~/Downloads/tcga_brca_380355cg"))
repo <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                 path.expand("~/Desktop/bioinformatics-research"))
out_dir <- file.path(repo, "Results")

inputs <- list(
  Normal = file.path(src_dir, "Sorted.BRCA.53Alive.Normal.380355cg.75col.May28.2026.txt"),
  Tumor  = file.path(src_dir, "Sorted.BRCA.53Alive.Tumor.380355cg.75col.May28.2026.txt"))
flags <- list(
  Normal = list(self=file.path(out_dir,"Task13_Normal_selfref_flags.csv"),
                ext =file.path(out_dir,"Task13_Normal_extref_flags.csv")),
  Tumor  = list(self=file.path(out_dir,"Task13_Tumor_selfref_flags.csv"),
                ext =file.path(out_dir,"Task13_Tumor_extref_flags.csv")))

N_SITES <- 100
ST <- c("L","LM","M","HM","H","R")
P_LEVEL <- 0.01
PAL <- c(L="#4E79A7", LM="#76B7B2", M="#B0B0B0", HM="#F28E2B", H="#E15759", R="#59A14F")

rd <- function(p, tab=FALSE) as.data.frame(data.table::fread(
  p, sep=if(tab) "\t" else ",", header=TRUE, check.names=FALSE, showProgress=FALSE))

build_bio <- function(meth, states) {
  o <- matrix(0, nrow(meth), ncol(meth), dimnames=dimnames(meth))
  for (i in seq_len(nrow(meth))) {
    st <- states[i]; v <- as.numeric(meth[i,])
    if (all(is.na(v))) { o[i,] <- NA; next }
    if (!st %in% c("L","LM","H","HM")) next
    hi <- st %in% c("L","LM"); p <- if (hi) 1-P_LEVEL else P_LEVEL
    thr <- quantile(v, p, na.rm=TRUE, names=FALSE)
    o[i,] <- if (hi) as.integer(v > thr) else -as.integer(v < thr)
    o[i, is.na(v)] <- NA
  }
  o
}

cat("=== TASK 27 - full 100-CpG sheet ===\n\n")

normal <- rd(inputs$Normal, TRUE); tumor <- rd(inputs$Tumor, TRUE)
n22 <- normal[normal$Chromosome=="chr22",]; n22 <- n22[order(n22$Start),]
sel <- n22[seq_len(N_SITES),]
ids <- sel$Composite.Element.REF

sheet <- data.frame(
  site.index = seq_len(N_SITES),
  cgID = ids, chr = sel$Chromosome, pos = sel$Start,
  gap.to.prev.bp = c(NA, diff(sel$Start)))

for (tissue in c("Normal","Tumor")) {
  src  <- if (tissue=="Normal") normal else tumor
  cols <- paste0(substr(tissue,1,1),1:53)
  meta <- src[match(ids, src$Composite.Element.REF),]
  B <- as.matrix(meta[,cols,drop=FALSE]); rownames(B) <- ids; storage.mode(B) <- "numeric"
  states <- meta$methy.state

  bio  <- build_bio(B, states)
  pick <- function(f){ m <- rd(f); mm <- as.matrix(m[match(ids,m$cgID),cols,drop=FALSE])
    storage.mode(mm) <- "numeric"; mm }
  self <- pick(flags[[tissue]]$self); ext <- pick(flags[[tissue]]$ext)

  pre <- substr(tissue,1,1)
  sheet[[paste0("state.",tolower(tissue))]] <- states
  sheet[[paste0("median.beta.",pre)]] <- round(apply(B,1,median,na.rm=TRUE),4)
  sheet[[paste0("q25.beta.",pre)]]    <- round(apply(B,1,quantile,.25,na.rm=TRUE,names=FALSE),4)
  sheet[[paste0("q75.beta.",pre)]]    <- round(apply(B,1,quantile,.75,na.rm=TRUE,names=FALSE),4)
  for (nm in c("bio","self","ext")) {
    M <- get(nm)
    sheet[[paste0(nm,".neg1.",pre)]] <- rowSums(M==-1, na.rm=TRUE)
    sheet[[paste0(nm,".pos1.",pre)]] <- rowSums(M== 1, na.rm=TRUE)
  }
  sheet[[paste0("ext.evaluable.",pre)]] <- rowSums(!is.na(ext))
  assign(paste0("keep.",tissue), list(B=B, states=states, bio=bio, self=self, ext=ext))
}

path <- file.path(out_dir,"Task27_Chr22_Full100CpG_SiteTable.csv")
write.csv(sheet, path, row.names=FALSE)
cat("written:", path, " (", nrow(sheet), " rows x ", ncol(sheet), " cols )\n\n", sep="")

cat("state distribution over the 100 sites:\n")
print(rbind(Normal=table(factor(sheet$state.normal,levels=ST)),
            Tumor =table(factor(sheet$state.tumor, levels=ST))))
cat("\nsites carrying at least one ext flag: Normal ",
    sum(sheet$ext.neg1.N+sheet$ext.pos1.N > 0), " / 100,  Tumor ",
    sum(sheet$ext.neg1.T+sheet$ext.pos1.T > 0), " / 100\n", sep="")
cat("first 12 rows:\n")
print(head(sheet[,c("site.index","cgID","pos","gap.to.prev.bp","state.normal",
                    "median.beta.N","bio.neg1.N","bio.pos1.N",
                    "self.neg1.N","self.pos1.N","ext.neg1.N","ext.pos1.N")],12),
      row.names=FALSE)

# ------------------------------------------------
# Fig 6 - the window as a map
# ------------------------------------------------
# x is SITE INDEX, not position: the 100 sites span 5.17 Mb with a
# 3.41 Mb gap, so a position axis would collapse 99 of them into a
# sliver. Position is carried on the axis labels instead.
png(file.path(out_dir,"Fig6_100CpG_WindowMap.png"), 2200, 1560, res=190)
layout(matrix(1:4,4,1), heights=c(3.1,0.66,1.75,1.75))
par(mgp=c(3.1,0.7,0), las=1, cex.axis=0.8, cex.lab=0.9)

x <- sheet$site.index

# --- panel 1: methylation landscape
par(mar=c(0.5,5.4,3.4,1.2))
plot(NA, xlim=c(0.5,100.5), ylim=c(0,1.02), xaxt="n",
     ylab=expression("methylation "*beta), xlab="")
title(main="The 100-CpG window, site by site (Normal tissue)", line=2.1, cex.main=1.05)
mtext("vertical bar = middle 50% of the 53 samples   |   dot = cohort median, coloured by state",
      3, line=0.55, cex=0.66, col="grey40")
segments(x, sheet$q25.beta.N, x, sheet$q75.beta.N, col="grey72", lwd=3.2)
points(x, sheet$median.beta.N, pch=19, cex=0.72, col=PAL[sheet$state.normal])
legend(1, 0.34, names(PAL), col=PAL, pch=19, bty="n", ncol=3, cex=0.74,
       pt.cex=1.0, title="methylation state", title.adj=0)

# --- panel 2: state tracks
par(mar=c(0.3,5.4,0.3,1.2))
plot(NA, xlim=c(0.5,100.5), ylim=c(0,2), axes=FALSE, xlab="", ylab="")
rect(x-0.5, 1.05, x+0.5, 1.95, col=PAL[sheet$state.normal], border=NA)
rect(x-0.5, 0.05, x+0.5, 0.95, col=PAL[sheet$state.tumor],  border=NA)
mtext("state, Normal", 2, at=1.5, las=1, cex=0.62, line=0.4)
mtext("state, Tumour", 2, at=0.5, las=1, cex=0.62, line=0.4)

# --- panels 3 and 4: flags per site, shared scale
ymax <- max(sheet$ext.neg1.N+sheet$ext.pos1.N, sheet$ext.neg1.T+sheet$ext.pos1.T)
par(mar=c(0.5,5.4,2.1,1.2))
barplot(rbind(sheet$ext.neg1.N, sheet$ext.pos1.N), beside=FALSE, space=0,
        col=c("#8FAADC","#E8A0A0"), border=NA, ylim=c(0,ymax*1.08),
        ylab="samples flagged", xaxt="n")
title(main="External reference \u2014 Normal: 65 of 100 sites carry a flag",
      line=0.7, cex.main=0.92, font.main=1)
legend("topleft", c("hypo (-1)","hyper (+1)"), fill=c("#8FAADC","#E8A0A0"),
       border=NA, bty="n", cex=0.74, horiz=TRUE)

par(mar=c(4.0,5.4,2.1,1.2))
bp <- barplot(rbind(sheet$ext.neg1.T, sheet$ext.pos1.T), beside=FALSE, space=0,
        col=c("#8FAADC","#E8A0A0"), border=NA, ylim=c(0,ymax*1.08),
        ylab="samples flagged", xaxt="n")
title(main="External reference \u2014 Tumour: 97 of 100 sites carry a flag, and far more per site",
      line=0.7, cex.main=0.92, font.main=1)
at <- seq(1,100,by=11)
axis(1, at=bp[at], labels=format(sheet$pos[at], big.mark=","), cex.axis=0.62)
mtext("chr22 position (bp)", 1, line=2.0, cex=0.72)
mtext("sites are drawn equally spaced; on the chromosome they are not \u2014 gaps run from 13 bp to 3.41 Mb",
      1, line=3.0, cex=0.64, col="grey40")
dev.off()
cat("\nwritten: Fig6_100CpG_WindowMap.png\n")

# ------------------------------------------------
# Fig 7 - flags per site, per method
# ------------------------------------------------
png(file.path(out_dir,"Fig7_100CpG_MethodBySite.png"), 2200, 1150, res=190)
par(mfrow=c(3,1), mar=c(0.5,5.4,2.0,1.2), mgp=c(3.2,0.7,0), las=1,
    cex.axis=0.8, cex.lab=0.9)
tot <- function(nm) sheet[[paste0(nm,".neg1.N")]] + sheet[[paste0(nm,".pos1.N")]]
ttl <- c(bio="bio - state-aware rule (0 or 1 per site, by construction)",
         self="self - own 53 samples as reference (exactly 2 per site, always)",
         ext="ext - external 747-sample tcga panel (varies: 0 to 8)")
ymx <- max(sapply(c("bio","self","ext"), function(n) max(tot(n))))
for (k in seq_along(ttl)) {
  nm <- names(ttl)[k]
  if (k==3) par(mar=c(3.4,5.4,2.0,1.2))
  barplot(rbind(sheet[[paste0(nm,".neg1.N")]], sheet[[paste0(nm,".pos1.N")]]),
          beside=FALSE, space=0, col=c("#8FAADC","#E8A0A0"), border=NA,
          ylim=c(0,ymx*1.15), ylab="samples flagged", xaxt="n",
          main=ttl[[nm]], cex.main=0.92, font.main=1)
  if (k==3) {
    axis(1, at=seq(0.5,99.5,by=11), labels=seq(1,100,by=11), cex.axis=0.72)
    mtext("site index within the window (1-100, ordered by position)", 1,
          line=2.1, cex=0.68, col="grey40")
  }
}
dev.off()
cat("written: Fig7_100CpG_MethodBySite.png\n")

cat("\nper-site flag totals, Normal (distinct values per method):\n")
for (nm in c("bio","self","ext"))
  cat(sprintf("  %-5s min %d  max %d  distinct %d\n", nm,
      min(tot(nm)), max(tot(nm)), length(unique(tot(nm)))))
cat("\nTASK 27 COMPLETE\n")
