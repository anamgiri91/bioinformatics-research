# ================================================================
# Task 45 - RRBS item 4 by methylation state, and nearby CG sites
#           in the same style as the earlier TCGA tasks
# Date: Sep 29, 2026
#
# THE REQUEST
#
#   Item 4, RRBS 18-sample lung data:
#     (1) 6-number summary + outlier.IQR2 + outlier.IQR3, NA count;
#         how many CG sites have 0 NA.
#     (2) methylation states for those sites; check the patterns of
#         CG sites close to each other.
#   Clarified 2026-09-29: report (1) by STATE, not site by site, and
#   check nearby sites the way earlier tasks did for the TCGA data.
#
# WHAT EARLIER TASKS DID, AND WHAT IS COPIED HERE
#
#   Task 9      boxplots of outliers.coef2 / coef3 per site, by state
#               -> Part A, Fig30 (b)
#   Tasks 17-18 a window of 100 consecutive CG sites; for each state,
#               the share of cells flagged -1 / 0 / +1
#               -> Part B, with IQR outliers low / normal / high in
#                  place of the -1 / 0 / +1 flags
#   Task 20     a second, TIGHT window: the 100 consecutive sites with
#               the smallest genomic span; both windows side by side;
#               is a flagged sample also flagged at neighbouring sites?
#               -> Part B windows index100 and tight100; Part C
#   Fig6        the window site by site: middle 50% and median per
#               site, coloured by state; a state strip; samples
#               flagged per site -> Fig31
#
# INPUT
#   Data/NN.hg38.18P.forw.chr22.w.header.txt   (gitignored)
#   Results/Task39_RRBS_CompleteSiteStates.csv.gz  (states, TSG rule)
#   Only the 38,564 complete sites (no NA) carry a state, so the NA
#   count is 0 in every state; the NA table stays in Task 39.
#
# DEFINITIONS
#   IQR outlier (per site, over its 18 values, type-7 quartiles):
#     low  = below Q1 - k * IQR,  high = above Q3 + k * IQR,
#     k = 2 (IQR2) or 3 (IQR3). Recomputed here and checked against
#     Task 39's outliers.coef2 / coef3.
#   Six-number summary by state: over all values (sites x 18) of the
#     sites in that state, pooled.
#   Neighbour: the previous or next complete site, if within 1,000
#     bases. An outlier "has a neighbour" when the same sample is an
#     outlier in the same direction at a neighbour.
#   Event: 3 or more consecutive complete sites, each within 1,000
#     bases of the next, where the same sample is an outlier in the
#     same direction (the Task 21 / epimutacions unit).
#   Chance: shuffle the 18 samples at every site independently (the
#     same shuffle for low and high), 20 times, seed 20260929.
#
# OUTPUTS (Results/)
#   Task45_RRBS_ByState_SixNumber.csv    Part A
#   Task45_RRBS_ByState_IQROutliers.csv  Part A
#   Task45_RRBS_WindowComparison.csv     Part B
#   Task45_RRBS_WindowStateSummary.csv   Part B
#   Task45_RRBS_OutlierNeighbours.csv    Part C
#   Task45_RRBS_OutlierRunLengths.csv    Part C
#   Fig30_RRBS_ByState.png               Part A
#   Fig31_RRBS_TwoWindows.png            Part B
#   Fig32_RRBS_OutlierRuns.png           Part C
#
# PRIVACY
#   Nothing here names a sample. Window positions are chromosome
#   coordinates only.
# ================================================================

suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(ggplot2)
  library(gridExtra); library(grid); library(scales)
})

pick_dir <- function(...) { for (d in c(...)) if (dir.exists(d)) return(d)
  stop("no candidate directory exists") }
repo_dir <- pick_dir("/mmfs1/home/wln26/Experiments.Outlier.July31.2026",
                     path.expand("~/Desktop/bioinformatics-research"))
res  <- file.path(repo_dir, "Results")
rrbs <- file.path(repo_dir, "Data", "NN.hg38.18P.forw.chr22.w.header.txt")
r4 <- function(x) round(x, 4)
pct <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
set.seed(20260929)
states <- c("L", "LM", "M", "HM", "H", "R")
N_SHUFFLE <- 20
NEAR <- 1000

# ----------------------------------------------------------------
# 0. data: the complete sites and their states
# ----------------------------------------------------------------
d <- fread(rrbs, na.strings = "NA", showProgress = FALSE)
samp <- grep("^P[0-9]+NF$", names(d), value = TRUE)
cs <- fread(file.path(res, "Task39_RRBS_CompleteSiteStates.csv.gz"))
row_in_file <- match(cs$pos, d$pos)
stopifnot(!anyNA(row_in_file), !is.unsorted(cs$pos, strictly = TRUE))
cm <- as.matrix(d[row_in_file, ..samp])
stopifnot(!anyNA(cm), ncol(cm) == 18)
rm(d); invisible(gc())
cs[, state := factor(state, levels = states)]
n <- nrow(cs); pos <- cs$pos

q <- rowQuantiles(cm, probs = c(.25, .75), type = 7L)
iqr <- q[, 2] - q[, 1]
low2  <- cm < q[, 1] - 2 * iqr; high2 <- cm > q[, 2] + 2 * iqr
low3  <- cm < q[, 1] - 3 * iqr; high3 <- cm > q[, 2] + 3 * iqr
stopifnot(rowSums(low2 | high2) == cs$outliers.coef2,
          rowSums(low3 | high3) == cs$outliers.coef3)
cat("complete sites:", n, "  IQR2 outlier cells:", sum(low2 | high2),
    "  IQR3 outlier cells:", sum(low3 | high3), "\n")

rows_of <- function(s) if (s == "All") seq_len(n) else which(cs$state == s)

# ----------------------------------------------------------------
# Part A. by state: six-number summary and IQR outliers
# ----------------------------------------------------------------
six <- rbindlist(lapply(c(states, "All"), function(s) {
  r <- rows_of(s); v <- as.vector(cm[r, , drop = FALSE])
  qq <- quantile(v, c(0, .25, .5, .75, 1), type = 7, names = FALSE)
  data.table(state = s, sites = length(r), values = length(v),
             Min = qq[1], Q1 = qq[2], Median = qq[3], Mean = mean(v), Q3 = qq[4], Max = qq[5],
             median_site_IQR = median(iqr[r]),
             median_site_SD = median(rowSds(cm[r, , drop = FALSE])))
}))
num <- c("Min", "Q1", "Median", "Mean", "Q3", "Max", "median_site_IQR", "median_site_SD")
six[, (num) := lapply(.SD, r4), .SDcols = num]
fwrite(six, file.path(res, "Task45_RRBS_ByState_SixNumber.csv"))
cat("\nSix-number summary of all values, by state:\n"); print(six)

iq <- rbindlist(lapply(c(states, "All"), function(s) {
  r <- rows_of(s); cells <- length(r) * 18
  one <- function(lo, hi, k) {
    nl <- sum(lo[r, ]); nh <- sum(hi[r, ]); cnt <- rowSums(lo[r, , drop = FALSE] | hi[r, , drop = FALSE])
    setNames(list(r4(mean(cnt > 0)), r4(mean(cnt)), max(cnt), nl, nh,
                  r4(nl / cells), r4(1 - (nl + nh) / cells), r4(nh / cells)),
             paste0(c("share_sites_any_", "mean_per_site_", "max_per_site_", "low_cells_",
                      "high_cells_", "pct_low_", "pct_normal_", "pct_high_"), k))
  }
  c(list(state = s, sites = length(r), cells = cells),
    one(low2, high2, "IQR2"), one(low3, high3, "IQR3"))
}))
fwrite(iq, file.path(res, "Task45_RRBS_ByState_IQROutliers.csv"))
cat("\nIQR outliers by state:\n"); print(iq)

# ----------------------------------------------------------------
# Part B. two windows of 100 consecutive complete sites (Task 20)
# ----------------------------------------------------------------
span <- pos[100:n] - pos[1:(n - 99)]
t0 <- which.min(span)
win <- list(index100 = 1:100, tight100 = t0:(t0 + 99))
wcmp <- rbindlist(lapply(names(win), function(w) {
  r <- win[[w]]; g <- diff(pos[r]); mix <- table(cs$state[r])
  data.table(window = w, first_pos = pos[r[1]], last_pos = pos[r[100]],
             span_bp = pos[r[100]] - pos[r[1]], median_gap_bp = median(g), largest_gap_bp = max(g),
             pairs_within_1kb = sum(g <= NEAR), truly_consecutive_pairs = sum(diff(row_in_file[r]) == 1),
             state_mix = paste(names(mix), as.integer(mix), collapse = " "))
}))
fwrite(wcmp, file.path(res, "Task45_RRBS_WindowComparison.csv"))
cat("\nThe two windows:\n"); print(wcmp)

wss <- rbindlist(lapply(names(win), function(w) rbindlist(lapply(states, function(s) {
  r <- win[[w]][cs$state[win[[w]]] == s]
  if (!length(r)) return(NULL)
  cells <- length(r) * 18
  nl2 <- sum(low2[r, ]); nh2 <- sum(high2[r, ]); nl3 <- sum(low3[r, ]); nh3 <- sum(high3[r, ])
  data.table(window = w, state = s, sites = length(r), cells = cells,
             IQR2_low = nl2, IQR2_normal = cells - nl2 - nh2, IQR2_high = nh2,
             pct_IQR2_low = r4(nl2 / cells), pct_IQR2_normal = r4(1 - (nl2 + nh2) / cells),
             pct_IQR2_high = r4(nh2 / cells),
             IQR3_low = nl3, IQR3_high = nh3,
             pct_IQR3_low = r4(nl3 / cells), pct_IQR3_high = r4(nh3 / cells))
}))))
fwrite(wss, file.path(res, "Task45_RRBS_WindowStateSummary.csv"))
cat("\nState summary in each window (IQR outlier cells, low / normal / high):\n"); print(wss)

# ----------------------------------------------------------------
# Part C. is an outlier sample also an outlier next door? (Task 20)
# ----------------------------------------------------------------
link <- c(FALSE, diff(pos) <= NEAR)          # site i within 1 kb of site i-1
has_nb <- function(M) {                        # same sample, same direction, next door
  prevM <- rbind(FALSE, M[-n, , drop = FALSE]) & link
  nextM <- rbind(M[-1, , drop = FALSE], FALSE) & c(link[-1], FALSE)
  M & (prevM | nextM)
}
events <- function(M) {                        # runs of >= 3 linked sites, per sample
  ev <- 0L; cells <- 0L
  for (j in seq_len(ncol(M))) {
    x <- M[, j]
    grp <- cumsum(!(x & c(FALSE, x[-n]) & link))
    rl <- tabulate(grp[x])
    ev <- ev + sum(rl >= 3); cells <- cells + sum(rl[rl >= 3])
  }
  c(events = ev, cells_in_events = cells)
}
RUN_BINS <- c("1", "2", "3", "4", "5", "6 or more")
run_hist <- function(M) {                      # how many runs of each length, per sample
  h <- integer(6)
  for (j in seq_len(ncol(M))) {
    x <- M[, j]
    grp <- cumsum(!(x & c(FALSE, x[-n]) & link))
    rl <- tabulate(grp[x]); rl <- rl[rl > 0]
    h <- h + tabulate(pmin(rl, 6L), nbins = 6)
  }
  setNames(h, RUN_BINS)
}
summarise_nb <- function(lo, hi, rows) {
  nb <- has_nb(lo) | has_nb(hi)
  out <- (lo | hi)[rows, , drop = FALSE]
  c(outlier_cells = sum(out), with_neighbour = sum(nb[rows, ]),
    share_with_neighbour = sum(nb[rows, ]) / sum(out))
}
shuffle_rows <- function() {                   # one shuffle of samples at every site
  P <- t(vapply(seq_len(n), function(i) sample.int(18L), integer(18)))
  idx <- cbind(rep(seq_len(n), 18), as.vector(P))
  list(lo = matrix(low2[idx], n, 18), hi = matrix(high2[idx], n, 18))
}
scopes <- c(list(`all complete sites` = seq_len(n)), win)
obs <- lapply(scopes, function(r) summarise_nb(low2, high2, r))
obs_ev <- events(low2) + events(high2)
obs_hist <- run_hist(low2) + run_hist(high2)
shuf <- replicate(N_SHUFFLE, {
  s <- shuffle_rows()
  c(sapply(scopes, function(r) summarise_nb(s$lo, s$hi, r)["share_with_neighbour"]),
    events(s$lo) + events(s$hi),
    setNames(run_hist(s$lo) + run_hist(s$hi), paste0("h.", RUN_BINS)))
})
nbt <- rbindlist(lapply(names(scopes), function(sc) data.table(
  scope = sc, outlier_cells = obs[[sc]][["outlier_cells"]],
  with_neighbour = obs[[sc]][["with_neighbour"]],
  share_with_neighbour = r4(obs[[sc]][["share_with_neighbour"]]),
  share_with_neighbour_shuffled = r4(mean(shuf[paste0(sc, ".share_with_neighbour"), ])))))
nbt[scope == "all complete sites", `:=`(
  events_3plus = obs_ev[["events"]], events_3plus_shuffled = round(mean(shuf["events", ]), 1),
  cells_in_events = obs_ev[["cells_in_events"]],
  cells_in_events_shuffled = round(mean(shuf["cells_in_events", ]), 1))]
fwrite(nbt, file.path(res, "Task45_RRBS_OutlierNeighbours.csv"))
cat("\nIQR2 outliers: same sample also an outlier next door (within 1,000 bases):\n"); print(nbt)
rlt <- data.table(run_length = RUN_BINS, runs_observed = as.integer(obs_hist),
                  runs_shuffled = round(rowMeans(shuf[paste0("h.", RUN_BINS), , drop = FALSE]), 1))
stopifnot(sum(rlt$runs_observed[3:6]) == obs_ev[["events"]])
fwrite(rlt, file.path(res, "Task45_RRBS_OutlierRunLengths.csv"))
cat("\nRuns of IQR2 outliers (same sample, same direction, neighbouring sites):\n"); print(rlt)

# ----------------------------------------------------------------
# Figures (same style as Task 44)
# ----------------------------------------------------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; GRAY <- "#c3c2b7"; BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
STATE_COL <- c(L = "#86b6ef", LM = "#5598e7", M = "#2a78d6", HM = "#1c5cab", H = "#104281", R = ORANGE)
theme_set(theme_minimal(base_size = 11, base_family = "sans") +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        panel.background = element_rect(fill = SURF, colour = NA),
        panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = INK2, size = 9.5), axis.title = element_text(colour = INK2, size = 10),
        plot.title = element_text(colour = INK, face = "bold", size = 12.5), plot.title.position = "plot",
        plot.subtitle = element_text(colour = INK2, size = 10, lineheight = 1.1),
        plot.caption = element_text(colour = MUTED, size = 8.5, hjust = 0), plot.caption.position = "plot",
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        legend.text = element_text(colour = INK2, size = 9.5),
        plot.margin = margin(10, 16, 8, 10)))
save_fig <- function(name, plot, w, h) {
  ggsave(file.path(res, name), plot, width = w, height = h, dpi = 200, bg = SURF)
  cat("  written:", file.path(res, name), "\n")
}
titled <- function(title, ..., widths = NULL, heights = NULL, ncol = 2) arrangeGrob(..., ncol = ncol,
  widths = widths, heights = heights, top = textGrob(title, x = unit(12, "pt"), hjust = 0,
  gp = gpar(fontsize = 14, fontface = "bold", col = INK)))

# Fig30 (a) the six-number summary of each state, drawn directly
b30 <- six[state %in% states][, state := factor(state, levels = states)]
p30a <- ggplot(b30, aes(x = state)) +
  geom_segment(aes(xend = state, y = Min, yend = Max), colour = MUTED, linewidth = 0.6) +
  geom_rect(aes(xmin = as.integer(state) - 0.28, xmax = as.integer(state) + 0.28, ymin = Q1, ymax = Q3),
            fill = BLUE, alpha = 0.18, colour = BLUE, linewidth = 0.4) +
  geom_segment(aes(x = as.integer(state) - 0.28, xend = as.integer(state) + 0.28, y = Median, yend = Median),
               colour = BLUE, linewidth = 1) +
  geom_point(aes(y = Mean), shape = 21, fill = ORANGE, colour = SURF, stroke = 0.7, size = 3) +
  geom_text(aes(y = 1.06, label = paste0(comma(sites), "\nsites")), colour = INK2, size = 2.9, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.2)) +
  labs(title = "(a) Six-number summary of all values, by state",
       subtitle = "Line: Min to Max. Box: Q1 to Q3, bar at the median. Orange dot: mean.",
       x = "Methylation state", y = "Methylation level (beta)",
       caption = "Source: Results/Task45_RRBS_ByState_SixNumber.csv (Task 45)")
b30b <- melt(cs[, .(state, outliers.coef2, outliers.coef3)], id.vars = "state")
b30b[, variable := factor(fifelse(variable == "outliers.coef2", "IQR2 outliers per site", "IQR3 outliers per site"),
                          levels = c("IQR2 outliers per site", "IQR3 outliers per site"))]
p30b <- ggplot(b30b, aes(state, value, fill = variable)) +
  geom_boxplot(position = position_dodge(width = 0.78), width = 0.66, colour = INK2,
               outlier.size = 0.7, outlier.colour = MUTED, linewidth = 0.35) +
  scale_fill_manual(values = c(BLUE, ORANGE)) +
  scale_y_continuous(breaks = 0:10) +
  labs(title = "(b) IQR outliers per site, by state (like Task 9)",
       subtitle = "Each site has 18 samples. Boxes: middle half of the sites.",
       x = "Methylation state", y = "Outliers at one site",
       caption = "Source: Results/Task39_RRBS_CompleteSiteStates.csv.gz (Task 39)")
save_fig("Fig30_RRBS_ByState.png",
         titled("RRBS chr22, complete sites by methylation state", p30a, p30b), 13, 5.8)

# Fig31 the two windows, site by site (like Fig6)
ymax31 <- max(sapply(win, function(r) max(rowSums(low2[r, ] | high2[r, ]))))
win_plots <- function(wname) {
  r <- win[[wname]]
  w <- data.table(i = seq_along(r), pos = pos[r], state = factor(cs$state[r], levels = states),
                  q1 = q[r, 1], q3 = q[r, 2], med = rowMedians(cm[r, ]),
                  low = rowSums(low2[r, ]), high = rowSums(high2[r, ]))
  present <- states[states %in% as.character(w$state)]
  brks <- round(seq(1, 100, length.out = 5))
  xs <- scale_x_continuous(breaks = brks, labels = comma(w$pos[brks]), expand = expansion(add = 4))
  wc <- wcmp[window == wname]
  lead <- c(index100 = "The first 100 complete sites", tight100 = "The tightest 100 complete sites")[[wname]]
  pa <- ggplot(w, aes(i)) +
    geom_segment(aes(xend = i, y = q1, yend = q3), colour = GRAY, linewidth = 1.3) +
    geom_point(aes(y = med, fill = state), shape = 21, colour = SURF, stroke = 0.6, size = 2.4) +
    scale_fill_manual(values = STATE_COL[present], breaks = present) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) + xs +
    labs(title = sprintf("%s: %s bases", lead, comma(wc$span_bp)),
         subtitle = sprintf("chr22:%s to %s. Grey bar: middle half of the 18 samples. Dot: median, by state.",
                            comma(wc$first_pos), comma(wc$last_pos)),
         x = NULL, y = "beta") +
    theme(axis.text.x = element_blank())
  pb <- ggplot(w, aes(i, 1, fill = state)) +
    geom_tile(height = 1, width = 0.95) +
    scale_fill_manual(values = STATE_COL, guide = "none") + xs +
    scale_y_continuous(breaks = 1, labels = "state") +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), plot.margin = margin(0, 16, 0, 10))
  wl <- melt(w[, .(i, low, high)], id.vars = "i")
  wl[, variable := factor(fifelse(variable == "high", "IQR2 outlier, too high", "IQR2 outlier, too low"),
                          levels = c("IQR2 outlier, too high", "IQR2 outlier, too low"))]
  pc <- ggplot(wl, aes(i, value, fill = variable)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = c(INK2, GRAY), drop = FALSE) + xs +
    scale_y_continuous(limits = c(0, ymax31), breaks = seq(0, ymax31, 2)) +
    labs(subtitle = "Samples that are IQR2 outliers at each site",
         x = "chr22 position (sites drawn equally spaced)", y = "Samples")
  gl <- lapply(list(pa, pb, pc), ggplotGrob)
  wmax <- do.call(unit.pmax, lapply(gl, function(g) g$widths))
  gl <- lapply(gl, function(g) { g$widths <- wmax; g })
  gl <- lapply(gl, function(g) { g$widths[length(g$widths)] <- unit(28, "pt"); g })
  arrangeGrob(grobs = gl, ncol = 1, heights = c(3.2, 0.55, 2.1))
}
g31 <- arrangeGrob(win_plots("index100"), win_plots("tight100"), ncol = 2,
  top = textGrob("RRBS chr22: two windows of 100 complete sites, as in Tasks 18 and 20",
                 x = unit(12, "pt"), hjust = 0, gp = gpar(fontsize = 14, fontface = "bold", col = INK)),
  bottom = textGrob(paste0("Sources: Data/NN.hg38.18P.forw.chr22.w.header.txt, ",
                           "Results/Task39_RRBS_CompleteSiteStates.csv.gz and Task45_RRBS_WindowComparison.csv (Task 45)"),
                    x = unit(12, "pt"), hjust = 0, gp = gpar(fontsize = 8.5, col = MUTED)))
save_fig("Fig31_RRBS_TwoWindows.png", g31, 16, 7.6)

# Fig32 runs: does the same sample stay an outlier at neighbouring sites?
r32 <- melt(rlt, id.vars = "run_length", variable.name = "series", value.name = "runs")
r32[, run_length := factor(run_length, levels = RUN_BINS)]
r32[, series := factor(fifelse(series == "runs_observed", "Real data", "After shuffling the samples at each site"),
                       levels = c("Real data", "After shuffling the samples at each site"))]
nb_all <- nbt[scope == "all complete sites"]
p32 <- ggplot(r32, aes(run_length, runs, colour = series, group = series)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(fill = series), shape = 21, colour = SURF, stroke = 0.7, size = 3) +
  geom_text(data = r32[series == "Real data"], aes(label = comma(runs), vjust = fifelse(run_length == "1", 2, -1)),
            colour = INK2, size = 3.2) +
  scale_colour_manual(values = c(BLUE, MUTED)) + scale_fill_manual(values = c(BLUE, MUTED)) +
  scale_y_log10(labels = comma, breaks = 10^(0:5), expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "The same sample tends to stay an outlier across neighbouring CG sites",
       subtitle = paste0("A run: consecutive complete sites, each within 1,000 bases of the next, where one sample is ",
                         "an IQR2 outlier in the same direction.\n", pct(nb_all$share_with_neighbour),
                         " of outlier cells have such a neighbour, against ",
                         pct(nb_all$share_with_neighbour_shuffled), " after shuffling (mean of ", N_SHUFFLE,
                         " shuffles). This shows that outliers cluster, not why."),
       x = "Run length (number of neighbouring sites)", y = "Number of runs (log scale)",
       caption = paste0("Sources: Results/Task45_RRBS_OutlierRunLengths.csv and ",
                        "Task45_RRBS_OutlierNeighbours.csv (Task 45)"))
save_fig("Fig32_RRBS_OutlierRuns.png", p32, 10, 5.8)

cat("\ndone\n")
