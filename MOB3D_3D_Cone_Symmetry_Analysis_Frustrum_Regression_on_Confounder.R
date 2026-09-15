################################################################################
## Sister-glomeruli regression, manuscript scenario: MitralMarker + RowOut, Rb = 200 um
##
## Self-contained: everything is read from ./data, everything is written to ./output.
##   Step 1  geometric centroid of each frustum
##   Step 2  pairwise centroid distances between all glomerulus pairs
##   Step 3  regression of pairwise correlation on sister status + distance + domain + class
##   Step 4  save coefficients and model summary
##   Step 5  boxplots with the fitted mean effects removed (3 panels, per Alex's request)
##
## Run:  module load r/4.5.1 ; Rscript run_pipeline.R
################################################################################

suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(grid)})
HERE <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) getwd())
if (is.na(HERE) || !nzchar(HERE)) HERE <- getwd()
DATA <- file.path(HERE, "data"); OUT <- file.path(HERE, "output")
dir.create(OUT, showWarnings = FALSE)

RT <- 50    # upper base radius (um), at the glomerulus
RB <- 200   # lower base radius (um), at the mitral end

################################################################################
## STEP 1 -- geometric centroid of each frustum
################################################################################
## frustum_top_bottom_center.txt is written by the frustum-generation script and holds,
## per glomerulus: G = glomerulus centroid, S = nearest mitral-layer spot (Leiden domain 6),
## B = S + 100 um along the G->S axis (the lower base centre).
tb <- read.table(file.path(DATA, "frustum_top_bottom_center.txt"))
stopifnot(all(c("glomerulus_id","g_x","g_y","g_z","b_x","b_y","b_z") %in% names(tb)))

## Centroid of a conical frustum, measured from the WIDE (mitral, B) end:
##   zbar / h = (Rb^2 + 2*Rb*Rt + 3*Rt^2) / (4 * (Rb^2 + Rb*Rt + Rt^2))
frac <- (RB^2 + 2*RB*RT + 3*RT^2) / (4 * (RB^2 + RB*RT + RT^2))     # = 0.32143 for 50/200
cat(sprintf("STEP 1  centroid sits %.5f of the height from the mitral end\n", frac))

tb <- tb %>% mutate(
  c_x = b_x + frac * (g_x - b_x),
  c_y = b_y + frac * (g_y - b_y),
  c_z = b_z + frac * (g_z - b_z))
rownames(tb) <- tb$glomerulus_id

################################################################################
## STEP 2 -- pairwise distance between frustum centroids
################################################################################
sg <- as.data.frame(readxl::read_xls(file.path(DATA, "glom_for_correlation_480.xls")))
colnames(sg) <- c("G1","G2")
sg <- sg[-104, ]; sg <- sg[-171, ]          # two duplicated-receptor pairs removed
gl428 <- unique(c(sg$G1, sg$G2))            # 428 glomeruli = 214 sister pairs
t4 <- tb[gl428, ]
cat(sprintf("STEP 2  %d glomeruli, %d pairs\n", nrow(t4), choose(nrow(t4), 2)))

pr <- t(combn(nrow(t4), 2))
A <- t4[pr[,1], ]; B <- t4[pr[,2], ]
dist_df <- data.frame(
  Glomerulus_1 = t4$glomerulus_id[pr[,1]],
  Glomerulus_2 = t4$glomerulus_id[pr[,2]],
  dist_centroid = sqrt((A$c_x-B$c_x)^2 + (A$c_y-B$c_y)^2 + (A$c_z-B$c_z)^2))
dist_df$key <- paste0(pmin(dist_df$Glomerulus_1, dist_df$Glomerulus_2), "|",
                      pmax(dist_df$Glomerulus_1, dist_df$Glomerulus_2))
saveRDS(dist_df, file.path(OUT, "frustum_pair_distances.rds"))

################################################################################
## STEP 3 -- assemble the model frame and fit
################################################################################
## The M-L domain and receptor-class labels are already stored per glomerulus in the
## correlation table (Glomerulus_*_ML, Glomerulus_*.Domain), so they are read from there.
## They are only canonicalised: the stored Domain_pair is order-dependent (classI-classII and
## classII-classI appear as separate levels), which would split one category into two dummies.
canon <- function(a, b) paste0(pmin(a, b), "-", pmax(a, b))
dist_map <- setNames(dist_df$dist_centroid, dist_df$key)

d <- readRDS(file.path(DATA, "MitralMarker_RowOut_Correlation_Results.rds")) %>%
  filter(!is.na(corr)) %>%
  mutate(key = paste0(pmin(Glomerulus_1, Glomerulus_2), "|",
                      pmax(Glomerulus_1, Glomerulus_2))) %>%
  distinct(key, .keep_all = TRUE) %>%                       # unordered pairs only
  mutate(dist_centroid = dist_map[key],
         ML_pair    = canon(Glomerulus_1_ML,     Glomerulus_2_ML),
         class_pair = canon(Glomerulus_1.Domain, Glomerulus_2.Domain)) %>%
  filter(!is.na(dist_centroid), !grepl("NA", class_pair)) %>%
  mutate(dist_c     = dist_centroid - mean(dist_centroid),  # mean-centered distance
         ML_pair    = relevel(factor(ML_pair),    ref = "lateral-medial"),
         class_pair = relevel(factor(class_pair), ref = "classI-classI"))
cat(sprintf("STEP 3  model frame: %d pairs (%d sister)\n", nrow(d), sum(d$is_sister)))

m <- lm(corr ~ is_sister + dist_c + ML_pair + class_pair, data = d)

################################################################################
## STEP 4 -- save the outcome
################################################################################
cf <- as.data.frame(summary(m)$coefficients)
names(cf) <- c("estimate","std_error","t_value","p_value")
cf$term <- rownames(cf); cf <- cf[, c("term","estimate","std_error","t_value","p_value")]
write.csv(cf, file.path(OUT, "regression_coefficients.csv"), row.names = FALSE)

## partial R^2 per term: (SSE_reduced - SSE_full) / SSE_reduced, factors dropped as a block
terms <- c("is_sister","dist_c","ML_pair","class_pair")
sse   <- function(f) sum(resid(lm(f, data = d))^2)
s_full <- sse(corr ~ is_sister + dist_c + ML_pair + class_pair)
pr2 <- do.call(rbind, lapply(terms, function(t) {
  s_red <- sse(as.formula(paste("corr ~", paste(setdiff(terms, t), collapse = " + "))))
  data.frame(term = t, partial_r2 = (s_red - s_full) / s_red)
}))
pr2$model_r2 <- summary(m)$r.squared
write.csv(pr2, file.path(OUT, "partial_r2.csv"), row.names = FALSE)
saveRDS(list(model = m, data = d), file.path(OUT, "fitted_model.rds"))
cat("STEP 4  wrote regression_coefficients.csv, partial_r2.csv, fitted_model.rds\n")
## rendered coefficient table for the response letter
suppressMessages(library(gridExtra))
pretty_p <- function(p) ifelse(p < 2e-16, "< 2e-16",
                        ifelse(p < 0.001, sprintf("%.1e", p), sprintf("%.3f", p)))
lab <- c("(Intercept)"               = "Intercept",
         "is_sisterTRUE"             = "Sister pair",
         "dist_c"                    = "Frustum distance (per um)",
         "ML_pairlateral-lateral"    = "Domain: lateral-lateral",
         "ML_pairmedial-medial"      = "Domain: medial-medial",
         "class_pairclassI-classII"  = "Class: I-II",
         "class_pairclassII-classII" = "Class: II-II")
tab <- data.frame(Term = lab[cf$term],
                  Estimate = ifelse(abs(cf$estimate) < 1e-3,
                                    sprintf("%.2e", cf$estimate), sprintf("%.4f", cf$estimate)),
                  `Std. error` = ifelse(cf$std_error < 1e-3,
                                    sprintf("%.2e", cf$std_error), sprintf("%.4f", cf$std_error)),
                  `t` = sprintf("%.1f", cf$t_value),
                  `p value` = pretty_p(cf$p_value),
                  check.names = FALSE)
## partial R2 is a per-term quantity: shown once per block, blank on the second dummy
blk <- c(NA, "is_sister", "dist_c", "ML_pair", NA, "class_pair", NA)
tab$`Partial R2` <- ifelse(is.na(blk), "",
                           sprintf("%.4f", pr2$partial_r2[match(blk, pr2$term)]))
write.csv(tab, file.path(OUT, "regression_table.csv"), row.names = FALSE)

th <- ttheme_minimal(base_size = 12,
        core    = list(fg_params = list(hjust = 0, x = 0.03,
                       fontface = c(rep(1, 1), 2, rep(1, 5)))),
        colhead = list(fg_params = list(fontface = 2)))
g <- tableGrob(tab, rows = NULL, theme = th)
cap <- textGrob(sprintf("MitralMarker + RowOut, Rb = 200 um   |   n = %s pairs (%d sister)   |   model R2 = %.3f",
                        format(nrow(d), big.mark = ","), sum(d$is_sister), summary(m)$r.squared),
                gp = gpar(fontsize = 10), hjust = 0, x = 0.02)
pdf(file.path(OUT, "regression_table.pdf"), width = 9, height = 3.2)
grid.arrange(g, cap, ncol = 1, heights = c(4, 0.5))
dev.off()
cat("STEP 4  wrote regression_table.pdf / .csv\n")

print(cf, row.names = FALSE, digits = 4)
print(pr2, row.names = FALSE, digits = 3)

################################################################################
## STEP 5 -- boxplots with the fitted mean effects removed
################################################################################
## Contributions are centered, so every panel stays on the correlation scale.
b_dist <- coef(m)["dist_c"]
ml_lev <- c("lateral-medial" = 0,
            "lateral-lateral" = unname(coef(m)["ML_pairlateral-lateral"]),
            "medial-medial"   = unname(coef(m)["ML_pairmedial-medial"]))
mlc <- ml_lev[as.character(d$ML_pair)]; mlc <- mlc - mean(mlc)
dc  <- b_dist * d$dist_c

PAN <- c("Unadjusted Correlation", "Distance effect removed", "Distance + M-L effects removed")
pd <- data.frame(group = factor(ifelse(d$is_sister, "sister", "non-sister"),
                                levels = c("non-sister","sister")),
                 `Unadjusted Correlation`         = d$corr,
                 `Distance effect removed`        = d$corr - dc,
                 `Distance + M-L effects removed` = d$corr - dc - mlc,
                 check.names = FALSE) %>%
  pivot_longer(all_of(PAN), names_to = "panel", values_to = "value") %>%
  mutate(panel = factor(panel, levels = PAN))

wt <- pd %>% group_by(panel) %>%
  summarise(p = suppressWarnings(wilcox.test(value ~ group)$p.value),
            n_sister = sum(group == "sister"), n_nonsister = sum(group == "non-sister"),
            .groups = "drop") %>%
  mutate(label = ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p)))
write.csv(wt, file.path(OUT, "wilcoxon_tests.csv"), row.names = FALSE)
print(as.data.frame(wt[, c("panel","p","n_sister","n_nonsister")]), row.names = FALSE, digits = 3)

## Fig. 4c styling: letter-value boxes, no outline, solid median, gray / sister red.
## Level i spans the central 1 - 2^-i of the data; the outermost box is extended to the
## true min/max so the tails are complete. Widths taper linearly -> thin tails.
lv_boxes <- function(v, K = 6, wmax = 0.72) {
  v <- sort(v); n <- length(v)
  dd <- (n + 1)/2; depths <- numeric(0)
  for (j in seq_len(K)) { depths <- c(depths, dd); dd <- (floor(dd) + 1)/2 }
  lo <- sapply(depths, function(z) (v[max(1, floor(z))] + v[max(1, ceiling(z))]) / 2)
  hi <- sapply(depths, function(z) (v[n + 1 - max(1, floor(z))] + v[n + 1 - max(1, ceiling(z))]) / 2)
  lo[K] <- v[1]; hi[K] <- v[n]
  data.frame(j = seq_len(K), ylo = lo, yhi = hi, w = wmax * (K - seq_len(K) + 1)/K)[-1, ]
}
lv <- pd %>% group_by(panel, group) %>% reframe(lv_boxes(value)) %>%
  mutate(x = as.numeric(group), xmin = x - w/2, xmax = x + w/2) %>% arrange(desc(j))
med <- pd %>% group_by(panel, group) %>% summarise(m = median(value), .groups = "drop") %>%
  mutate(x = as.numeric(group))
rng <- pd %>% group_by(panel) %>% summarise(lo = min(value), hi = max(value), .groups = "drop")
ylo_ax <- min(rng$lo) - 0.03 * (max(rng$hi) - min(rng$lo))
wt$y   <- max(rng$hi) + 0.07 * (max(rng$hi) - min(rng$lo))
yhi_ax <- max(wt$y)   + 0.06 * (max(rng$hi) - min(rng$lo))

p <- ggplot() +
  geom_rect(data = lv, aes(xmin = xmin, xmax = xmax, ymin = ylo, ymax = yhi, fill = group),
            colour = NA) +
  geom_segment(data = med, aes(x = x - 0.36, xend = x + 0.36, y = m, yend = m),
               colour = "black", linewidth = 0.7) +
  geom_text(data = wt, aes(x = 1.5, y = y, label = label), size = 10) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = c("non-sister" = "#A9AAAC", "sister" = "#BE202F")) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Non-sisters","Sisters"), limits = c(0.45, 2.55)) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.25)) +
  coord_cartesian(ylim = c(ylo_ax, yhi_ax)) +
  labs(x = NULL, y = "Correlation") +
  theme_classic(base_size = 25) +
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text  = element_text(size = 30),
        axis.text   = element_text(size = 25, colour = "black"),
        axis.title  = element_text(size = 25),
        axis.line   = element_line(colour = "black", linewidth = 0.5),
        axis.ticks  = element_line(colour = "black"))
ggsave(file.path(OUT, "boxplot_sister_vs_nonsister.pdf"), p, width = 20, height = 7)
cat("STEP 5  wrote boxplot_sister_vs_nonsister.pdf\n")
