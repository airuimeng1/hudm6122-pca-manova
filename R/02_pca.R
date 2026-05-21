## ============================================================
## 02_pca.R — PCA on 22 indicators (panel, n = 310)
## - log-transform heavily right-skewed indicators (|skew| > 1.5)
##   so that the PCA correlation structure is not dominated by outliers
## - z-score standardize
## - KMO + Bartlett's test of sphericity
## - PCA, scree plot, loadings, retain PCs by Kaiser + Cumulative >= 80%
## - export PC scores to be reused as MANOVA DVs in Step 5
## ============================================================
source("R/00_setup.R")

df <- load_panel()
xs <- paste0("x", 1:22)
X  <- as.data.frame(df[, xs])

## ---- 1. log1p transform for heavy-skew variables ----
skew0 <- sapply(X, psych::skew)
heavy <- names(skew0)[abs(skew0) > 1.5]
cat("Variables log1p-transformed (|skew|>1.5):\n"); print(heavy)

X_t <- X
X_t[heavy] <- lapply(X_t[heavy], function(v) log1p(v - min(v) + 1e-6))

skew1 <- sapply(X_t, psych::skew)
cat("\nSkewness before / after transform:\n")
print(round(cbind(before = skew0, after = skew1), 2))

## ---- 2. z-score standardize ----
Z <- scale(X_t)  # column-wise z-score
attr(Z, "scaled:center") <- NULL  # tidy

## ---- 3. KMO + Bartlett's test of sphericity ----
R <- cor(Z)
kmo <- psych::KMO(R)
bart <- psych::cortest.bartlett(R, n = nrow(Z))

cat("\n## KMO measure of sampling adequacy ##\n")
cat("Overall MSA =", round(kmo$MSA, 3),
    " (>=.6 acceptable, >=.8 meritorious)\n")
cat("\nMSA per indicator:\n"); print(round(kmo$MSAi, 3))

cat("\n## Bartlett's test of sphericity ##\n")
cat("Chi-sq =", round(bart$chisq, 2),
    " df =", bart$df,
    " p =", format.pval(bart$p.value), "\n")

## ---- 4. PCA ----
pca <- prcomp(Z, center = FALSE, scale. = FALSE)  # already standardized
eig  <- pca$sdev^2
prop <- eig / sum(eig)
cum  <- cumsum(prop)

eig_tab <- data.frame(
  PC = seq_along(eig),
  Eigenvalue = round(eig, 4),
  PropVar = round(prop, 4),
  CumVar = round(cum, 4)
)
cat("\n## Eigenvalues / variance explained ##\n")
print(eig_tab[1:10, ])
write.csv(eig_tab, file.path(TBL_DIR, "02_pca_eigen.csv"), row.names = FALSE)

## Retention rule: Kaiser (eig >= 1) + cumulative >= 80%
n_kaiser <- sum(eig >= 1)
n_cum80  <- which(cum >= 0.80)[1]
n_keep   <- max(n_kaiser, n_cum80)
cat(sprintf("\nKaiser (eig>=1): %d  |  Cumulative>=80%%: %d  |  Retained k = %d\n",
            n_kaiser, n_cum80, n_keep))

## ---- 5. Loadings (rotated for interpretability via Varimax on retained PCs) ----
loadings_unrot <- pca$rotation[, 1:n_keep] %*% diag(pca$sdev[1:n_keep])
colnames(loadings_unrot) <- paste0("PC", 1:n_keep)
rownames(loadings_unrot) <- xs

vmx <- varimax(loadings_unrot)
loadings_rot <- unclass(vmx$loadings)
colnames(loadings_rot) <- paste0("RC", 1:n_keep)

load_tab <- data.frame(
  variable  = xs,
  label     = unname(IND_LABELS[xs]),
  dimension = unname(IND_DIM[xs]),
  round(loadings_unrot, 3),
  round(loadings_rot, 3)
)
cat("\n## Loadings (unrotated PC1..PCk and Varimax-rotated RC1..RCk) ##\n")
print(load_tab, row.names = FALSE)
write.csv(load_tab, file.path(TBL_DIR, "02_pca_loadings.csv"), row.names = FALSE)

## ---- 6. PC scores ----
scores <- as.data.frame(pca$x[, 1:n_keep])
colnames(scores) <- paste0("PC", 1:n_keep)
out <- cbind(df[, c("year","id","province","prov_en","region","period")], scores)
write.csv(out, file.path(OUT_DIR, "pc_scores.csv"), row.names = FALSE)
saveRDS(list(pca = pca, n_keep = n_keep, heavy = heavy,
             eig_tab = eig_tab, load_tab = load_tab,
             vmx = vmx, scores = out, X = X, X_t = X_t, Z = Z, df = df),
        file.path(OUT_DIR, "pca_objects.rds"))

## ---- 7. Plots: scree, cumulative, loading heatmap ----
p_scree <- ggplot(eig_tab[1:10, ], aes(PC, Eigenvalue)) +
  geom_col(fill = "steelblue") + geom_line() + geom_point() +
  geom_hline(yintercept = 1, linetype = 2, colour = "red") +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Scree plot (first 10 PCs)",
       subtitle = "Red dashed line: Kaiser eigenvalue = 1") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG_DIR, "02_scree.png"), p_scree, width = 7, height = 4.5, dpi = 300)

p_cum <- ggplot(eig_tab[1:10, ], aes(PC, CumVar)) +
  geom_col(fill = "darkorange") +
  geom_hline(yintercept = 0.80, linetype = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Cumulative variance explained",
       subtitle = "Dashed: 80% target",
       y = "Cumulative variance") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG_DIR, "02_cumvar.png"), p_cum, width = 7, height = 4.5, dpi = 300)

heat_df <- as.data.frame(loadings_rot)
heat_df$variable <- factor(rownames(heat_df), levels = rev(xs))
heat_df$label    <- factor(IND_LABELS[as.character(heat_df$variable)],
                           levels = rev(IND_LABELS[xs]))
heat_long <- pivot_longer(heat_df, cols = starts_with("RC"),
                          names_to = "PC", values_to = "loading")
p_load <- ggplot(heat_long, aes(PC, label, fill = loading)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", loading)), size = 3) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = sprintf("Varimax-rotated loadings (k = %d retained PCs)", n_keep),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
ggsave(file.path(FIG_DIR, "02_loadings_heatmap.png"), p_load,
       width = 8, height = 8, dpi = 300)

cat("\n[DONE 02_pca]  Retained k =", n_keep,
    "PCs explaining", round(cum[n_keep] * 100, 1), "% variance.\n")
