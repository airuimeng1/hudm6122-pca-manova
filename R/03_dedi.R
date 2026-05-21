## ============================================================
## 03_dedi.R — Construct DEDI composite index from retained PCs
## ============================================================
source("R/00_setup.R")

obj <- readRDS(file.path(OUT_DIR, "pca_objects.rds"))
df  <- obj$df
pca <- obj$pca
k   <- obj$n_keep
eig <- pca$sdev^2

## ---- 1. Variance-share weights for retained PCs ----
prop  <- eig[1:k] / sum(eig)        # share of total variance
wts   <- prop / sum(prop)           # renormalised within retained
cat("Retained k =", k, "PCs\n")
cat("Variance shares:\n"); print(round(prop, 4))
cat("Renormalised weights (sum to 1):\n"); print(round(wts, 4))

## NOTE: PC2 in the unrotated solution has *negative* loadings on intensity vars
## (x10, x12, x16, x18), so a province with HIGH per-capita penetration ends
## up with a LOW PC2 score. To make all PCs point in the "more digital = higher"
## direction, we flip any PC whose total signed correlation with the raw scaled
## indicators is negative.
Z <- obj$Z
sign_flip <- numeric(k)
for (j in 1:k) {
  s <- sign(sum(cor(pca$x[, j], Z)))
  sign_flip[j] <- ifelse(s == 0, 1, s)
}
cat("Sign flips (1 = keep, -1 = reverse):", sign_flip, "\n")

scores_signed <- sweep(pca$x[, 1:k], 2, sign_flip, "*")
colnames(scores_signed) <- paste0("PC", 1:k)

## ---- 2. DEDI = weighted sum of (sign-aligned) PC scores ----
dedi_raw <- as.numeric(scores_signed %*% wts)

## min-max rescale to 0-100 for interpretability (only DEDI, not PCs)
rescale01 <- function(v) (v - min(v)) / (max(v) - min(v))
dedi_100 <- rescale01(dedi_raw) * 100

panel_out <- df[, c("year","id","province","prov_en","region","period")]
panel_out <- cbind(panel_out, scores_signed, DEDI_raw = dedi_raw, DEDI = dedi_100)
write.csv(panel_out, file.path(OUT_DIR, "panel_dedi.csv"), row.names = FALSE)

## ---- 3. National DEDI trajectory by region (yearly mean) ----
region_year <- panel_out %>%
  group_by(region, year) %>%
  summarise(DEDI = mean(DEDI), .groups = "drop")

p1 <- ggplot(region_year, aes(year, DEDI, colour = region)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2) +
  scale_x_continuous(breaks = 2013:2022) +
  scale_colour_manual(values = c(East="#d62728", Central="#1f77b4",
                                 West="#2ca02c", Northeast="#9467bd")) +
  labs(title = "Regional Digital Economy Development Index (DEDI), 2013-2022",
       subtitle = "Yearly mean of provinces within region (rescaled 0-100)",
       y = "DEDI (0-100)", x = NULL, colour = "Region") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(FIG_DIR, "03_region_trend.png"), p1, width = 8, height = 5, dpi = 300)

## ---- 4. Province-level heatmap (province x year) ----
prov_order <- panel_out %>%
  filter(year == 2022) %>% arrange(desc(DEDI)) %>% pull(prov_en)

heat_df <- panel_out
heat_df$prov_en <- factor(heat_df$prov_en, levels = rev(prov_order))

p2 <- ggplot(heat_df, aes(year, prov_en, fill = DEDI)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_x_continuous(breaks = 2013:2022, expand = c(0, 0)) +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(title = "DEDI by province and year",
       subtitle = "Provinces ordered by 2022 DEDI",
       y = NULL, x = NULL, fill = "DEDI") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0))
ggsave(file.path(FIG_DIR, "03_dedi_heatmap.png"), p2, width = 9, height = 9, dpi = 300)

## ---- 5. 2019 cross-section ranking ----
top2019 <- panel_out %>% filter(year == 2019) %>%
  arrange(desc(DEDI)) %>%
  select(province = prov_en, region, DEDI, PC1, PC2, PC3, PC4)
top2019[, 3:7] <- lapply(top2019[, 3:7], function(v) round(v, 2))
write.csv(top2019, file.path(TBL_DIR, "03_dedi_2019_ranking.csv"), row.names = FALSE)
cat("\n## 2019 DEDI ranking (top 10 / bottom 5) ##\n")
print(rbind(head(top2019, 10), tail(top2019, 5)), row.names = FALSE)

## ---- 6. Region x period mean cell table (preview for Step 5) ----
cell_means <- panel_out %>%
  group_by(region, period) %>%
  summarise(n = n(), DEDI = round(mean(DEDI), 2),
            PC1 = round(mean(PC1), 2), PC2 = round(mean(PC2), 2),
            PC3 = round(mean(PC3), 2), PC4 = round(mean(PC4), 2),
            .groups = "drop")
cat("\n## Region x Period cell means (preview) ##\n")
print(cell_means)
write.csv(cell_means, file.path(TBL_DIR, "03_cell_means.csv"), row.names = FALSE)

## ---- 7. Convergence/divergence diagnostic: SD across provinces over time ----
sd_year <- panel_out %>%
  group_by(year) %>%
  summarise(sd_DEDI = sd(DEDI),
            cv_DEDI = sd(DEDI) / mean(DEDI),
            .groups = "drop")
cat("\n## Inter-provincial dispersion of DEDI ##\n")
print(sd_year)
write.csv(sd_year, file.path(TBL_DIR, "03_dispersion.csv"), row.names = FALSE)

p3 <- ggplot(sd_year, aes(year, sd_DEDI)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_x_continuous(breaks = 2013:2022) +
  labs(title = "Inter-provincial dispersion of DEDI",
       subtitle = "Standard deviation across 31 provinces, by year",
       y = "SD of DEDI", x = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG_DIR, "03_dispersion.png"), p3, width = 7, height = 4.5, dpi = 300)

cat("\n[DONE 03_dedi]\n")
