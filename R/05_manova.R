## ============================================================
## 05_manova.R — Factorial MANOVA: region (4) x period (2)
## DV: PC1..PC4 (orthogonal by construction)
## Primary statistic: Pillai's trace (robust to Box's M violation)
## Follow-up: univariate ANOVAs (Welch where Levene fails) + Games-Howell post hoc
## ============================================================
source("R/00_setup.R")

dat <- read.csv(file.path(OUT_DIR, "panel_dedi.csv"), stringsAsFactors = FALSE)
dat$region <- factor(dat$region, levels = c("East","Central","West","Northeast"))
dat$period <- factor(dat$period, levels = c("Early(2013-2017)","Late(2018-2022)"))
DV <- c("PC1","PC2","PC3","PC4")
Y  <- as.matrix(dat[, DV])

## Use Type II SS (orthogonal interpretation; default car::Manova)
## with balanced/near-balanced design. heplots::etasq() gives partial eta^2.

## ---- 1. Factorial MANOVA ----
fit <- lm(Y ~ region * period, data = dat)

cat("\n## Factorial MANOVA: region * period — All 4 statistics ##\n")
mvn_pillai  <- car::Manova(fit, test.statistic = "Pillai",            type = "II")
mvn_wilks   <- car::Manova(fit, test.statistic = "Wilks",             type = "II")
mvn_hl      <- car::Manova(fit, test.statistic = "Hotelling-Lawley",  type = "II")
mvn_roy     <- car::Manova(fit, test.statistic = "Roy",               type = "II")

print(mvn_pillai)
cat("\n--- Wilks ---\n");            print(mvn_wilks)
cat("\n--- Hotelling-Lawley ---\n"); print(mvn_hl)
cat("\n--- Roy ---\n");              print(mvn_roy)

## ---- 2. Multivariate effect sizes (partial eta-squared) ----
cat("\n## Multivariate effect sizes (partial eta^2 via Pillai) ##\n")
es <- heplots::etasq(fit, test.statistic = "Pillai", anova = TRUE)
print(es)
write.csv(as.data.frame(unclass(es)),
          file.path(TBL_DIR, "05_manova_effectsizes.csv"))

## Save MANOVA tables
manova_summary <- function(m, name) {
  s <- summary(m)
  s$multivariate.tests
}

## ---- 3. Follow-up univariate ANOVAs per DV ----
cat("\n## Follow-up: univariate ANOVAs (Type II) per DV ##\n")
uni_list <- lapply(DV, function(v) {
  fr <- as.formula(paste(v, "~ region * period"))
  aov_fit <- lm(fr, data = dat)
  aov_tab <- car::Anova(aov_fit, type = "II")
  list(dv = v, table = aov_tab, fit = aov_fit)
})
for (u in uni_list) {
  cat("\n### DV =", u$dv, "###\n"); print(u$table)
}

## Compile univariate ANOVA into one table
uni_compiled <- do.call(rbind, lapply(uni_list, function(u) {
  tab <- as.data.frame(u$table)
  tab$Term <- rownames(tab)
  tab$DV <- u$dv
  tab
}))
uni_compiled <- uni_compiled[, c("DV","Term","Sum Sq","Df","F value","Pr(>F)")]
names(uni_compiled) <- c("DV","Term","SS","df","F","p")
uni_compiled$F <- round(uni_compiled$F, 3)
uni_compiled$p <- signif(uni_compiled$p, 3)
uni_compiled$SS <- round(uni_compiled$SS, 2)
print(uni_compiled, row.names = FALSE)
write.csv(uni_compiled, file.path(TBL_DIR, "05_univariate_anova.csv"), row.names = FALSE)

## ---- 4. Welch ANOVA + Games-Howell post hoc on REGION (DV by DV) ----
##      We use Games-Howell because Levene flagged unequal variances on PC1, PC2, PC4.
cat("\n## Region post hoc: Games-Howell (DV by DV, collapsing period) ##\n")
gh_all <- list()
for (v in DV) {
  cat("\n### DV =", v, "###\n")
  fr <- as.formula(paste(v, "~ region"))
  ## Welch's ANOVA
  w <- oneway.test(fr, data = dat, var.equal = FALSE)
  cat(sprintf("Welch's F(%g, %g) = %.3f, p = %s\n",
              w$parameter[1], w$parameter[2], w$statistic, format.pval(w$p.value, 3)))
  gh <- rstatix::games_howell_test(dat, formula = fr)
  gh_all[[v]] <- gh
  print(gh)
}
gh_compiled <- do.call(rbind, lapply(names(gh_all), function(v) {
  d <- as.data.frame(gh_all[[v]]); d$DV <- v; d
}))
write.csv(gh_compiled, file.path(TBL_DIR, "05_games_howell_region.csv"), row.names = FALSE)

## ---- 5. Period contrast per DV (just t-test since only 2 levels) ----
cat("\n## Period contrast: Welch t-test (Late - Early) per DV ##\n")
period_tab <- do.call(rbind, lapply(DV, function(v) {
  tt <- t.test(dat[[v]] ~ dat$period, var.equal = FALSE)
  data.frame(DV = v,
             mean_Early = round(tt$estimate[1], 3),
             mean_Late  = round(tt$estimate[2], 3),
             diff = round(diff(tt$estimate), 3),
             t = round(tt$statistic, 3),
             df = round(tt$parameter, 1),
             p = signif(tt$p.value, 3),
             ci_low = round(tt$conf.int[1], 3),
             ci_high = round(tt$conf.int[2], 3))
}))
print(period_tab, row.names = FALSE)
write.csv(period_tab, file.path(TBL_DIR, "05_period_ttest.csv"), row.names = FALSE)

## ---- 6. Visualisations ----
library(scales)
plot_df <- dat %>%
  pivot_longer(all_of(DV), names_to = "PC", values_to = "score")
pc_labels <- c(PC1 = "PC1: Industry & Innovation Scale",
               PC2 = "PC2: Per-capita Penetration",
               PC3 = "PC3: Infrastructure Connectivity",
               PC4 = "PC4: Enterprise Digitalization")
plot_df$PC <- factor(plot_df$PC, levels = names(pc_labels), labels = pc_labels)

p_box <- ggplot(plot_df, aes(region, score, fill = period)) +
  geom_boxplot(outlier.size = 0.7, position = position_dodge(0.8)) +
  facet_wrap(~ PC, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Early(2013-2017)" = "#a6cee3",
                               "Late(2018-2022)"  = "#1f78b4")) +
  labs(title = "PC scores by region and period",
       x = NULL, y = "PC score", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"))
ggsave(file.path(FIG_DIR, "05_boxplots.png"), p_box, width = 10, height = 7, dpi = 300)

## Interaction plot (cell means)
cell_summary <- dat %>%
  pivot_longer(all_of(DV), names_to = "PC", values_to = "score") %>%
  group_by(region, period, PC) %>%
  summarise(mean = mean(score), se = sd(score)/sqrt(n()), .groups = "drop")
cell_summary$PC <- factor(cell_summary$PC, levels = names(pc_labels), labels = pc_labels)

p_int <- ggplot(cell_summary, aes(period, mean, colour = region, group = region)) +
  geom_point(size = 2.5) + geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  facet_wrap(~ PC, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(East="#d62728", Central="#1f77b4",
                                 West="#2ca02c", Northeast="#9467bd")) +
  labs(title = "Region x Period interaction on PC scores",
       subtitle = "Mean +/- 1 SE; non-parallel lines indicate interaction",
       y = "Mean PC score", x = NULL, colour = "Region") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"))
ggsave(file.path(FIG_DIR, "05_interaction.png"), p_int, width = 10, height = 7, dpi = 300)

## Save all results
saveRDS(list(fit = fit, pillai = mvn_pillai, wilks = mvn_wilks,
             hl = mvn_hl, roy = mvn_roy, etasq = es,
             univariate = uni_compiled, games_howell = gh_compiled,
             period = period_tab),
        file.path(OUT_DIR, "manova_results.rds"))

cat("\n[DONE 05_manova]\n")
