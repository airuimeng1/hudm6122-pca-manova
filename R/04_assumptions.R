## ============================================================
## 04_assumptions.R — MANOVA assumption diagnostics
## DV: PC1..PC4 scores (sign-aligned)
## IVs: region (4 levels) x period (2 levels)
## Check: multivariate normality, Box's M (cov homogeneity),
##        univariate Levene tests, multicollinearity, outliers
## ============================================================
source("R/00_setup.R")

dat <- read.csv(file.path(OUT_DIR, "panel_dedi.csv"), stringsAsFactors = FALSE)
dat$region <- factor(dat$region, levels = c("East","Central","West","Northeast"))
dat$period <- factor(dat$period, levels = c("Early(2013-2017)","Late(2018-2022)"))
dat$cell   <- interaction(dat$region, dat$period, drop = TRUE)
DV <- c("PC1","PC2","PC3","PC4")
Y  <- as.matrix(dat[, DV])

cat("Cell sample sizes (region x period):\n")
print(table(dat$region, dat$period))

## ---- 1. Multivariate normality (Mardia) ----
cat("\n## Multivariate normality — Mardia's test (overall, n = 310) ##\n")
mvn_overall <- MVN::mvn(as.data.frame(Y), mvn_test = "mardia",
                        descriptives = FALSE, tidy = TRUE)
print(mvn_overall$multivariate_normality)

cat("\nUnivariate Shapiro-Wilk by DV:\n")
sw_tab <- t(sapply(DV, function(v) {
  s <- shapiro.test(dat[[v]])
  c(W = round(s$statistic, 4), p = signif(s$p.value, 3))
}))
print(sw_tab)

## Per-cell Mardia (because MANOVA assumes MVN within each cell)
cat("\n## Multivariate normality within each region x period cell ##\n")
mvn_summary <- do.call(rbind, lapply(split(as.data.frame(Y), dat$cell), function(blk) {
  res <- try(MVN::mvn(blk, mvn_test = "mardia", descriptives = FALSE, tidy = TRUE),
             silent = TRUE)
  if (inherits(res, "try-error")) return(NULL)
  res$multivariate_normality
}))
mvn_summary$cell <- rep(levels(dat$cell), each = 2)
mvn_summary <- mvn_summary[, c("cell", "Test", "Statistic", "p.value", "MVN")]
print(mvn_summary, row.names = FALSE)
write.csv(mvn_summary, file.path(TBL_DIR, "04_mvn_per_cell.csv"), row.names = FALSE)

## ---- 2. Homogeneity of covariance matrices — Box's M ----
cat("\n## Box's M test of covariance homogeneity ##\n")
cat("Across 8 cells (region x period):\n")
boxM_cell <- heplots::boxM(Y, dat$cell)
print(boxM_cell)
cat("Across 4 regions (collapsing period):\n")
boxM_reg <- heplots::boxM(Y, dat$region)
print(boxM_reg)
cat("Across 2 periods (collapsing region):\n")
boxM_per <- heplots::boxM(Y, dat$period)
print(boxM_per)

## Box's M is notoriously sensitive; report at alpha = .001 (Tabachnick & Fidell).

## ---- 3. Univariate Levene's tests per DV (cell-level) ----
cat("\n## Univariate Levene's test per DV (across 8 cells) ##\n")
lev_tab <- do.call(rbind, lapply(DV, function(v) {
  fr <- as.formula(paste(v, "~ region * period"))
  l  <- car::leveneTest(fr, data = dat)
  data.frame(DV = v, F = round(l[1, "F value"], 3),
             df1 = l[1, "Df"], df2 = l[2, "Df"],
             p = signif(l[1, "Pr(>F)"], 3))
}))
print(lev_tab, row.names = FALSE)
write.csv(lev_tab, file.path(TBL_DIR, "04_levene.csv"), row.names = FALSE)

## ---- 4. Multicollinearity among DVs ----
cat("\n## Correlation among PC scores (should be ~ 0 by construction) ##\n")
print(round(cor(Y), 3))

## ---- 5. Multivariate outliers — Mahalanobis distance ----
md  <- mahalanobis(Y, colMeans(Y), cov(Y))
crit <- qchisq(0.999, df = ncol(Y))
n_out <- sum(md > crit)
cat(sprintf("\nMultivariate outliers (Mahalanobis^2 > chi-sq(.999, df=%d)=%.2f): %d obs\n",
            ncol(Y), crit, n_out))
if (n_out > 0) {
  print(dat[md > crit, c("year","prov_en","region","period","PC1","PC2","PC3","PC4")])
}

## ---- 6. Save assumption-check summary ----
saveRDS(list(mvn_overall = mvn_overall,
             boxM_cell = boxM_cell,
             boxM_reg = boxM_reg,
             boxM_per = boxM_per,
             levene = lev_tab,
             mahalanobis = md,
             dat = dat),
        file.path(OUT_DIR, "assumption_checks.rds"))

cat("\n[DONE 04_assumptions]\n")
