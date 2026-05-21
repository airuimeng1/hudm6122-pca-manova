## ============================================================
## 06_sensitivity.R — Sensitivity analyses
## A. Single-year MANOVA (2018, 2019, 2020) using first 3 PCs
##    (the cumulative >= 80% set) — addresses N=3 Northeast issue
## B. Robustness of PCA loading structure across years
##    (do PCs from full panel match year-specific PCs?)
## C. Refit factorial MANOVA dropping Tibet 2015-2016 outliers
## ============================================================
source("R/00_setup.R")

dat <- read.csv(file.path(OUT_DIR, "panel_dedi.csv"), stringsAsFactors = FALSE)
dat$region <- factor(dat$region, levels = c("East","Central","West","Northeast"))
dat$period <- factor(dat$period, levels = c("Early(2013-2017)","Late(2018-2022)"))

##  =========================================================
##  A. Single-year MANOVA (one-way, factor = region only)
##  =========================================================
cat("\n## A. Single-year one-way MANOVA on REGION ##\n")
DV3 <- c("PC1","PC2","PC3")  # use 3 PCs (cum var >= 80%) — works with n=3 NE cell

single_year_results <- list()
for (yr in c(2018, 2019, 2020)) {
  cat(sprintf("\n=== Year %d (n = 31) ===\n", yr))
  d_yr <- dat[dat$year == yr, ]
  Y3 <- as.matrix(d_yr[, DV3])
  fit_yr <- lm(Y3 ~ region, data = d_yr)
  m <- car::Manova(fit_yr, test.statistic = "Pillai", type = "II")
  print(m)
  ## Box's M (warn about singular cov in Northeast n=3, k=3)
  bm <- try(heplots::boxM(Y3, d_yr$region), silent = TRUE)
  if (!inherits(bm, "try-error")) print(bm) else
    cat("Box's M failed (likely singular covariance in Northeast; n=3, k=3).\n")

  single_year_results[[as.character(yr)]] <-
    list(manova = m, boxM = if (!inherits(bm, "try-error")) bm else NULL)
}

## Summary table — capture via printed text and parse numbers from the
## first row (region) of the test table
extract_pillai <- function(manova_obj) {
  txt <- capture.output(print(manova_obj))
  ## find row that starts with "region"
  reg_line <- txt[grepl("^region", txt)]
  toks <- strsplit(trimws(reg_line), "\\s+")[[1]]
  ## tokens: "region" df  test_stat  approx_F  num_Df  den_Df  p
  data.frame(
    Pillai   = as.numeric(toks[3]),
    F_approx = as.numeric(toks[4]),
    df1      = as.numeric(toks[5]),
    df2      = as.numeric(toks[6]),
    p        = as.numeric(sub("\\*+$", "", toks[7]))
  )
}

sy_tab <- do.call(rbind, lapply(names(single_year_results), function(yr) {
  e <- extract_pillai(single_year_results[[yr]]$manova)
  cbind(year = yr, e)
}))
sy_tab$Pillai   <- round(sy_tab$Pillai, 4)
sy_tab$F_approx <- round(sy_tab$F_approx, 3)
sy_tab$p        <- signif(sy_tab$p, 3)
cat("\n## Summary: single-year MANOVA on region (Pillai) ##\n")
print(sy_tab, row.names = FALSE)
write.csv(sy_tab, file.path(TBL_DIR, "06A_singleyear_manova.csv"), row.names = FALSE)

##  =========================================================
##  B. Year-specific PCA loading stability
##  =========================================================
obj <- readRDS(file.path(OUT_DIR, "pca_objects.rds"))
xs <- paste0("x", 1:22)
heavy <- obj$heavy

cat("\n## B. Year-specific PCA: stability of top-2 PC loadings ##\n")

cor_load <- function(L_full, L_yr) {
  ## match each PC of yr to PC of full by max |corr|
  k <- ncol(L_full)
  out <- numeric(k)
  for (j in 1:k) {
    cors <- sapply(1:k, function(i) cor(L_full[, j], L_yr[, i]))
    out[j] <- max(abs(cors))
  }
  out
}

L_full <- obj$pca$rotation[, 1:obj$n_keep] %*% diag(obj$pca$sdev[1:obj$n_keep])

stab_tab <- do.call(rbind, lapply(c(2013, 2018, 2022), function(yr) {
  d_yr <- obj$df[obj$df$year == yr, ]
  X_yr <- as.data.frame(d_yr[, xs])
  X_yr[heavy] <- lapply(X_yr[heavy], function(v) log1p(v - min(v) + 1e-6))
  Z_yr <- scale(X_yr)
  pca_yr <- prcomp(Z_yr, center = FALSE, scale. = FALSE)
  L_yr <- pca_yr$rotation[, 1:obj$n_keep] %*% diag(pca_yr$sdev[1:obj$n_keep])
  data.frame(year = yr,
             PC = paste0("PC", 1:obj$n_keep),
             best_abs_corr_with_full = round(cor_load(L_full, L_yr), 3))
}))
print(stab_tab, row.names = FALSE)
write.csv(stab_tab, file.path(TBL_DIR, "06B_pca_stability.csv"), row.names = FALSE)

##  =========================================================
##  C. Refit factorial MANOVA dropping Tibet 2015-2016
##  =========================================================
cat("\n## C. Factorial MANOVA without Tibet 2015-2016 outliers ##\n")
DV4 <- c("PC1","PC2","PC3","PC4")
## key off Tibet's GB id (540000) to avoid prov_en encoding/lookup issues
drop_idx <- which(dat$id == 540000 & dat$year %in% c(2015, 2016))
d_no <- dat[-drop_idx, ]
Y4 <- as.matrix(d_no[, DV4])
fit_no <- lm(Y4 ~ region * period, data = d_no)
m_no <- car::Manova(fit_no, test.statistic = "Pillai", type = "II")
print(m_no)
cat(sprintf("\nDropped %d outlier obs (Tibet 2015-2016). New n = %d\n",
            length(drop_idx), nrow(d_no)))

cat("\n[DONE 06_sensitivity]\n")
