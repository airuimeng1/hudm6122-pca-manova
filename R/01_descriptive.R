## ============================================================
## 01_descriptive.R  —  Descriptive statistics + missing audit
## ============================================================
source("R/00_setup.R")

df <- load_panel()
xs <- paste0("x", 1:22)

cat("\n## Panel structure ##\n")
cat("Years:", paste(sort(unique(df$year)), collapse=", "), "\n")
cat("# provinces:", length(unique(df$province)), "\n")
cat("# obs:", nrow(df), "\n")
cat("Region counts (per year):\n"); print(table(df$region) / length(unique(df$year)))
cat("\n## Missing-value audit (panel data) ##\n")
print(colSums(is.na(df[, xs])))

## Compare original vs interpolated: how many cells were imputed?
orig <- read_excel(file.path(DATA_DIR, "original_data.xlsx"), skip = 2)
interp <- read_excel(file.path(DATA_DIR, "interpolation_processing.xlsx"), skip = 2)
names(orig) <- trimws(names(orig)); names(interp) <- trimws(names(interp))
## keep only 2013-2022 to align with panel
orig <- orig[orig$year %in% 2013:2022, ]
interp <- interp[interp$year %in% 2013:2022, ]

## the columns that became x1..x22 — match by approximate name
ind_orig_cols <- c(
  "Long-Distance Optical Cable Length (km)",
  "Internet Broadband Access Ports (10,000 units)",
  "Mobile Phone Base Stations (10,000 units)",
  "Number of Internet Domain Names (10,000 units)",
  "Number of IPv4 Addresses (10,000 units)",
  "Number of Internet Websites (10,000 units)",
  "Software Business Revenue (100 million yuan)",
  "Telecommunications Business Revenue (100 million yuan)",
  "Number of Websites per 100 Enterprises (units)",
  "Percentage of Enterprises with E-commerce Activities (%)",
  "E-commerce Sales Revenue (100 million yuan)",
  "Number of Computers per 100 People (units)",
  "Mobile Internet Users (10,000 households)",
  "Mobile Phone Subscribers (Year-End, 10,000 households)",
  "Digital TV Subscribers (10,000 households)",
  "Digital Inclusive Finance Index",
  "Digital-Real Economy Integration Index",
  "Percentage of IT-related Employment (%)",
  "Number of Bachelor's Degree Graduates (persons)",
  "Full-Time Equivalent (FTE) of R&D Personnel (person-years)",
  "Number of R&D Institutions (units)",
  "Region-Specific Patent Grants (units)"
)

orig_x18 <- orig[["Percentage of IT-related Employment (%)"]]
if (is.null(orig_x18)) {
  ## x18 is computed in original from numerator/denominator
  num <- orig[["Employment in Telecommunications, Software, and IT Services (10,000 persons, Non-private Sector)"]]
  den <- orig[["Total Employment in Non-private Sectors (10,000 persons)"]]
  orig_x18 <- num / den * 100
}

orig_mat <- as.data.frame(matrix(NA, nrow = nrow(orig), ncol = 22))
colnames(orig_mat) <- xs
for (i in seq_along(ind_orig_cols)) {
  nm <- ind_orig_cols[i]
  if (nm %in% names(orig)) orig_mat[[i]] <- orig[[nm]]
}
orig_mat$x18 <- orig_x18

miss_per_var <- colSums(is.na(orig_mat))
miss_pct <- miss_per_var / nrow(orig_mat) * 100
miss_tab <- data.frame(
  variable = xs,
  label    = unname(IND_LABELS[xs]),
  missing  = miss_per_var,
  missing_pct = round(miss_pct, 2)
)
miss_tab <- miss_tab[order(-miss_tab$missing), ]
cat("\n## Imputation rate per indicator (original 2013-2022, n=310) ##\n")
print(miss_tab, row.names = FALSE)
write.csv(miss_tab, file.path(TBL_DIR, "01_imputation_rate.csv"), row.names = FALSE)

## ---- Descriptive statistics table ----
desc_stats <- data.frame(
  variable = xs,
  label = unname(IND_LABELS[xs]),
  dimension = unname(IND_DIM[xs]),
  mean   = sapply(df[, xs], mean),
  sd     = sapply(df[, xs], sd),
  min    = sapply(df[, xs], min),
  median = sapply(df[, xs], median),
  max    = sapply(df[, xs], max),
  skew   = sapply(df[, xs], psych::skew),
  kurt   = sapply(df[, xs], psych::kurtosi)
)
desc_stats[, 4:10] <- lapply(desc_stats[, 4:10], function(v) signif(v, 4))
cat("\n## Descriptive statistics (panel, n=310) ##\n")
print(desc_stats, row.names = FALSE)
write.csv(desc_stats, file.path(TBL_DIR, "01_descriptive_stats.csv"), row.names = FALSE)

## ---- Region balance table ----
region_tbl <- df %>%
  filter(year == 2019) %>%
  count(region, name = "n_provinces")
cat("\n## Region balance (2019 cross-section) ##\n")
print(region_tbl)
write.csv(region_tbl, file.path(TBL_DIR, "01_region_balance.csv"), row.names = FALSE)

cat("\n[DONE 01_descriptive]\n")
