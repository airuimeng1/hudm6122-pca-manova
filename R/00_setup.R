## ============================================================
## 00_setup.R  —  Common setup for all analysis scripts
## HUDM 6122 Final Project: PCA + Factorial MANOVA on
## China Provincial Digital Economy Development (2013–2022)
## ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(psych)
  library(car)
  library(MVN)
  library(rstatix)
  library(heplots)
  library(knitr)
  library(scales)
})

## All paths are relative to the repository root.
## Scripts assume getwd() is the repo root: setwd("path/to/hudm6122-pca-manova").
PROJ_ROOT <- normalizePath(".", mustWork = TRUE)
DATA_DIR  <- file.path(PROJ_ROOT, "data")
OUT_DIR   <- file.path(PROJ_ROOT, "outputs")
FIG_DIR   <- file.path(OUT_DIR, "figures")
TBL_DIR   <- file.path(OUT_DIR, "tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

## NBS official 4-region classification, keyed by GB province codes
## (avoids any UTF-8 encoding pitfalls between Excel-read text and source literals)
REGION_BY_ID <- list(
  East      = c(110000, 120000, 130000, 310000, 320000, 330000,
                350000, 370000, 440000, 460000),
  Central   = c(140000, 340000, 360000, 410000, 420000, 430000),
  West      = c(150000, 450000, 500000, 510000, 520000, 530000,
                540000, 610000, 620000, 630000, 640000, 650000),
  Northeast = c(210000, 220000, 230000)
)

assign_region <- function(id) {
  for (r in names(REGION_BY_ID)) if (id %in% REGION_BY_ID[[r]]) return(r)
  NA_character_
}

## English province labels for plotting (keyed by GB id)
PROV_EN <- c(
  "110000"="Beijing","120000"="Tianjin","130000"="Hebei",
  "140000"="Shanxi","150000"="Inner Mongolia",
  "210000"="Liaoning","220000"="Jilin","230000"="Heilongjiang",
  "310000"="Shanghai","320000"="Jiangsu","330000"="Zhejiang",
  "340000"="Anhui","350000"="Fujian","360000"="Jiangxi",
  "370000"="Shandong","410000"="Henan","420000"="Hubei",
  "430000"="Hunan","440000"="Guangdong","450000"="Guangxi",
  "460000"="Hainan","500000"="Chongqing","510000"="Sichuan",
  "520000"="Guizhou","530000"="Yunnan","540000"="Tibet",
  "610000"="Shaanxi","620000"="Gansu","630000"="Qinghai",
  "640000"="Ningxia","650000"="Xinjiang"
)

## Period split (factorial MANOVA second factor)
PERIOD_SPLIT <- function(year) ifelse(year <= 2017, "Early(2013-2017)", "Late(2018-2022)")

## Indicator labels (x1..x22)
IND_LABELS <- c(
  x1  = "Long-distance optical cable",
  x2  = "Internet broadband ports",
  x3  = "Mobile phone base stations",
  x4  = "Internet domain names",
  x5  = "IPv4 addresses",
  x6  = "Internet websites",
  x7  = "Software business revenue",
  x8  = "Telecom business revenue",
  x9  = "Websites per 100 enterprises",
  x10 = "% enterprises w/ e-commerce",
  x11 = "E-commerce sales revenue",
  x12 = "Computers per 100 people",
  x13 = "Mobile internet users",
  x14 = "Mobile phone subscribers",
  x15 = "Digital TV subscribers",
  x16 = "Digital inclusive finance index",
  x17 = "Digital-real economy integration",
  x18 = "% IT-related employment",
  x19 = "Bachelor's graduates",
  x20 = "R&D personnel (FTE)",
  x21 = "R&D institutions",
  x22 = "Patent grants"
)

IND_DIM <- c(
  x1="Infrastructure", x2="Infrastructure", x3="Infrastructure",
  x4="Infrastructure", x5="Infrastructure", x6="Infrastructure",
  x7="Industry", x8="Industry", x9="Industry", x10="Industry", x11="Industry",
  x12="Environment", x13="Environment", x14="Environment", x15="Environment",
  x16="Environment", x17="Environment", x18="Environment", x19="Environment",
  x20="Environment", x21="Environment", x22="Environment"
)

## Loader: returns a tidy data frame with year, id, province, region, period, x1..x22
load_panel <- function() {
  df <- read_excel(file.path(DATA_DIR, "panel_data.xlsx"))
  names(df) <- trimws(names(df))
  df$region   <- vapply(df$id, assign_region, character(1))
  df$period   <- PERIOD_SPLIT(df$year)
  df$prov_en  <- unname(PROV_EN[format(df$id, scientific = FALSE, trim = TRUE)])
  df$region   <- factor(df$region, levels = c("East","Central","West","Northeast"))
  df$period   <- factor(df$period, levels = c("Early(2013-2017)","Late(2018-2022)"))
  df
}
