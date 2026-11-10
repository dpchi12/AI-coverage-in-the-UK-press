# =============================================================================
# Load data 
# =============================================================================

library(tseries)    # adf.test, kpss.test, pp.test
library(urca)       # ur.df for manual ADF with augmentations

# Load data
df <- read.csv("weekly_series_full.csv")
df$date <- as.Date(df$date)

cat("Loaded:", nrow(df), "rows\n")
cat("Outlets:", unique(df$outlet), "\n")
cat("Topics:", unique(df$topic), "\n")
cat("Total series:", length(unique(df$outlet)) * length(unique(df$topic)), "\n")

# =============================================================================
# Define outlets and topics
# =============================================================================

outlets <- c("Guardian", "Sun")

topics <- c(
  "ai, technology, robots",
  "brexit, economy, emissions",
  "disease, health, covid",
  "driverless, vehicles, uber",
  "bank, payments, crypto",
  "phone, device, apps",
  "digital rights, privacy, concerns"
)

# =============================================================================
# Stationarity test functions
# =============================================================================

# ADF test 
# H0: unit root (non-stationary)
# Reject H0 (p < 0.05) → stationary
run_adf <- function(y) {
  result <- tryCatch(
    adf.test(y, alternative="stationary"),
    error = function(e) NULL
  )
  if (is.null(result)) return(list(stat=NA, pval=NA, decision="failed"))
  decision <- ifelse(result$p.value < 0.05, "Stationary", "Non-stationary")
  list(stat=round(result$statistic, 4),
       pval=round(result$p.value, 4),
       decision=decision)
}

# Phillips-Perron test
# H0: unit root (non-stationary)
# Reject H0 (p < 0.05) → stationary
run_pp <- function(y) {
  result <- tryCatch(
    pp.test(y, alternative="stationary"),
    error = function(e) NULL
  )
  if (is.null(result)) return(list(stat=NA, pval=NA, decision="failed"))
  decision <- ifelse(result$p.value < 0.05, "Stationary", "Non-stationary")
  list(stat=round(result$statistic, 4),
       pval=round(result$p.value, 4),
       decision=decision)
}

# KPSS test
# H0: stationary (REVERSED null)
# Reject H0 (p < 0.05) → NON-stationary
run_kpss <- function(y) {
  result <- tryCatch(
    kpss.test(y, null="Level"),
    error = function(e) NULL
  )
  if (is.null(result)) return(list(stat=NA, pval=NA, decision="failed"))
  decision <- ifelse(result$p.value < 0.05, "Non-stationary", "Stationary")
  list(stat=round(result$statistic, 4),
       pval=round(result$p.value, 4),
       decision=decision)
}

# =============================================================================
# Run all tests for all 14 series
# =============================================================================

results <- list()

for (outlet in outlets) {
  for (topic in topics) {
   
    sub <- df[df$outlet == outlet & df$topic == topic, ]
    sub <- sub[order(sub$date), ]
    y   <- sub$share[!is.na(sub$share)]
    
    cat("\n", rep("=", 60), "\n")
    cat(outlet, "|", topic, "\n")
    cat("Observations:", length(y), "\n")
    cat(rep("=", 60), "\n")
  
    adf  <- run_adf(y)
    pp   <- run_pp(y)
    kpss <- run_kpss(y)
    
    cat(sprintf("  ADF   stat = %8.4f  p = %.4f  → %s\n",
                adf$stat,  adf$pval,  adf$decision))
    cat(sprintf("  PP    stat = %8.4f  p = %.4f  → %s\n",
                pp$stat,   pp$pval,   pp$decision))
    cat(sprintf("  KPSS  stat = %8.4f  p = %.4f  → %s\n",
                kpss$stat, kpss$pval, kpss$decision))
    
    n_stationary <- sum(c(adf$decision, pp$decision, kpss$decision) == "Stationary",
                        na.rm = TRUE)
  
    
    results[[paste(outlet, topic, sep=" | ")]] <- data.frame(
      outlet       = outlet,
      topic        = topic,
      n_obs        = length(y),
      adf_stat     = adf$stat,
      adf_pval     = adf$pval,
      adf_decision = adf$decision,
      pp_stat      = pp$stat,
      pp_pval      = pp$pval,
      pp_decision  = pp$decision,
      kpss_stat    = kpss$stat,
      kpss_pval    = kpss$pval,
      kpss_decision= kpss$decision,
      stringsAsFactors = FALSE
    )
  }
}

# =============================================================================
# Summary table
# =============================================================================

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL

cat("\n\n", rep("=", 60), "\n")
cat("STATIONARITY TEST SUMMARY — ALL 14 SERIES\n")
cat(rep("=", 60), "\n\n")

print(results_df[, c("outlet","topic","adf_decision",
                     "pp_decision","kpss_decision")])

write.csv(results_df, "stationarity_results.csv", row.names=FALSE)
cat("\nSaved to stationarity_results.csv\n")