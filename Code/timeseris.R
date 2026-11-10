# =============================================================================
# SETUP
# =============================================================================

library(forecast)
library(tseries)
library(lmtest)
library(sandwich)
library(zoo)

df           <- read.csv("weekly_series_full.csv")
df$date      <- as.Date(df$date)
CHATGPT_DATE <- as.Date("2022-11-30")


outlets <- c("Guardian", "Sun")
topics  <- c(
  "ai, technology, robots",
  "brexit, economy, emissions",
  "disease, health, covid",
  "driverless, vehicles, uber",
  "bank, payments, crypto",
  "phone, device, apps",
  "digital rights, privacy, concerns"
)

sig_stars <- function(p) {
  if (is.na(p)) return("n/a")
  ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*","n.s.")))
}

safe_coeftest <- function(m) {
  tryCatch(suppressWarnings(coeftest(m)), error=function(e) NULL)
}


# =============================================================================
# ITS FUNCTIONS
# =============================================================================

run_its_models <- function(y, n, TB) {
  trend <- 1:n
  DU    <- c(rep(0, TB), rep(1, n - TB))
  DT    <- c(rep(0, TB), 1:(n - TB))
  mA    <- lm(y ~ trend + DU)
  mB    <- lm(y ~ trend + DT)
  mC    <- lm(y ~ trend + DU + DT)
  aics  <- c(A=AIC(mA), B=AIC(mB), C=AIC(mC))
  bics  <- c(A=BIC(mA), B=BIC(mB), C=BIC(mC))
  list(mA=mA, mB=mB, mC=mC,
       aic=aics, bic=bics,
       best=c("A","B","C")[which.min(aics)],
       trend=trend, DU=DU, DT=DT)
}

plot_its <- function(y, dates, its, outlet_name, topic_name) {
  smooth    <- as.numeric(rollmean(y, k=3, fill=NA, align="center"))
  yr_ticks  <- seq(min(dates), max(dates), by="year")
  title_col <- ifelse(outlet_name == "Guardian", "#185FA5", "#C0392B")
  
  col_obs <- "grey80"
  col_mA  <- "#E74C3C"
  col_mB  <- "#3498DB"
  col_mC  <- "#2ECC71"
  
  # Fixed y-axis 0 to 100 for all topics
  plot(dates, smooth, type="l", col=col_obs, lwd=1,
       ylim=c(0, 100),
       main="", xlab="Year", ylab="Topic share (%)", xaxt="n")
  
  title(main     = paste0("ITS - ", outlet_name, " | ", topic_name),
        col.main = title_col, cex.main = 1.5, font.main = 2)
  
  axis(1, at=as.numeric(yr_ticks),
       labels=format(yr_ticks, "%Y"), las=0, cex.axis=0.9)
  
  abline(v=as.numeric(CHATGPT_DATE), col="black", lty=2, lwd=1.5)
  
  lines(dates, fitted(its$mA), col=col_mA, lwd=2, lty=1)
  lines(dates, fitted(its$mB), col=col_mB, lwd=2, lty=2)
  lines(dates, fitted(its$mC), col=col_mC, lwd=2, lty=3)
}

plot_arimax <- function(y, dates, arimax, outlet_name, topic_name) {
  smooth    <- as.numeric(rollmean(y, k=3, fill=NA, align="center"))
  yr_ticks  <- seq(min(dates), max(dates), by="year")
  title_col <- ifelse(outlet_name == "Guardian", "#185FA5", "#C0392B")
  
  col_obs <- "grey80"
  col_c1  <- "#E74C3C"
  col_c2  <- "#3498DB"
  col_c3  <- "#F39C12"
  col_c4  <- "#9B59B6"
  
  case_info <- list(
    list(key="case1", col=col_c1, lty=1),
    list(key="case2", col=col_c2, lty=2),
    list(key="case3", col=col_c3, lty=3),
    list(key="case4", col=col_c4, lty=4)
  )
  
  # Fixed y-axis 0 to 100 for all topics
  plot(dates, smooth, type="l", col=col_obs, lwd=1,
       ylim=c(0, 100),
       main="", xlab="Year", ylab="Topic share (%)", xaxt="n")
  
  title(main     = paste0("ARIMAX - ", outlet_name, " | ", topic_name),
        col.main = title_col, cex.main = 2, font.main = 2)
  
  axis(1, at=as.numeric(yr_ticks),
       labels=format(yr_ticks, "%Y"), las=0, cex.axis=1.5)
  
  abline(v=as.numeric(CHATGPT_DATE), col="black", lty=2, lwd=1.5)
  
  for (ci in case_info) {
    m <- arimax[[ci$key]]
    if (!is.null(m)) {
      lines(dates, fitted(m), col=ci$col, lwd=1.8, lty=ci$lty)
    }
  }
}

get_its_best_coefs <- function(its) {
  best_m <- list(mA=its$mA, mB=its$mB, mC=its$mC)[[paste0("m", its$best)]]
  ct     <- coef(summary(best_m))
  
  get_c <- function(name) {
    if (name %in% rownames(ct))
      list(est=round(ct[name,"Estimate"],4),
           sig=sig_stars(ct[name,"Pr(>|t|)"]))
    else list(est=NA, sig="")
  }
  
  list(model  = its$best,
       aic    = round(its$aic[its$best], 2),
       bic    = round(its$bic[its$best], 2),
       r2     = round(summary(best_m)$r.squared, 4),
       trend  = get_c("trend"),
       DU     = get_c("DU"),
       DT     = get_c("DT"))
}

get_arimax_best_coefs <- function(arimax) {
  best_m <- arimax[[arimax$best_case]]
  if (is.null(best_m)) return(NULL)
  
  int_var  <- ifelse(grepl("case2|case4", arimax$best_case), "PULSE", "SHIFT")
  ct       <- safe_coeftest(best_m)
  valid    <- !is.null(ct) && !any(is.nan(ct[,"Std. Error"]))
  
  n_obs    <- length(best_m$residuals)
  k_params <- length(coef(best_m))
  bic_val  <- round(-2 * best_m$loglik + k_params * log(n_obs), 2)
  
  get_c <- function(name) {
    if (valid && name %in% rownames(ct))
      list(est=round(ct[name,"Estimate"],4),
           sig=sig_stars(ct[name,"Pr(>|z|)"]))
    else if (name %in% names(coef(best_m)))
      list(est=round(coef(best_m)[name],4), sig="conv.")
    else
      list(est=NA, sig="")
  }
  
  list(case      = arimax$best_case,
       aic       = round(min(arimax$aics, na.rm=TRUE), 2),
       bic       = bic_val,
       int_var   = int_var,
       intercept = get_c("intercept"),
       trend     = get_c("trend"),
       int_coef  = get_c(int_var))
}

run_arimax_cases <- function(y, n, TB) {
  trend <- 1:n
  SHIFT <- c(rep(0, TB), rep(1, n - TB))
  PULSE <- c(rep(0, TB-1), 1, rep(0, n - TB))
  
  results <- list()
  results$case1 <- tryCatch(arima(y, order=c(0,0,1), include.mean=TRUE,
                                  xreg=cbind(trend=trend, SHIFT=SHIFT)), error=function(e) NULL)
  results$case2 <- tryCatch(arima(y, order=c(0,0,1), include.mean=TRUE,
                                  xreg=cbind(trend=trend, PULSE=PULSE)), error=function(e) NULL)
  results$case3 <- tryCatch(arima(y, order=c(1,0,1), include.mean=TRUE,
                                  xreg=cbind(trend=trend, SHIFT=SHIFT)), error=function(e) NULL)
  results$case4 <- tryCatch(arima(y, order=c(1,0,1), include.mean=TRUE,
                                  xreg=cbind(trend=trend, PULSE=PULSE)), error=function(e) NULL)
  
  aics <- sapply(results[c("case1","case2","case3","case4")],
                 function(m) if (!is.null(m)) AIC(m) else NA)
  
  results$aics      <- aics
  results$best_case <- names(which.min(aics))
  results$SHIFT     <- SHIFT
  results$PULSE     <- PULSE
  results$trend     <- trend
  return(results)
}


# =============================================================================
# STANDALONE LEGEND PLOTS
# =============================================================================

dir.create("plots_its",    showWarnings=FALSE)
dir.create("plots_arimax", showWarnings=FALSE)

png("plots_its/LEGEND_its.png", width=700, height=100, res=130)
par(mar=c(0,0,0,0))
plot.new()
legend("center",
       legend=c("Model A - Level only (DU)",
                "Model B - Slope only (DT)",
                "Model C - Level + Slope",
                "ChatGPT launch (Nov 2022)",
                "Observed (3-week avg)"),
       col=c("#E74C3C","#3498DB","#2ECC71","black","grey80"),
       lty=c(1,2,3,2,1), lwd=c(2,2,2,1.5,1),
       ncol=3, bty="n", cex=0.8)
dev.off()
cat("Saved: plots_its/LEGEND_its.png\n")

png("plots_arimax/LEGEND_arimax.png", width=800, height=100, res=130)
par(mar=c(0,0,0,0))
plot.new()
legend("center",
       legend=c("Case 1: ARIMA(0,0,1)+SHIFT - no memory, permanent",
                "Case 2: ARIMA(0,0,1)+PULSE - no memory, transient",
                "Case 3: ARIMA(1,0,1)+SHIFT - memory, permanent",
                "Case 4: ARIMA(1,0,1)+PULSE - memory, transient",
                "ChatGPT launch (Nov 2022)",
                "Observed (3-week avg)"),
       col=c("#E74C3C","#3498DB","#F39C12","#9B59B6","black","grey80"),
       lty=c(1,2,3,4,2,1), lwd=c(1.8,1.8,1.8,1.8,1.5,1),
       ncol=2, bty="n", cex=0.8)
dev.off()
cat("Saved: plots_arimax/LEGEND_arimax.png\n")


# =============================================================================
# MAIN LOOP — ITS
# =============================================================================

its_summary <- list()

for (outlet in outlets) {
  for (topic in topics) {
    
    cat("\n", rep("=",60), "\n")
    cat("ITS |", outlet, "|", topic, "\n")
    cat(rep("=",60), "\n")
    
    sub   <- df[df$outlet==outlet & df$topic==topic, ]
    sub   <- sub[order(sub$date), ]
    sub   <- sub[!is.na(sub$share), ]
    y     <- sub$share
    dates <- sub$date
    n     <- length(y)
    TB <- which(dates >= CHATGPT_DATE)[1] - 1
    cat("N observations:", n, "| TB:", TB, "\n")
    if (n < 50) next
    
    its <- run_its_models(y, n, TB)
    
    cat("\n--- Model A ---\n"); print(summary(its$mA))
    cat("AIC:", round(its$aic["A"],2), " BIC:", round(its$bic["A"],2), "\n")
    cat("\n--- Model B ---\n"); print(summary(its$mB))
    cat("AIC:", round(its$aic["B"],2), " BIC:", round(its$bic["B"],2), "\n")
    cat("\n--- Model C ---\n"); print(summary(its$mC))
    cat("AIC:", round(its$aic["C"],2), " BIC:", round(its$bic["C"],2), "\n")
    cat("\nBest by AIC: Model", its$best, "\n")
    
    best_m <- list(mA=its$mA, mB=its$mB, mC=its$mC)[[paste0("m",its$best)]]
    cat("\nBG tests (best model", its$best, "):\n")
    for (ord in 1:5) {
      bg <- bgtest(best_m, order=ord)
      cat("  Order", ord, ": p=", round(bg$p.value,4),
          ifelse(bg$p.value>0.05,"✓","⚠"), "\n")
    }
    cat("\nHAC corrected coefficients:\n")
    print(coeftest(best_m, vcov=NeweyWest(best_m)))
    
    fname <- gsub("[,/ ]", "_", paste0(outlet, "_", topic))
    png(paste0("plots_its/its_", fname, ".png"),
        width=1400, height=700, res=130)
    plot_its(y, dates, its, outlet, topic)
    dev.off()
    cat("Saved: plots_its/its_", fname, ".png\n")
    
    bc <- get_its_best_coefs(its)
    its_summary[[paste(outlet,topic)]] <- data.frame(
      outlet     = outlet,
      topic      = topic,
      best_model = bc$model,
      aic        = bc$aic,
      bic        = bc$bic,
      r2         = bc$r2,
      trend_est  = bc$trend$est,
      trend_sig  = bc$trend$sig,
      DU_est     = bc$DU$est,
      DU_sig     = bc$DU$sig,
      DT_est     = bc$DT$est,
      DT_sig     = bc$DT$sig,
      stringsAsFactors = FALSE
    )
  }
}
cat("\nAll ITS plots saved in: plots_its/\n")


# =============================================================================
# MAIN LOOP — ARIMAX
# =============================================================================

arimax_summary <- list()

for (outlet in outlets) {
  for (topic in topics) {
    
    cat("\n", rep("=",60), "\n")
    cat("ARIMAX |", outlet, "|", topic, "\n")
    cat(rep("=",60), "\n")
    
    sub   <- df[df$outlet==outlet & df$topic==topic, ]
    sub   <- sub[order(sub$date), ]
    sub   <- sub[!is.na(sub$share), ]
    y     <- sub$share
    dates <- sub$date
    n     <- length(y)
    if (n < 50) next
    
    arimax <- run_arimax_cases(y, n, TB)
    
    case_labels <- list(
      case1 = "Case 1: ARIMA(0,0,1)+SHIFT",
      case2 = "Case 2: ARIMA(0,0,1)+PULSE",
      case3 = "Case 3: ARIMA(1,0,1)+SHIFT",
      case4 = "Case 4: ARIMA(1,0,1)+PULSE"
    )
    
    for (cn in names(case_labels)) {
      m <- arimax[[cn]]
      if (!is.null(m)) {
        cat("\n---", case_labels[[cn]], "---\n")
        ct <- safe_coeftest(m)
        if (!is.null(ct) && !any(is.nan(ct[,"Std. Error"]))) {
          print(ct)
        } else {
          cat("Note: convergence issue — coefficients only\n")
          print(round(coef(m), 5))
        }
        n_obs    <- length(m$residuals)
        k_params <- length(coef(m))
        bic_m    <- round(-2 * m$loglik + k_params * log(n_obs), 2)
        cat("AIC:", round(AIC(m),2), " BIC:", bic_m, "\n")
      }
    }
    cat("\nBest case:", arimax$best_case,
        "| AIC =", round(min(arimax$aics,na.rm=TRUE),1), "\n")
    
    fname <- gsub("[,/ ]", "_", paste0(outlet, "_", topic))
    png(paste0("plots_arimax/arimax_", fname, ".png"),
        width=1400, height=700, res=130)
    plot_arimax(y, dates, arimax, outlet, topic)
    dev.off()
    cat("Saved: plots_arimax/arimax_", fname, ".png\n")
    
    bc <- get_arimax_best_coefs(arimax)
    if (!is.null(bc)) {
      arimax_summary[[paste(outlet,topic)]] <- data.frame(
        outlet      = outlet,
        topic       = topic,
        best_case   = bc$case,
        aic         = bc$aic,
        bic         = bc$bic,
        int_type    = bc$int_var,
        int_est     = bc$int_coef$est,
        int_sig     = bc$int_coef$sig,
        trend_est   = bc$trend$est,
        trend_sig   = bc$trend$sig,
        stringsAsFactors = FALSE
      )
    }
  }
}
cat("\nAll ARIMAX plots saved in: plots_arimax/\n")


# =============================================================================
# SUMMARY TABLES
# =============================================================================

its_df <- do.call(rbind, its_summary)
rownames(its_df) <- NULL
cat("\n\n=== ITS BEST MODEL RESULTS ===\n")
print(its_df)
write.csv(its_df, "its_best_results.csv", row.names=FALSE)

arimax_df <- do.call(rbind, arimax_summary)
rownames(arimax_df) <- NULL
cat("\n\n=== ARIMAX BEST CASE RESULTS ===\n")
print(arimax_df)
write.csv(arimax_df, "arimax_best_results.csv", row.names=FALSE)

cat("\nSaved: its_best_results.csv | arimax_best_results.csv\n")

