# validate/generate-reference.R
# Generates reference values from R metafor for JS validation
library(metafor)

# BCG vaccine data
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

cat("=== BCG Data (yi, vi) ===\n")
cat("yi:", paste(sprintf("%.8f", dat$yi), collapse=", "), "\n")
cat("vi:", paste(sprintf("%.8f", dat$vi), collapse=", "), "\n\n")

# REML
res_reml <- rma(yi, vi, data=dat, method="REML")
cat("=== REML ===\n")
cat("tau2:", sprintf("%.8f", res_reml$tau2), "\n")
cat("theta:", sprintf("%.8f", res_reml$beta), "\n")
cat("se:", sprintf("%.8f", res_reml$se), "\n")
cat("ci:", sprintf("%.8f", res_reml$ci.lb), sprintf("%.8f", res_reml$ci.ub), "\n")
cat("I2:", sprintf("%.4f", res_reml$I2), "\n")
cat("H2:", sprintf("%.4f", res_reml$H2), "\n")
cat("Q:", sprintf("%.4f", res_reml$QE), "\n")
cat("Qp:", sprintf("%.8f", res_reml$QEp), "\n\n")

# DL
res_dl <- rma(yi, vi, data=dat, method="DL")
cat("=== DL ===\n")
cat("tau2:", sprintf("%.8f", res_dl$tau2), "\n")
cat("theta:", sprintf("%.8f", res_dl$beta), "\n\n")

# PM
res_pm <- rma(yi, vi, data=dat, method="PM")
cat("=== PM ===\n")
cat("tau2:", sprintf("%.8f", res_pm$tau2), "\n")
cat("theta:", sprintf("%.8f", res_pm$beta), "\n\n")

# HKSJ
res_hksj <- rma(yi, vi, data=dat, method="REML", test="knha")
cat("=== HKSJ (REML + knha) ===\n")
cat("theta:", sprintf("%.8f", res_hksj$beta), "\n")
cat("ci:", sprintf("%.8f", res_hksj$ci.lb), sprintf("%.8f", res_hksj$ci.ub), "\n\n")

# Prediction interval
pred <- predict(res_reml)
cat("=== Prediction Interval ===\n")
cat("pi:", sprintf("%.8f", pred$pi.lb), sprintf("%.8f", pred$pi.ub), "\n\n")

# Fixed-effect
res_fe <- rma(yi, vi, data=dat, method="FE")
cat("=== Fixed-Effect ===\n")
cat("theta:", sprintf("%.8f", res_fe$beta), "\n")
cat("Q:", sprintf("%.8f", res_fe$QE), "\n\n")

# Write JSON for automated comparison
ref <- list(
  bcg = list(
    yi = as.numeric(dat$yi),
    vi = as.numeric(dat$vi),
    reml = list(tau2=res_reml$tau2, theta=as.numeric(res_reml$beta),
                se=as.numeric(res_reml$se),
                ci=c(res_reml$ci.lb, res_reml$ci.ub), I2=res_reml$I2,
                H2=res_reml$H2, Q=res_reml$QE, Qp=res_reml$QEp),
    dl = list(tau2=res_dl$tau2, theta=as.numeric(res_dl$beta)),
    pm = list(tau2=res_pm$tau2, theta=as.numeric(res_pm$beta)),
    hksj = list(theta=as.numeric(res_hksj$beta),
                ci=c(res_hksj$ci.lb, res_hksj$ci.ub)),
    pi = c(pred$pi.lb, pred$pi.ub),
    fe = list(theta=as.numeric(res_fe$beta), Q=res_fe$QE)
  )
)

jsonlite::write_json(ref, "test/reference-values.json", pretty=TRUE, digits=10)
cat("Reference values written to test/reference-values.json\n")
