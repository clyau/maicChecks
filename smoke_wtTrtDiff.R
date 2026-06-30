## IPD vs AD: 3 scenarios x 3 variance methods -> one CSV
suppressMessages(devtools::load_all('.', quiet = TRUE))

rows <- list()
for (k in seq_len(nrow(eAD))) {

  ad   <- eAD[k, ]
  lab  <- paste0("scenario_", k)

  ## try to compute weights; on failure -> NA rows
  w <- tryCatch(
         suppressWarnings(maicWt(eIPD[, c("y1","y2")], ad[, c("y1","y2")])$maic.wt.rs),
         error = function(e) NULL)

  if (is.null(w) || any(!is.finite(w))) {
    warning(lab, ": maicWt failed (AD likely outside convex hull) -> NA")
    s <- data.frame(var.method = c("paper","weighted_ess","sq_residual"),
                    wt.y1 = NA, wt.y2 = NA, wt.diff = NA,
                    se1 = NA, se2 = NA, se.diff = NA,
                    ci.lower = NA, ci.upper = NA,
                    conf.level = 0.95, ess1 = NA, ess2 = NA)
    wm.y1 <- NA; wm.y2 <- NA; status <- "WEIGHTS_FAILED"
  } else {
    s <- wtTrtDiff(ipd1.te = eIPD$r.cont, w1 = w,
                   ad.mean = ad$r.cont.mean,
                   ad.sd   = ad$r.cont.sd,
                   ad.n    = ad$r.cont.n,
                   var.method = "all")$summary
    wm.y1 <- sum(w * eIPD$y1) / sum(w)
    wm.y2 <- sum(w * eIPD$y2) / sum(w)
    status <- "OK"
  }

  rows[[k]] <- data.frame(
    scenario      = lab,
    weight.method = "maicWt",
    response      = "r.cont",
    n.ipd         = nrow(eIPD),
    n.ad          = ad$r.cont.n,
    ad.y1         = ad$y1,
    ad.y2         = ad$y2,
    ad.r.mean     = ad$r.cont.mean,
    ad.r.sd       = ad$r.cont.sd,
    wm.y1.ipd     = wm.y1,
    wm.y2.ipd     = wm.y2,
    s,
    status        = status,
    stringsAsFactors = FALSE
  )
}

ad.full <- do.call(rbind, rows)

## write CSV and show ABSOLUTE path
path <- file.path(getwd(), "var_compare_IPDvsAD.csv")
write.csv(ad.full, path, row.names = FALSE)

cat("\n=====================================================\n")
cat("CSV WRITTEN TO:\n  ", path, "\n", sep = "")
cat("rows =", nrow(ad.full), "  cols =", ncol(ad.full), "\n")
cat("file size on disk:", file.info(path)$size, "bytes\n")
cat("=====================================================\n")