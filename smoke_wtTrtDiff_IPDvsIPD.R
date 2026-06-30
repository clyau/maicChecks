## IPD vs IPD: 2 responses (continuous Y, binary Y.bin) x 3 variance methods -> one CSV
suppressMessages(devtools::load_all('.', quiet = TRUE))

## --- the two IPD datasets shipped in `sim110` ------------------------
ipd1 <- sim110[sim110$study == "IPD A", ]
ipd2 <- sim110[sim110$study == "IPD B", ]

## --- symmetric exact matching: weights for BOTH studies --------------
w.out <- tryCatch(
  exmWt.2ipd(ipd1, ipd2,
             vars_to_match  = paste0("X", 1:5),
             cat_vars_to_01 = paste0("X", 1:3)),
  error = function(e) { warning("exmWt.2ipd failed: ", conditionMessage(e)); NULL }
)

if (is.null(w.out)) {
  w1 <- NULL; w2 <- NULL
} else {
  w1 <- w.out$ipd1$exm.wts
  w2 <- w.out$ipd2$exm.wts
}

## --- helper: one response -> 3-row data.frame ------------------------
one_response <- function(label, y1, y2) {

  base <- data.frame(
    response       = label,
    weight.method  = "exmWt.2ipd (symmetric QP)",
    n.ipd1         = length(y1),
    n.ipd2         = length(y2),
    ess.ipd1       = if (is.null(w1)) NA else sum(w1)^2 / sum(w1^2),
    ess.ipd2       = if (is.null(w2)) NA else sum(w2)^2 / sum(w2^2),
    stringsAsFactors = FALSE
  )

  if (is.null(w1) || is.null(w2)) {
    s <- data.frame(var.method = c("paper","weighted_ess","sq_residual"),
                    wt.y1 = NA, wt.y2 = NA, wt.diff = NA,
                    se1 = NA, se2 = NA, se.diff = NA,
                    ci.lower = NA, ci.upper = NA,
                    conf.level = 0.95, ess1 = NA, ess2 = NA)
    return(cbind(base, s, status = "WEIGHTS_FAILED", stringsAsFactors = FALSE))
  }

  s <- tryCatch(
    wtTrtDiff(ipd1.te = y1, w1 = w1,
              ipd2.te = y2, w2 = w2,
              var.method = "all")$summary,
    error = function(e) { warning("wtTrtDiff failed (", label, "): ",
                                  conditionMessage(e)); NULL })

  if (is.null(s)) {
    s <- data.frame(var.method = c("paper","weighted_ess","sq_residual"),
                    wt.y1 = NA, wt.y2 = NA, wt.diff = NA,
                    se1 = NA, se2 = NA, se.diff = NA,
                    ci.lower = NA, ci.upper = NA,
                    conf.level = 0.95, ess1 = NA, ess2 = NA)
    return(cbind(base, s, status = "ESTIMATION_FAILED", stringsAsFactors = FALSE))
  }

  cbind(base, s, status = "OK", stringsAsFactors = FALSE)
}

## --- assemble both responses ----------------------------------------
ipd.full <- rbind(
  one_response("Y (continuous)", ipd1$Y,     ipd2$Y),
  one_response("Y.bin (binary)", ipd1$Y.bin, ipd2$Y.bin)
)
rownames(ipd.full) <- NULL

## --- write CSV and report absolute path -----------------------------
path <- file.path(getwd(), "var_compare_IPDvsIPD.csv")
write.csv(ipd.full, path, row.names = FALSE)

cat("\n=====================================================\n")
cat("CSV WRITTEN TO:\n  ", path, "\n", sep = "")
cat("rows =", nrow(ipd.full), "  cols =", ncol(ipd.full), "\n")
cat("file size on disk:", file.info(path)$size, "bytes\n")
cat("=====================================================\n")