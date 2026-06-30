## Strategy 1: simulate 3 IPD pairs with controlled overlap -> one CSV
if (!file.exists("DESCRIPTION") && file.exists("maicChecks/DESCRIPTION"))
  setwd("maicChecks")
suppressMessages(devtools::load_all('.', quiet = TRUE))

## ---------------------------------------------------------------------
## (1) simulator -- block 1 of Glimm & Yau (2026) supplement, simplified
## ---------------------------------------------------------------------
simulate_pair <- function(shift, n_per_arm = 300, rho = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- matrix(rho, 5, 5); diag(Sigma) <- 1

  A <- mvtnorm::rmvnorm(n_per_arm, mean = rep(0, 5),               sigma = Sigma)
  B <- mvtnorm::rmvnorm(n_per_arm, mean = c(shift, 0, 0, 0, 0),    sigma = Sigma)
  colnames(A) <- colnames(B) <- paste0("X", 1:5)
  A <- as.data.frame(A); B <- as.data.frame(B)

  ## continuous response (no treatment effect in this DGP)
  make_Y <- function(d) 0.3 * d$X1 + 0.2 * d$X3 + 0.1 * d$X4 +
                        rnorm(nrow(d), 0, 1)
  A$Y <- make_Y(A); B$Y <- make_Y(B)

  ## dichotomize X1, X2, X3 (paper thresholds for X1, X2; X3 here also binary)
  thr <- list(X1 = qnorm(0.238), X2 = qnorm(0.312), X3 = qnorm(0.50))
  for (v in names(thr)) {
    A[[v]] <- as.integer(A[[v]] > thr[[v]])
    B[[v]] <- as.integer(B[[v]] > thr[[v]])
  }

  ## binary response: dichotomize Y at pooled median
  ymed  <- median(c(A$Y, B$Y))
  A$Y.bin <- as.integer(A$Y > ymed)
  B$Y.bin <- as.integer(B$Y > ymed)

  A$study <- "IPD A"; B$study <- "IPD B"
  rbind(A, B)[, c("study", paste0("X", 1:5), "Y", "Y.bin")]
}

## ---------------------------------------------------------------------
## (2) three pairs with increasing mean-shift on X1
## ---------------------------------------------------------------------
pairs <- list(
  list(label = "sim_shift_0.2 (high overlap)",     d = simulate_pair(0.2, seed = 101)),
  list(label = "sim_shift_0.6 (moderate overlap)", d = simulate_pair(0.6, seed = 102)),
  list(label = "sim_shift_1.2 (low overlap)",      d = simulate_pair(1.2, seed = 103))
)

cat("\nOverlap diagnostic (SMD on X1, the shifted covariate):\n")
for (p in pairs) {
  a <- p$d[p$d$study == "IPD A", "X1"]
  b <- p$d[p$d$study == "IPD B", "X1"]
  pooled <- sqrt((var(a) + var(b)) / 2)
  cat(sprintf("  %-40s  mean(A)=%.3f  mean(B)=%.3f  SMD=%.3f\n",
              p$label, mean(a), mean(b), abs(mean(a) - mean(b)) / pooled))
}

## ---------------------------------------------------------------------
## (3) analysis identical to smoke_wtTrtDiff_IPDvsIPD_3pairs.R
## ---------------------------------------------------------------------
run_one <- function(d, dataset.name) {

  ipd1 <- d[d$study == "IPD A", ]
  ipd2 <- d[d$study == "IPD B", ]

  w.out <- tryCatch(
    exmWt.2ipd(ipd1, ipd2,
               vars_to_match  = paste0("X", 1:5),
               cat_vars_to_01 = paste0("X", 1:3)),
    error = function(e) {
      warning(dataset.name, ": exmWt.2ipd failed: ", conditionMessage(e)); NULL })

  if (is.null(w.out)) { w1 <- NULL; w2 <- NULL }
  else { w1 <- w.out$ipd1$exm.wts; w2 <- w.out$ipd2$exm.wts }

  one_response <- function(label, y1, y2) {
    base <- data.frame(
      dataset       = dataset.name,
      response      = label,
      weight.method = "exmWt.2ipd (symmetric QP)",
      n.ipd1        = length(y1),
      n.ipd2        = length(y2),
      ess.ipd1      = if (is.null(w1)) NA else sum(w1)^2 / sum(w1^2),
      ess.ipd2      = if (is.null(w2)) NA else sum(w2)^2 / sum(w2^2),
      stringsAsFactors = FALSE)

    if (is.null(w1)) {
      s <- data.frame(var.method = c("paper","weighted_ess","sq_residual"),
                      wt.y1 = NA, wt.y2 = NA, wt.diff = NA,
                      se1 = NA, se2 = NA, se.diff = NA,
                      ci.lower = NA, ci.upper = NA,
                      conf.level = 0.95, ess1 = NA, ess2 = NA)
      return(cbind(base, s, status = "WEIGHTS_FAILED", stringsAsFactors = FALSE))
    }

    s <- wtTrtDiff(ipd1.te = y1, w1 = w1,
                   ipd2.te = y2, w2 = w2,
                   var.method = "all")$summary

    cbind(base, s, status = "OK", stringsAsFactors = FALSE)
  }

  rbind(
    one_response("Y (continuous)", ipd1$Y,     ipd2$Y),
    one_response("Y.bin (binary)", ipd1$Y.bin, ipd2$Y.bin)
  )
}

ipd.full <- do.call(rbind, lapply(pairs, function(p) run_one(p$d, p$label)))
rownames(ipd.full) <- NULL

path <- file.path(getwd(), "var_compare_IPDvsIPD_3pairs_sim.csv")
write.csv(ipd.full, path, row.names = FALSE)

cat("\n=====================================================\n")
cat("CSV WRITTEN TO:\n  ", path, "\n", sep = "")
cat("rows =", nrow(ipd.full), "  cols =", ncol(ipd.full), "\n")
cat("file size on disk:", file.info(path)$size, "bytes\n")
cat("=====================================================\n")