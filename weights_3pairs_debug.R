## Strategy 1: simulate 3 IPD pairs with controlled overlap -> one CSV
## Weights via exmWt.2ipd(); variance via wtTrtDiff() with var.method = "all".

## ---- preamble ------------------------------------------------------
if (basename(getwd()) == "R" ||
    grepl("[/\\\\]R$", normalizePath(getwd(), mustWork = FALSE))) {
  stop("This script must NOT be placed inside the package's R/ folder.")
}
.tries <- 0
while (!file.exists("DESCRIPTION") && .tries < 5) {
  if (file.exists("maicChecks/DESCRIPTION")) {
    setwd("maicChecks")
    break
  }
  setwd("..")
  .tries <- .tries + 1
}
if (!file.exists("DESCRIPTION")) {
  stop("Could not locate the maicChecks package root.")
}
cat("Working directory:", getwd(), "\n")
suppressMessages(devtools::load_all('.', quiet = TRUE))
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("install.packages('mvtnorm')")
}
## ---- end preamble --------------------------------------------------

## ---------------------------------------------------------------------
## simulator -- accepts a 5-vector of mean shifts (one per covariate)
## ---------------------------------------------------------------------
simulate_pair <- function(shift_vec, n_per_arm = 300, rho = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(shift_vec) == 5)
  Sigma <- matrix(rho, 5, 5)
  diag(Sigma) <- 1
  A <- mvtnorm::rmvnorm(n_per_arm, mean = rep(0, 5), sigma = Sigma)
  B <- mvtnorm::rmvnorm(n_per_arm, mean = shift_vec, sigma = Sigma)
  colnames(A) <- paste0("X", 1:5)
  colnames(B) <- paste0("X", 1:5)
  A <- as.data.frame(A)
  B <- as.data.frame(B)

  make_Y <- function(d) {
    0.3 * d$X1 + 0.2 * d$X3 + 0.1 * d$X4 + rnorm(nrow(d), 0, 1)
  }
  A$Y <- make_Y(A)
  B$Y <- make_Y(B)

  thr <- list(X1 = qnorm(0.238), X2 = qnorm(0.312), X3 = qnorm(0.5))
  for (v in names(thr)) {
    A[[v]] <- as.integer(A[[v]] > thr[[v]])
    B[[v]] <- as.integer(B[[v]] > thr[[v]])
  }

  ymed <- median(c(A$Y, B$Y))
  A$Y.bin <- as.integer(A$Y > ymed)
  B$Y.bin <- as.integer(B$Y > ymed)

  A$study <- "IPD A"
  B$study <- "IPD B"
  rbind(A, B)[, c("study", paste0("X", 1:5), "Y", "Y.bin")]
}

## ---------------------------------------------------------------------
## three pairs with increasing overlap loss
## ---------------------------------------------------------------------
pairs <- list(
  list(label = "high_overlap     (shifts 0.1/0/0/0.1/0)",
       d = simulate_pair(c(0.1, 0, 0, 0.1, 0),  seed = 101)),
  list(label = "moderate_overlap (shifts 0.8/0/0/1.0/0)",
       d = simulate_pair(c(0.8, 0, 0, 1.0, 0),  seed = 102)),
  list(label = "low_overlap      (shifts 1.5/0/0/1.8/1.2)",
       d = simulate_pair(c(1.5, 0, 0, 1.8, 1.2), seed = 103))
)

cat("\nOverlap diagnostic (per covariate, mean A / mean B / SMD):\n")
for (p in pairs) {
  cat(sprintf("  %s\n", p$label))
  a <- p$d[p$d$study == "IPD A", paste0("X", 1:5)]
  b <- p$d[p$d$study == "IPD B", paste0("X", 1:5)]
  for (v in paste0("X", 1:5)) {
    pooled <- sqrt((var(a[[v]]) + var(b[[v]])) / 2)
    smd    <- abs(mean(a[[v]]) - mean(b[[v]])) / pooled
    cat(sprintf("    %s  A=%6.3f  B=%6.3f  SMD=%5.3f\n",
                v, mean(a[[v]]), mean(b[[v]]), smd))
  }
}

## ---------------------------------------------------------------------
## numerical safeguard for QP boundary noise
## ---------------------------------------------------------------------
clip_qp_noise <- function(w, tol = 1e-8) {
  if (any(w < -tol)) {
    warning("QP returned weights more negative than tol; min = ",
            signif(min(w), 3))
  }
  pmax(w, 0)
}

## ---------------------------------------------------------------------
## analysis: weights via exmWt.2ipd(), SEs via wtTrtDiff()
## ---------------------------------------------------------------------
run_one <- function(d, dataset.name) {
  ipd1 <- d[d$study == "IPD A", ]
  ipd2 <- d[d$study == "IPD B", ]

  w.out <- tryCatch(
    exmWt.2ipd(ipd1, ipd2,
               vars_to_match  = paste0("X", 1:5),
               cat_vars_to_01 = paste0("X", 1:3)),
    error = function(e) {
      warning(dataset.name, ": exmWt.2ipd failed: ", conditionMessage(e))
      NULL
    })

  if (is.null(w.out)) {
    w1 <- NULL
    w2 <- NULL
  } else {
    w1 <- clip_qp_noise(w.out$ipd1$exm.wts)
    w2 <- clip_qp_noise(w.out$ipd2$exm.wts)
  }

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
      s <- data.frame(var.method = c("paper", "weighted_ess", "sq_residual"),
                      wt.y1 = NA, wt.y2 = NA, wt.diff = NA,
                      se1 = NA, se2 = NA, se.diff = NA,
                      ci.lower = NA, ci.upper = NA,
                      conf.level = 0.95, ess1 = NA, ess2 = NA)
      return(cbind(base, s, status = "WEIGHTS_FAILED",
                   stringsAsFactors = FALSE))
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
cat("=====================================================\n")