## IPD vs IPD for 3 datasets (sim110, sim120, sim130) -> one CSV
suppressMessages(devtools::load_all('.', quiet = TRUE))

run_one <- function(d, dataset.name) {

  ipd1 <- d[d$study == "IPD A", ]
  ipd2 <- d[d$study == "IPD B", ]

  w.out <- exmWt.2ipd(ipd1, ipd2,
                      vars_to_match  = paste0("X", 1:5),
                      cat_vars_to_01 = paste0("X", 1:3))
  w1 <- w.out$ipd1$exm.wts
  w2 <- w.out$ipd2$exm.wts

  one_response <- function(label, y1, y2) {
    base <- data.frame(
      dataset       = dataset.name,
      response      = label,
      weight.method = "exmWt.2ipd (symmetric QP)",
      n.ipd1        = length(y1),
      n.ipd2        = length(y2),
      ess.ipd1      = sum(w1)^2 / sum(w1^2),
      ess.ipd2      = sum(w2)^2 / sum(w2^2),
      stringsAsFactors = FALSE)

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

ipd.full <- rbind(
  run_one(sim110, "sim110"),
  run_one(sim120, "sim120"),
  run_one(sim130, "sim130")
)
rownames(ipd.full) <- NULL

path <- file.path(getwd(), "var_compare_IPDvsIPD_3pairs.csv")
write.csv(ipd.full, path, row.names = FALSE)

cat("\n=====================================================\n")
cat("CSV WRITTEN TO:\n  ", path, "\n", sep = "")
cat("rows =", nrow(ipd.full), "  cols =", ncol(ipd.full), "\n")
cat("file size on disk:", file.info(path)$size, "bytes\n")
cat("=====================================================\n")