## smoke test for wtTrtDiff()
suppressMessages(devtools::load_all('.', quiet = TRUE))

cat('\n========== (1) IPD vs AD, with ad.sd / ad.n -- continuous ==========\n')
w.out <- maicWt(eIPD[, c('y1', 'y2')], eAD[1, c('y1', 'y2')])
r1    <- wtTrtDiff(ipd1.te = eIPD$r.cont, w1 = w.out$maic.wt.rs,
                   ad.mean = eAD$r.cont.mean[1],
                   ad.sd   = eAD$r.cont.sd[1],
                   ad.n    = eAD$r.cont.n[1])
cat('wt.y1     :', round(r1$wt.y1, 4), '\n')
cat('wt.y2     :', round(r1$wt.y2, 4), '\n')
cat('wt.diff   :', round(r1$wt.diff, 4), '\n')
cat('se1       :', round(r1$se1, 4), '\n')
cat('se2       :', round(r1$se2, 4), '\n')
cat('se.diff   :', round(r1$se.diff, 4), '\n')
cat('CI (95)   : [', round(r1$ci.lower, 4), ',', round(r1$ci.upper, 4), ']\n')
cat('ess1      :', round(r1$ess1, 2), '\n')
cat('ess2      :', r1$ess2, '\n')

cat('\n========== (2) IPD vs AD, AD as constant -- binary ==========\n')
r2 <- wtTrtDiff(ipd1.te = eIPD$r.bin, w1 = w.out$maic.wt.rs,
                ad.mean = eAD$r.bin.p[1])
cat('wt.y2     :', round(r2$wt.y2, 4), '  (should equal eAD$r.bin.p[1])\n')
cat('se2       :', r2$se2, '  (should be 0)\n')
cat('ess2      :', r2$ess2, '  (should be NA)\n')
cat('se.diff   :', round(r2$se.diff, 4), '  (should equal r2$se1)\n')
stopifnot(r2$se2 == 0)
stopifnot(is.na(r2$ess2))
stopifnot(abs(r2$se.diff - r2$se1) < 1e-12)

cat('\n========== (3) IPD vs AD, binary with ad.sd = sqrt(p(1-p)) ==========\n')
r3 <- wtTrtDiff(ipd1.te = eIPD$r.bin, w1 = w.out$maic.wt.rs,
                ad.mean = eAD$r.bin.p[1],
                ad.sd   = sqrt(eAD$r.bin.p[1] * (1 - eAD$r.bin.p[1])),
                ad.n    = eAD$r.bin.n[1])
cat('wt.y1   :', round(r3$wt.y1, 4), '\n')
cat('wt.y2   :', round(r3$wt.y2, 4), '\n')
cat('wt.diff :', round(r3$wt.diff, 4), '\n')
cat('se.diff :', round(r3$se.diff, 4), '\n')
cat('CI (95) : [', round(r3$ci.lower, 4), ',', round(r3$ci.upper, 4), ']\n')

cat('\n========== (4) IPD vs IPD -- continuous (sim110$Y) ==========\n')
ipd1  <- sim110[sim110$study == 'IPD A', ]
ipd2  <- sim110[sim110$study == 'IPD B', ]
w.out <- exmWt.2ipd(ipd1, ipd2,
                    vars_to_match  = paste0('X', 1:5),
                    cat_vars_to_01 = paste0('X', 1:3))
r4 <- wtTrtDiff(ipd1.te = ipd1$Y, w1 = w.out$ipd1$exm.wts,
                ipd2.te = ipd2$Y, w2 = w.out$ipd2$exm.wts)
cat('wt.y1 (A) :', round(r4$wt.y1, 4), '\n')
cat('wt.y2 (B) :', round(r4$wt.y2, 4), '\n')
cat('wt.diff   :', round(r4$wt.diff, 4), '\n')
cat('se.diff   :', round(r4$se.diff, 4), '\n')
cat('CI (95)   : [', round(r4$ci.lower, 4), ',', round(r4$ci.upper, 4), ']\n')
cat('ess1 / ess2:', round(r4$ess1, 1), '/', round(r4$ess2, 1), '\n')

cat('\n========== (5) IPD vs IPD -- binary (sim110$Y.bin) ==========\n')
r5 <- wtTrtDiff(ipd1.te = ipd1$Y.bin, w1 = w.out$ipd1$exm.wts,
                ipd2.te = ipd2$Y.bin, w2 = w.out$ipd2$exm.wts)
cat('wt.y1 (A) :', round(r5$wt.y1, 4), '\n')
cat('wt.y2 (B) :', round(r5$wt.y2, 4), '\n')
cat('wt.diff   :', round(r5$wt.diff, 4), '\n')
cat('se.diff   :', round(r5$se.diff, 4), '\n')
cat('CI (95)   : [', round(r5$ci.lower, 4), ',', round(r5$ci.upper, 4), ']\n')

cat('\n========== (6) algebraic sanity check ==========\n')
## var = sum(w^2) / sum(w)^2 * mean((y - mean(y))^2) -- recompute manually
w        <- w.out$ipd1$exm.wts
y        <- ipd1$Y
manual.v <- (sum(w^2) / sum(w)^2) * mean((y - mean(y))^2)
cat('manual var1 :', signif(manual.v, 6), '\n')
cat('r4$se1^2    :', signif(r4$se1^2, 6), '\n')
stopifnot(abs(manual.v - r4$se1^2) < 1e-12)

cat('\n========== (7) input validation errors ==========\n')
tryCatch(wtTrtDiff(ipd1.te = 1:5, w1 = 1:4),
         error = function(e) cat('length mismatch caught:', conditionMessage(e), '\n'))
tryCatch(wtTrtDiff(ipd1.te = 1:5, w1 = c(1, 2, -1, 1, 1)),
         error = function(e) cat('negative weight caught:', conditionMessage(e), '\n'))
tryCatch(wtTrtDiff(ipd1.te = 1:5, w1 = rep(1, 5)),
         error = function(e) cat('no IPD 2 / AD caught:', conditionMessage(e), '\n'))
tryCatch(wtTrtDiff(ipd1.te = 1:5, w1 = rep(1, 5),
                   ipd2.te = 1:3, w2 = rep(1, 3),
                   ad.mean = 1),
         error = function(e) cat('both IPD 2 & AD caught:', conditionMessage(e), '\n'))
tryCatch(wtTrtDiff(ipd1.te = 1:5, w1 = rep(1, 5),
                   ad.mean = 1, ad.sd = 1),     ## ad.n missing
         error = function(e) cat('ad.sd w/o ad.n caught:', conditionMessage(e), '\n'))

cat('\n========== WTTRTDIFF SMOKE TESTS PASSED ==========\n')
