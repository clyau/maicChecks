## smoke test for the modernized maicChecks R/ functions
## NOTE: this file lives at repo root, gitignored / Rbuildignored.
##
## In v0.3.0 eIPD gained response columns (r.cont, r.bin). The matching
## columns are y1, y2 (named historically -- they play the role of X).
## All MAIC checks/weights require IPD and AD to have the *same* set of
## matching columns in the same order, hence the explicit
## `eIPD[, c('y1', 'y2')]` subsetting below.
suppressMessages(devtools::load_all('.', quiet = TRUE))

ipd.mx <- eIPD[, c('y1', 'y2')]
ad1.mx <- eAD[1, c('y1', 'y2')]
ad3.mx <- eAD[3, c('y1', 'y2')]

cat('\n========== maicLP: scenario A (inside hull) ==========\n')
print(maicLP(ipd.mx, ad1.mx))     # expect lp.check = 0

cat('\n========== maicLP: scenario C (outside hull) ==========\n')
print(maicLP(ipd.mx, ad3.mx))     # expect lp.check = 2

cat('\n========== maicWt: scenario A ==========\n')
m1 <- maicWt(ipd.mx, ad1.mx)
cat('ESS:', m1$ipd.ess, '\n')
cat('weighted means (should equal ad1.mx):\n')
print(m1$ipd.wtsumm)
print(ad1.mx)

cat('\n========== maicWt: scenario C (should stop via maicLP gate) ==========\n')
res <- tryCatch(maicWt(ipd.mx, ad3.mx),
                error = function(e) e)
if (inherits(res, 'error')) {
  cat('Got expected error:\n  ', conditionMessage(res), '\n')
} else {
  cat('!!! UNEXPECTED: maicWt did not stop on infeasible AD\n')
}

cat('\n========== maxessWt: scenario A ==========\n')
m0 <- maxessWt(ipd.mx, ad1.mx)
cat('ESS:', m0$ipd.ess, '\n')
print(m0$ipd.wtsumm)

cat('\n========== maicT2Test: scenario A (returns list now) ==========\n')
t2 <- maicT2Test(ipd.mx, ad1.mx)
cat('returned object:\n')
str(t2)

cat('\n========== maicMD: scenario A ==========\n')
md <- maicMD(ipd.mx, ad1.mx)
cat('md.check =', md$md.check, ' (expect 2 -- A is close to center)\n')
cat('plot object class:', paste(class(md$md.plot), collapse = '/'), '\n')

cat('\n========== maicPCA: scenario A ==========\n')
pca <- maicPCA(ipd.mx, ad1.mx)
cat('pc.check =', pca$pc.check, '\n')
cat('plot object class:', paste(class(pca$pc.dplot), collapse = '/'), '\n')

cat('\n========== exmLP.2ipd: sim110 X1..X5 ==========\n')
ipd1 <- sim110[sim110$study == 'IPD A', ]
ipd2 <- sim110[sim110$study == 'IPD B', ]
x <- exmLP.2ipd(ipd1, ipd2,
                vars_to_match  = paste0('X', 1:5),
                cat_vars_to_01 = paste0('X', 1:3))
cat('lp.check =', x$lp.check, '\n')

cat('\n========== SMOKE TESTS PASSED ==========\n')
