## smoke test for exmWt.2ipd new modes
suppressMessages(devtools::load_all('.', quiet = TRUE))

ipd1 <- sim110[sim110$study == 'IPD A', ]
ipd2 <- sim110[sim110$study == 'IPD B', ]

cat('\n========== (1) symmetric QP, default ==========\n')
x0 <- exmWt.2ipd(ipd1, ipd2,
                 vars_to_match  = paste0('X', 1:5),
                 cat_vars_to_01 = paste0('X', 1:3))
cat('lp.check =', x0$lp.check, '\n')
cat('method   =', x0$method, '\n')
cat('target.ipd is NULL:', is.null(x0$target.ipd), '\n')
print(x0$wtd.summ)

cat('\n========== (2) target.ipd = "ipd1", method = "maicWt" ==========\n')
x1 <- exmWt.2ipd(ipd1, ipd2,
                 vars_to_match  = paste0('X', 1:5),
                 cat_vars_to_01 = paste0('X', 1:3),
                 target.ipd     = 'ipd1',
                 method         = 'maicWt')
cat('lp.check =', x1$lp.check, '\n')
cat('method   =', x1$method, '\n')
cat('target.ipd =', x1$target.ipd, '\n')
print(x1$wtd.summ)
cat('IPD A weights all 1?', all(x1$ipd1$exm.wts == 1), '\n')
cat('weighted means of IPD B should match colMeans of expanded IPD A...\n')

cat('\n========== (3) target.ipd = "ipd2", method = "maxessWt" ==========\n')
x2 <- exmWt.2ipd(ipd1, ipd2,
                 vars_to_match  = paste0('X', 1:5),
                 cat_vars_to_01 = paste0('X', 1:3),
                 target.ipd     = 'ipd2',
                 method         = 'maxessWt')
cat('lp.check =', x2$lp.check, '\n')
cat('method   =', x2$method, '\n')
cat('target.ipd =', x2$target.ipd, '\n')
print(x2$wtd.summ)
cat('IPD B weights all 1?', all(x2$ipd2$exm.wts == 1), '\n')

cat('\n========== (4) bad target.ipd value should error ==========\n')
res <- tryCatch(exmWt.2ipd(ipd1, ipd2,
                           vars_to_match  = paste0('X', 1:5),
                           cat_vars_to_01 = paste0('X', 1:3),
                           target.ipd     = 'nonsense'),
                error = function(e) e)
cat('Got expected error:', conditionMessage(res), '\n')

cat('\n========== EXMWT.2IPD SMOKE TESTS PASSED ==========\n')
