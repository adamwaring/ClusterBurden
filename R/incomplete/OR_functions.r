# calculate ORs and standard errors
# different intervals possible
calc_or = function(n1, n2, s1, s2) (n1*s2) / (n2*s1)

calc_sd = function(n1, n2, s1, s2) sqrt(rowSums(1/cbind(n1, n2, s1, s2)))

calc_se = function(n1, n2, s1, s2, siglevel) qnorm(1 - siglevel/2) * calc_sd(n1, n2, s1, s2)

calc_lwr = function(or, se) exp(log(or) - se)
calc_upr = function(or, se) exp(log(or) + se)

calc_over_x = function(or, sd_logor, x) 1 - pnorm(x, or, sd_logor)
