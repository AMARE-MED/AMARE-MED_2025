
# First run Exercise_E_part1.R

library(dsem)

tsdata = t( n_at )
colnames(tsdata) = paste0("a",ages)
tsdata = ts(log(tsdata))
plot( tsdata )

time_term = "
  a1 -> a2, 1, rho1
  a2 -> a3, 1, rho1
  a3 -> a4, 1, rho1
  a4 -> a5, 1, rho1
  a5 -> a6, 1, rho1
  a6 -> a7, 1, rho1
  a7 -> a8, 1, rho1
  a8 -> a9, 1, rho1
  a9 -> a10, 1, rho1

  a1 <-> a1, 0, sd1
  a2 <-> a2, 0, sd2
  a3 <-> a3, 0, sd2
  a4 <-> a4, 0, sd2
  a5 <-> a5, 0, sd2
  a6 <-> a6, 0, sd2
  a7 <-> a7, 0, sd2
  a8 <-> a8, 0, sd2
  a9 <-> a9, 0, sd2
  a10 <-> a10, 0, sd2
"

# Fit full model
fit = dsem(
  tsdata = tsdata,
  sem = time_term
)

# Drop last 5 years
tsdata2 = tsdata
tsdata2[T-0:4,] = NA
fit2 = dsem(
  tsdata = tsdata2,
  sem = time_term
)

# Make into data frame
tspred = fit2$internal$parhat$x_tj
DF = expand.grid( t = seq_len(T), a = seq_len(A) )
DF = cbind(DF, obs = as.vector(tsdata), pred = as.vector(tspred) )

# Plot observed vs. predicted
library(ggplot2)
ggplot(DF) +
  geom_point( aes(x=t, y=obs), col = "black" ) +
  geom_point( aes(x=t, y=pred), col = "red" ) +
  facet_wrap( vars(a) )

