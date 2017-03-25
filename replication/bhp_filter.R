# R code to take BHP data and filter slightly
#
library("foreign")
bhp <- read.dta("gasoline_demand_BHP2012.dta")
#
Q <- exp(bhp$log_q)
P <- exp(bhp$log_p)
Y <- exp(bhp$log_y)
dist <- exp(bhp$distance_oil1000)
hhsize <- floor(exp(data_new$log_hhsize))
driver <- floor(exp(bhp$log_driver))
#
bhp_new <- data.frame(Q,P,Y,dist,hhsize,driver)
ix <- sum(is.na(dist)==FALSE)
bhp_new <- bhp_new[(1:ix),]
#
write.csv(bhp_new, file = "bhp_new.csv")
