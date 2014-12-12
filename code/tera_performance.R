###########
## This script will evaluate the performance of tera() versus its parallel counterpart tera_parallel() 

source("tera_parallel.R")
source("tera.R")

M <- matrix(rnorm((1E+2)^2), ncol = 1E+2 ); diag(M) <- 0
s_min <- seq(0, 1, length.out = 100)
country.vec <- as.factor(rep(LETTERS[1:10], 10))

library(microbenchmark)
compare <- microbenchmark(tera(M, s_min), tera.parallel(M, s_min), times = 100)

library(ggplot2)
autoplot(compare)
 
# Rprof(tmp <- tempfile())
# tera(M.NT, s_min, export = FALSE, path = paste(path.g,"NT/", sep = ""))
# Rprof()
# summaryRprof(tmp)
# system.time(tera(M.NT, s_min, export = FALSE, path = paste(path.g,"NT/", sep = "")))
# system.time(tera.parallel(M.NT, s_min, export = FALSE, path = paste(path.g,"NT/", sep = "")))
