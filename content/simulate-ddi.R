library(tidyverse)
library(PKPDmisc)
library(mrgsolve)

mod <- mread("model/yoshikado.cpp", end = 12, delta  = 0.1)

ddi <- c(
  ev(amt = 2000, cmt = 2, time = 0), 
  ev(amt = 30, cmt = 1, time = 1)
)

n <- 2000
idata <- tibble(ikiu = rlnorm(n, log(mod$ikiu),sqrt(0.09)))

out <- mrgsim_ei(mod, events = ddi, idata = idata)

head(out)

summ <- 
  out %>% 
  group_by(ID) %>% 
  summarise(auc = auc_partial(time,CP), .groups = "drop")

ggplot(summ, aes(x = auc)) + geom_histogram(col = "white")


