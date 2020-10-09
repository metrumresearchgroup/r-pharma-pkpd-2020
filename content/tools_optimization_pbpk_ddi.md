Estimate parameters in a PBPK model
================
Metrum Research Group

  - [Packages and setup](#packages-and-setup)
  - [Reference](#reference)
  - [Data](#data)
  - [PBPK model: pitavastatin / CsA
    DDI](#pbpk-model-pitavastatin-csa-ddi)
  - [Objective function](#objective-function)
      - [Prediction function](#prediction-function)
  - [Data grooming](#data-grooming)
  - [Optimize](#optimize)
      - [`nloptr::newuoa`: minimization without
        derivatives](#nloptrnewuoa-minimization-without-derivatives)
          - [Get some predictions to look at how the fit
            went](#get-some-predictions-to-look-at-how-the-fit-went)
          - [Make some plots](#make-some-plots)
          - [A nicer plot](#a-nicer-plot)
          - [The final objective function value and
            estimates](#the-final-objective-function-value-and-estimates)

# Packages and setup

``` r
library(tidyverse)
library(PKPDmisc)
library(mrgsolve)
source("script/functions.R")
source("script/global.R")
```

``` r
set.seed(10101)
```

``` r
theme_set(theme_bw() + theme(legend.position = "top"))
scale_colour_discrete <- function(...) scale_color_brewer(palette="Set2")
```

Models are located here:

``` r
model_dir <- "model"
```

# Reference

**Quantitative Analyses of Hepatic OATP-Mediated Interactions Between
Statins and Inhibitors Using PBPK Modeling With a Parameter Optimization
Method**

  - T Yoshikado, K Yoshida, N Kotani, T Nakada, R Asaumi, K Toshimoto, K
    Maeda, H Kusuhara and Y Sugiyama

  - CLINICAL PHARMACOLOGY & THERAPEUTICS | VOLUME 100 NUMBER 5 |
    NOVEMBER 2016

  - <https://www.ncbi.nlm.nih.gov/pubmed/27170342>

# Data

  - Example taken from figure 4a from the publication
  - Using this as example data to fit

<!-- end list -->

``` r
data.file <- "data/fig4a.csv"

data <-
  data.file %>% 
  read_csv() %>% 
  mutate(
    profile = NULL, 
    type=ID, 
    typef=factor(ID, labels = c("Statin", "Statin+CsA")), 
    DV = ifelse(DV==-1, NA_real_, DV)
  )


data
```

    . # A tibble: 23 x 8
    .       ID  time     DV  evid   amt   cmt  type typef     
    .    <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <fct>     
    .  1     2  0     NA        1  2000     2     2 Statin+CsA
    .  2     2  1     NA        1    30     1     2 Statin+CsA
    .  3     2  1.49  73.7      0     0     0     2 Statin+CsA
    .  4     2  1.99 102.       0     0     0     2 Statin+CsA
    .  5     2  2.49  59.9      0     0     0     2 Statin+CsA
    .  6     2  3.00  37.6      0     0     0     2 Statin+CsA
    .  7     2  3.97  15.7      0     0     0     2 Statin+CsA
    .  8     2  5.01   9.24     0     0     0     2 Statin+CsA
    .  9     2  6.99   3.54     0     0     0     2 Statin+CsA
    . 10     2  9.01   2.22     0     0     0     2 Statin+CsA
    . # … with 13 more rows

``` r
data %>% filter(evid==1)
```

    . # A tibble: 3 x 8
    .      ID  time    DV  evid   amt   cmt  type typef     
    .   <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <fct>     
    . 1     2     0    NA     1  2000     2     2 Statin+CsA
    . 2     2     1    NA     1    30     1     2 Statin+CsA
    . 3     1     1    NA     1    30     1     1 Statin

  - The goal is to fit the pitavastatin data either alone or in
    combination with cyclosporin administered 1 hour before the
    pitavastatin

<!-- end list -->

``` r
ggplot(data=data,aes(time,DV)) + 
  geom_point(aes(col = typef), size = 3) + 
  geom_line(col = "darkgrey", aes(group = typef)) + 
  scale_y_continuous(trans="log", limits=c(0.1,300), breaks=logbr()) 
```

<img src="figures/pbpk-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

# PBPK model: pitavastatin / CsA DDI

  - Check out the model / data with a quick simulation

<!-- end list -->

``` r
mod <- mread_cache("yoshikado", model_dir)
```

Make some persistent updates to the model

  - Simulate out to 14 hours
  - Only interested in `CP`, the pitavastatin concentration

<!-- end list -->

``` r
mod <- mod %>% update(end=14, delta=0.1) %>% Req(CP) 
```

A practice simulation

``` r
dose <- filter(data, evid==1) %>% mutate(typef=NULL)

sims <- 
  mod %>% 
  mrgsim_d(dose, obsaug=TRUE) %>% 
  mutate(type = typef(ID))

ggplot(sims, aes(time,CP,col=type)) + 
  geom_line(lwd = 1) + 
  scale_x_continuous(breaks = seq(0,12,2)) + 
  scale_y_log10(name = "Pitavastatin concentration")
```

<img src="figures/pbpk-unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

``` r
sims %>% 
  group_by(type) %>% 
  summarise(auc = auc_partial(time,CP)) %>% 
  mutate(fold_increase = auc /first(auc))
```

    . # A tibble: 2 x 3
    .   type                 auc fold_increase
    .   <fct>              <dbl>         <dbl>
    . 1 Pitavastatin alone  44.1          1   
    . 2 Pitavastatin + CsA 161.           3.65

# Objective function

  - Least squares objective function
  - Weighted by the observations

Arguments:

  - `dv` the observed data
  - `pred` the predicted data

<!-- end list -->

``` r
wss <- function(dv, pred, weight = 1/dv) {
  sum(((dv-pred)*weight)^2,na.rm=TRUE) 
}
```

### Prediction function

  - Let’s go through step by step what each line is doing for us

Arguments:

  - `p` the parameters proposed by the optimizer
  - `.data` the simulation template (doses and observation records)
  - `yobs` a vector of observed data which matches observations in
    `.data`
  - `pred` logical; if `TRUE`, just return predicted data

<!-- end list -->

``` r
sim_ofv <- function(p, data, pred = FALSE) {
  
  names(p) <- names(theta)
  
  p <- lapply(p,exp)
  
  mod <- param(mod, p)
  
  out <- mrgsim_q(mod, data = data, output="df")
  
  if(pred) return(out)
  
  ofv <- wss(data[["DV"]], out[["CP"]])
  
  return(ofv)
  
  #return(-1*sum(dnorm(log(yobs),log(out$CP),.par$sigma,log=TRUE)))
  
}
```

What this function does:

1.  Take in arguments; focus is on a new set of parameters `p` proposed
    by the optimizer; other arguments are just fixed data that we need
2.  Get the parameters out of log scale
3.  Also, put names on the list of parameters; this is crutial
4.  Update the model object with the new parameters
5.  (optionally simulate and return)
6.  Simulate from the data set, taking only observed values
7.  Calculate and return the objective function value

# Data grooming

  - Pick out the observations
  - Drop the non-numeric columns

<!-- end list -->

``` r
data <-  dplyr::select(data, -typef)
```

# Optimize

First, set up the initial estimates

``` r
theta <- c(
  fbCLintall = 1.2, 
  ikiu = 1.2, 
  fbile = 0.9, 
  ka = 0.1, 
  ktr = 0.1
) %>% log()
```

## `nloptr::newuoa`: minimization without derivatives

``` r
fit <- nloptr::newuoa(x0 = theta, fn = sim_ofv, data = data)
```

``` r
fit
```

    . $par
    . [1] -0.20421286 -4.51432097 -1.06749446 -0.01109318 -0.37133662
    . 
    . $value
    . [1] 0.6860763
    . 
    . $iter
    . [1] 351
    . 
    . $convergence
    . [1] 4
    . 
    . $message
    . [1] "NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached."

### Get some predictions to look at how the fit went

Recall that our (transformed) parameters are

``` r
fit$par
```

    . [1] -0.20421286 -4.51432097 -1.06749446 -0.01109318 -0.37133662

We can generate a prediction that matches our data like this

``` r
sim_ofv(fit$par, data = dose, pred = TRUE) %>% filter(time >= 1) %>% head
```

    .   ID time       CP
    . 1  2  1.0  0.00000
    . 2  2  1.0  0.00000
    . 3  2  1.1 18.19345
    . 4  2  1.2 28.55603
    . 5  2  1.3 36.22961
    . 6  2  1.4 42.13621

We can also get the predictions under the initial conditions by passing
in `theta` rather than `fit$par`

In the next block, generate

1.  Predictions with the final estimates
2.  Predications with the initial estimates
3.  Observed data to overlay

<!-- end list -->

``` r
df_pred <- sim_ofv(fit$par, dose, pred=TRUE) %>% mutate(type = typef(ID))
df_init <- sim_ofv(theta,   dose, pred=TRUE) %>% mutate(type = typef(ID))
df_obs <-  mutate(data, type=typef(ID))
```

### Make some plots

``` r
ggplot(df_pred, aes(time,CP)) + 
  geom_line(lwd=1) + 
  geom_point(data = df_obs, aes(time,DV),col="firebrick",size=2) + 
  facet_wrap(~type) + scale_y_log10() 
```

<img src="figures/pbpk-unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

### A nicer plot

``` r
ggplot(data=df_pred) + 
  geom_line(data=df_init,aes(time,CP,lty="A"), col="black", lwd=0.7) +
  geom_line(aes(time,CP,lty="B"),col="black",lwd=0.7) + 
  geom_point(data=df_obs,aes(time,DV,col=type),size=3) + 
  facet_wrap(~type) + 
  scale_y_continuous(trans="log",breaks=10^seq(-4,4), 
                     limits=c(0.1,100),
                     "Pitavastatin concentration (ng/mL)") +
  scale_x_continuous(name="Time (hours)", breaks=seq(0,14,2)) +
  scale_linetype_manual(values= c(2,1), guide = FALSE,
                        labels=c("Initial estimates", "Final estimates"), name="") +
  theme_bw() + theme(legend.position="top") 
```

<img src="figures/pbpk-unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

### The final objective function value and estimates

``` r
sim_ofv(fit$par,data=data)
```

    . [1] 0.6860763

``` r
exp(fit$par)
```

    . [1] 0.81528881 0.01095104 0.34386902 0.98896812 0.68981170
