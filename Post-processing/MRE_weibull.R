library("tidyverse")
library("brms")
library("tidybayes")

#Simulate data
n_obs = 30 
x = seq(1, 5)
shape = 25
intercept = 2
beta = 1.2
mu = intercept + beta*x
scale = mu / gamma(1 + 1 /shape)

set.seed(1)

data_gen <- tibble(x = x, mu = mu, shape = shape, scale = scale)
data <- data_gen %>% 
  expand(nesting(x, mu, shape, scale), trial = 1:n_obs) %>% 
  mutate(y = rweibull(n(), shape = shape, scale = scale))


m1 <- brm(data = data, formula = y ~ x,
          family = weibull(link = "identity"), #inits = my_inits,
          iter = 5000, cores = 4, refresh = 500, chains = 1, warmup = 1000)


m2 <- brm(data = data, formula = y ~ factor(x), 
          family = weibull(link = "identity"), #inits = my_inits,
          iter = 5000, cores = 4, refresh = 500, chains = 1, warmup = 1000)
