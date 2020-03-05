library("brms")


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


plot(data$x, data$y)

mypriors <-  c(set_prior("gamma(1,1)", coef = "Intercept"), set_prior("gamma(1,1)", coef = "x"))
mypriors <-  c(set_prior("gamma(1,1)", class = "Intercept"), set_prior("gamma(1,1)", coef = "x"))
myform <- bf(y ~ 1 + x, center = FALSE)
init_f <- function () list(Intercept = rnorm(1, 15, 2.), shape = max(.01,rnorm(1, 3, 1)), scale =max(.01,rnorm(1, 3, 1))  )
my_inits = list(list(b = c(1,1), shape = 1))
my_inits_func = function() list(b = c(max(.5, rnorm(1,3,1)),max(.5, rnorm(1,3,1))), shape = max(.5, rnorm(1,3,1)))
my_inits = list(list(x = 1, shape = 1, Intercept = 1))
m1 <- brm(data = data, formula = myform,
          family = weibull(link = "identity"), prior = mypriors, inits = my_inits_func,
          iter = 5000, cores = 4, refresh = 500, chains = 1, warmup = 1000)

summary(m1)
#setting priors over intercept and beta don't seem to get read of error.
get_prior(formula = myform, data = data, family = weibull(link = "identity"))

#prior_summary(m1)
#m2 <- brm(data = data, formula = y ~ factor(x), 
#          family = weibull(link = "identity"), #inits = my_inits,
#          iter = 5000, cores = 4, refresh = 500, chains = 1, warmup = 1000)
