library(tidyverse)


V <- seq(-70, 50, by = 2)
# t <- seq()
# Vt <- 

ass_fn <- function(V) {return(1/(1+exp(-(V+22.5)/7.7)))}
ass <- ass_fn(V)
qplot(V, ass)

iss_fn <- function(V) {return(1/(1+exp((V+4.52)/5.7)))}
tau_iur_fn <- function(V) {return(4995.4-170/(1+exp((V+76.1)/5.7)))}
iur_fn <- function(t, V) {
  iss <- iss_fn(V)
  tau_iur <- tau_iur_fn(V)
  return(iss - (iss - 0.998543)*exp(-t/tau_iur))
}
iss <- iss_fn(V)
qplot(V, iss)
tau_iur <- tau_iur_fn(V)
qplot(V, tau_iur)
iur <- 

tau_tis_fn <- function(V) {return(270+1050/(1+exp((V+45.2)/5.7)))}
tau_tis <- tau_tis_fn(V)
qplot(V, tau_tis)

tau_tas_fn <- function(V) {return(0.493*exp(-0.0628*V)+2.058)}
tau_tas <- tau_tas_fn(V)
qplot(V, tau_tas)

alpha_a_fn <- function(V) {return(0.18064*exp(0.03577*(V+70)))}
alpha_a <- alpha_a_fn(V)
qplot(V, alpha_a)

beta_a_fn <- function(V) {return(0.3956*exp(-0.06237*(V+30)))}
beta_a <- beta_a_fn(V)
qplot(V, beta_a)

alpha_i_fn <- function(V) {return((0.000152*exp(-(V+13.5)/7.0))/(0.067083*exp(-(V+33.5)/7.0)+1))}
alpha_i <- alpha_i_fn(V)
qplot(V, alpha_i)

beta_i_fn <- function(V) {return((0.00095*exp((V+70.5)/7.0))/(0.051335*exp((V+70.5)/7.0)+1))}
beta_i <- beta_i_fn(V)
plot(V, beta_i)

## Dongping -----
std_ss <- function(V, x1, x2) {return(1/(1+exp(-(V-x1)/x2)))}
ext_ss <- function(V, x1, x2, x3, x4) {
  oper1 <- (0.21)/(1+exp(-(V-x1)/x2))
  oper2 <- (0.79)/(1+exp(-(V-x3)/x3))
  return(oper1 + oper2)
}

# WT
ass <- std_ss(V, -0.3, 12)
iss <- std_ss(V, -37.40, -4.10)
ars <- std_ss(V, -14.70, 8.20)
irs <- ext_ss(V, -81.46, -5.98, -34.31, -4.31)

ass <- std_ss(V, 7.40, 14)
iss <- std_ss(V, -36.20, -4.90)

qplot(V, ass)
qplot(V, iss)
qplot(V, ars)
qplot(V, irs)

# ST3Gal4
