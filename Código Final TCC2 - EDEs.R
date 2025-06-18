# Load packages
library(crypto2)
library(dplyr)
library(pracma)
library(ggplot2)
library(forecast) 
library(gridExtra)
library(scales)
library(tidyr)
library(DEoptim)
library(lmtest)


#  Data
btc_list <- crypto_list(only_active = TRUE) %>% 
  filter(name == "Bitcoin")
btc <- crypto_history(coin_list=btc_list,
                      start_date="20190101", end_date="20241231")


# Log-prices and log-returns
price <- as.numeric(btc$close)
log_price <- log(price)
ret       <- diff(log_price)   
dt        <- 1              
data <-as.Date(btc$timestamp)

# Plot of closing prices

ggplot(btc, aes(x = data, y = price)) +
  geom_line(linewidth = 0.6, color = "#a11d21") +
  scale_y_continuous(labels = label_comma(big.mark = ".", decimal.mark = ","))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "1 year")+
  labs(x = "Data", y = "Preço (US$)") +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 10),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 9.5),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "top")

ggsave("Preço_tempo.pdf", width = 158, height = 93, units = "mm")

# Gamma function for GOU-FE
Gamma_theta <- function(t, theta1, theta3, theta2) {
  nu0 <- theta1 * theta3 / 2
  nu  <- -nu0^2 + theta2^2 * (1 - theta3)
  if (nu > 0) {
    k1 <- theta2^2*theta3 - 2*theta1*(1-theta3)*nu0
    k2 <- (1/nu)*(theta1*theta2^2 - nu0*(theta2^2*theta3 - theta1*(1-theta3)*2*nu0) - theta1*(1-theta3)*(nu0^2+nu))
    exp((-t*theta1*theta3)/2)*(k1*cos(nu*t)+k2*sin(nu*t))
  } else if (nu < 0) {
    k1 <- theta2^2*theta3 - 2*theta1*(1-theta3)*nu0
    k2 <- (1/nu)*(theta1*theta2^2 - nu0*(theta2^2*theta3 - theta1*(1-theta3)*2*nu0) - theta1*(1-theta3)*(nu0^2-nu))
    exp((-t*theta1*theta3)/2)*(k1*cosh(nu*t)+k2*sinh(nu*t))
  } else {
    k1 <- theta2^2*theta3 - 2*theta1*(1-theta3)*nu0
    k2 <- theta1*theta2^2 - nu0*(theta2^2*theta3 - theta1*(1-theta3)*2*nu0) - theta1*(1-theta3)*(nu0^2)
    exp((-t*theta1*theta3)/2)*(k1 + k2*t)
  }
}

#  Exact transition sampling from CIR process
rCIR <- function(y0, kappa, omega, xi, dt, npart) {
  c0     <- xi^2 * (1 - exp(-kappa*dt)) / (4*kappa)
  d      <- 4*kappa*omega / xi^2
  lambda <- 4*kappa*y0*exp(-kappa*dt) /
    (xi^2*(1 - exp(-kappa*dt)))
  # y/c0 ~ χ²(df=d, ncp=lambda)
  c0 * rchisq(npart, df = d, ncp = lambda)
}


# Simulating FBM via Riemann–Liouville approach
simulate_fbm_mvn <- function(n, H, dt) {
  eps <- rnorm(n, mean = 0, sd = sqrt(dt))
  w  <- sapply(1:n, function(k) k^(H - 0.5) - (k - 1)^(H - 0.5))
  fgn <- numeric(n)
  for (i in 1:n) {
    fgn[i] <- sum( w[1:i] * rev(eps[1:i]) )
  }
  fbm <- c(0, cumsum(fgn))
  dBH <- c(NA, diff(fbm))
  list(fbm = fbm, dBH = dBH)
}

# —————————————————————————————————————————————————————————————
# Log-likelihood of GOUFE-CIR via Particle Filter (Monte Carlo)
# —————————————————————————————————————————————————————————————

loglik_GOUFE_CIR_pf <- function(params, x, dt, npart) {
  θ1 <- params[1]; θ2 <- params[2]; θ3 <- params[3]
  κ  <- params[4]; ω  <- params[5]; ξ  <- params[6]
  H  <- params[7]
  n    <- length(x)
  incs <- diff(x)            
  
  fbm_out <- simulate_fbm_mvn(n, H, dt)
  dBH     <- fbm_out$dBH     
  
 
  Yp   <- rep(ω, npart)   
  w    <- rep(1/npart, npart)
  logL <- 0
  
  for (i in 2:n) {
    j   <- seq_len(i - 1)
    mem <- sum(Gamma_theta(j * dt, θ1, θ3, θ2) * x[i - j]) * dt
    
    # drift
    muX <- ( -θ1*(1-θ3)*x[i-1] - mem ) * dt
    
    # propagate Y through the exact CIR transition
    Yp <- rCIR(Yp, κ, ω, ξ, dt, npart)
    
    # weight update using the density of ΔX_i | Yp
    scale <- x[i - 1] * sqrt(Yp) * abs(dBH[i]) 
    
    dens  <- dnorm(incs[i - 1], mean = muX, sd = scale)
    
    w      <- w * dens
    sw     <- sum(w)
    if (sw == 0) return(1e+200)
    
    # accumulate log-weight and resample
    logL <- logL + log(sw)
    w    <- w / sw
    idx  <- sample.int(npart, npart, TRUE, prob = w)
    Yp   <- Yp[idx]
    w    <- rep(1/npart, npart)
  }
  
  return(-logL)
}

# —————————————————————————————————————————————————————————————
# Log-likelihood of the GOU-FE model with constant sigma
# —————————————————————————————————————————————————————————————

loglik_GOUFE_CONST <- function(params, x, dt) {
  θ1    <- params[1]
  θ2    <- params[2]
  θ3    <- params[3]
  sigma <- params[4]    
  H     <- params[5]
  
  n    <- length(x)
  incs <- diff(x) 
  
  fbm_out <- simulate_fbm_mvn(n, H, dt)
  dBH     <- fbm_out$dBH
  
  logL <- 0
  
  for (i in 2:n) {
    # memory
    j   <- seq_len(i - 1)
    mem <- sum( Gamma_theta(j * dt, θ1, θ3, θ2) * x[i - j] ) * dt
    
    # drift
    muX <- ( -θ1 * (1 - θ3) * price[i - 1] - mem ) * dt
    
    # scale
    scale <- price[i - 1] * sigma * abs(dBH[i])
    
    logL <- logL +
      dnorm(incs[i - 1], mean = muX, sd = scale, log = TRUE)
  }
  
  -logL
}
# ————————————————————————————————————————————————————————
# Log-likelihood of the GFBM-CIR model via Particle Filter
# ————————————————————————————————————————————————————————

loglik_GFBM_CIR_pf <- function(params, x, dt, npart) {
  μ     <- params[1]
  κ     <- params[2]
  ω     <- params[3]
  ξ     <- params[4]
  H     <- params[5]
  incs  <- diff(x)        
  n     <- length(x)
  
  # initial particles at ω
  Yp    <- rep(ω, npart)
  w     <- rep(1/npart, npart)
  logL  <- 0
  
  # simulate FBM
  fbm_out <- simulate_fbm_mvn(n, H, dt)
  dBH     <- fbm_out$dBH
  
  for (i in 2:n) {
    # latent state prediction via CIR
    Ynew <- rCIR(Yp, κ, ω, ξ, dt, npart)
    
    # conditionally normal price increments:
    # ΔS_i | Ynew ~ N(μ * S_{i-1} * dt,  (sqrt(Ynew)*S_{i-1} * |dBH[i]|)^2 )
    sd_inc <- sqrt(Ynew) * x[i-1] * abs(dBH[i])
    wi     <- dnorm(incs[i-1], mean = μ * x[i-1] * dt, sd = sd_inc) 
    
    # weight update
    w    <- w * wi
    sw   <- sum(w)
    if (sw == 0) return(1e+200)
    logL <- logL + log(sw)
    
    # normalize and resample
    w  <- w / sw
    idx <- sample.int(npart, npart, replace = TRUE, prob = w)
    Yp <- Ynew[idx]
    w  <- rep(1/npart, npart)
  }
  
  -logL
}

# ————————————————————————————————————————————————————————
# Log-likelihood of the GFBM model with constant sigma
# ————————————————————————————————————————————————————————

loglik_GFBM_CONST<- function(params, x, dt) {
  μ     <- params[1]             
  sigma <- params[2]
  H     <- params[3]             
  n     <- length(x)
  incs  <- diff(x)              
  Sprev <- x[-n]                
  
  # simulate FBM
  fbm_out <- simulate_fbm_mvn(n, H, dt)
  dBH     <- fbm_out$dBH[2:n]     
  
  # conditional standard deviation:
  # ΔS_i | params ~ N( μ·Sprev·dt, (σ_fixed·Sprev·|ΔB^H_i|)^2 )
  sd_inc  <- sigma * Sprev * abs(dBH)
  # sum of the log-densities
  ll <- sum(dnorm(incs, mean = μ * Sprev * dt, sd   = sd_inc, log  = TRUE))

  -ll
}

# ——————————————————
# Initial estimates
# ——————————————————

mu_init    <- mean(ret)
sq_ret    <- ret^2 
ar1       <- ar(sq_ret, aic = FALSE, order.max = 1)
phi       <- ar1$ar             
kappa_init<- min(max(-log(phi) / dt, 0.1), 5)
omega_init<- mean(sq_ret)        

# simple returns R_i
R      <- diff(price) / price[-1]

# mean of returns
R_bar  <- mean(R)

n         <- length(price)
sigma_S <- sqrt( sum((R - R_bar)^2) / ((n-1)*dt) )
sigma_log <- sqrt(sum(ret^2) / ((n-1)*dt))
sigma_init <- sd(ret)
xi_init   <- min(max(sqrt(2 * kappa_init * var(sq_ret) / omega_init), 0.1), 5)

#GOUFE-CIR
H0    <- hurstexp(ret)[[1]]
init_GOUFE_CIR <- c(θ1 = 0.01, θ2 = 1e-4,θ3 = 0.9, κ=kappa_init, ω = omega_init, ξ=xi_init,H= H0)

set.seed(124)

res_GOUFE_CIR <- optim(
  par    = init_GOUFE_CIR,
  fn     = loglik_GOUFE_CIR_pf,
  x      = price,
  dt     = dt,
  npart  = 200,
  method = "Nelder-Mead")

res_GOUFE_CIR

#GOUFE-CONST

init_GOUFE_CONST <- c(θ1 = 0.01, θ2 = 1e-4,θ3 = 0.9, ξ=xi_init,H= H0)

lower <- c(
  theta1 = 1e-2,
  theta2 = 1e-3,
  theta3 = 1e-1,
  sigma  =  1e-1,
  H      = 0.54
)
upper <- c(
  theta1 = 1,
  theta2 = 1,
  theta3 = 1,
  sigma  = 0.8,
  H      = 0.9
)

set.seed(124)

res_GOUFE_CONST <- optim(
  par    = init_GOUFE_CONST,
  fn     = loglik_GOUFE_CONST,
  x      = price,
  dt     = dt,
  lower = lower,
  upper = upper,
  method = "L-BFGS-B")

res_GOUFE_CONST

#GFBM-CIR

init_pars_GFBM_CIR <- c(μ    = 0.001, κ = kappa_init,
                        ω = omega_init, ξ = xi_init, H = H0)

lower <- c(
  μ = 1e-4,
  κ = 1,
  ω = 1e-1,
  ξ  = 1e-1,
  H      = 0.53
)
upper <- c(
  μ = 1,
  κ = 3,
  ω = 1,
  ξ  =  0.8,
  H      = 0.6
)

set.seed(124)

res_GFBM_CIR <- optim(
  par    = init_pars_GFBM_CIR,
  fn     = loglik_GFBM_CIR_pf,
  x      = price,
  dt     = dt,
  npart  = 500,
  lower = lower,
  upper = upper,
  method = "L-BFGS-B")

res_GFBM_CIR

# GFBM-CONST

init_GFBM_CONST <- c(
  μ    = unname(res_GFBM_CIR$par["μ"]),
  ξ     = unname(res_GFBM_CIR$par["ξ"]),
  H     = H0)

mu_lower <- 1e-4
mu_upper <- 1
sigma_lower <- 1e-1
sigma_upper <- 0.7
lower_GFBM_CONST <- c(mu = mu_lower, sigma_lower, H = 0.53)
upper_GFBM_CONST <- c(mu = mu_upper, sigma_upper, H = 0.9)

set.seed(124)

res_GFBM_CONST<- optim(
  par    = init_GFBM_CONST,
  fn     = loglik_GFBM_CONST,
  x      = price,
  dt     = dt,
  lower = lower_GFBM_CONST,
  upper = upper_GFBM_CONST,
  method = "L-BFGS-B")

res_GFBM_CONST

# ——————————————————————————————
# Results of the Model Estimates
# ——————————————————————————————

logLik_GOUFE_CIR   <- res_GOUFE_CIR$value
logLik_GOUFE_CONST <- res_GOUFE_CONST$value
logLik_GFBM_CIR    <- res_GFBM_CIR$value
logLik_GFBM_CONST  <- res_GFBM_CONST$value

pars_GOUFE_CIR   <- res_GOUFE_CIR$par
pars_GOUFE_CONST <- res_GOUFE_CONST$par
pars_GFBM_CIR    <- res_GFBM_CIR$par
pars_GFBM_CONST  <- res_GFBM_CONST$par

all_param_names <- unique(c(
  names(pars_GOUFE_CIR),
  names(pars_GOUFE_CONST),
  names(pars_GFBM_CIR),
  names(pars_GFBM_CONST)
))

# Function to create the results in a single data frame
make_result_row <- function(model_name, pars, logLik) {
  full_pars <- setNames(rep(NA_real_, length(all_param_names)), all_param_names)
  full_pars[names(pars)] <- pars
  data.frame(
    Modelo = model_name,
    LogLik = logLik,
    as.list(full_pars),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Results
df_estimates <- bind_rows(
  make_result_row("GOUFE_CIR",   pars_GOUFE_CIR,   logLik_GOUFE_CIR),
  make_result_row("GOUFE_CONST", pars_GOUFE_CONST, logLik_GOUFE_CONST),
  make_result_row("GFBM_CIR",    pars_GFBM_CIR,    logLik_GFBM_CIR),
  make_result_row("GFBM_CONST",  pars_GFBM_CONST,  logLik_GFBM_CONST)
)

df_estimates

# ——————————————————————————————————————————————————————
# One Step Ahead  forecast of the GOUFE-CIR model
# ——————————————————————————————————————————————————————

# Extract estimated parameters
pars_GOUFE_CIR   <- res_GOUFE_CIR$par   # θ1, θ2, θ3, κ, ω, ξ, H
pars_GOUFE_CONST <- res_GOUFE_CONST$par # θ1, θ2, θ3, σ, H


# Initialize prediction vectors
n <- length(price)
pred_GOUFE_CIR   <- numeric(n);   pred_GOUFE_CIR[1]   <- price[1]
pred_GOUFE_CONST <- numeric(n);   pred_GOUFE_CONST[1] <- price[1]

# One‐step‐ahead para GOUFE-CIR

npart = 200

Yp <- rep(pars_GOUFE_CIR[5],npart)  # Y₀ = ω
for (i in 2:n) {
  # memory
  j   <- seq_len(i - 1)
  mem <- sum(
    Gamma_theta(j * dt,
                theta1 = pars_GOUFE_CIR[1],
                theta3 = pars_GOUFE_CIR[3],
                theta2 = pars_GOUFE_CIR[2]) *
      price[i - j]
  ) * dt
  
  # drift
  muX <- (-pars_GOUFE_CIR[1] * (1 - pars_GOUFE_CIR[3]) * price[i - 1] - mem) * dt
  
  # propagate the CIR
  Yp <- rCIR(Yp,
             kappa = pars_GOUFE_CIR[4],
             omega = pars_GOUFE_CIR[5],
             xi    = pars_GOUFE_CIR[6],
             dt    = dt,
             npart = npart)
  
  # prediction of X_i
  pred_GOUFE_CIR[i] <- price[i - 1] + muX
}

# ——————————————————————————————————————————————————————
# One Step Ahead forecast of the GOUFE-CONST model
# ——————————————————————————————————————————————————————

for (i in 2:n) {
  # memory
  j   <- seq_len(i - 1)
  mem <- sum(
    Gamma_theta(j * dt,
                theta1 = pars_GOUFE_CONST[1],
                theta3 = pars_GOUFE_CONST[3],
                theta2 = pars_GOUFE_CONST[2]) *
      price[i - j]
  ) * dt
  
  # drift
  muX <- (-pars_GOUFE_CONST[1] * (1 - pars_GOUFE_CONST[3]) * price[i - 1] - mem) * dt
  
  # prediction of X_i
  pred_GOUFE_CONST[i] <- price[i - 1] + muX
}


# ——————————————————————————————————————————————————————————————
# One-Step Ahead forecast of the GFBM-CIR and GFBM-CONST models
# ——————————————————————————————————————————————————————————————

# Extract estimated parameters
mu_GFBM_CIR     <- res_GFBM_CIR$par["μ"]
mu_GFBM_CONST  <- res_GFBM_CONST$par["μ"]

# Initialize prediction vectors
pred_GFBM_CIR     <- numeric(n)
pred_GFBM_CONST  <- numeric(n)

# first forecast equal to the first observed value
pred_GFBM_CIR [1]    <- price[1]
pred_GFBM_CONST[1] <- price[1]


# one-step-ahead using the last observed price
for (i in 2:n) {
  # GFBM-CIR
  pred_GFBM_CIR [i]    <- price[i-1] + mu_GFBM_CIR  * price[i-1]    * dt
  
  # GFBM-CONST
  pred_GFBM_CONST [i] <- price[i-1] + mu_GFBM_CONST * price[i-1] * dt
}


# ——————————————————————————————————————————————————————
# Plot of the GOUFE and GFBM models
# ——————————————————————————————————————————————————————

## Plot of observed vs. predicted prices
df <- data.frame(
  Date     = data,
  Observed = price,
  Pred_GFBM_cir     = pred_GFBM_CIR,
  Pred_GFBM_const=  pred_GFBM_CONST,
  Pred_GOUFE_cir = pred_GOUFE_CIR,
  Pred_GOUFE_const = pred_GOUFE_CONST
)


model_labels <- c(
  Pred_GFBM_cir        = "GFBM–CIR",
  Pred_GFBM_const  = "GFBM–CONST",
  Pred_GOUFE_cir       = "GOUFE–CIR",
  Pred_GOUFE_const = "GOUFE–CONST"
)

df_long <- df %>%
  pivot_longer(cols = starts_with("Pred"),
               names_to = "Model", values_to = "Pred") %>%
  mutate(Model = model_labels[Model])


ggplot(df_long, aes(Date)) +
  geom_line(aes(y = Observed,    color = "Observado"),
            linewidth = 0.6, linetype = "solid") +
  geom_line(aes(y = Pred,         color = Model),
            linewidth = 0.8, linetype = "dashed") +
  scale_y_continuous(
    labels = label_number(
      scale = 1e-3,         
      suffix = "",            
      big.mark = ".", 
      decimal.mark = ","
    ))+
  scale_x_date(date_labels = "%b\n%Y")+
  scale_color_manual(values = c(
    `Observado`        = "#000000",
    `GFBM–CIR`        = "#E41A1C",
    `GOUFE–CIR`       = "#377EB8",
    `GFBM–CONST`  = "#4DAF4A",
    `GOUFE–CONST` = "#FF7F00"
  )) +
  facet_wrap(~ Model, ncol = 2,labeller = as_labeller(model_labels)) +
  labs(x = "Data", y = "Preço (mil US$)", color = "") +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 9.5),
    axis.title.x = element_text(colour = "black", size = 9.5),
    legend.text = element_text(size = 9.5),
    axis.text = element_text(colour = "black", size = 9.5),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position = "top")

ggsave("Observado_vs_Predito.pdf", width = 170, height = 93, units = "mm")

# Calculation of errors
e_GOUFE_CIR   <- price - pred_GOUFE_CIR
e_GOUFE_CONST  <- price - pred_GOUFE_CONST
e_GFBM_CIR   <- price - pred_GFBM_CIR
e_GFBM_CONST  <- price - pred_GFBM_CONST

std_GOUFE_CIR   <- as.numeric(scale(e_GOUFE_CIR[-1]))
std_GOUFE_CONST  <- as.numeric(scale(e_GOUFE_CONST[-1]))
std_GFBM_CIR   <- as.numeric(scale(e_GFBM_CIR[-1]))
std_GFBM_CONST  <- as.numeric(scale(e_GFBM_CONST[-1]))

# Function to calculate metrics
calc_metrics <- function(res, k, logL){
  n <- length(res)
  AIC  <- 2*k - 2*logL
  BIC  <- log(n)*k - 2*logL
  EDC  <- log(log(n))*k - 2*logL
  RMSE <- sqrt(mean(res^2))
  MAE  <- mean(abs(res))
  MAPE <- mean(abs(res / lag(price)[-1]))
  BIAS <- mean(res)
  y    <- price[-1]
  SSR  <- sum(res^2)
  SST  <- sum((y - mean(y))^2)
  R2   <- 1 - SSR/SST
  data.frame(AIC, BIC, EDC, RMSE, MAE, MAPE, BIAS, R2)
}

# extract log-likelihoods and number of parameters
k_GOUFE_CIR   <- length(pred_GOUFE_CIR);    logL_GOUFE_CIR   <- -res_GOUFE_CIR$value
k_GOUFE_CONST  <- length(pred_GOUFE_CONST);   logL_GOUFE_CONST  <- -res_GOUFE_CONST$value
k_GFBM_CIR   <- length(pred_GFBM_CIR);    logL_GFBM_CIR   <- -res_GFBM_CIR$value
k_GFBM_CONST  <- length(pred_GFBM_CONST);   logL_GFBM_CONST  <- -res_GFBM_CONST$value

met_GOUFE_CIR   <- calc_metrics(e_GOUFE_CIR[-1],   k_GOUFE_CIR,  logL_GOUFE_CIR)
met_GOUFE_CONST  <- calc_metrics(e_GOUFE_CONST[-1],  k_GOUFE_CONST, logL_GOUFE_CONST)
met_GFBM_CIR   <- calc_metrics(e_GFBM_CIR[-1],   k_GFBM_CIR ,  logL_GFBM_CIR)
met_GFBM_CONST  <- calc_metrics(e_GFBM_CONST[-1],  k_GFBM_CONST, logL_GFBM_CONST)

# Results

all_metrics <- bind_rows(
  GOUFE_CIR   = met_GOUFE_CIR,
  GOUFE_CONST = met_GOUFE_CONST,
  GFBM_CIR    = met_GFBM_CIR,
  GFBM_CONST  = met_GFBM_CONST,
  .id = "Model"
)
all_metrics


# Descriptive Quantile Analysis
resid_stats <- df_long %>%
  mutate(Residuals = Pred - Observed) %>%
  group_by(Model) %>%
  summarise(
    Min   = min(Residuals),
    Q1    = quantile(Residuals, 0.25),
    Med   = median(Residuals),
    Q3    = quantile(Residuals, 0.75),
    Max   = max(Residuals),
    SD    = sd(Residuals)
  )

resid_stats

# Boxplot of Residuals

ggplot(df_long, aes(x = Model, y = Pred - Observed)) +
  geom_boxplot(fill = "#a11d21", width = 0.5) +
  stat_summary(fun = "mean", geom = "point", shape = 22, size = 2, fill = "white") +
  labs(
    x     = "Modelo",
    y     = "Resíduo"
  ) +
  theme_bw()+
  theme(axis.title.y=element_text(colour="black", size=12),
        axis.title.x = element_text(colour="black", size=12),
        axis.text = element_text(colour = "black", size=9.5),
        axis.text.x =element_text(colour = "black", size=9.5),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave("Boxplot_Residuos.pdf", width = 158, height = 93, units = "mm")

# Function to extract ACF into a data frame
get_acf_df <- function(residuals, model, lag.max = 60) {
  acf_obj <- acf(residuals, plot = FALSE, lag.max = lag.max)
  data.frame(
    lag   = acf_obj$lag[-1],
    acf   = acf_obj$acf[-1],
    model = model
  )
}

# Create a data frame with the ACFs of the four models
df_acf <- bind_rows(
  get_acf_df(std_GFBM_CIR,  "GFBM-CIR"),
  get_acf_df(std_GOUFE_CIR, "GOUFE-CIR"),
  get_acf_df(std_GOUFE_CONST,"GFBM-CONST"),
  get_acf_df(std_GFBM_CONST,"GOUFE-CONST")
)

# Confidence level for the bands (±1.96/√n)
n_obs <- length(std_GFBM_CIR) 
ci    <- qnorm(0.975) / sqrt(n_obs)

# ACF plot
ggplot(df_acf, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_segment(aes(xend = lag, yend = 0)) +
  geom_hline(yintercept =  c(-ci, ci),
             linetype   = "dashed",
             color      = "blue") +
  facet_wrap(~ model, ncol = 2) +
  labs(
    x     = "Lag",
    y     = "ACF"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 10),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 9.5),
    axis.line = element_line(colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("ACF.pdf", width = 158, height = 93, units = "mm")

# Residuals Scatter Plot
df_scatter <- data.frame(
  price    = rep(price, times = 4),
  residual = c(e_GFBM_CIR, e_GFBM_CONST, e_GOUFE_CIR, e_GOUFE_CONST),
  model    = rep(
    c("GFBM-CIR", "GFBM-CONST", "GOUFE-CIR", "GOUFE-CONST"),
    each = length(price)
  )
)

ggplot(df_scatter, aes(x = price, y = residual)) +
  geom_point(alpha = 0.5, size = 1.5,color = "black") +
  facet_wrap(~ model, ncol = 2, scales = "free") +
  labs(
    x = "Preço Observado",
    y = "Resíduo"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 10),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 9.5),
    axis.line = element_line(colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("Dispersao_res.pdf", width = 158, height = 93, units = "mm")

# QQnorm Plot of the Residuals
df_qq <- data.frame(
  residual = c(std_GFBM_CIR, std_GFBM_CONST, std_GOUFE_CIR, std_GOUFE_CONST),
  model = factor(rep(
    c("GFBM-CIR", "GFBM-CONST", "GOUFE-CIR", "GOUFE-CONST"),
    times = c(
      length(std_GFBM_CIR),
      length(std_GFBM_CONST),
      length(std_GOUFE_CIR),
      length(std_GOUFE_CONST)
    )
  ),
  levels = c("GFBM-CIR", "GFBM-CONST", "GOUFE-CIR", "GOUFE-CONST"))
)

ggplot(df_qq, aes(sample = residual)) +
  stat_qq(shape = 1) +  
  stat_qq_line(color = "red", linewidth = 0.7) +
  facet_wrap(~ model, ncol = 2) +
  labs(
    x = "Quantis Teóricos N(0,1)",
    y = "Quantis dos Resíduos Padronizados"
  ) +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 10),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 9.5),
    axis.line = element_line(colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
ggsave("Residuos_Padronizados.pdf", width = 158, height = 93, units = "mm")


# Ljung-Box test for autocorrelation in the residuals
n       <- length(std_GFBM_CIR)
lag_opt <- 60

tests <- list(
  GFBM_CIR    = Box.test(std_GFBM_CIR,    lag = lag_opt, type = "Ljung-Box"),
  GFBM_CONST  = Box.test(std_GFBM_CONST,  lag = lag_opt, type = "Ljung-Box"),
  GOUFE_CIR   = Box.test(std_GOUFE_CIR,   lag = lag_opt, type = "Ljung-Box"),
  GOUFE_CONST = Box.test(std_GOUFE_CONST, lag = lag_opt, type = "Ljung-Box")
)

# Creating the data frame to extract the results
df_box <- data.frame(
  Modelo      = names(tests),
  Estatística = sapply(tests, function(x) unname(x$statistic)),
  DF          = sapply(tests, function(x) unname(x$parameter)),
  P_Value     = sapply(tests, function(x) x$p.value),
  row.names   = NULL,
  stringsAsFactors = FALSE
)

df_box

# Normality test (Shapiro-Wilk)
res_list <- list(
  GFBM_CIR    = std_GFBM_CIR,
  GFBM_CONST  = std_GFBM_CONST,
  GOUFE_CIR   = std_GOUFE_CIR,
  GOUFE_CONST = std_GOUFE_CONST
)
tests_shapiro <- lapply(res_list, shapiro.test)

# Generating the results in the data frame
df_shapiro <- data.frame(
  Modelo   = names(tests_shapiro),
  W        = sapply(tests_shapiro, function(t) unname(t$statistic)),
  p_value  = sapply(tests_shapiro, function(t) t$p.value),
  row.names = NULL,
  stringsAsFactors = FALSE
)

df_shapiro

# Expected Shortfall
res_list <- list(
  "GOUFE–CIR"   = e_GOUFE_CIR[-1],
  "GOUFE–CONST" = e_GOUFE_CONST[-1],
  "GFBM–CIR"    = e_GFBM_CIR[-1],
  "GFBM–CONST"  = e_GFBM_CONST[-1]
)

# Create a data frame of standard deviations
sd_df <- tibble(
  Modelo = names(res_list),
  sd_res = vapply(res_list, sd, numeric(1), na.rm = TRUE)
)

# Expected Shortfall (95%)
es_95 <- lapply(res_list, function(err) {
  alpha <- 0.05
  var_a <- quantile(err, probs = alpha, na.rm = TRUE)
  mean(err[err <= var_a], na.rm = TRUE)
})

es_df <- data.frame(
  Modelo = names(es_95),
  ES_95  = unlist(es_95)
)

es_df
