
###########################################################
## Packages
###########################################################

library(ACDm)
library(zipfR)
library(tidyverse)
library(dplyr)
library(zoo)
library(ggplot2)
###########################################################
## Dados
###########################################################

btc ="BTCUSD_Ticks_20.01.2025-20.01.2025.csv"

## data
mydata = read.csv(file=btc, header=T)
mydata <- mydata
head(mydata)


## Rename and create variables
mydata$time = mydata$Local.time
mydata$price = (mydata$Ask+mydata$Bid)/2
mydata$volume = max(mydata$AskVolume, mydata$BidVolume)

## Remove un-necessary variables
mydata$Local.time = NULL
mydata$Ask = NULL
mydata$Bid = NULL
mydata$AskVolume = NULL
mydata$BidVolume = NULL



## date/time has GMT so remove
mydata$time = as.factor(substr(strptime(mydata$time, "%d.%m.%Y %H:%M:%S"),1,25))
mydata <- mydata[-c(1, 2), ]
head(mydata)
tail(mydata)

# ================================================
# Sensitivity Tests
# ================================================

#Testing different values of priceDiff (0.05%, 0.075%, 0.10%, 0.15%)


mydata2 = mydata

priceDiff_values <- c(52.37, 78.56213, 104.75, 157.1243)
sensitivity_results <- list()
durations_list <- list()

# Loop to compute durations

for (th in priceDiff_values) {
  durations <- computeDurations(mydata2, 
                                open = "0:30:00", close = "23:30:00",
                                rm0dur = TRUE, type = "price", 
                                priceDiff = th)
  
  key <- paste0("priceDiff_", round(th, 2))
  sensitivity_results[[key]] <- summary(durations)
  # Store durations in a list
  durations_list[[key]] <- durations
}

# Convert the list of durations into a data frame
durations_df <- do.call(rbind, lapply(names(durations_list), function(name) {
  data.frame(priceDiff = name, duration = durations_list[[name]])
}))

priceDiff_labels <- paste0("priceDiff_", round(sort(priceDiff_values), 2))

# Convert priceDiff to an ordered factor with the correct levels
durations_df$priceDiff <- factor(durations_df$priceDiff, levels = priceDiff_labels)

# Compute standard deviation for each priceDiff
std_dev_results <- durations_df %>%
  group_by(priceDiff) %>%
  summarise(std_dev = sd(duration.durations))

std_dev_results

# Boxplot
ggplot(durations_df, aes(x = priceDiff, y = duration.durations)) +
  geom_boxplot(fill = "#a11d21", width = 0.5) +
  stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white")+
  labs( x = "Diferença de Preço Relativo (%)",y = "Duração (segundos)") +
  scale_x_discrete(labels = c("0,05%", "0,075%", "0,10%", "0,15%"))+
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 9.5),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "top")

ggsave("boxplot_durations.pdf", width = 158, height = 93, units = "mm")


## Duration data (cutting of first and last 30 mins
## These are price duration with a min threshold of 0.01
## only if type = "price". Price durtions are (here) defind as the 
## duration until the price has changed by at least 'priceDiff' 
## in absolute value.


relative_pct <- 0.1 / 100  # 0,1%
avg_price <- mean(mydata$price, na.rm = TRUE)
priceDiff_relative <- avg_price * relative_pct


# Compute durations using a relative threshold

durDataRelative <- computeDurations(mydata, 
                                    open = "0:30:00", close = "23:30:00",
                                    rm0dur = TRUE, type = "price", 
                                    priceDiff = priceDiff_relative)


durDataShort <- computeDurations(mydata, 
                                 open = "0:30:00", close = "23:30:00",
                                 rm0dur = TRUE, type = "price", 
                                 priceDiff = priceDiff_relative)

## Results
str(durDataShort)
head(durDataShort)

ggplot(durDataShort, aes(x = time, y = durations)) +
  geom_linerange(aes(ymin = 0, ymax = durations), color = "black") +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 9.5),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "top") +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "2 hours")+
  labs(x = "Tempo", y = "Duração")

ggsave("duration_short.pdf", width = 158, height = 93, units = "mm")

## Perform diurnal transform

durDataShortAdj = diurnalAdj(durDataShort, aggregation = "none", method = "supsmu")

ggsave("peform_diurnal.pdf", width = 158, height = 93, units = "mm")

## Results

ggplot(durDataShortAdj, aes(x = adjDur)) +
  geom_histogram(binwidth = 1, fill = "#a11d21", color = "black") + 
  labs(x = "Duração Ajustada",y = "Frequência") +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 10),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 9.5),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "top")

ggsave("duration_adj.pdf", width = 158, height = 93, units = "mm")

mydata <- as.vector(durDataShortAdj$adjDur)

length(mydata)

####################################################
## Descriptive statistics
####################################################

descriptiveSummary <- function(x){
  n            <- length(x)
  skewness     <- (1 / n) * sum(((x - mean(x)) / sd(x)) ^ 3)
  kurtosis     <- ((1 / n) * sum(((x - mean(x)) / sd(x)) ^ 4) - 3 ) 
  cVariation   <- (sd(x) / mean(x))
  statistics   <- list(minimum              = round(min(x), 3),
                       median               = round(median(x), 3),
                       mean                 = round(mean(x), 3),
                       maximum              = round(max(x), 3),
                       standardDeviation    = round(sd(x), 3),
                       coefficientVariation = round(cVariation * 100, 3),
                       coefficientSkewness  = round(skewness, 3),
                       coefficientKurtosis  = round(kurtosis, 3),
                       range                = round(max(x) - min(x), 3)
  )
  return(statistics)
}


descriptiveSummary(mydata)
sd(mydata)

####################################################
## Estimation - Exponential model
####################################################

## ACDm fit (initial values)
fitModelEXP <- acdFit(durations = mydata, model = "LACD1",
                      dist = "exponential", order = c(1,1))

omegahatEXP  <-  as.numeric(fitModelEXP$mPar[1])
alpha1hatEXP <-  as.numeric(fitModelEXP$mPar[3])
beta1hatEXO  <-  as.numeric(fitModelEXP$mPar[2])


residualsEXP <- fitModelEXP$residuals
muhatEXP     <- mydata/residualsEXP

fitModelEXP

fitted_survivalEXP <- function(x, muhat){
  
  survival   <- 1.0 - pexp(x, rate = muhat)
  
  return(survival)
}

## Estimated SF
fittedSurvivalEXP <- fitted_survivalEXP(x=mydata,mu=1/muhatEXP)

## Cox-Snell residuals
CSresidualModelEXP <- -log(fittedSurvivalEXP)

## QQ plot
a <- ppoints(3114)
QE <- qexp(a)  

df_exp <- data.frame(QE = QE, CSresidualModelEXP = CSresidualModelEXP)
ggplot(df_exp, aes(sample = CSresidualModelEXP)) +
  stat_qq(distribution = qexp) + 
  stat_qq_line(distribution = qexp, color = "red") + 
  theme_minimal() +
  labs(x = "Quantis Exp(1)", y = "Quantis Amostrais") +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 9.5),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "top")

ggsave("Exponential_model.pdf", width = 158, height = 93, units = "mm")

####################################################
## Estimation - Weibull model
####################################################

## ACDm fit (initial values)
fitModelWEI <- acdFit(durations = mydata, model = "LACD1",
                      dist = "weibull", order = c(1,1))

omegahatWEI  <- as.numeric(fitModelWEI$mPar[1])
alpha1hatWEI <- as.numeric(fitModelWEI$mPar[2])
beta1hatWEI  <- as.numeric(fitModelWEI$mPar[3])
gammahatWEI <- as.numeric(fitModelWEI$dPara)


residualsWEI <- fitModelWEI$residuals
muhatWEI     <- mydata / residualsWEI  

fitModelWEI

fitted_survivalWEI <- function(x, muhat, gamma) {
  scale <- muhat / gamma(1 + 1/gamma) 
  survival <- 1 - pweibull(x, shape = gamma, scale = scale)
  return(survival)
}

## Estimated SF
fittedSurvivalWEI <- fitted_survivalWEI(x = mydata, muhat = muhatWEI, gamma = gammahatWEI)

## Cox-Snell residuals
CSresidualModelWEI <- -log(fittedSurvivalWEI)

## QQ plot
a <- ppoints(3114)
QW <- qexp(a)  

df_wei <- data.frame(QW = QW, CSresidualModelWEI = CSresidualModelWEI)

ggplot(df_wei, aes(sample = CSresidualModelWEI)) +
  stat_qq(distribution = qexp) + 
  stat_qq_line(distribution = qexp, color = "red") +
  theme_minimal() +
  labs(x = "Quantis Exp(1)", y = "Quantis Amostrais") +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 9.5),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "top")

ggsave("Weibull_model.pdf", width = 158, height = 93, units = "mm")


####################################################
## Estimation - GG model
####################################################

## ACDm fit (initial values)
fitModelGG <- acdFit(durations = mydata, model = "LACD1",dist = "gengamma", 
                     order = c(1,1))

omegahatGG  <-  as.numeric(fitModelGG$mPar[1])
alpha1hatGG <-  as.numeric(fitModelGG$mPar[3])
beta1hatGG  <-  as.numeric(fitModelGG$mPar[2])

kappahatGG <- as.numeric(fitModelGG$dPar[1])
zetahatGG  <- as.numeric(fitModelGG$dPar[2])  

fitModelGG

residualsGG <- fitModelGG$residuals
muhatGG     <- mydata/residualsGG

fitted_survivalGG <- function(x,mu,kappa,zeta){
  
  phikt <- (gamma(kappa)/gamma(kappa+zeta^(-1)))
  
  at <- (mydata/(phikt*mu))^(zeta)
  
  survival   <- 1.0 - (Igamma(kappa,at,lower=TRUE)/gamma(kappa))
  
  return(survival)
}

## Estimated SF
fittedSurvivalGG <- fitted_survivalGG(x=mydata,mu=muhatGG,
                                      kappa=kappahatGG,zeta=zetahatGG)

## Cox-Snell residuals
CSresidualModelGG <- -log(fittedSurvivalGG)

## QQ plot
a <- ppoints(3114)
QGG <- qexp(a)  #quantis da Exp(1)

df_gg <- data.frame(QGG = QGG, CSresidualModelGG = CSresidualModelGG)
ggplot(df_gg, aes(sample = CSresidualModelGG)) +
  stat_qq(distribution = qexp) +  # QQ-Plot para distribuição exponencial
  stat_qq_line(distribution = qexp, color = "red") +  # Linha de referência
  theme_minimal() +
  labs(x = "Quantis Exp(1)", y = "Quantis Amostrais") +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 10),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 9.5),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "top")

ggsave("GG_model.pdf", width = 158, height = 93, units = "mm")



# Combining the residuals from all models into a single data frame

df_qq <- data.frame(
  residual = c(CSresidualModelWEI,CSresidualModelEXP, CSresidualModelGG),
  model = rep(
    c("Weibull", "Exponencial", "Generalized Gamma"),
    times = c(length(CSresidualModelWEI),length(CSresidualModelEXP), length(CSresidualModelGG))
  )
)

# # Ordering the models
df_qq$model <- factor(df_qq$model,
                      levels = c("Weibull", "Exponencial", "Generalized Gamma"))

# QQ plot for all models

ggplot(df_qq, aes(sample = residual)) +
  stat_qq(distribution = qexp) +
  stat_qq_line(distribution = qexp, color = "red") +
  facet_wrap(~ model) +
  theme_bw() +
  labs(
    x = "Quantis Exp(1)",
    y = "Quantis Amostrais"
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

ggsave("GGEXPWEI_model.pdf", width = 158, height = 93, units = "mm")


# 1) Function Expected Shortfall
expected_shortfall <- function(x, p) {
  var_p <- quantile(x, probs = p, na.rm = TRUE)
  es   <- mean(x[x >= var_p], na.rm = TRUE)
  return(es)
}

# 2) ES 95%

es_EXP_95 <- expected_shortfall(CSresidualModelEXP, p = 0.95)
es_WEI_95  <- expected_shortfall(CSresidualModelWEI,  p = 0.95)
es_GG_95  <- expected_shortfall(CSresidualModelGG,  p = 0.95)


# 4) Results

es_results <- data.frame(
  Modelo      = c("Exponencial", "Weibull", "Generalized Gamma"),
  ES_95       = c(es_EXP_95, es_WEI_95,  es_GG_95)
)

es_results

# Function to extract ACF into a data frame
get_acf_df <- function(residuals, model, lag.max = 35) {
  acf_obj <- acf(residuals, plot = FALSE, lag.max = lag.max)
  data.frame(
    lag   = acf_obj$lag[-1],
    acf   = acf_obj$acf[-1],
    model = model
  )
}

# Create a data frame with the ACFs for all models
df_acf <- bind_rows(
  get_acf_df(CSresidualModelEXP,  "Exponencial"),
  get_acf_df(CSresidualModelWEI,  "Weibull"),
  get_acf_df(CSresidualModelGG,"Generalized Gamma")
)

# Confidence level for the bands (±1.96/√n)
n_obs <- length(CSresidualModelEXP)
ci    <- qnorm(0.975) / sqrt(n_obs)

#ACF plot
ggplot(df_acf, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_segment(aes(xend = lag, yend = 0)) +
  geom_hline(yintercept =  c(-ci, ci),
             linetype   = "dashed",
             color      = "blue") +
  facet_wrap(~ model, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 35, by = 7),limits = c(0, 35)) +
  labs( x     = "Lag",y     = "ACF") +
  theme_bw() +
  theme(
    axis.title.y = element_text(colour = "black", size = 10),
    axis.title.x = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 9.5),
    axis.line = element_line(colour = "black"),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("ACF_acd.pdf", width = 158, height = 93, units = "mm")


# Ljung-Box Test
residuals_list <- list(
  "Exponencial"        = CSresidualModelEXP,
  "Weibull"            = CSresidualModelWEI,
  "Generalized Gamma" = CSresidualModelGG
)

run_ljung_box <- function(res, p, q) {
  lag  <- 35
  tb   <- Box.test(res, lag = lag, type = "Ljung-Box", fitdf = p + q)
  data.frame(
    Modelo      = NA,                # será preenchido abaixo
    Lag         = lag,
    Estatística = unname(tb$statistic),
    GrausLib    = unname(tb$parameter),
    p_value     = tb$p.value
  )
}

results <- do.call(rbind, lapply(names(residuals_list), function(model_name) {
  df <- run_ljung_box(residuals_list[[model_name]], p = 1, q = 1)
  df$Modelo <- model_name
  df
}))

results


