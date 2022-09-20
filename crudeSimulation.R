
library(tmle)
library(SuperLearner)
library(kableExtra)
library(tidyverse)
library(ggplot2)

set.seed(7) # so results are reproducible

generate_data <- function(n, adjSetType) {
  #n = 1000
  W1 <- rbinom(n, size = 1, prob=0.2) # binary confounder
  W2 <- rbinom(n, size = 1, prob=0.5) #binary confounder
  W3 <- rbinom(n, size = 1, prob=0.3)
  W4 <- rbinom(n, size = 1, prob=0.2)
  A <- rbinom(n, size=1, prob = plogis(-2+0.2*W1+0.1*W2+0.1*W3+0.1*W4)) # binary treatment
  P1 <- rbinom(n, size=1, prob=0.4)
  P2 <- rbinom(n, size=1, prob=0.3)
  P3 <- rbinom(n, size=1, prob=0.1)
  Y <- rbinom(n, size =1, prob=plogis(-1+A-0.1*W1+0.2*W2+0.1*W3+0.1*W4+0.1*P1+0.3*P2+0.1*P3))
  M1 <- rbinom(n, size=1, prob=plogis(1+0.4*A))
  C1 <- rbinom(n, size=1, prob=plogis(-1+0.3*A+0.4*Y))
  M2 <- rbinom(n, size = 1, prob = plogis(-1+0.8*A))
  WM <- rbinom(n, size=1, prob=plogis(1+0.4*W3+0.7*M2))
  C2 <- rbinom(n, size=1, prob=plogis(-1+0.3*A+0.4*Y))
  M3 <- rbinom(n, size = 1, prob = plogis(-1+0.3*A))
  Chimera <- rbinom(n, size=1, prob=plogis(2+0.3*W4+0.1*M3+.4*C2))
  if (adjSetType == "all") {
     dat <- tibble(W1, W2, M1, C1, W3, M2, WM, C2, W4, M3, P1, P2, P3, Chimera, A, Y)
  } else if (adjSetType == "allPure") {
    dat <- tibble(W1, W2, M1, C1, WM, P1, P2, P3, Chimera, A, Y)
  } else if (adjSetType == "onlyConfounders") { dat <- tibble(W1, W2, W3, W4, A, Y) 
  } else if (adjSetType == "onlyConfoundersAndPrecisionVars") {
    dat <- tibble(W1, W2, P1, P2, P3, A, Y)
  } else if (adjSetType == "Confounders_C1_WM_Chimera") {
     dat <- tibble(W1, W2, W3, W4, C1, WM, Chimera, A, Y)
  } else if (adjSetType == "onlyBadCovariates") {
    dat <- tibble(C1, M1, M2, C2, WM, Chimera, A, Y)
  }
  return(dat)
}


#quantifyEffectWithSubset()

n = 20000
#dat_obs <- generate_data(n, "onlyBadCovariates")
dat_obs <- generate_data(n, "onlyConfounders")
dat_obs <- generate_data(n, "onlyConfoundersAndPrecisionVars")
#dat_obs <- generate_data(n, "all")
#dat_obs <- generate_data(n, "allPure")

dat <- data.frame(dat_obs)
forbid <- list(c('A', 'W1'), c('Y', 'W1'), c('A', 'W2'), c('Y', 'W2'), c('Y', 'P1'), c('Y', 'P2'), c('Y', 'P3'))
require <- list(c('A', 'Y'), c('W1', 'A'), c('W1', 'Y'), c('W2', 'A'), c('W2', 'Y'), c('W3', 'A'), c('W3', 'Y'))
prior <- priorKnowledge(requiredirect = require, forbiddirect = forbid)
##

dd <- tetradrunner(algoId = 'fges', df = dat, dataType = 'discrete', scoreId = 'bdeu', 
                   maxDegree=-1,
                   faithfulnessAssumed=TRUE, verbose = TRUE, priorKnowledge = prior)

dd <- tetradrunner(algoId = 'gfci', df = dat, dataType = 'discrete', scoreId = 'bdeu', verbose = TRUE)
dd$graph
graphDot <- tetradrunner.tetradGraphToDot(dd$graph)
dot(graphDot)
##te_data(n, "allPure")
dat <- data.frame(dat_obs)


dd <- tetradrunner(algoId = 'fges', df = dat, dataType = 'discrete', scoreId = 'bdeu', 
                   #maxDegree=-1,
                   faithfulnessAssumed=TRUE, verbose = TRUE)
dd$graph
graphDot <- tetradrunner.tetradGraphToDot(dd$graph)
dot(graphDot)




kable(head(dat_obs), digits=2, caption = "Simulated dataset.")

sl_libs <- c('SL.glmnet', 'SL.ranger', 'SL.nnet', 'SL.earth') # (penalized regression, random forests, and multivariate adaptive regression splines)

Y <- dat_obs$Y
#W_A <- dat_obs %>% select(-Y) # remove the outcome to make a matrix of predictors (A, W1, W2, W3, W4) for SuperLearner

A <- dat_obs$A
W <- dat_obs %>% select(-Y, -A) # matrix of predictors that only contains the confounders W1, W2, W3, and W4

tmle_fit <-
  tmle::tmle(Y = Y, # outcome vector
             A = A, # treatment vector
             W = W, # matrix of confounders W1, W2, W3, W4
             Q.SL.library = sl_libs, # superlearning libraries from earlier for outcome regression Q(A,W)
             g.SL.library = sl_libs) # superlearning libraries from earlier for treatment regression g(W)

tmle_fit

ATE <- c(tmle_fit$estimates$ATE$psi, tmle_fit$estimates$ATE$CI[1], tmle_fit$estimates$ATE$CI[2], tmle_fit$estimates$ATE$var.psi, tmle_fit$estimates$ATE$pvalue)
ATT <- c(tmle_fit$estimates$ATT$psi, tmle_fit$estimates$ATT$CI[1], tmle_fit$estimates$ATT$CI[2], tmle_fit$estimates$ATT$var.psi, tmle_fit$estimates$ATT$pvalue)
ATC <- c(tmle_fit$estimates$ATC$psi, tmle_fit$estimates$ATC$CI[1], tmle_fit$estimates$ATC$CI[2], tmle_fit$estimates$ATC$var.psi, tmle_fit$estimates$ATC$pvalue)

#return(ATE, ATT, ATC)
print(ATE)
print(ATT)
print(ATC)



}

################
################
# Collider example


race <- rbinom(n, 1, 0.5) 
disease_severity <- rbinom(n, 1, 0.5)
ses <- rbinom(n, 1, plogis(-race))
hospitalized <- rbinom(n, 1, plogis(-3 + 2 * disease_severity + 2 * ses + 3 * race))
death <- rbinom(n, 1, plogis(-(2 * disease_severity + hospitalized)))

covid_df <-
  data.frame(
    race = factor(race, levels=c(1,0), labels=c("White","Black")),
    disease_severity = factor(disease_severity, levels=c(1,0), labels=c("Mild","Severe")),
    ses = factor(ses, levels=c(0,1), labels=c("Above poverty line", "Below poverty line")),
    hospitalized = factor(hospitalized, levels=0:1, labels=c("No","Yes")),
    death = factor(death, levels=0:1, labels=c("No","Yes"))
  )

hospitalized_df <-
  covid_df %>%
  filter(hospitalized == "Yes")

race_death_count <-
  hospitalized_df %>%
  group_by(race) %>%
  count(death) %>%
  mutate(sum_n = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/sum_n)

race_death_count %>%
  filter(death == "Yes") %>%
  ggplot(aes(x=race, y=prop)) +
  geom_bar(stat="identity") +
  theme_classic() +
  labs(x="Race",y="In-hospital Mortality",title="Race and Mortality in COVID-19",subtitle="Among Hospitalized Patients") +
  scale_y_continuous(labels = scales::percent_format(), expand=c(0,0), limits=c(0,.25)) +
  geom_text(aes(label=paste(n, "/", sum_n), y=prop+.01))


pretty_logistic_table <- function(model_fit) {
  model_fit %>%
    broom::tidy(exponentiate=T, conf.int=T) %>%
    filter(term != "(Intercept)") %>%
    mutate( term = case_when(term == "raceBlack" ~ "Race: Black",
                             term == "disease_severitySevere" ~ "Disease severity: Severe",
                             term == "sesBelow poverty line" ~ "SES: Below poverty level",
                             TRUE ~ term),
            odds_ratio = paste0(round(estimate,1)," (",round(conf.low,1),", ",round(conf.high,1),")"),
            p_value = case_when(p.value < .001 ~ "<.001",
                                TRUE ~ as.character(round(p.value, 3)))) %>%
    select(term, odds_ratio, p_value) %>%
    gt::gt() %>%
    gt::cols_label(
      term = "Coefficient",
      odds_ratio = "Odds Ratio (95% CI)",
      p_value = "P-value"
    )
}

glm(death ~ race, data = hospitalized_df, family = binomial()) %>%
  pretty_logistic_table()




