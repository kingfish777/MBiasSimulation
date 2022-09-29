
# install all of these libraries if you haven't yet
# first: install.packages("BiocManager")
# then: 
# BiocManager::install(c("tmle", "SuperLearner", "kableExtra", "tidyverse", "ggplot2", "earth", "ranger", "dagitty"), dependencies = TRUE, force = TRUE)
#
library(tmle) 
library(SuperLearner)
library(kableExtra)
library(tidyverse)
library(ggplot2)
library(earth) 
library(ranger)
library(dagitty)
library(gt)
library(simcausal)

set.seed(7) # so results are reproducible

generate_data <- function(n, adjSetType) {
  #n = 1000BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
  W1 <- rbinom(n, size = 1, prob=0.2) # binary confounder
  W2 <- rbinom(n, size = 1, prob=0.5) #binary confounder
  W3 <- rbinom(n, size = 1, prob=0.3)
  W4 <- rbinom(n, size = 1, prob=0.2)
  A <- rbinom(n, size=1, prob = plogis(-2+0.2*W1+0.1*W2+0.1*W3+0.1*W4)) # binary treatment
  P1 <- rbinom(n, size=1, prob=0.4)
  P2 <- rbinom(n, size=1, prob=0.3)
  P3 <- rbinom(n, size=1, prob=0.1)
#  Y <- rbinom(n, size =1, prob=plogis(-1+A-0.1*W1+0.2*W2+0.1*W3+0.1*W4+0.1*P1+0.3*P2+0.1*P3))
  # counterfactual
  Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  # Observed outcome
  Y <- Y.1*A + Y.0*(1 - A)
   M1 <- rbinom(n, size=1, prob=plogis(1+0.4*A))
  C1 <- rbinom(n, size=1, prob=plogis(-1+0.3*A+0.4*Y))
  M2 <- rbinom(n, size = 1, prob = plogis(-1+0.8*A))
  WM <- rbinom(n, size=1, prob=plogis(1+0.4*W3+0.7*M2))
  C2 <- rbinom(n, size=1, prob=plogis(-1+0.3*A+0.4*Y))
  M3 <- rbinom(n, size = 1, prob = plogis(-1+0.3*A))
  Chimera <- rbinom(n, size=1, prob=plogis(2+0.3*W4+0.1*M3+.4*C2))
  if (adjSetType == "all") {
    dat <- data.frame(W1, W2, M1, C1, W3, M2, WM, C2, W4, M3, P1, P2, P3, Chimera, A, Y, Y.1, Y.0)
  } else if (adjSetType == "allPure") {
    dat <- data.frame(W1, W2, M1, C1, WM, P1, P2, P3, Chimera, A, Y, Y.1, Y.0)
  } else if (adjSetType == "onlyConfounders") { dat <- data.frame(W1, W2, W3, W4, A, Y, Y.1, Y.0) 
  } else if (adjSetType == "onlyConfoundersAndPrecisionVars") {
    dat <- data.frame(W1, W2, P1, P2, P3, A, Y, Y.1, Y.0)
  } else if (adjSetType == "Confounders_C1_WM_Chimera") {
    dat <- data.frame(W1, W2, W3, W4, C1, WM, Chimera, A, Y, Y.1, Y.0)
  } else if (adjSetType == "onlyBadCovariates") {
    dat <- data.frame(C1, M1, M2, C2, WM, Chimera, A, Y, Y.1, Y.0)
  }
  return(dat)
}

generateData<- function(n){
  w1 <- rbinom(n, size=1, prob=0.5)
  w2 <- rbinom(n, size=1, prob=0.65)
  w3 <- round(runif(n, min=0, max=4), digits=0)
  w4 <- round(runif(n, min=0, max=5), digits=0)
  A <- rbinom(n, size=1, prob= plogis(-5 + 0.05*w2 + 0.25*w3 + 0.6*w4 + 0.4*w2*w4))
  # the counterfactuals aka potential outcomes
  Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  # Observed outcome
  Y <- Y.1*A + Y.0*(1 - A)
  # return data.frame
  data.frame(w1, w2, w3, w4, A, Y, Y.1, Y.0)
}
#quantifyEffectWithSubset()

n = 20000
#dat_obs <- generate_data(n, "onlyBadCovariates")
dat_obs <- generate_data(n, "onlyConfounders")
#dat_obs <- generate_data(n, "onlyConfoundersAndPrecisionVars")
#dat_obs <- generate_data(n, "all")
#dat_obs <- generate_data(n, "allPure")

#dat <- data.frame(dat_obs)
#
#kable(head(dat_obs), digits=2, caption = "Simulated dataset.")

# So we will start out with a simulation containing only
# four "well-behaved" / "pure" confounders (W1, W2, W3, and W4)
# plus the exposure (A) and outcome variables (Y)
#############

# Now starting with where the original tutorial left off
#  with a more solid data generating process

# Source: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fsim.7628&file=sim7628-sup-0003-Appendix.pdf
#   also at: https://github.com/migariane/SIM-TMLE-tutorial/blob/master/Appendix_SIM_2018.pdf

# the gold standard: 

# True ATE in the population
set.seed(7777)

#dat_obs <- generate_data(n, "onlyBadCovariates")
#dat_obs <- generate_data(n, "onlyConfounders")

ObsDataTrueATE <- generate_data(n = 500000, "onlyConfounders")
True_EY.1 <- mean(ObsDataTrueATE$Y.1)
True_EY.0 <- mean(ObsDataTrueATE$Y.0)
True_ATE <- True_EY.1-True_EY.0 ;True_ATE
# Below is the OR = ad/bc
True_MOR <- (True_EY.1*(1-True_EY.0))/((1-True_EY.1)*True_EY.0);True_MOR
cat("\n True_ATE:", abs(True_ATE))
#

# Data for analysis
set.seed(7722)
ObsData <- generate_data(n = 20000, "onlyConfounders")
write.csv(ObsData, "obsData.csv")


# naive estimation - conditional odds ratio

naive = glm(data = ObsData, Y ~ A + W1 + W2 + W3 + W4, family = binomial)
summary(naive)
exp(naive$coef[2])
exp(confint(naive))
# exponentiated odds ratio = 2.608
# True mOR = 2.52 (!!!)

## TMLE implementation by hand
# Step 1 estimation and prediction of the model for the outcome (G-computation)
gm <- glm(data = ObsData, Y ~ A + W1 + W2 + W3 + W4, family = binomial)
# Prediction for A, A = 1 and, A = 0
QAW_0 <- predict(gm, type = "response")
Q1W_0 = predict(gm, newdata=data.frame(A = 1, ObsData [,c("W1","W2","W3","W4")]), type = "response")
Q0W_0 = predict(gm, newdata=data.frame(A = 0, ObsData [,c("W1","W2","W3","W4")]), type = "response")
# Estimated mortality risk difference
mean(Q1W_0 - Q0W_0)
# Estimated MOR
mean(Q1W_0)*(1-mean(Q0W_0))/((1-mean(Q1W_0))*mean(Q0W_0))
# Step 2 estimation and prediction of the propensity score (ps)
psm <- glm(A ~ W1 + W2 + W3 + W4, family = binomial, data = ObsData)
gW = predict(psm, type = "response")
summary(gW)
# Step 3 computation of H (clever covariates) and estimation of epsilon
H1W = ObsData$A / gW
H0W = (1-ObsData$A) / (1 - gW)
epsilon <- coef(glm(ObsData$Y ~ -1 + H0W + H1W + offset(qlogis(QAW_0)), family = binomial))
# Step 4 Targeted estimate of the ATE and Marginal Odds Ratio
Q1W_1 <- plogis(qlogis(Q1W_0) + epsilon[2] / gW)
Q0W_1 <- plogis(qlogis(Q0W_0) + epsilon[1] / (1-gW))
# ATE
ATEtmle1 <- mean(Q1W_1 - Q0W_1); ATEtmle1
cat("\n ATEtmle1_bias:", abs(True_ATE - ATEtmle1))
cat("\n ATEtmle1_rel_bias:",abs(True_ATE - ATEtmle1)/True_ATE,"%")
# Marginal OR
tmle1.MOR <- mean(Q1W_1) * (1 - mean(Q0W_1)) / ((1 - mean(Q1W_1)) * mean(Q0W_1)); tmle1.MOR
# Table to visualize the data
psi <- Q1W_1 - Q0W_1

library(DT)
df <- round(cbind(Q1W_0, Q0W_0, gW, eps1=epsilon[1], eps2=epsilon[2], psi), digits = 4)
renderDataTable({datatable(head(df, n = nrow(df)), options = list(pageLength = 5, digits = 3))})
# Step 5 statistical inference (efficient influence curve)
# Efficient influence curve ATE
EY1tmle<-mean(Q1W_1)
EY0tmle<-mean(Q0W_1)
d1 <- ((ObsData$A) * (ObsData$Y - Q1W_1)/gW) + Q1W_1 - EY1tmle
d0 <- ((1 - ObsData$A) * (ObsData$Y - Q0W_1))/(1 - gW) + Q0W_1 - EY0tmle
IC <- d1 - d0
n <- nrow(ObsData)
varHat.IC <- var(IC) / n
ATEtmle1CI <- c(ATEtmle1 - 1.96 * sqrt(varHat.IC), ATEtmle1 + 1.96 * sqrt(varHat.IC)); ATEtmle1;
ATEtmle1CI


# Efficient influence curve MOR
ICmor_tmle <- (1 - EY0tmle) / EY0tmle / (1 - EY1tmle)^2 * d1 - EY1tmle / (1 - EY1tmle) / EY0tmle^2 *
  d0
varHat2.IC <- var(ICmor_tmle) / n
tmle1_ORCI <- tmle1.MOR + c(-1.96,1.96)*sqrt(varHat2.IC); tmle1.MOR; tmle1_ORCI
# Augmented inverse probability treatment weighting (AIPTW) estimator
EY1aiptw <- mean((ObsData$A) * (ObsData$Y - Q1W_0) / gW + Q1W_0)
EY0aiptw <- mean((1 - ObsData$A) * (ObsData$Y - Q0W_0) / (1 - gW) + Q0W_0)
AIPTW_ATE <- EY1aiptw - EY0aiptw; AIPTW_ATE
cat("\n AIPTW_bias:", abs(True_ATE - AIPTW_ATE))
cat("\n AIPTW_rel_bias:",abs(True_ATE - AIPTW_ATE) / True_ATE,"%")
D1 <- (ObsData$A) * (ObsData$Y - Q1W_0) / gW + Q1W_0 - EY1aiptw
D0 <- (1 - ObsData$A) * (ObsData$Y - Q0W_0) / (1 - gW) + Q0W_0 - EY0aiptw
varHat_AIPTW <- var(D1 - D0) / n
# AIPTW ATE 95%CI
ATEaiptw_CI <- c(AIPTW_ATE - 1.96 * sqrt(varHat_AIPTW), AIPTW_ATE + 1.96 *
                   sqrt(varHat_AIPTW)); AIPTW_ATE; ATEaiptw_CI
# AIPTW MOR 95%CI
AIPTW_MOR <- (EY1aiptw * (1 - EY0aiptw))/((1 - EY1aiptw) * EY0aiptw);AIPTW_MOR
ICmor_aiptw <- (1 - EY0aiptw) / EY0aiptw / (1 - EY1aiptw)^2 * D1 - EY1aiptw / (1 - EY1aiptw) /
  EY0aiptw^2 * D0
varHat_AIPTW2 <- var(ICmor_aiptw) / n
MORaiptw_CI <-c(AIPTW_MOR - 1.96*sqrt(varHat_AIPTW2), AIPTW_MOR +
                  1.96*sqrt(varHat_AIPTW2)); AIPTW_MOR; MORaiptw_CI
# R-package tmle (base implementation includes SL.step, SL.glm and SL.glm.interaction)
library(tmle)
library(SuperLearner)
TMLE2 <- tmle(Y = ObsData$Y, A = ObsData$A, W = ObsData[,c("W1", "W2", "W3", "W4")], family =
                "binomial")
#NOTE:
#Note that the tmle function default bounds the probabilities in the clever covariate denominators at 0.025.
#You can remove this bound by specifying gbound=0
ATEtmle2 <- TMLE2$estimates$ATE$psi;ATEtmle2
TMLE2$estimates$ATE$CI
MORtmle2 <- TMLE2$estimates$OR$psi;MORtmle2
TMLE2$estimates$OR$CI
cat("\n ATEtmle2_bias:", abs(True_ATE - ATEtmle2))
cat("\n ATEtmle2_Rel_bias:",abs(True_ATE - ATEtmle2) / True_ATE,"%")

# R-package tmle with user-selected Super learner libraries
library(tmle)
library(SuperLearner)
SL.library <- c("SL.glm","SL.step","SL.step.interaction", "SL.glm.interaction","SL.gam",
                "SL.randomForest", "SL.rpart")
TMLE3 <- tmle(Y = ObsData$Y,A = ObsData$A,W = ObsData [,c("W1", "W2", "W3", "W4")],
              family = "binomial", Q.SL.library = SL.library,g.SL.library = SL.library)
ATEtmle3 <- TMLE3$estimates$ATE$psi;ATEtmle3
TMLE3$estimates$ATE$CI
MORtmle3 <- TMLE3$estimates$OR$psi;MORtmle3
TMLE3$estimates$OR$CI
cat("\n ATEtmle3_bias:", abs(True_ATE - ATEtmle3))
cat("\n ATEtmle3_rel_bias:",abs(True_ATE - ATEtmle3) / True_ATE,"%")




#More complex data-generating process
#Readers interested in simulating more complex dependence structures among the covariates W1-W4 could
#potentially use the R-package simcausal (Sofrygin O, van der Laan MJ, Neugebauer R (2015). simcausal:
#                                           Simulating Longitudinal Data with Causal Inference Applications. R package version 0.5).
#See the example here below:
devtools::install_github('osofr/simcausal', build_vignettes = FALSE, force = TRUE, dependencies = TRUE)

library(simcausal)
M <- DAG.empty()
M <- M +
  node("W1", # age (0/1); 1 -> high age
       distr = "rbern",
       prob = .5) +
  node("W2", # ses (1/2/3/4/5); higher age, higher probability of belonging to upper class
       distr = "rcat.b1",
       probs = {
         plogis(-3.1 + 0.05*W1); # upper middle class, 4%
         plogis(-1.25 + 0.04*W1); # middle class, 22%
         plogis(-0.05 + 0.03*W1); # lower middle 49%
         plogis(-1.65 + 0.02*W1); # working class 16%
         plogis(-2.3 + 0.01*W1)}) +
  node("W3", #comorbidities (1/2/3/4/5);
       distr = "rcat.b1",
       probs = {
         plogis(-0.8 + 0.005*W1 + 0.1*W2);
         plogis(-0.1 + 0.010*W1 + 0.12*W2);
         plogis(-1.2 + 0.015*W1 + 0.15*W2);
         plogis(-1.6 + 0.020*W1 + 0.2*W2);
         plogis(-2.5 + 0.025*W1 + 0.25*W2)}) +
  node("W4", # stage (1/2/3/4); # the higher the worse
       distr = "rcat.b1",
       probs = {
         plogis(-1 + 0.01*W1 - 0.04*W2);
         plogis(-0.2 + 0.02*W1 - 0.05*W2);
         plogis(-0.8 + 0.03*W1 - 0.055*W2);
         plogis(-2 + 0.04*W1 - 0.1*W2)}) + 
node("A", # a = 0 mono therapy; a = 1 dualtherapy
     distr = "rbern",
     prob = plogis(-1.4 + 0.05*W1 + 0.25*W3 + 0.1*exp(W4))) +
node("Y", #y = 0 -> death; y = 1 -> alive
       distr = "rbern",
       prob = plogis(-3.4 + 0.75*A - 0.1*W1 + 0.35*W2 + 0.25*W3 + 0.20*sqrt(1/W4) - 0.9*A*W2 + 1.1*A*W3))
Mset <- set.DAG(M)
# simulate observed data
Odat <- simcausal::sim(DAG = Mset, n = 10000, rndseed = 7693)
# specify the two interventions
a1 <- node("A", distr = "rbern", prob = 1)
Mset <- Mset + action("a1", nodes = a1)
a0 <- node("A", distr = "rbern", prob = 0)
Mset <- Mset + action("a0", nodes = a0)
# counterfactual data
dat <- simcausal::sim(DAG = Mset, actions = c("a1", "a0"), n = 1000000, rndseed = 7693)
head(dat[["a1"]]); head(dat[["a0"]])
# E(y) under a=1 (chemo)
Mset <- set.targetE(Mset, outcome = "Y", param = "a1")
eval.target(Mset, data = dat)$res
# E(y) unter a=0 (chemo and radio)
Mset <- set.targetE(Mset, outcome = "Y", param = "a0")
eval.target(Mset, data = dat)$res
# ATE (additive scale)
Mset <- set.targetE(Mset, outcome = "Y", param = "a1-a0")
eval.target(Mset, data = dat)$res
# multiplicative scale
Mset <- set.targetE(Mset, outcome = "Y", param = "a1/a0")
eval.target(Mset, data = dat)$res
#DAG
plotDAG(Mset, xjitter = 0.3, yjitter = 0.04,edge_attrs = list(width = 0.5, arrow.width = 0.2, arrow.size = 0.3),
        vertex_attrs = list(size = 12, label.cex = 0.8))

### pay no attention beneath here!



sl_libs <- c('SL.glmnet', 'SL.ranger', 'SL.nnet', 'SL.earth') # (penalized regression, random forests, and multivariate adaptive regression splines)

Y <- ObsData$Y
#W_A <- dat_obs %>% select(-Y) # remove the outcome to make a matrix of predictors (A, W1, W2, W3, W4) for SuperLearner

A <- ObsData$A
W <- ObsData %>% select(-Y, -A) # matrix of predictors that only contains the confounders W1, W2, W3, and W4

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

#######

# Trying to get rcausal running


library(rcausal)
dd <- tetradrunner(algoId = 'gfci', df = dat, dataType = 'discrete', scoreId = 'bdeu', 
                   #maxDegree=-1,
                   faithfulnessAssumed=TRUE, verbose = TRUE)
dd$graph
graphDot <- tetradrunner.tetradGraphToDot(dd$graph)
dot(graphDot)
data("charity")    #Load the charity dataset
tetradrunner <- tetradrunner(algoId = 'fges',df = charity,scoreId = 'sem-bic', dataType = 'continuous',faithfulnessAssumed=TRUE,maxDegree=-1,verbose=TRUE)    #Compute FGES search
tetradrunner$nodes #Show the result's nodes
tetradrunner$edges #Show the result's edges

graph <- tetradrunner$graph
graph$getAttribute('BIC')

nodes <- graph$getNodes()
for(i in 0:as.integer(nodes$size()-1)){
  node <- nodes$get(i)
  cat(node$getName(),": ",node$getAttribute('BIC'),"\n")
}




