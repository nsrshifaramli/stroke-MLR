# =====================================================================
# Logistic Regression Practice — stroke-fatality-cooked2.dta
# Follows the chapter "Logistic Regression Analysis in Health and Medicine"
# (Kamarul Imran Musa). Outcome: 1 = dead at discharge, 0 = alive.
# Run one section at a time in RStudio (Ctrl/Cmd + Enter).
# =====================================================================

# ---- 0. Packages ----------------------------------------------------
# Run once if you don't have them:
# install.packages(c("ResourceSelection", "pROC", "car", "parameters"))
library(haven)        # read Stata .dta
library(tidyverse)    # data wrangling + ggplot2
library(broom)        # tidy model output
library(gtsummary)    # publication tables
library(mgcv)         # generalized additive models
library(ResourceSelection)  # Hosmer-Lemeshow test
library(pROC)               # ROC curve / AUC

# ---- 1. Read and prepare the data ----------------------------------
setwd("~/Desktop/MLR/practice")   # this script's folder
stroke <- read_dta("../data/stroke-fatality-cooked2.dta")

glimpse(stroke)

# Turn the Stata labelled variables into R factors.
# The first level is the reference category.
stroke <- stroke %>%
  mutate(across(where(is.labelled), as_factor)) %>%
  mutate(
    stroketype   = fct_recode(stroketype,
                              Ischaemic    = "CI (Ischaemic Stroke)",
                              Haemorrhagic = "Haemorrhagic Stroke"),
    referralcase = fct_recode(referralcase, Hospital = "hospital",
                              `GP/Home` = "GP/Home/Missing"),
    sex          = fct_relevel(sex, "male"),
    outcome_num  = if_else(outcome == "dead", 1, 0)   # numeric 0/1 copy
  )

# ---- 2. Describe the data ------------------------------------------
summary(stroke)
count(stroke, outcome) %>% mutate(pct = n / sum(n) * 100)  # 62 dead / 164 alive

stroke %>%
  select(-rn, -outcome_num) %>%
  tbl_summary(by = outcome) %>%
  add_p() %>%
  add_overall()

# ---- 3. Simple (univariable) logistic regression -------------------
# One predictor at a time
slr_gcs <- glm(outcome ~ gcs, data = stroke, family = binomial)
summary(slr_gcs)
tidy(slr_gcs, exponentiate = TRUE, conf.int = TRUE)   # odds ratio + 95% CI

slr_stroke <- glm(outcome ~ stroketype, data = stroke, family = binomial)
tidy(slr_stroke, exponentiate = TRUE, conf.int = TRUE)

# Likelihood ratio test vs the intercept-only model
slr_null <- glm(outcome ~ 1, data = stroke, family = binomial)
anova(slr_null, slr_gcs, test = "LRT")

# Screen ALL candidate predictors at once (Table 1 of the chapter)
stroke %>%
  select(-rn, -outcome_num) %>%
  tbl_uvregression(method = glm, y = outcome,
                   method.args = list(family = binomial),
                   exponentiate = TRUE) %>%
  add_global_p()       # LRT p-value for each variable

# ---- 4. Assess the simple model: residuals -------------------------
stroke_aug <- augment(slr_gcs, type.predict = "response",
                      type.residuals = "pearson")
stroke_aug$std_pearson <- rstandard(slr_gcs, type = "pearson")

# Raw residuals vs fitted: always two bands for a binary outcome
stroke_aug$raw_resid <- stroke$outcome_num - stroke_aug$.fitted
ggplot(stroke_aug, aes(.fitted, raw_resid)) +
  geom_point() + labs(x = "Fitted probability", y = "Raw residual")

# Standardised Pearson residuals vs the predictor (the useful plot)
ggplot(stroke_aug, aes(gcs, std_pearson)) +
  geom_jitter(width = 0.2) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed")
sum(abs(stroke_aug$std_pearson) > 2)   # how many exceed |2|

# ---- 5. Multiple logistic regression: main-effect model ------------
# GCS is NOT included: it is a mediator (stroke type -> GCS -> death)
mlr_main <- glm(outcome ~ stroketype + wbc + referralcase + sex,
                data = stroke, family = binomial)
summary(mlr_main)
tidy(mlr_main, conf.int = TRUE)                       # log-odds (beta)
tidy(mlr_main, exponentiate = TRUE, conf.int = TRUE)  # adjusted OR
drop1(mlr_main, test = "LRT")                         # LRT for each term

# Overall model test (G statistic)
anova(slr_null, mlr_main, test = "LRT")

# Crude vs adjusted ORs side by side
crude <- stroke %>%
  select(outcome, stroketype, wbc, referralcase, sex) %>%
  tbl_uvregression(method = glm, y = outcome,
                   method.args = list(family = binomial), exponentiate = TRUE)
adj <- tbl_regression(mlr_main, exponentiate = TRUE)
tbl_merge(list(crude, adj), tab_spanner = c("**Crude**", "**Adjusted**"))

# ---- 6. Linearity in the logit for wbc -----------------------------
# 6a. Empirical logits by quartile
stroke %>%
  mutate(wbc_q = cut_number(wbc, 4)) %>%
  group_by(wbc_q) %>%
  summarise(n = n(), died = sum(outcome_num),
            prop = died / n, emp_logit = log(prop / (1 - prop)))

# 6b. Quadratic term and log transformation
mlr_quad <- update(mlr_main, . ~ . + I(wbc^2))
anova(mlr_main, mlr_quad, test = "LRT")
mlr_log  <- glm(outcome ~ stroketype + log(wbc) + referralcase + sex,
                data = stroke, family = binomial)
AIC(mlr_main, mlr_log)

# 6c. Generalized additive model (edf close to 1 = linear is fine)
gam_wbc <- gam(outcome_num ~ stroketype + s(wbc) + referralcase + sex,
               data = stroke, family = binomial, method = "REML")
summary(gam_wbc)
plot(gam_wbc, shade = TRUE, rug = TRUE)

# ---- 7. Interaction ------------------------------------------------
# 7a. Screen each two-way interaction (Table 5)
pairs <- c("stroketype:referralcase", "stroketype:wbc", "referralcase:sex",
           "stroketype:sex", "wbc:sex", "wbc:referralcase")
map_dfr(pairs, function(term) {
  m <- update(mlr_main, as.formula(paste(". ~ . +", term)))
  lrt <- anova(mlr_main, m, test = "LRT")
  tibble(interaction = term, chisq = lrt$Deviance[2],
         p = lrt$`Pr(>Chi)`[2], AIC = AIC(m))
}) %>% arrange(p)

# 7b. Final model with stroke type x referral source
mlr_final <- glm(outcome ~ stroketype * referralcase + wbc + sex,
                 data = stroke, family = binomial)
summary(mlr_final)
tidy(mlr_final, exponentiate = TRUE, conf.int = TRUE)
anova(mlr_main, mlr_final, test = "LRT")

# 7c. Stratum-specific ORs (needed because of the interaction)
# OR haemorrhagic vs ischaemic among GP/Home = exp(b_haem + b_interaction)
b <- coef(mlr_final); V <- vcov(mlr_final)
haem <- "stroketypeHaemorrhagic"; gp <- "referralcaseGP/Home"
int  <- "stroketypeHaemorrhagic:referralcaseGP/Home"
strat_or <- function(terms) {
  est <- sum(b[terms]); se <- sqrt(sum(V[terms, terms]))
  round(exp(c(OR = est, lower = est - 1.96 * se, upper = est + 1.96 * se)), 2)
}
strat_or(haem)          # haemorrhagic vs ischaemic, hospital-referred
strat_or(c(haem, int))  # haemorrhagic vs ischaemic, GP/Home-referred
strat_or(gp)            # GP/Home vs hospital, ischaemic
strat_or(c(gp, int))    # GP/Home vs hospital, haemorrhagic

# 7d. Illustrative continuous x categorical interaction (GCS x sex)
m_gcs <- glm(outcome ~ gcs * sex + wbc, data = stroke, family = binomial)
summary(m_gcs)
confint.default(m_gcs)  # Wald CI
confint(m_gcs)          # profile-likelihood CI (better under near-separation)
table(stroke$sex, cut(stroke$gcs, c(2, 8, 12, 15)), stroke$outcome)

# ---- 8. Assess the final model -------------------------------------
# 8a. Goodness of fit: Hosmer-Lemeshow with different groupings
for (g in c(10, 8, 6)) print(hoslem.test(stroke$outcome_num, fitted(mlr_final), g = g))
hoslem.test(stroke$outcome_num, fitted(mlr_main), g = 10)  # compare

# 8b. Discrimination: ROC curve and AUC
roc_final <- roc(stroke$outcome_num, fitted(mlr_final))
auc(roc_final); ci.auc(roc_final)
plot(roc_final, print.auc = TRUE)

# Tjur's R2: mean fitted prob in dead minus in alive
with(stroke, mean(fitted(mlr_final)[outcome_num == 1]) -
             mean(fitted(mlr_final)[outcome_num == 0]))

# 8c. Diagnostics: leverage, Cook's distance, delta chi-square
diag_final <- augment(mlr_final, type.predict = "response") %>%
  mutate(rn = stroke$rn,
         std_pearson = rstandard(mlr_final, type = "pearson"),
         delta_chisq = std_pearson^2)
ggplot(diag_final, aes(.fitted, delta_chisq)) + geom_point() +
  labs(x = "Fitted probability", y = "Delta chi-square")
diag_final %>% arrange(desc(.cooksd)) %>%
  select(rn, outcome, stroketype, referralcase, wbc, .fitted, .cooksd) %>%
  head(5)

# ---- 9. Present and predict ----------------------------------------
tbl_regression(mlr_final, exponentiate = TRUE)

# Predicted probability of death by wbc, stroke type and referral (males)
newdat <- expand_grid(wbc = seq(4, 27, by = 0.5),
                      stroketype = levels(stroke$stroketype),
                      referralcase = levels(stroke$referralcase),
                      sex = "male")
pred <- predict(mlr_final, newdata = newdat, type = "link", se.fit = TRUE)
newdat <- newdat %>%
  mutate(prob  = plogis(pred$fit),
         lower = plogis(pred$fit - 1.96 * pred$se.fit),
         upper = plogis(pred$fit + 1.96 * pred$se.fit))
ggplot(newdat, aes(wbc, prob, colour = stroketype, fill = stroketype)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line() + facet_wrap(~ referralcase) +
  labs(x = "White cell count", y = "Predicted probability of death")

# Predict for one new patient
predict(mlr_final, type = "response",
        newdata = data.frame(stroketype = "Haemorrhagic", referralcase = "GP/Home",
                             wbc = 15, sex = "female"))
