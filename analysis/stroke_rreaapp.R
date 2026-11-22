# =====================================================================
# Stroke fatality at discharge: binary logistic regression (RREAAPP)
# Plain R script version of analysis/stroke_rreaapp.qmd
# (generated with knitr::purl(); the .qmd stays the main source)
#
# 1. REASON
#   Question: which admission characteristics are associated with death
#   before discharge after acute stroke?
#   Outcome:  outcome (alive = 0, dead = 1); the model estimates odds of death.
#   GCS is a mediator (stroke type -> GCS -> death), so it is used in the
#   univariable step only and left out of the multivariable model.
#   No interaction terms.
#
# How to run:
#   - Put stroke-fatality-cooked2.dta in ../data/
#   - Set the working directory to this analysis/ folder, e.g.
#       setwd("~/Desktop/MLR/analysis")
#   - Run section by section (Cmd/Ctrl + Enter), or the whole file with
#     Cmd/Ctrl + Shift + Enter (source with echo) to see every table and plot.
#   - Teaching dataset: results illustrate method only.
# =====================================================================

# ---- setup ----
invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))
library(haven)
library(tidyverse)
library(gtsummary)
library(broom)
library(ResourceSelection)
library(pROC)
library(car)
set.seed(2026)
theme_set(theme_minimal(base_size = 12))
theme_gtsummary_compact()

# =====================================================================
# 2. READ: import with haven, fix names and codes
# =====================================================================

# ---- read-import ----
raw <- read_dta("../data/stroke-fatality-cooked2.dta")
dim(raw)

# ---- read-labels ----
tibble(variable = names(raw),
       label    = map_chr(raw, ~ attr(.x, "label") %||% ""),
       class    = map_chr(raw, ~ class(.x)[1]))

# ---- read-clean ----
stroke <- raw %>%
  zap_label() %>%
  rename(id = rn, heart_rate = hr, haemoglobin = hb, platelet = plat,
         sodium = na, potassium = potas, blood_sugar = cbs,
         married = maritalstat, referral = referralcase,
         stroke_type = stroketype) %>%
  mutate(across(where(is.labelled), as_factor)) %>%
  mutate(
    outcome     = fct_relevel(outcome, "alive"),
    sex         = fct_relevel(sex, "male"),
    stroke_type = fct_recode(stroke_type,
                             Ischaemic    = "CI (Ischaemic Stroke)",
                             Haemorrhagic = "Haemorrhagic Stroke") %>%
                  fct_relevel("Ischaemic"),
    referral    = fct_recode(referral, Hospital = "hospital",
                             `GP/Home` = "GP/Home/Missing") %>%
                  fct_relevel("Hospital"),
    across(c(hypertension, married, malay, diabetes, dyslipidaemia),
           ~ fct_relevel(.x, "no")),
    dead        = as.integer(outcome == "dead")      # numeric 0/1 copy
  ) %>%
  mutate(across(where(is.factor), fct_drop))

# variable labels for tables
var_labels <- list(
  sex = "Sex", age = "Age (years)", malay = "Malay ethnicity",
  married = "Married", gcs = "Glasgow Coma Scale", sbp = "Systolic BP (mmHg)",
  dbp = "Diastolic BP (mmHg)", heart_rate = "Heart rate (bpm)",
  haemoglobin = "Haemoglobin (g/dL)", platelet = "Platelet (×10⁹/L)",
  wbc = "White cell count (×10⁹/L)", sodium = "Sodium (mmol/L)",
  potassium = "Potassium (mmol/L)", blood_sugar = "Capillary blood sugar (mmol/L)",
  urea = "Urea (mmol/L)", hypertension = "Hypertension", diabetes = "Diabetes",
  dyslipidaemia = "Dyslipidaemia", referral = "Referral source",
  stroke_type = "Stroke type", outcome = "Outcome at discharge")

glimpse(stroke)
count(stroke, outcome) %>% mutate(percent = round(100 * n / sum(n), 1))

# ---- read-groups ----
num_vars <- c("age", "gcs", "sbp", "dbp", "heart_rate", "haemoglobin", "platelet",
              "wbc", "sodium", "potassium", "blood_sugar", "urea")
cat_vars <- c("sex", "malay", "married", "hypertension", "diabetes",
              "dyslipidaemia", "referral", "stroke_type")

# =====================================================================
# 3. EXPLORE: missing data, errors, outliers, plots
# =====================================================================

# ---- explore-missing ----
stroke %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  arrange(desc(n_missing))

# ---- explore-errors ----
ranges <- tribble(
  ~variable,     ~lower, ~upper,
  "age",             18,    110,
  "gcs",              3,     15,
  "sbp",             60,    300,
  "dbp",             30,    200,
  "heart_rate",      30,    220,
  "haemoglobin",      4,     22,
  "platelet",        10,   1000,
  "wbc",              1,     50,
  "sodium",         110,    170,
  "potassium",        2,      8,
  "blood_sugar",      1,     40,
  "urea",             1,     60)

range_check <- stroke %>%
  select(id, all_of(num_vars)) %>%
  pivot_longer(-id, names_to = "variable") %>%
  left_join(ranges, by = "variable") %>%
  filter(value < lower | value > upper)
range_check

tibble(
  check  = c("Duplicate IDs", "DBP >= SBP", "Non-integer GCS",
             "Values outside plausible range"),
  n_flagged = c(sum(duplicated(stroke$id)), sum(stroke$dbp >= stroke$sbp),
                sum(stroke$gcs != round(stroke$gcs)), nrow(range_check)))

map(stroke[cat_vars], table)   # check factor levels

# ---- explore-outliers ----
outliers <- stroke %>%
  select(id, all_of(num_vars)) %>%
  pivot_longer(-id, names_to = "variable") %>%
  group_by(variable) %>%
  mutate(q1 = quantile(value, 0.25), q3 = quantile(value, 0.75),
         iqr = q3 - q1,
         outlier = value < q1 - 1.5 * iqr | value > q3 + 1.5 * iqr) %>%
  ungroup()

outliers %>% group_by(variable) %>%
  summarise(min = min(value), q1 = first(q1), median = median(value),
            q3 = first(q3), max = max(value), n_outliers = sum(outlier)) %>%
  arrange(desc(n_outliers))

outliers %>% filter(outlier) %>% select(id, variable, value) %>%
  arrange(variable, value) %>% print(n = Inf)

# ---- explore-hist ----
stroke_long <- stroke %>%
  select(outcome, all_of(num_vars)) %>%
  pivot_longer(-outcome, names_to = "variable")

ggplot(stroke_long, aes(value)) +
  geom_histogram(bins = 20, fill = "steelblue", colour = "white") +
  facet_wrap(~ variable, scales = "free", ncol = 4) +
  labs(title = "Histograms of numeric variables", x = NULL, y = "Count")

# ---- explore-box ----
ggplot(stroke_long, aes(x = "", y = value)) +
  geom_boxplot(fill = "lightblue", outlier.colour = "red") +
  facet_wrap(~ variable, scales = "free", ncol = 4) +
  labs(title = "Boxplots of numeric variables (red = IQR outliers)", x = NULL, y = NULL)

# ---- explore-box-outcome ----
ggplot(stroke_long, aes(outcome, value, fill = outcome)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(alive = "#4C9F70", dead = "#C0504D")) +
  labs(title = "Numeric variables by outcome", x = NULL, y = NULL) +
  theme(legend.position = "none")

# ---- explore-bar ----
cat_long <- stroke %>%
  select(outcome, all_of(cat_vars)) %>%
  pivot_longer(-outcome, names_to = "variable", values_to = "level")

ggplot(cat_long, aes(level)) +
  geom_bar(fill = "steelblue") +
  facet_wrap(~ variable, scales = "free", ncol = 4) +
  labs(title = "Bar charts of categorical variables", x = NULL, y = "Count")

# ---- explore-bar-outcome ----
ggplot(cat_long, aes(level, fill = outcome)) +
  geom_bar(position = "fill") +
  facet_wrap(~ variable, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = c(alive = "#4C9F70", dead = "#C0504D")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion who died, by category", x = NULL, y = "Percent")

# =====================================================================
# 3. EXPLORE (cont.): descriptive table
# =====================================================================

# ---- tbl-descriptive ----
stroke %>%
  select(-id, -dead) %>%
  tbl_summary(by = outcome, label = var_labels,
              statistic = all_continuous() ~ "{median} ({p25}, {p75})") %>%
  add_overall() %>%
  add_p() %>%
  bold_p() %>%
  modify_spanning_header(all_stat_cols(FALSE) ~ "**Outcome at discharge**")

# =====================================================================
# 4. ASSUME: binomial outcome, linearity, specification
# =====================================================================

# ---- assume-binary ----
table(stroke$outcome, stroke$dead)   # confirms 0/1 coding of the event

# =====================================================================
# 5. APPLY: univariable, selection, multivariable (no interaction)
# =====================================================================

# ---- apply-univariable ----
candidates <- c(num_vars, cat_vars)
m_null <- glm(outcome ~ 1, data = stroke, family = binomial)

uni <- map_dfr(candidates, function(v) {
  m   <- glm(reformulate(v, "outcome"), data = stroke, family = binomial)
  lrt <- anova(m_null, m, test = "LRT")
  tidy(m, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    mutate(variable = v, lrt_chisq = lrt$Deviance[2], lrt_p = lrt$`Pr(>Chi)`[2])
}) %>%
  transmute(variable, term,
            beta = estimate, beta_low = conf.low, beta_high = conf.high,
            OR = exp(estimate), OR_low = exp(conf.low), OR_high = exp(conf.high),
            wald_p = p.value, lrt_chisq, lrt_p) %>%
  arrange(lrt_p)

uni %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% print(n = Inf)

# ---- apply-select ----
screen <- uni %>% distinct(variable, lrt_p) %>%
  mutate(p_lt_0.25 = lrt_p < 0.25,
         decision = case_when(
           variable == "gcs"           ~ "Excluded: mediator",
           variable %in% c("sex", "age") ~ "Included: demographic adjustment",
           p_lt_0.25                   ~ "Included: p < 0.25",
           TRUE                        ~ "Not carried forward"))
screen %>% mutate(lrt_p = signif(lrt_p, 3)) %>% print(n = Inf)

selected <- screen %>% filter(str_starts(decision, "Included")) %>% pull(variable)
selected

# ---- apply-full ----
m_full <- glm(reformulate(selected, "outcome"), data = stroke, family = binomial)
tidy(m_full, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# ---- apply-backward ----
forced <- c("sex", "age")
m_cur  <- m_full
log_steps <- list()
repeat {
  dr <- drop1(m_cur, test = "LRT") %>% as.data.frame() %>%
    rownames_to_column("variable") %>% filter(variable != "<none>")
  cand <- dr %>% filter(!variable %in% c(forced, unlist(map(log_steps, "kept")))) %>%
    arrange(desc(`Pr(>Chi)`)) %>% slice(1)
  if (nrow(cand) == 0 || cand$`Pr(>Chi)` < 0.05) break
  m_try <- update(m_cur, as.formula(paste(". ~ . -", cand$variable)))
  sig_terms <- tidy(m_cur) %>% filter(term != "(Intercept)", p.value < 0.05) %>% pull(term)
  sig_terms <- intersect(sig_terms, names(coef(m_try)))
  pct <- abs((coef(m_try)[sig_terms] - coef(m_cur)[sig_terms]) / coef(m_cur)[sig_terms]) * 100
  confounder <- length(pct) > 0 && max(pct) > 20
  log_steps[[length(log_steps) + 1]] <- tibble(
    step = length(log_steps) + 1, variable = cand$variable,
    lrt_p = round(cand$`Pr(>Chi)`, 3),
    max_change_pct = round(if (length(pct)) max(pct) else 0, 1),
    action = if (confounder) "kept (confounder)" else "removed",
    kept = if (confounder) cand$variable else NA_character_)
  if (!confounder) m_cur <- m_try
}
bind_rows(log_steps) %>% select(-kept)
m_prelim <- m_cur
formula(m_prelim)

# ---- apply-addback ----
excluded <- setdiff(candidates, c(all.vars(formula(m_prelim))[-1], "gcs"))
key_terms <- tidy(m_prelim) %>% filter(term != "(Intercept)", p.value < 0.05) %>% pull(term)
map_dfr(excluded, function(v) {
  m_add <- update(m_prelim, as.formula(paste(". ~ . +", v)))
  lrt <- anova(m_prelim, m_add, test = "LRT")
  pct <- abs((coef(m_add)[key_terms] - coef(m_prelim)[key_terms]) / coef(m_prelim)[key_terms]) * 100
  tibble(variable = v, lrt_p = round(lrt$`Pr(>Chi)`[2], 3),
         max_change_pct = round(max(pct), 1))
}) %>% arrange(lrt_p)

# ---- apply-linearity ----
num_in_model <- intersect(num_vars, all.vars(formula(m_prelim)))

# (a) empirical logits by quartile
map_dfr(num_in_model, function(v) {
  stroke %>% mutate(q = cut_number(.data[[v]], 4)) %>%
    group_by(q) %>%
    summarise(n = n(), died = sum(dead), .groups = "drop") %>%
    mutate(variable = v, emp_logit = log((died + 0.5) / (n - died + 0.5)))
}) %>% mutate(emp_logit = round(emp_logit, 3)) %>% select(variable, everything())

# (b) LRT for an added quadratic term
map_dfr(num_in_model, function(v) {
  m_q <- update(m_prelim, as.formula(paste0(". ~ . + I(", v, "^2)")))
  tibble(variable = v, quad_lrt_p = round(anova(m_prelim, m_q, test = "LRT")$`Pr(>Chi)`[2], 3))
})

# ---- apply-linearity-plot ----
map_dfr(num_in_model, function(v) {
  stroke %>% mutate(q = cut_number(.data[[v]], 4)) %>% group_by(q) %>%
    summarise(mid = median(.data[[v]]), n = n(), died = sum(dead), .groups = "drop") %>%
    mutate(variable = v, emp_logit = log((died + 0.5) / (n - died + 0.5)))
}) %>%
  ggplot(aes(mid, emp_logit)) + geom_point(size = 3) + geom_line() +
  facet_wrap(~ variable, scales = "free_x") +
  labs(x = "Quartile median", y = "Empirical logit", title = "Linearity in the logit")

# ---- apply-vif ----
vif(m_prelim)

# ---- apply-final ----
m_final <- m_prelim
summary(m_final)
glance(m_final)
anova(m_null, m_final, test = "LRT")   # overall model test

# ---- helper ----
g_lrt <- anova(m_null, m_final, test = "LRT")

# =====================================================================
# 6. ASSESS: Hosmer-Lemeshow, diagnostics, influence, ROC, accuracy
# =====================================================================

# ---- assess-hl ----
hl <- map_dfr(c(10, 8, 6), function(g) {
  h <- hoslem.test(stroke$dead, fitted(m_final), g = g)
  tibble(groups = g, C_hat = round(unname(h$statistic), 3),
         df = unname(h$parameter), p_value = round(h$p.value, 4))
})
hl

h10 <- hoslem.test(stroke$dead, fitted(m_final), g = 10)
cbind(h10$observed, round(h10$expected, 1))   # observed vs expected by decile

# ---- assess-diag ----
p_par <- length(coef(m_final))          # parameters incl. intercept
n     <- nrow(stroke)

diag <- augment(m_final, type.predict = "response", type.residuals = "pearson") %>%
  mutate(id          = stroke$id,
         std_pearson = .resid / sqrt(1 - .hat),
         dchisq      = std_pearson^2,
         dbeta       = std_pearson^2 * .hat / (1 - .hat),
         leverage    = .hat)

cutoffs <- c(leverage = 2 * p_par / n, dchisq = 4, dbeta = 1, cooksd = 4 / n)
round(cutoffs, 3)

diag %>% select(.fitted, leverage, std_pearson, dchisq, dbeta, .cooksd) %>%
  summary()

# ---- assess-diag-plots ----
diag_long <- diag %>%
  select(id, .fitted, leverage, dchisq, dbeta) %>%
  pivot_longer(c(leverage, dchisq, dbeta), names_to = "statistic") %>%
  mutate(cut = cutoffs[statistic],
         statistic = factor(statistic, c("leverage", "dchisq", "dbeta"),
                            c("Leverage (h)", "Δχ² (dchisquare)", "Δβ (dbeta)")))

ggplot(diag_long, aes(.fitted, value)) +
  geom_point(alpha = 0.6) +
  geom_hline(aes(yintercept = cut), linetype = "dashed", colour = "red") +
  geom_text(data = ~ filter(.x, value > cut), aes(label = id),
            size = 2.6, vjust = -0.6, check_overlap = TRUE) +
  facet_wrap(~ statistic, scales = "free_y", ncol = 1) +
  labs(x = "Fitted probability of death", y = NULL,
       title = "Diagnostic statistics vs fitted probability (dashed = cut-off)")

# ---- assess-diag-bubble ----
ggplot(diag, aes(.fitted, dchisq, size = dbeta, colour = factor(outcome))) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 4, linetype = "dashed") +
  scale_colour_manual(values = c(alive = "#4C9F70", dead = "#C0504D"), name = "Outcome") +
  labs(x = "Fitted probability", y = "Δχ²", size = "Δβ",
       title = "Δχ² vs fitted probability, sized by Δβ")

# ---- assess-influential ----
influential <- diag %>%
  mutate(high_leverage = leverage > cutoffs["leverage"],
         high_dchisq   = dchisq   > cutoffs["dchisq"],
         high_dbeta    = dbeta    > cutoffs["dbeta"],
         high_cooks    = .cooksd  > cutoffs["cooksd"],
         n_flags = high_leverage + high_dchisq + high_dbeta + high_cooks) %>%
  filter(n_flags > 0) %>%
  arrange(desc(dbeta))

influential %>%
  select(id, outcome, all_of(all.vars(formula(m_final))[-1]),
         fitted = .fitted, leverage, dchisq, dbeta, cooksd = .cooksd, n_flags) %>%
  mutate(across(where(is.double), ~ round(.x, 3))) %>%
  print(n = Inf)

tibble(criterion = c("Leverage > 2p/n", "Δχ² > 4", "Δβ > 1", "Cook's D > 4/n"),
       n_flagged = c(sum(influential$high_leverage), sum(influential$high_dchisq),
                     sum(influential$high_dbeta), sum(influential$high_cooks)))

# ---- assess-sensitivity ----
drop_ids <- if (any(diag$dbeta > 1)) diag$id[diag$dbeta > 1] else
  diag %>% slice_max(dbeta, n = 5) %>% pull(id)
drop_ids

m_sens <- update(m_final, data = filter(stroke, !id %in% drop_ids))

full_join(
  tidy(m_final, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(term, OR_all = estimate, CI_all = sprintf("%.2f–%.2f", conf.low, conf.high)),
  tidy(m_sens, exponentiate = TRUE, conf.int = TRUE) %>%
    transmute(term, OR_sens = estimate, CI_sens = sprintf("%.2f–%.2f", conf.low, conf.high)),
  by = "term") %>%
  filter(term != "(Intercept)") %>%
  mutate(pct_change = round(100 * (OR_sens - OR_all) / OR_all, 1),
         across(c(OR_all, OR_sens), ~ round(.x, 3)))

# ---- assess-roc ----
roc_final <- roc(stroke$dead, fitted(m_final), quiet = TRUE)
auc_ci <- ci.auc(roc_final)
auc_ci

ggroc(roc_final, legacy.axes = TRUE, colour = "steelblue", linewidth = 1) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.65, y = 0.2, size = 5,
           label = sprintf("AUC = %.3f (95%% CI %.3f–%.3f)", auc_ci[2], auc_ci[1], auc_ci[3])) +
  labs(title = "ROC curve: final model", x = "1 − Specificity", y = "Sensitivity")

# ---- assess-accuracy ----
youden <- coords(roc_final, "best", best.method = "youden",
                 ret = c("threshold", "sensitivity", "specificity"))

classify <- function(cut) {
  pred <- factor(if_else(fitted(m_final) >= cut, "dead", "alive"), c("alive", "dead"))
  tab  <- table(Predicted = pred, Observed = stroke$outcome)
  TP <- tab["dead", "dead"]; TN <- tab["alive", "alive"]
  FP <- tab["dead", "alive"]; FN <- tab["alive", "dead"]
  list(table = tab,
       metrics = tibble(cut_off = round(cut, 3),
                        accuracy = (TP + TN) / sum(tab),
                        sensitivity = TP / (TP + FN), specificity = TN / (TN + FP),
                        PPV = TP / (TP + FP), NPV = TN / (TN + FN)))
}
c05 <- classify(0.5)
cyj <- classify(youden$threshold[1])

c05$table
cyj$table
acc_tbl <- bind_rows(`0.5` = c05$metrics, Youden = cyj$metrics, .id = "rule") %>%
  mutate(across(accuracy:NPV, ~ round(.x, 3)))
acc_tbl

# =====================================================================
# 7. PRESENT: tables of log-odds and odds ratios with 95% CI
# =====================================================================

# ---- tbl-univariable ----
uni_data <- stroke %>% select(outcome, all_of(candidates))
uni_beta <- tbl_uvregression(uni_data, method = glm, y = outcome,
                             method.args = list(family = binomial),
                             exponentiate = FALSE, label = var_labels,
                             hide_n = TRUE) %>%
  modify_column_hide(p.value)
uni_or <- tbl_uvregression(uni_data, method = glm, y = outcome,
                           method.args = list(family = binomial),
                           exponentiate = TRUE, label = var_labels) %>%
  add_global_p() %>% bold_p(t = 0.25)
tbl_merge(list(uni_beta, uni_or),
          tab_spanner = c("**Log-odds (β)**", "**Crude odds ratio**")) %>%
  modify_footnote_header("p-value: likelihood ratio test; bold = p < 0.25 (screening threshold)",
                         columns = p.value_2, replace = TRUE)

# ---- tbl-multivariable ----
mv_beta <- tbl_regression(m_final, exponentiate = FALSE, label = var_labels,
                          intercept = TRUE) %>%
  modify_column_hide(p.value)
mv_or <- tbl_regression(m_final, exponentiate = TRUE, label = var_labels,
                        intercept = TRUE) %>%
  add_global_p(keep = TRUE) %>% bold_p()
tbl_merge(list(mv_beta, mv_or),
          tab_spanner = c("**Log-odds (β)**", "**Adjusted odds ratio**"))

# ---- tbl-multivariable-broom ----
lrt_terms <- drop1(m_final, test = "LRT") %>% as.data.frame() %>%
  rownames_to_column("variable") %>% filter(variable != "<none>") %>%
  select(variable, LRT_chisq = LRT, LRT_p = `Pr(>Chi)`)
tidy(m_final, conf.int = TRUE) %>%
  transmute(term,
            `β (95% CI)`  = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
            SE = round(std.error, 3),
            `OR (95% CI)` = if_else(term == "(Intercept)", "—",
                                    sprintf("%.3f (%.3f, %.3f)", exp(estimate),
                                            exp(conf.low), exp(conf.high))),
            `Wald p` = format.pval(p.value, digits = 3, eps = 0.001)) %>%
  knitr::kable()
lrt_terms %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% knitr::kable()

# ---- tbl-crude-adjusted ----
final_vars <- all.vars(formula(m_final))[-1]
crude_tbl <- stroke %>% select(outcome, all_of(final_vars)) %>%
  tbl_uvregression(method = glm, y = outcome,
                   method.args = list(family = binomial),
                   exponentiate = TRUE, label = var_labels, hide_n = TRUE)
adj_tbl <- tbl_regression(m_final, exponentiate = TRUE, label = var_labels)
tbl_merge(list(crude_tbl, adj_tbl), tab_spanner = c("**Crude OR**", "**Adjusted OR**"))

# ---- tbl-performance ----
tribble(
  ~Measure, ~Value,
  "Deviance (df)", sprintf("%.2f (%d)", deviance(m_final), df.residual(m_final)),
  "AIC", sprintf("%.2f", AIC(m_final)),
  "Model LRT vs null", sprintf("G = %.2f, df = %d, p %s", g_lrt$Deviance[2], g_lrt$Df[2],
                               format.pval(g_lrt$`Pr(>Chi)`[2], eps = 0.001)),
  "Hosmer–Lemeshow, g = 10", sprintf("Ĉ = %.3f, df = %d, p = %.4f", hl$C_hat[1], hl$df[1], hl$p_value[1]),
  "Hosmer–Lemeshow, g = 8",  sprintf("Ĉ = %.3f, df = %d, p = %.4f", hl$C_hat[2], hl$df[2], hl$p_value[2]),
  "Hosmer–Lemeshow, g = 6",  sprintf("Ĉ = %.3f, df = %d, p = %.4f", hl$C_hat[3], hl$df[3], hl$p_value[3]),
  "AUC (95% CI)", sprintf("%.3f (%.3f–%.3f)", auc_ci[2], auc_ci[1], auc_ci[3]),
  "Accuracy / Sens / Spec at 0.5",
    sprintf("%.1f%% / %.1f%% / %.1f%%", 100*acc_tbl$accuracy[1], 100*acc_tbl$sensitivity[1], 100*acc_tbl$specificity[1]),
  sprintf("Accuracy / Sens / Spec at Youden (%.3f)", acc_tbl$cut_off[2]),
    sprintf("%.1f%% / %.1f%% / %.1f%%", 100*acc_tbl$accuracy[2], 100*acc_tbl$sensitivity[2], 100*acc_tbl$specificity[2]),
  "Influential observations (Δβ > 1)", as.character(sum(diag$dbeta > 1))
) %>% knitr::kable()

# ---- helper ----
or_tab <- tidy(m_final, exponentiate = TRUE, conf.int = TRUE)
fmt_or <- function(pattern) {
  r <- filter(or_tab, str_detect(term, pattern))
  p_txt <- if (r$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", r$p.value)
  sprintf("%.2f (95%% CI %.2f, %.2f; %s)", r$estimate, r$conf.low, r$conf.high, p_txt)
}

# =====================================================================
# 8. PUBLISH: session information
# =====================================================================

# ---- session-info ----
sessionInfo()
