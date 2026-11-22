# Implementation Plan: Stroke Fatality Logistic Regression (RREAAPP)

**Status:** Approved and executed on 2026-11-22. Report: `reports/stroke_rreaapp.html`.
**Deliverable:** `analysis/stroke_rreaapp.qmd`, rendered to a self-contained `reports/stroke_rreaapp.html`.
**Rules:** see `AGENTS.md`.

## Decisions agreed with the user

| # | Topic | Decision |
|---|---|---|
| 1 | Variable selection | Purposeful selection (univariable p < 0.25 + clinical judgement, then backward elimination by LRT). **GCS is excluded** because it is a mediator. |
| 2 | Outliers / errors | Flag and report in a table; **keep all records**. |
| 3 | Publish | Quarto → self-contained HTML. |
| 4 | Accuracy cut-off | Report at **0.5** and at the **Youden-optimal** cut-off from the ROC curve. |
| 5 | Variable names | Rename to descriptive snake_case **and** convert 0/1 codes to labelled factors. |
| 6 | Influential obs. | Flag with standard cut-offs; refit without them only as a **sensitivity analysis**. |
| 7 | AGENTS.md | Rules for any agent working in this folder (created). |

---

## Step 1: Reason

- **Research question:** which admission characteristics are associated with death before discharge in patients with acute stroke?
- **Outcome:** `outcome` (dead = 1 / alive = 0). The model estimates the odds of **death**.
- **Candidate predictors (20):** sex, age, ethnicity (Malay), marital status, GCS, SBP, DBP, heart rate, haemoglobin, platelets, WBC, sodium, potassium, blood sugar, urea, hypertension, diabetes, dyslipidaemia, referral source, stroke type.
- **Analysis goal:** estimate associations (adjusted ORs), not build a prediction model. This is why the mediator GCS is left out.
- **Design note:** one record per admission, so observations are independent.

## Step 2: Read

```r
library(haven); library(tidyverse)
raw <- read_dta("../data/stroke-fatality-cooked2.dta")
```

1. Inspect with `glimpse()`, `names()` and the variable labels (`attr(x, "label")`).
2. **Rename** to descriptive names:

| Original | New | Original | New |
|---|---|---|---|
| rn | id | plat | platelet |
| sex | sex | wbc | wbc |
| gcs | gcs | na | sodium |
| sbp | sbp | potas | potassium |
| dbp | dbp | cbs | blood_sugar |
| hr | heart_rate | urea | urea |
| hb | haemoglobin | age | age |
| outcome | outcome | maritalstat | married |
| hypertension | hypertension | malay | malay |
| diabetes | diabetes | dyslipidaemia | dyslipidaemia |
| referralcase | referral | stroketype | stroke_type |

3. **Convert** labelled variables with `as_factor()`, then set clear levels and reference categories:
   - `outcome`: alive (ref), dead
   - `sex`: male (ref), female
   - `stroke_type`: Ischaemic (ref), Haemorrhagic
   - `referral`: Hospital (ref), GP/Home
   - yes/no variables: no (ref), yes
4. Keep a numeric 0/1 copy (`dead`), which the Hosmer–Lemeshow and ROC functions need.
5. Check that the outcome counts are 62 dead and 164 alive.

## Step 3: Explore

**3a. Missing data:** count missing values per variable (`summarise(across(everything(), ~ sum(is.na(.x))))`). Expected: none.

**3b. Possible errors:** compare each numeric variable with plausible clinical ranges and flag values outside them:

| Variable | Plausible range |
|---|---|
| age | 18–110 years |
| gcs | 3–15 (integer) |
| sbp / dbp | 60–300 / 30–200 mmHg; check that dbp < sbp |
| heart_rate | 30–220 bpm |
| haemoglobin | 4–22 g/dL |
| platelet | 10–1000 × 10⁹/L |
| wbc | 1–50 × 10⁹/L |
| sodium | 110–170 mmol/L |
| potassium | 2–8 mmol/L |
| blood_sugar | 1–40 mmol/L |
| urea | 1–60 mmol/L |

Also check for duplicate `id`s and for factor levels with no or very few patients, e.g. a "2 = Other" stroke type mentioned in the label.

**3c. Outliers:** use the IQR rule (below Q1 − 1.5·IQR or above Q3 + 1.5·IQR) for each numeric variable. Present a table of flagged patients (id, variable, value). **Report them; no records are removed.**

**3d. Plots (tidyverse/ggplot2)**
- Numeric variables: reshape to long format with `pivot_longer()`, then draw a **histogram** and a **boxplot** for each, using `facet_wrap(scales = "free")`. Also draw boxplots by outcome.
- Categorical variables: draw **bar charts** of counts, and a second set filled by outcome (`geom_bar(position = "fill")`).

**3e. Descriptive table (gtsummary):** `tbl_summary(by = outcome)` with `add_overall()` and `add_p()`. Numbers are shown as median (IQR) and categories as n (%).

## Step 4: Assume

Each assumption is stated together with how it is checked:

| Assumption | Check |
|---|---|
| Binomial outcome | Outcome is binary 0/1; `family = binomial(link = "logit")` |
| Independent observations | One admission per patient; no duplicate ids (Step 3b) |
| Linearity in the logit (numeric predictors) | Empirical logits by quartile; LRT for an added quadratic term; optional GAM smooth (edf ≈ 1 means linear) |
| Model correctly specified | Purposeful selection; no important omitted confounders; Hosmer–Lemeshow test and residuals (Step 6) |
| No multicollinearity | `car::vif()` < 5 |
| Enough events per variable | 62 deaths; aim for ≥ 10 events per parameter, so at most about 6 parameters |

If a numeric predictor is clearly non-linear, report it and discuss categorising or transforming it. Do not change it without telling the user.

## Step 5: Apply

**5a. Univariable logistic regression:** fit `glm(outcome ~ x, family = binomial)` for each of the 20 candidates.
- `gtsummary::tbl_uvregression(exponentiate = TRUE) %>% add_global_p()` gives crude ORs, 95% CIs and LRT p-values.
- Use `broom::tidy()` to get β (log-odds) for each model.

**5b. Identify possible independent variables:**
- Keep variables with univariable LRT p < 0.25, plus clinically important variables (sex and age as demographic adjustments).
- Exclude GCS (mediator) and document the reason.
- Avoid putting closely related variables in together (e.g. SBP and DBP); keep the more clinically relevant one.

**5c. Multivariable logistic regression (main effects only, no interaction):**
1. Fit the full model with all selected variables.
2. Backward elimination: remove the least significant variable (`drop1(test = "LRT")`, p > 0.05). Before dropping it, check whether any remaining β changes by more than 20%; if so, keep it as a confounder.
3. Add back each excluded variable one at a time, to check for significance or confounding.
4. Check linearity of the numeric variables in the final model (Step 4 methods).
5. Report the overall model LRT against the null model.

Expected (from the chapter): the preliminary model will probably be close to stroke type + WBC + referral + sex, but the final set will be whatever the data and the selection steps produce.

## Step 6: Assess

**6a. Goodness of fit**
- **Hosmer–Lemeshow** test (`ResourceSelection::hoslem.test`) with g = 10, plus g = 8 and 6 to check stability.
- Deviance, AIC (`broom::glance()`).

**6b. Diagnostic statistics (broom)**
`augment(final_model, type.predict = "response")` gives `.fitted`, `.hat` (leverage), `.std.resid` and `.cooksd`. Add:
- Standardised Pearson residual `r_s = rstandard(model, type = "pearson")`
- **Δχ² (dchisquare)** = r_s²
- **Δβ (dbeta)** = r_s² × h / (1 − h)
- **Leverage** h = `.hat`

**Plots (Hosmer–Lemeshow style):** Δχ² vs fitted probability, Δβ vs fitted probability, leverage vs fitted probability, and Δχ² vs fitted probability with points sized by Δβ.

**6c. Influential observations**
- Flag patients with: leverage > 2(p+1)/n; |r_s| > 2 (Δχ² > 4); Δβ > 1; Cook's D > 4/n.
- Table of the flagged patients: id, covariates, observed outcome, fitted probability, statistics.
- **Sensitivity analysis:** refit without the most influential patients and compare the ORs side by side. The main model is not changed.

**6d. Discrimination and classification**
- **ROC curve and AUC** with 95% CI (`pROC::roc`, `ci.auc`).
- **Accuracy**, sensitivity, specificity, PPV and NPV from a 2×2 classification table at:
  - cut-off 0.5, and
  - the Youden-optimal cut-off (`coords(roc, "best")`).

## Step 7: Present

1. **Table 1:** descriptive characteristics by outcome (gtsummary).
2. **Table 2:** univariable results: **β (log-odds)**, 95% CI of β, **crude OR**, 95% CI, LRT p.
3. **Table 3:** multivariable results: **β (log-odds)**, 95% CI, **adjusted OR**, 95% CI, Wald p and LRT p.
4. **Table 4:** crude vs adjusted ORs side by side (`tbl_merge`).
5. **Table 5:** model fit and performance: deviance, AIC, Hosmer–Lemeshow (each g), AUC (95% CI), accuracy, sensitivity and specificity at both cut-offs.
6. **Table 6:** influential observations, plus the sensitivity-analysis comparison.
7. Figures: exploration plots, ROC curve, diagnostic plots.
8. A written interpretation of each adjusted OR, in the chapter's style ("holding other variables constant…").

β and its CI come from `broom::tidy(conf.int = TRUE)`, ORs from `exponentiate = TRUE`, and formatted tables from `gtsummary::tbl_regression()` with `exponentiate = FALSE` and `TRUE`, merged together.

## Step 8: Publish

- Write everything in `analysis/stroke_rreaapp.qmd`: one top-level section per RREAAPP step, with table of contents, numbered sections, foldable code and `embed-resources: true`.
- Render with the RStudio-bundled Quarto: `quarto render analysis/stroke_rreaapp.qmd`, then move `stroke_rreaapp.html` into `reports/`.
- Check that it renders without errors and that the outcome counts are 62/164.
- Include the teaching-dataset disclaimer and the session info (`sessionInfo()`) at the end.

## Files

| File | Purpose | Status |
|---|---|---|
| `data/stroke-fatality-cooked2.dta` | Raw data (read-only) | exists |
| `AGENTS.md` | Agent rules | created |
| `implementation.md` | This plan | created |
| `analysis/stroke_rreaapp.qmd` | Analysis + report source | created |
| `reports/stroke_rreaapp.html` | Rendered report | created |
| `practice/logistic_practice.R`, `practice/logistic_practice_report.*` | Earlier practice files (with interaction) | moved; data path updated |

## Packages to confirm before running

Installed: haven, tidyverse, broom, gtsummary, mgcv, ResourceSelection, pROC, car, parameters. Nothing else needs installing.
