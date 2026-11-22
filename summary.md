# Summary: Stroke Fatality at Discharge (Binary Logistic Regression)

**Online report:** <https://nsrshifaramli.github.io/stroke-MLR/>
**Code:** <https://github.com/nsrshifaramli/stroke-MLR>
**Method reference:** *Logistic Regression Analysis in Health and Medicine*, Kamarul Imran Musa (2026)

> Teaching dataset with values modified by its author. The results illustrate method only. They are not evidence about stroke care.

## 1. Reason

- **Question:** which admission characteristics are associated with death before discharge after acute stroke?
- **Outcome:** `outcome`, status at discharge (alive = 0, dead = 1).
- **Goal:** estimate adjusted associations (odds ratios), not build a prediction tool.
- **Key decision:** the Glasgow Coma Scale (GCS) is a **mediator** (stroke type → depressed consciousness → death). It is examined on its own but **left out of the multivariable model**.

## 2. Read

- Imported `stroke-fatality-cooked2.dta` with `haven::read_dta()`: **226 patients, 22 variables**.
- Renamed variables to descriptive names (e.g. `na` → `sodium`, `cbs` → `blood_sugar`, `referralcase` → `referral`).
- Converted the Stata labels to factors and set reference groups: alive, male, ischaemic, hospital referral, "no".
- **62 deaths (27.4%)** and 164 survivors.

## 3. Explore

| Check | Result |
|---|---|
| Missing data | None |
| Duplicate patient IDs | None |
| DBP ≥ SBP, non-integer GCS | None |
| Values outside plausible clinical range | 2 (sodium 109 and 177 mmol/L): flagged, kept |
| IQR outliers | 72 values in 52 patients, mostly sodium and urea: plausible extremes, kept |

- Histograms and boxplots for the 12 numeric variables, overall and by outcome.
- Bar charts for the 8 categorical variables, with the proportion who died in each category.
- Table 1 (descriptive statistics by outcome) made with `gtsummary::tbl_summary()`.

## 4. Assume

| Assumption | Check | Result |
|---|---|---|
| Binomial outcome | Binary 0/1, `family = binomial` | Met |
| Independent observations | One admission per patient | Met |
| Linear in the logit | Empirical logits by quartile, quadratic LRT | Met (WBC p = 0.22, age p = 0.16) |
| Correct specification | Purposeful selection, fit tests, diagnostics | Borderline (see Hosmer–Lemeshow) |
| No multicollinearity | VIF | Met (all about 1.0–1.2) |
| Enough events | 62 deaths, 5 predictors | About 12 events per predictor |

## 5. Apply

**Univariable screening.** The table shows the variables with LRT p < 0.25:

| Variable | Crude OR (95% CI) |
|---|---|
| GCS (per point) | 0.53 (0.44, 0.61): strongest, excluded as a mediator |
| WBC (per 1 × 10⁹/L) | 1.24 (1.14, 1.36) |
| Haemorrhagic vs ischaemic | 5.08 (2.74, 9.60) |
| GP/Home vs hospital referral | 0.24 (0.13, 0.45) |
| Sodium (per mmol/L) | 1.08 (1.01, 1.15) |
| Dyslipidaemia | 0.28 (0.07, 0.85) |
| Heart rate, haemoglobin, sex, blood sugar | p between 0.06 and 0.21 |

**Variable selection: purposeful selection**
- Carry forward variables with p < 0.25. Exclude GCS. Keep sex and age as demographic adjustments.
- Backward elimination by LRT removed haemoglobin, sodium, blood sugar, dyslipidaemia and heart rate. None changed the remaining coefficients by more than 20%, so none was a confounder.
- Adding the excluded variables back one at a time: none was significant or a confounder.
- No interaction terms (by design).

**Final model:** `outcome ~ wbc + stroke_type + referral + sex + age`

## 6. Assess

| Measure | Result |
|---|---|
| Model LRT vs null | G = 51.39, df = 5, p < 0.001 |
| Deviance / AIC | 214.17 on 220 df / 226.17 |
| Hosmer–Lemeshow, g = 10 | Ĉ = 13.15, p = 0.107 |
| Hosmer–Lemeshow, g = 8 / g = 6 | p = 0.029 / p = 0.009 |
| AUC (95% CI) | **0.801 (0.743–0.859)**: excellent discrimination |
| Accuracy at 0.5 | 73.9% (sensitivity 35.5%, specificity 88.4%) |
| Accuracy at Youden cut-off 0.198 | 70.8% (sensitivity 88.7%, specificity 64.0%) |

**Goodness of fit is borderline.** It passes with 10 groups but fails with 8 and 6. The likely reason is the stroke type × referral interaction, which this main-effect model leaves out by design. In the chapter, adding that interaction raised the Hosmer–Lemeshow p to 0.75.

**Diagnostics** (`broom::augment()` plus leverage, Δχ² and Δβ):
- 13 patients have high leverage (> 2p/n), 10 have Δχ² > 4, and 15 have Cook's D > 4/n.
- No patient has Δβ > 1.
- Sensitivity analysis without the 5 most influential patients: the ORs for WBC, stroke type and referral changed by less than 10%, and the conclusions held.

## 7. Present: final model

| Variable | β (log-odds) | Adjusted OR (95% CI) | p |
|---|---|---|---|
| Intercept | −4.053 | — | < 0.001 |
| WBC (per 1 × 10⁹/L) | 0.183 | **1.20 (1.10, 1.33)** | < 0.001 |
| Haemorrhagic vs ischaemic | 0.942 | **2.57 (1.24, 5.34)** | 0.011 |
| GP/Home vs hospital referral | −0.935 | **0.39 (0.19, 0.80)** | 0.010 |
| Female vs male | 0.209 | 1.23 (0.62, 2.48) | 0.55 |
| Age (per year) | 0.017 | 1.02 (0.99, 1.04) | 0.17 |

**Interpretation** (each OR holds the other variables constant):
- **WBC:** each 1 × 10⁹/L increase is associated with **20% higher odds** of death.
- **Stroke type:** haemorrhagic stroke has about **2.6 times the odds** of death compared with ischaemic stroke.
- **Referral:** patients referred by a GP or from home have **61% lower odds** of death than those transferred from another hospital. Transfers are probably more severe cases.
- **Sex and age:** no significant association after adjustment.
- **Crude vs adjusted:** adjustment pulled the stroke-type OR from 5.08 down to 2.57 and the referral OR from 0.24 up to 0.39. Part of each crude effect was due to the other variables.

## 8. Publish

- Analysis: `analysis/stroke_rreaapp.qmd`, rendered with Quarto to `reports/stroke_rreaapp.html`.
- Public GitHub repository with GitHub Pages: <https://nsrshifaramli.github.io/stroke-MLR/>
- The raw `.dta` file is **not** included in the repository (see `data/README.md`).
- Project rules are in `AGENTS.md` and the full plan in `implementation.md`.

## Key lessons

1. Choosing covariates is a causal question: the strongest predictor (GCS) was left out because it is a mediator.
2. Crude and adjusted ORs can differ a lot, so report both.
3. A borderline Hosmer–Lemeshow result often points to a missing interaction, and the result depends on the number of groups.
4. With an unbalanced outcome, accuracy at a 0.5 cut-off hides poor sensitivity. Report sensitivity, specificity and AUC with it.
5. Flag and report outliers and influential patients. Never delete them to improve the fit.
