# AGENTS.md — Stroke Fatality Logistic Regression Project

Rules for any AI agent (or human) working in this folder.

## Project

- **Goal:** binary logistic regression of stroke fatality at discharge, following the RREAAPP workflow (see `implementation.md`).
- **Reference:** *Logistic Regression Analysis in Health and Medicine*, Kamarul Imran Musa (2026), https://drkamarul.github.io/logistic_regression_2026/logistic-regression-chapter.html
- **Data:** `data/stroke-fatality-cooked2.dta`, a teaching dataset: 226 patients, 22 variables.
- **Outcome:** `outcome` (0 = alive, 1 = dead at discharge). The modelled event is **dead**.

## Folder layout

```
MLR/
├── AGENTS.md                 # this file: rules
├── implementation.md         # RREAAPP plan
├── index.html                # GitHub Pages landing page (https://nsrshifaramli.github.io/stroke-MLR/)
├── data/                     # raw data (read-only)
│   └── stroke-fatality-cooked2.dta
├── analysis/                 # Quarto/R source code
│   └── stroke_rreaapp.qmd
├── reports/                  # rendered outputs (HTML)
│   └── stroke_rreaapp.html
└── practice/                 # earlier practice files (chapter model with interaction)
    ├── logistic_practice.R
    ├── logistic_practice_report.qmd
    └── logistic_practice_report.html
```

- Code in `analysis/` and `practice/` reads data with the relative path `../data/stroke-fatality-cooked2.dta`.
- New analysis code goes in `analysis/` and rendered outputs go in `reports/`. Never write into `data/`.
- To render the main report, run this from the `MLR/` folder:

  ```bash
  quarto render analysis/stroke_rreaapp.qmd && mv analysis/stroke_rreaapp.html reports/
  ```

## Workflow: RREAAPP (follow in order)

1. **Reason:** state the question, the outcome and the candidate predictors.
2. **Read:** import the data with `haven`, then rename and label the variables.
3. **Explore:** look for missing data, outliers and impossible values. Plot numeric variables with histograms and boxplots, and categorical variables with bar charts.
4. **Assume:** state the assumptions and how each is checked: binomial outcome, linearity in the logit, a correctly specified model, independence.
5. **Apply:** univariable models, then variable selection, then a multivariable model with **no interaction terms**.
6. **Assess:** goodness of fit, diagnostics, influential observations, ROC/AUC, Hosmer–Lemeshow test and accuracy.
7. **Present:** tables of log-odds (β) and odds ratios with 95% CIs.
8. **Publish:** render the Quarto report to self-contained HTML.

Do not skip a step or merge steps. Each step is its own section in the report.

## Tools and packages

- R (≥ 4.3) and Quarto. Quarto is bundled with RStudio: `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto`.
- Required packages: `haven`, `tidyverse`, `gtsummary`, `broom`, `ResourceSelection` (Hosmer–Lemeshow), `pROC` (ROC/AUC), `car`, `parameters` (needed by `gtsummary::add_global_p`), `mgcv` (optional linearity check).
- Use **tidyverse** for data handling and plotting (`ggplot2`), **gtsummary** for tables, and **broom** (`tidy()`, `glance()`, `augment()`) for model output and diagnostics.
- Fit models with `glm(..., family = binomial(link = "logit"))`.

## Data rules

- **Never modify, overwrite or delete `data/stroke-fatality-cooked2.dta`.** All cleaning happens in code.
- Keep every record. Flag outliers and suspicious values in a table but **do not drop or change them**. Do not set a value to `NA` without the user's approval.
- Influential observations are **reported**. They are removed only in a separately labelled sensitivity analysis, never in the main model.
- Never change data to improve model fit.
- Treat the dataset as teaching data. Do not present results as epidemiological findings.

## Modelling rules

- Variable selection follows **purposeful selection**: univariable p < 0.25 (LRT) plus clinical judgement, then backward elimination by LRT with a check for confounding (a coefficient change of more than 20%).
- **Do not include GCS** in the multivariable model. It is a mediator (stroke type → GCS → death). Report it in the univariable step only.
- No interaction terms in this analysis.
- Report β, OR, 95% CI and p-values. Use likelihood ratio tests to compare models.

## Coding style

- Use one Quarto file for the report (`analysis/stroke_rreaapp.qmd`), with code chunks named after their RREAAPP step.
- Use the native pipe `|>` or `%>%` consistently (default `%>%`). Use snake_case names.
- Set `set.seed()` before anything random. Use relative file paths only.
- Every table or figure gets a caption and 1–3 sentences of interpretation.

## Before finishing any task

- Render the report, move the HTML to `reports/`, and confirm it builds without errors.
- The repository is **public** and GitHub Pages serves `main` from the root. Never commit the `.dta` data file.
- Check that the outcome count is 62 dead / 164 alive (n = 226).
- Do not run the analysis unless the user has approved the plan in `implementation.md`.
