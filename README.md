# Stroke Fatality at Discharge: Binary Logistic Regression

A binary logistic regression analysis in R of in-hospital death after acute stroke, following the **RREAAPP** workflow:

**R**eason → **R**ead → **E**xplore → **A**ssume → **A**pply → **A**ssess → **P**resent → **P**ublish

Based on *Logistic Regression Analysis in Health and Medicine* by Kamarul Imran Musa (2026):
<https://drkamarul.github.io/logistic_regression_2026/logistic-regression-chapter.html>

> **Teaching dataset.** The data were modified by the author so that the patterns being taught show clearly. The results illustrate method only. They are not evidence about stroke care.

## Main results (final main-effect model, n = 226, 62 deaths)

| Variable | Adjusted OR (95% CI) | p |
|---|---|---|
| White cell count (per 1 × 10⁹/L) | 1.20 (1.10, 1.33) | < 0.001 |
| Haemorrhagic vs ischaemic stroke | 2.57 (1.24, 5.34) | 0.011 |
| GP/Home vs hospital referral | 0.39 (0.19, 0.80) | 0.010 |
| Female vs male | 1.23 (0.62, 2.48) | 0.55 |
| Age (per year) | 1.02 (0.99, 1.04) | 0.17 |

- **AUC:** 0.80 (0.74–0.86).
- **Hosmer–Lemeshow:** p = 0.107 with 10 groups. The fit is borderline with 8 and 6 groups.
- **GCS:** excluded from the multivariable model because it is a mediator (stroke type → GCS → death).

## Folder layout

```
├── AGENTS.md            # rules for anyone or any AI agent working on the project
├── implementation.md    # the RREAAPP analysis plan and decisions
├── data/                # put stroke-fatality-cooked2.dta here (not included)
├── analysis/
│   └── stroke_rreaapp.qmd   # main analysis (Quarto)
├── reports/
│   ├── index.html           # GitHub Pages landing page
│   └── stroke_rreaapp.html  # rendered report
└── practice/            # earlier practice files (chapter model with interaction)
```

## How to reproduce

1. Put `stroke-fatality-cooked2.dta` in `data/` (from the course materials).
2. Install the R packages:

   ```r
   install.packages(c("haven", "tidyverse", "gtsummary", "broom",
                      "ResourceSelection", "pROC", "car", "parameters", "mgcv"))
   ```

3. Render the report from the project folder:

   ```bash
   quarto render analysis/stroke_rreaapp.qmd && mv analysis/stroke_rreaapp.html reports/
   ```
