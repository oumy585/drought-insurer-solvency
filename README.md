# When Does Drought Risk Matter for Insurer Solvency? A Portfolio-Structure Perspective

Replication code for the paper by Oumy Dione (Université Paris Nanterre, EconomiX, and Lendys), submitted to the *European Actuarial Journal*.

The code builds a commune-level drought hazard index from the Soil Wetness Index (SWI), estimates the two-stage CatNat procedure (request, then recognition), converts recognitions into an annual climate loss calibrated on ONRN costs, simulates financial losses with a random-vector functional link (RVFL) economic scenario generator, and assembles total losses (VaR and Expected Shortfall at 99.5 %) for a stylised insurer over a grid of market shares, balance-sheet sizes, hazard shifts and recognition regimes.

## Contents

| Step | File | Purpose |
|---|---|---|
| 0 | `build_swi_communal.py` | Monthly communal SWI (area-weighted) and event-level SWI minimum for CatNat drought events |
| 1 | `01_preparation_CatNat_ONRN_couts.R` | CatNat decrees, ONRN cost classes, two-part cost model, cost weights |
| 2 | `02_construction_alea_SWI.R` | Recognition-rate series, regime-break tests, renormalised communal SWI |
| 3 | `04_alea_SWI_et_transmission_CatNat.R` | Hazard index, commune-year request and recognition models, hazard figures |
| 4 | `05_perte_climatique.R` | Annual climate loss calibrated on the ONRN bridge (scale factor and residuals) |
| 5 | `06_macro_esg.R` | RVFL economic scenario generator and stylised-portfolio financial loss |
| 6 | `07_assemblage.R` | Total loss, VaR/ES grid, robustness checks, tables and figures |

Script numbering follows the order of the project; there is no script `03` in this repository.

## Order of execution

Run all scripts from the repository root, in this order:

1. `build_swi_communal.py`
2. `01_preparation_CatNat_ONRN_couts.R`
3. `02_construction_alea_SWI.R`
4. `04_alea_SWI_et_transmission_CatNat.R`
5. `05_perte_climatique.R`
6. `06_macro_esg.R`
7. `07_assemblage.R`

Intermediate files are written to `output/etape0` (scripts 01 to 05 and 07) and `data/macro` (script 06); `07_assemblage.R` reads both.

## Data layout

Input data are not included in this repository (size and licensing). Place them as follows:

```
data/
  onrn/      ONRN cost files (file names are set at the top of script 01)
  catnat/    Arretes-CatNat-commune-1.xlsx
             CatNat_drought_final_analysis.xlsx
             CatNat_secheresse_couvertes_clean.xlsx
             CatNat_drought_events_SWI.xlsx      (written by build_swi_communal.py)
  swi/       SWI_2000_2025.csv                    (monthly SWI by grid cell)
             mapping_SWI_communes_weighted.csv    (commune - grid cell area weights)
  macro/     ESG_data_2000_2024_CAC40corrige.rds  (monthly macro-financial series)
output/
  etape0/    created by the scripts
```

Sources:
- CatNat decrees of natural-disaster recognition and ONRN cost data (French public data).
- SWI: Météo-France, monthly SWI by grid cell, 2000-2025.
- Macro-financial series (CAC 40 returns, interest rates, inflation): proprietary data (Bloomberg and CCR). They are **not** included. Without `ESG_data_2000_2024_CAC40corrige.rds`, script 06 cannot run, and script 07 stops at its robustness section (Section 7 of the script). The grid, VaR/ES results and figures of script 07 before that section only need `simulation_financiere_reference.rds`, produced by script 06.

## Reproducibility

- Seeds are set in each script: 7 (script 02), 42 and 1 (script 05), 2028 (script 06), 5000, 7000+i, 9000 and 11000+i (script 07).
- Reference simulation (script 06): RVFL with ridge penalties on both layers, `ahead` R package, 5 lags, 20 Sobol nodes, ReLU activation, penalties lambda1 = 2.05 and lambda2 = 7.1 fixed after Bayesian optimisation, moving block bootstrap with blocks of 6 months, 10,000 simulations. Other block lengths (3, 9 and 12 months) were compared in sensitivity runs; 12 months was rejected because it produced too few distinct simulated losses.
- Reference configuration (script 07): market share Q = 5 %, insurer assets of EUR 20 billion, no hazard shift, recognition regime from 2017.
- Software: R (version: *fill in*), Python 3 with `pandas`, `numpy`, `pyarrow`, `openpyxl`. R packages: `dplyr`, `tidyr`, `stringr`, `lubridate`, `readxl`, `Matrix`, `ggplot2`, `patchwork`, `scales`, `ahead`; optional: `sandwich`, `lmtest`.

## Citation

If you use this code, please cite the paper and the archived version of this repository:

> Dione, O. (2026). *Code for "When Does Drought Risk Matter for Insurer Solvency? A Portfolio-Structure Perspective"* (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23505966

## Licence

MIT, see `LICENSE`.
