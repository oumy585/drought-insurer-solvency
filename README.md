# drought-insurer-solvency

Code for "When Does Drought Risk Matter for Insurer Solvency? A Portfolio-Structure Perspective" (European Actuarial Journal, submitted).

## Contents

| Script | Role |
|---|---|
| `01_preparation_CatNat_ONRN_couts.R` | Preparation of CatNat recognition and ONRN cost data |
| `02_construction_alea_SWI.R` | Construction of the SWI-based hazard index |
| `04_alea_SWI_et_transmission_CatNat.R` | Hazard-to-recognition transmission (demand and recognition) |
| `05_perte_climatique.R` | Climate loss model |
| `06_macro_esg.R` | Economic scenario generator (macro-financial module) |
| `07_assemblage.R` | Aggregation, VaR/ES, 72-configuration grid |

## Data

Climate, municipal and financial data are public. The CCR recognition database was built by the author from a dedicated query and cleaning; it is available on request.

## Reproducibility

Each configuration of the grid uses its own random seed. R version: [...]; main packages: [...].
