# -*- coding: utf-8 -*-
"""
build_swi_communal.py
Monthly communal SWI (area-weighted, no renormalisation) and event-level SWI
minimum for the CatNat drought events.
Run from the repository root (see README for the data/ layout).
Requires: pandas, numpy, pyarrow, openpyxl.
"""

import os
import numpy as np
import pandas as pd

dir_swi = os.path.join("data", "swi")
dir_catnat = os.path.join("data", "catnat")
dir_out = os.path.join("output", "etape0")
os.makedirs(dir_out, exist_ok=True)

# ----------------------------------------------------------
# Inputs
# ----------------------------------------------------------
swi = pd.read_csv(os.path.join(dir_swi, "SWI_2000_2025.csv"))                       # monthly SWI by grid cell
mapping = pd.read_csv(os.path.join(dir_swi, "mapping_SWI_communes_weighted.csv"))   # commune-grid cell weights
df_catnat = pd.read_excel(os.path.join(dir_catnat, "CatNat_secheresse_couvertes_clean.xlsx"))

print(swi.shape, mapping.shape, df_catnat.shape)
print(swi.columns.tolist())
print(mapping.columns.tolist())
print(df_catnat.columns.tolist())

# ----------------------------------------------------------
# 1. Harmonise variable names and types
# ----------------------------------------------------------
swi = swi.rename(columns={"NUMERO": "grid_id", "DATE": "month", "SWI_UNIF_MENS": "swi"})
swi["grid_id"] = swi["grid_id"].astype(int)

mapping = mapping.rename(columns={"num_maille": "grid_id"})
mapping["grid_id"] = mapping["grid_id"].astype(int)
mapping["code_insee"] = mapping["code_insee"].astype(str).str.zfill(5)

print("Period:", swi["month"].min(), "->", swi["month"].max())

# ----------------------------------------------------------
# 2. Monthly communal SWI: weighted sum over grid cells
# ----------------------------------------------------------
municipality_swi = swi.merge(mapping, on="grid_id", how="inner")
print(municipality_swi.shape)

municipality_swi["weighted_swi"] = municipality_swi["swi"] * municipality_swi["weight"]

monthly_swi = (
    municipality_swi
    .groupby(["code_insee", "month"], as_index=False)
    .agg(swi=("weighted_swi", "sum"))
)

print(monthly_swi.shape, monthly_swi["code_insee"].nunique())
print(monthly_swi["month"].min(), monthly_swi["month"].max())

monthly_swi.to_parquet(
    os.path.join(dir_out, "SWI_monthly_communes_2000_2025.parquet"),
    index=False,
    compression="snappy",
)

# ----------------------------------------------------------
# 3. Prepare CatNat events
# ----------------------------------------------------------
dfc = df_catnat.copy()
dfc["code_insee"] = dfc["code_insee"].astype(str).str.zfill(5)

dfc["start_month"] = pd.to_datetime(dfc["Début d'Evénement"]).values.astype("datetime64[M]")
dfc["end_month"] = pd.to_datetime(dfc["Fin d'Evénement"]).values.astype("datetime64[M]")

dfc = dfc.reset_index(drop=True)
dfc["event_id"] = np.arange(len(dfc))
print(dfc.shape)

# ----------------------------------------------------------
# 4. Expand events month by month
# ----------------------------------------------------------
def months_between(row):
    if pd.isna(row.start_month) or pd.isna(row.end_month):
        return []
    return pd.date_range(row.start_month, row.end_month, freq="MS")


expanded = (
    dfc[["event_id", "code_insee", "start_month", "end_month"]]
    .assign(month=dfc.apply(months_between, axis=1))
    .explode("month", ignore_index=True)
    .dropna(subset=["month"])
)
print(expanded.shape)

monthly_swi["month"] = pd.to_datetime(monthly_swi["month"].astype(str), format="%Y%m")

joined = expanded.merge(monthly_swi, on=["code_insee", "month"], how="left")
print(joined.shape)

# ----------------------------------------------------------
# 5. Minimum SWI per event
# ----------------------------------------------------------
idx = joined.groupby("event_id")["swi"].idxmin().dropna().astype(int)

swi_event = (
    joined.loc[idx, ["event_id", "swi", "month"]]
    .rename(columns={"swi": "swi_min", "month": "month_swi_min"})
    .set_index("event_id")
)

print("CatNat events:", len(dfc))
print("Events with SWI:", len(swi_event))
print("Without SWI:", len(dfc) - len(swi_event))

# ----------------------------------------------------------
# 6. Final CatNat + SWI table
# ----------------------------------------------------------
final = dfc.set_index("event_id").join(swi_event, how="left").reset_index(drop=True)
print(final.shape, "| share without SWI:", final["swi_min"].isna().mean())

final.to_excel(os.path.join(dir_catnat, "CatNat_drought_events_SWI.xlsx"), index=False)
