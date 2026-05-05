# Synoptic-DustEvents_CentralValley

Code and workflows for identifying Central Valley dust events from surface observations and linking them to synoptic-scale circulation patterns.
Data is available at https://doi.org/10.6084/m9.figshare.32049018

This repository supports analyses such as:
- dust observation → station event → regional “dusty day” summaries  
- seasonal/diurnal statistics 
- synoptic pattern classification (e.g., SOM/clustering on Z500)
- hydroclimate context (precipitation onset and drought categories)

---

## What this repo does

**Dust event detection**
- Ingests quality-controlled surface weather observations (ASOS/AWOS/METAR-derived)
- Applies dust criteria (e.g., present weather codes and visibility/wind thresholds as configured)
- Aggregates to:
  1) **dust observations** (individual reports meeting criteria)  
  2) **station-level dust events** (event grouping at each station)  
  3) **regional dusty days** (any local day with ≥1 station reporting dust)

**Synoptic classification**
- Builds daily synoptic fields (e.g., Z500, surfcae wind speed, soil moisture) for dusty days
- Classifies patterns using SOM / clustering
- Produces composites and significance tests

**Hydroclimate linkages**
- Rainfall onset metrics using GridMET precipitation (spatial coverage onset)
- Drought categories (e.g., SPEI-based)

---




