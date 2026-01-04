# figure 2b [Study area Map]

import numpy as np
import pandas as pd
import pygmt
import geopandas as gpd
from io import StringIO

# ── styling. We use PYGMT, install the library ───────────
pygmt.config(
    MAP_FRAME_TYPE="plain",
    MAP_TICK_LENGTH_PRIMARY="0.1c",
    FONT_LABEL="10p,Helvetica",
    FONT_ANNOT_PRIMARY="10p,Helvetica",
)

# ── Meteorlogical Stations (15) ASOS/AWOS ─────────────────────────────────────────
data = """station,lat,lon
BFL,35.43440,-119.0542
BAB,39.13609,-121.4366
FAT,36.78000,-119.7194
HJO,36.31139,-119.6232
MAE,36.98486,-120.1107
MYV,39.10203,-121.5688
MCE,37.28603,-120.5179
OVE,39.49000,-121.6200
PTV,36.02732,-119.0629
RBL,40.15190,-122.2536
RDD,40.50900,-122.2934
SAC,38.50690,-121.4950
SMF,38.69542,-121.5908
SCK,37.89417,-121.2383
VIS,36.31867,-119.3929
"""
stations = (
    pd.read_csv(StringIO(data))
      .sort_values("lat", ascending=False)  # legend order: N→S
      .reset_index(drop=True)
)

# ── Region Projection ───────────────────
region    = [-125, -112, 31, 42]
proj_main = "L-119.5/37.0/33/42/15c"

n = 200
lon0, lon1, lat0, lat1 = region
f = np.linspace(0, 1, n)
top    = np.column_stack((np.linspace(lon0, lon1, n), lat1 - 0.10*np.sin(np.pi*f)))
right  = np.column_stack(([lon1]*n,        np.linspace(lat1, lat0, n)))
bottom = np.column_stack((np.linspace(lon1, lon0, n), lat0 - 0.30*np.sin(np.pi*f)))
left   = np.column_stack(([lon0]*n,        np.linspace(lat0, lat1, n)))
pinch  = np.vstack([top, right, bottom, left, top[:1]])

# ──  turn Polygon or MultiPolygon into list of parts ─────────
def parts_list(gdf):
    geom = gdf.geometry.union_all().simplify(0.01)
    return list(getattr(geom, "geoms", [geom]))

# ── Central Valley & California outlines ──────────────────
cv = gpd.read_file("central_valley_alluvial_boundary_usgs.shp").to_crs("EPSG:4326")
cv_parts = parts_list(cv)

ca = gpd.read_file("california_shapefile.shp").to_crs("EPSG:4326")
ca_parts = parts_list(ca)

# ── Figure & grayscale topo (water white) ─────────────────
fig = pygmt.Figure()
fig.grdimage(
    grid="@earth_relief_01m", region=region, projection=proj_main,
    cmap="gray", shading=True, nan_transparent=True
)
fig.coast(
    region=region, projection=proj_main, water="white",
    shorelines="0.7p,black", borders=["1/0.6p,gray40","2/0.6p,gray40"],
    frame=["WSne+t"]
)

# Pinched frame + outlines
fig.plot(x=pinch[:,0], y=pinch[:,1], pen="1.4p,black")
for g in ca_parts:
    x, y = g.exterior.xy
    fig.plot(x=x, y=y, pen="0.9p,black")
for g in cv_parts:
    x, y = g.exterior.xy
    fig.plot(x=x, y=y, pen="1.2p,red")

# ── Colored symbols (STATIONS)──────────────────
color = {
    "RDD":"#1f77b4","RBL":"#ff7f0e","OVE":"#2ca02c","BAB":"#d62728",
    "MYV":"#9467bd","SMF":"#8c564b","SAC":"#e377c2","SCK":"#7f7f7f",
    "MCE":"#bcbd22","MAE":"#17becf","FAT":"#1b9e77",
    "VIS":"#7570b3","HJO":"#00b3b3","PTV":"#66a61e","BFL":"#ffb000"
}
shape = {
    "RDD":"c","RBL":"s","OVE":"d","BAB":"t",
    "MYV":"i","SMF":"a","SAC":"n","SCK":"h",
    "MCE":"p","MAE":"x","FAT":"+",
    "VIS":"g","HJO":"t","PTV":"r","BFL":"s"    
}
SIZE_NORMAL = "0.38c"
SIZE_BIG    = "0.48c"     # larger for BFL & HJO
EDGE        = "0.55p,white"

for _, r in stations.iterrows():
    sym = shape[r.station]
    sz  = SIZE_BIG if r.station in {"BFL","HJO"} else SIZE_NORMAL
    fig.plot(
        x=[r.lon], y=[r.lat],
        style=f"{sym}{sz}",
        fill=color[r.station],
        pen=EDGE,
        label=r.station
    )



# create Legend
fig.legend(position="JMR+jMR+o0.9c/0c", box="+gwhite+p0.9p", S="0.42c")

# ── USA inset (CA & CV highlighted) ───────────────────────
with fig.shift_origin(xshift="1.0c", yshift="1.0c"):
    fig.coast(
        region=[-125, -66, 24.39, 50], projection="L-100/40/33/45/5c",
        land="lightgray", water="white",
        shorelines="0.5p,gray40", borders=["1/0.4p,black","2/0.4p,black"],
        frame="n"
    )
    for g in ca_parts:
        x, y = g.exterior.xy
        fig.plot(x=x, y=y, pen="0.8p,black", fill="white")
    for g in cv_parts:
        x, y = g.exterior.xy
        fig.plot(x=x, y=y, pen="0.6p,red", fill="#9be3a8")

# ── Save ──────────────────────────────────────────────────
fig.savefig("fig2a.png", dpi=300)
print("Saved: fig2a-Study area map.png")
