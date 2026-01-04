# ------------ we use the from somoclu python for Self-organizing mapping of dust synoptic pattern.---------- Read the documentation of Somoclu for extra information

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import geopandas as gpd
from somoclu import Somoclu

# ---------------------- projection robinson ----------------------
lon0, lon1 = -130, -105
lat0, lat1 =   30,   50
n = 200
top    = np.column_stack((np.linspace(lon0, lon1, n), [lat1]*n))
right  = np.column_stack(([lon1]*n, np.linspace(lat1, lat0, n)))
bottom = np.column_stack((np.linspace(lon1, lon0, n), [lat0]*n))
left   = np.column_stack(([lon0]*n, np.linspace(lat0, lat1, n)))
pts    = np.vstack([top, right, bottom, left, top[:1]])
codes  = [Path.MOVETO] + [Path.LINETO]*(len(pts)-2) + [Path.CLOSEPOLY]
box_path = Path(pts, codes)

def clip_to_curve(ax, linewidth=4, edgecolor="black"):
    ax.set_boundary(box_path, transform=ccrs.PlateCarree())
    ax.patch.set_facecolor("none")
    ax.add_patch(PathPatch(box_path, transform=ccrs.PlateCarree(),
                           facecolor="none", edgecolor=edgecolor,
                           linewidth=linewidth, zorder=10))

# ------------------------------ CV shapefiles ---------------------------------
cv_gdf  = gpd.read_file("central_valley_alluvial_boundary_usgs.shp").to_crs(epsg=4326)
cv1_gdf = gpd.read_file("california_shapefile.shp").to_crs(epsg=4326)

# --------------------------- load geopotential z500 of composite dust events ----------------------------
ds = xr.open_dataset("son_z_dusty_m.nc")

# dateclass
time_name = "date" if "date" in ds.coords or "date" in ds.dims else \
            "valid_time" if "valid_time" in ds.coords or "valid_time" in ds.dims else None
if time_name is None:
    raise KeyError("No time coordinate named 'date' or 'time' found.")

# varname
z = ds["z"]
# select 500 hPa
if "pressure_level" in z.dims:
    z = z.sel(pressure_level=500)

# ----------------------------- subset septemner-november (SON) season -------------------------------
son = z.sel({time_name: ds[time_name].dt.month.isin([9,10,11])})
son = (
    son
    .where((son.latitude >= lat0) & (son.latitude <= lat1), drop=True)
    .where((son.longitude >= lon0) & (son.longitude <= lon1), drop=True)
)

# ---------------------------- flatten & train SOM ------------------------------
data = son.values               # (ntime, nlat, nlon)
ntime, nlat, nlon = data.shape
flat  = data.reshape(ntime, nlat*nlon)

som = Somoclu(n_columns=2, n_rows=2, maptype="planar", initialization="pca")
np.random.seed(777)
som.train(flat, epochs=1000)

# ------------------------- composites & frequencies ----------------------------
# Somoclu BMUs are (row, col). Convert to a single node index
labels = som.bmus[:,0] * 2 + som.bmus[:,1]

comps, freqs = {}, {}
for node in range(4):
    mask = (labels == node)
    freqs[node] = mask.sum() / ntime * 100.0
    comps[node] = np.nanmean(data[mask], axis=0) if mask.any() else np.full((nlat, nlon), np.nan)

# ------------------------------- contour levels --------------------------------
vmin = min(np.nanmin(c) for c in comps.values())
vmax = max(np.nanmax(c) for c in comps.values())
levels = np.arange((vmin//15)*15, (vmax//15+1)*15 + 1, 15)

# ---------------------------- grids & projection --------------------------------
lon2d, lat2d = np.meshgrid(ason.longitude, ason.latitude)
proj = ccrs.LambertConformal(central_longitude=-117, central_latitude=35,
                             standard_parallels=(33, 45))

# ------------------- ORDER NODES BY FREQUENCY (DESCENDING) ---------------------
sorted_nodes = sorted(freqs.keys(), key=lambda n: freqs[n], reverse=True)
# rank=1 is the most frequent, rank=4 is the least frequent.

# ------------------------------- plotting --------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(16, 6), subplot_kw={"projection": proj})
plt.tight_layout(rect=[0, 0, 1, 0.88])

for rank, node in enumerate(sorted_nodes, start=1):
    ax = axes[rank-1]
    clip_to_curve(ax)
    ax.set_extent([lon0, lon1, lat0, lat1], crs=ccrs.PlateCarree())

    cs = ax.contour(lon2d, lat2d, comps[node], levels=levels,
                    colors="red", linewidths=1.2, transform=ccrs.PlateCarree())
    ax.clabel(cs, fmt="%d", fontsize=8)

    ax.add_geometries(cv_gdf.geometry,  crs=ccrs.PlateCarree(),
                      edgecolor="green", facecolor="none", linewidth=2.4)
    ax.add_geometries(cv1_gdf.geometry, crs=ccrs.PlateCarree(),
                      edgecolor="black", facecolor="none", linewidth=2.2)

    ax.coastlines(resolution="50m", linewidth=1.2)
    ax.add_feature(cfeature.BORDERS, linewidth=1.2)
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=1, zorder=2)

    # label the Title: Node rank by frequency + its percentage
    ax.set_title(f"Z500mb-Type {rank}\n{freqs[node]:.1f}%", fontsize=20, pad=4)
    
    ax0=axes.flat[0]
    
    b = ax0.get_position()
    fig.text(0.01, b.y1,"A", ha="left",va="bottom", fontsize=34,fontweight="bold",
             fontfamily="sans-serif",transform=fig.transFigure)

fig.suptitle("", fontsize=20)
fig.savefig("SOm_pattern.png", dpi=500, bbox_inches="tight")
plt.show()
