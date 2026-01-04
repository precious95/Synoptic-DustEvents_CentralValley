import os

os.chdir('/Users/precious/Downloads')

import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader

from matplotlib.path import Path
from matplotlib.patches import PathPatch

import geopandas as gpd
from shapely.vectorized import contains

from goes2go import GOES

# -----------------------------------------------------------------------------
# 0) Load Central Valley alluvial boundary
# -----------------------------------------------------------------------------
cv = gpd.read_file("central_valley_alluvial_boundary_usgs.shp").to_crs(epsg=4326)

# -----------------------------------------------------------------------------
# 1) Load California boundary from Natural Earth
# -----------------------------------------------------------------------------
shpfile = shapereader.natural_earth(
    resolution="50m", category="cultural",
    name="admin_1_states_provinces"
)
reader = shapereader.Reader(shpfile)
CA_geom = None

for rec in reader.records():

    if rec.attributes["name"] == "California":
        CA_geom = rec.geometry
        break
if CA_geom is None:
    raise RuntimeError("Could not find California geometry")

# -----------------------------------------------------------------------------
# 2) projection (robinson)
# -----------------------------------------------------------------------------
lon0, lon1 = -125, -112
lat0, lat1 =   31,   42.1
n = 200
f          = np.linspace(0, 1, n)
amp_top    = 0.28
amp_bot    = 0.32

top    = np.column_stack((np.linspace(lon0, lon1, n), lat1 - amp_top * np.sin(np.pi * f)))
right  = np.column_stack(([lon1]*n,     np.linspace(lat1, lat0, n)))
bottom = np.column_stack((np.linspace(lon1, lon0, n), lat0 - amp_bot * np.sin(np.pi * f)))
left   = np.column_stack(([lon0]*n,     np.linspace(lat0, lat1, n)))
pts    = np.vstack([top, right, bottom, left, top[:1]])
codes  = [Path.MOVETO] + [Path.LINETO] * (len(pts)-2) + [Path.CLOSEPOLY]
frame_path = Path(pts, codes)

def clip_to_curve(ax, linewidth=4, edgecolor="black", zorder=10):
    ax.set_boundary(frame_path, transform=ccrs.PlateCarree())
    ax.patch.set_facecolor("none")
    ax.add_patch(PathPatch(
        frame_path, transform=ccrs.PlateCarree(),
        facecolor="none", edgecolor=edgecolor,
        linewidth=linewidth, zorder=zorder
    ))

# -----------------------------------------------------------------------------
# 3) #Change to your specific date u want----here we check the october 11 2021 dust event in central valley
# -----------------------------------------------------------------------------
times = ["2021-10-11 21:00:00", "2021-10-11 22:00:00", "2021-10-11 23:00:00"]
proj  = ccrs.LambertConformal(
    central_longitude=-117, central_latitude=35,
    standard_parallels=(33, 45)
)
fig, axs = plt.subplots(
    1, 3, figsize=(18, 7),
    subplot_kw={"projection": proj}
)
plt.subplots_adjust(wspace=0.03, left=0.02, right=0.98, top=0.90, bottom=0.08)

# -----------------------------------------------------------------------------
# 4) mask Dust RGB over california
# -----------------------------------------------------------------------------
for ax, ts in zip(axs, times):
    # a) Set extent & clip frame
    ax.set_extent([lon0, lon1, lat0, lat1], ccrs.PlateCarree())
    clip_to_curve(ax)

    # b) Fetch GOES Dust RGB
    sat  = GOES(satellite=18).nearesttime(ts)

    ds   = sat

    dust = ds.rgb.Dust().copy()

    # c) Invert GOES grid → lon2d, lat2d
    h       = ds["goes_imager_projection"].perspective_point_height
    x       = ds["x"].values * h
    y       = ds["y"].values * h
    X2d, Y2d = np.meshgrid(x, y)
    geo_pts  = ccrs.PlateCarree().transform_points(ds.rgb.crs, X2d, Y2d)
    lon2d    = geo_pts[...,0]
    lat2d    = geo_pts[...,1]

    # d) Mask outside California
    mask_ca = contains(CA_geom, lon2d, lat2d)
    for k in range(3):
        dust[..., k] = np.where(mask_ca, dust[..., k], np.nan)

    # e) Plot masked Dust RGB
    im_kwargs = ds.rgb.imshow_kwargs.copy()
    im_kwargs.pop("interpolation", None)
    im_kwargs.pop("resample", None)
    ax.imshow(
        dust,
        interpolation="bicubic",
        resample=True,
        **im_kwargs
    )

    # f) Add geographic features & outlines
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white", zorder=1)
    ax.add_feature(cfeature.LAKES.with_scale("50m"), facecolor="lightgrey", zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=1, zorder=2)
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=1, zorder=2)
    ax.add_geometries([CA_geom], crs=ccrs.PlateCarree(),
                      facecolor="none", edgecolor="black", linewidth=2, zorder=3)
    ax.add_geometries(cv.geometry, crs=ccrs.PlateCarree(),
                      facecolor="none", edgecolor="green", linewidth=1.5, zorder=4)

    # g) Title
    ax.set_title(f"GOES Dust RGB\n{ts[:-3]} UTC", fontsize=22)

# -----------------------------------------------------------------------------
# 5) ouput of GOES dust RGB (pink mangenta indicate dust)
# -----------------------------------------------------------------------------
plt.savefig("9geos_2021_three_timestamps.png", dpi=500, bbox_inches="tight")
plt.show()

