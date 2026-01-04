#----------here we check if there was any convective activities, such as thunderstorm during a dust event case study using radar-------------------------------
#-----------install the Herbie python library 
from herbie import Herbie
from herbie.toolbox import EasyMap, pc
import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader
import matplotlib.pyplot as plt
import pandas as pd

# 1) Define your map projection
proj = ccrs.PlateCarree()

# 2) Load the reflectivity field for 2021-10-11 20:00 UTC for example
hb = Herbie("2024-11-11 23:00")
ds = hb.xarray("REFC:entire")
ds["refc"] = ds.refc.where(ds.refc > 0)  # mask non-returns

# 
emap = EasyMap(crs=ds.herbie.crs, theme="dark")
emap.STATES(linewidth=3,edgecolor="g").BORDERS().OCEAN()
ax = emap.ax

# 4) Zoom into the western U.S. coverage area of interest
ax.set_extent([-125, -105, 30, 50], crs=proj)

# 5) Title with the actual timestamp from the data
time0 = pd.to_datetime(ds.time.values)
title_time = time0.strftime("%Y-%m-%d %H:%M UTC")
ax.set_title(f"", fontsize=19)

# 6) Overlay the Central Valley shapefile in red

import geopandas as gpd
cv = gpd.read_file("central_valley_alluvial_boundary_usgs.shp").to_crs(epsg=4326)
central_feat = ShapelyFeature(
    cv.geometry,
    ccrs.PlateCarree(),
    edgecolor="red",
    facecolor="none",
    linewidth=1.2,
)
ax.add_feature(central_feat)



# 7) Plot the reflectivity
ds.refc.plot(
    x="longitude",
    y="latitude",
    ax=ax,
    transform=pc,
    cbar_kwargs={"shrink": 0.5, "orientation": "horizontal", "pad": 0.01},
)
plt.title("2024-11-11 16:00 LST")
plt.savefig("B3.png", dpi=500, bbox_inches="tight")
plt.show()
