# 2021-10-11 CV dust event: Wind, Visibility, RH (5-min from METAR T/Td), and PM10 (hourly)


import os
os.chdir('/Users/precious/Downloads')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from matplotlib.dates import DateFormatter
from dateutil import tz
import requests, io, re

# ------------------- style -------------------
plt.rcParams.update({
    'font.size':        18,
    'axes.titleweight': 'bold',
    'axes.labelsize':   18,
    'xtick.labelsize':  16,
    'ytick.labelsize':  16
})

# -------------- METAR → T/Td → RH ------------- we check the METAR text [read ASOS manual for METAR CODES]
_re_pairC   = re.compile(r'\b(M?\d{2})/(M?\d{2})\b')
_re_Tgroup  = re.compile(r'\bT(\d{8})\b')

def _decode_pairC(s):
    return -float(s[1:]) if s.startswith('M') else float(s)

def temp_dew_from_metar(metar):
    if not isinstance(metar, str):
        return np.nan, np.nan
    m = _re_Tgroup.search(metar)
    if m:
        d = m.group(1)
        sign_t  = -1 if d[0] == '1' else 1
        sign_td = -1 if d[4] == '1' else 1
        T  = sign_t  * (int(d[1:4]) / 10.0)
        Td = sign_td * (int(d[5:8]) / 10.0)
        return T, Td
    m = _re_pairC.search(metar)
    if m:
        t_str, td_str = m.groups()
        return _decode_pairC(t_str), _decode_pairC(td_str)
    return np.nan, np.nan

def rh_from_t_td_c(Tc, Tdc):
    T  = np.asarray(Tc,  dtype='float64')
    Td = np.asarray(Tdc, dtype='float64')
    es  = np.exp((17.625 * T )/(243.04 + T ))
    esd = np.exp((17.625 * Td)/(243.04 + Td))
    rh  = 100.0 * (esd / es)
    return np.clip(rh, 0, 100)

# -------------- ASOS download ----------------WE got the asos extraction download request from the official mesonent website
def download_asos_data(station_id, year):
    base = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py/"
    params = dict(
        station=station_id,
        data='sped,drct,vsby,sknt,relh,metar,wxcodes',
        year1=year, month1='10', day1='11',
        year2=year, month2='10', day2='11',
        tz='UTC', format='comma', direct='yes'
    )
    r = requests.get(base, params=params)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), skiprows=5)

# ---------------- stations (panels only) ----------------
stations = ["SCK","MCE","MAE","FAT","VIS","PTV"]  

# -------------- precessing ------------------
dfs = []
for st in stations:
    df0 = download_asos_data(st, 2021)
    df0['station'] = st
    dfs.append(df0)
asos = pd.concat(dfs, ignore_index=True)

asos.rename(columns={'drct':'wind_direction','sknt':'wind_speed_knots','vsby':'visibility_miles'},
            inplace=True)

# Time handling (UTC → local)
asos['valid']       = pd.to_datetime(asos.valid, errors='coerce', utc=True)
asos['valid_local'] = asos['valid'].dt.tz_convert('America/Los_Angeles')

# Unit conversions
asos['wind_speed_mps'] = pd.to_numeric(asos.wind_speed_knots, errors='coerce') * 0.514444
asos['visibility_km']  = pd.to_numeric(asos.visibility_miles, errors='coerce') * 1.60934

# T/Td → RH
TT = asos['metar'].apply(lambda s: pd.Series(temp_dew_from_metar(s), index=['temp_c','dewpoint_c']))
asos[['temp_c','dewpoint_c']] = TT
asos['RH_calc'] = rh_from_t_td_c(asos['temp_c'], asos['dewpoint_c'])

#  plot using local times
mask_d = asos.valid.dt.date == datetime(2021,10,11).date()
mask_t = asos.valid.dt.hour.between(14,23)  # 07–16 LST
df = asos[mask_d & mask_t].dropna(subset=['wind_speed_mps','visibility_km'])

# -------------- PM10 hourly -------------------
# PM10 was extracted from EPA officialy website
pm10_raw = pd.DataFrame({
    "time":[14,15,16,17,18,19,20,21,22,23],
    "SCK":[127.1, np.nan, np.nan, np.nan, 900.6, 621.8, 507.9, 471.1, 503.8, 433.5],
    "MAE":[264.2, 567.1, 772.4, 581.0, 898.5, 802.4, 847.3, 982.4, 887.2, np.nan],
    "FAT":[123.2, 264.4, 329.8, 345.8, 474.2, 535.4, 529.4, 494.5, 541.6, 246.6],
    "VIS":[61.7, 133.8, 431.7, 551.7, 693.9, 601.5, 985.0, 542.4, 985.0, 928.9],
    "BFL":[156.5, 264.5, 468.2, 985.0, 985.0, 777.2, 806.0, 985.0, 985.0, 985.0],
})
pm10 = pm10_raw.set_index('time')
pm10['PTV'] = pm10['BFL']  # copy BFL → PTV
pm10 = pm10.drop(columns=['BFL'])

# Convert UTC hours to local time index for plotting
pm10.index = (pd.to_datetime([f"2021-10-11 {h:02d}:00" for h in pm10.index], utc=True)
              .tz_convert('America/Los_Angeles'))
PM_MAX = float(pm10.max().max()) * 1.10

# ---- 
PM_LABEL_STATIONS = {'SCK', 'FAT'}
SUPPRESS_LABEL_STATIONS = set(stations) - PM_LABEL_STATIONS

# -------------- thresholds & ranges -----------
wind_thresh = 6
vis_thresh  = 10
rh_thresh   = 70
wind_max    = df.wind_speed_mps.max() * 1.1
vis_max     = df.visibility_km.max()  * 1.1

# -------------- figure: 2 rows × 3 cols -------
fig, axes = plt.subplots(2, 3, figsize=(21, 12), sharex=True)
axes = axes.flatten()
fig.subplots_adjust(left=0.10, right=0.88, top=0.92, bottom=0.14)
fig.suptitle("(a) October 11 2021 Dust Event", y=0.90, x=0.10)

# Local tz 
pac_tz = tz.gettz('America/Los_Angeles')
hour_fmt = DateFormatter('%H:%M', tz=pac_tz)
idx5 = pd.date_range('2021-10-11 07:00', '2021-10-11 16:59',
                     freq='5min', tz='America/Los_Angeles')

for i, st in enumerate(stations):
    ax = axes[i]
    sd = df[df.station==st].copy()
    sd = sd.set_index('valid_local').sort_index()
    sd = sd[~sd.index.duplicated(keep='last')]

    ax.set_title(st)
    ax.grid(True, linestyle='--', alpha=0.5)

    # wind (left, blue)
    ax.plot(sd.index, sd['wind_speed_mps'], color='tab:blue')
    ax.axhline(wind_thresh, color='tab:blue', linestyle='--')
    ax.set_ylim(0, wind_max)
    if i % 3 == 0:   # left column
        ax.set_ylabel('Wind (m/s)', color='tab:blue')
    else:
        ax.yaxis.label.set_visible(False)
    ax.tick_params(axis='y', labelcolor='tab:blue')

    # PM10 (outer left, purple)
    if st in pm10.columns:
        ax_pm = ax.twinx()
        ax_pm.spines['left'].set_position(('outward', 56))
        ax_pm.set_ylim(0, PM_MAX)
        s_pm = pm10[st].dropna()
        ax_pm.plot(s_pm.index, s_pm.values, color='tab:purple',
                   linewidth=1.6, marker='o', drawstyle='steps-post')
        ax_pm.spines['right'].set_visible(False)
        ax_pm.yaxis.set_label_position('left')
        ax_pm.yaxis.set_ticks_position('left')

        if st in SUPPRESS_LABEL_STATIONS:
            ax_pm.spines['left'].set_visible(False)
            ax_pm.tick_params(left=False, labelleft=False, right=False, labelright=False)
            ax_pm.set_ylabel('')
        else:
            ax_pm.spines['left'].set_visible(True)
            ax_pm.tick_params(left=True, labelleft=True, right=False, labelright=False)
            ax_pm.set_ylabel(r'PM$_{10}$ ($\mu$g m$^{-3}$)', color='tab:purple', labelpad=12)

    # visibility (inner right, red)
    ax2 = ax.twinx()
    ax2.plot(sd.index, sd['visibility_km'], color='tab:red')
    ax2.axhline(vis_thresh, color='tab:red', linestyle='--')
    ax2.set_ylim(0, vis_max)
    if i % 3 == 2:   # right column
        ax2.set_ylabel('Visibility (km)', color='tab:red')
    else:
        ax2.tick_params(labelright=False, right=False)

    # RH (outer right, green)
    if sd['temp_c'].notna().sum() >= 2 and sd['dewpoint_c'].notna().sum() >= 2:
        T5  = sd['temp_c']    .reindex(sd.index.union(idx5)).sort_index().interpolate('time').reindex(idx5)
        Td5 = sd['dewpoint_c'].reindex(sd.index.union(idx5)).sort_index().interpolate('time').reindex(idx5)
        RH5 = rh_from_t_td_c(T5, Td5)

        ax3 = ax.twinx()
        ax3.spines['right'].set_position(('outward', 56))
        ax3.plot(idx5, RH5, color='tab:green', linewidth=1.6)
        ax3.axhline(rh_thresh, color='tab:green', linestyle='--', alpha=0.85)
        ax3.set_ylim(0, 100)
        if i % 3 == 2:
            ax3.set_ylabel('RH (%)', color='tab:green', labelpad=10)
        else:
            ax3.tick_params(labelright=False, right=False)

    # x-axis ticks — LOCAL time
    ax.xaxis.set_major_locator(mdates.HourLocator(tz=pac_tz))
    ax.xaxis.set_major_formatter(hour_fmt)
    ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=[0,30], tz=pac_tz))
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    # shading in LOCAL time: vis<16 window with ≥1 (wind>=6 & vis<=10)
    low16 = sd['visibility_km'] < 16
    runs  = (low16 != low16.shift()).cumsum()
    for _, seg in sd.groupby(runs):
        if not low16.loc[seg.index[0]]:
            continue
        t0, t1 = seg.index[0], seg.index[-1]
        has_dust = ((seg['visibility_km'] <= 10) & (seg['wind_speed_mps'] >= 6)).any()
        if has_dust:
            ax.axvspan(t0 - pd.Timedelta('5min'), t1 + pd.Timedelta('5min'),
                       facecolor='tab:red', alpha=0.18)

# Legend
lines = [
    plt.Line2D([0],[0], color='tab:blue'),
    plt.Line2D([0],[0], color='tab:red'),
    plt.Line2D([0],[0], color='tab:green'),
    plt.Line2D([0],[0], color='tab:purple'),
]
fig.legend(lines, ['Wind Speed (m/s)', 'Visibility (km)', 'RH (%)',
                   r'PM$_{10}$ ($\mu$g m$^{-3}$)'],
           loc='lower center', ncol=4, frameon=False, fontsize=16,
           bbox_to_anchor=(0.5, 0.06))

plt.tight_layout(rect=[0,0.12,1,0.90])
plt.savefig('WID-OCT.png',
            dpi=600, bbox_inches='tight')
plt.show()
