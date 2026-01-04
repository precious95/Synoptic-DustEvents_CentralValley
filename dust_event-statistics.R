library(ggplot2)
library(dplyr)
library(maps)
library(cowplot)
library(gridExtra)
library(grid) 
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(ggrepel)
library(ggspatial)

setwd('/Users/precious/Downloads/ASOS-PAPER')
# 1) LOAD central valley boundary shapefile
cv_boundary <- st_read("central_valley_alluvial_boundary_usgs.shp", quiet = TRUE) %>%
  st_zm(drop=TRUE) %>%
  st_transform(4326)

tables_total=read.csv('dust_events_summary.csv',header=T)

year_totals5=tables_total$yearly
colnames(year_totals5)=c('year','count')

month_totals5=tables_total$monthly
colnames(month_totals5)=c('month','count')

dur_counts5=tables_total$duration
colnames(dur_counts5)=c('duration_group','count')

diurnal_totals5=tables_total$diurnal
colnames(diurnal_totals5)=c('hour','count')

#-------------------------------------------------------------------
library(scales)

# make sure the factor order is what you want
dur_counts5 <- dur_counts5 %>%
  mutate(duration_group = factor(duration_group,
                                 levels = c("≤1h","2h","3h","4h","5h","6h","7h","8h","9h","≥10h"),
                                 ordered = TRUE)) %>%
  arrange(duration_group) %>%
  mutate(
    pct     = count / sum(count),
    pct_lab = percent(pct, accuracy = 0.1) 
  )

#-----------------

library(dplyr)
library(scales)

dur_counts5 <- dur_counts5 %>%
  mutate(
    duration_group = factor(
      as.character(duration_group),
      levels  = c(as.character(1:9), ">10"),
      ordered = TRUE
    )
  ) %>%
  arrange(duration_group) %>%
  mutate(
    pct     = count / sum(count),
    pct_lab = percent(pct, accuracy = 0.1)
  )

#---------

#––– Panel C: Duration + inset map –––
ymax_c <- max(dur_counts5$count) * 1.15  
#-----------------------------------------------------------------------


#––– 5) Prepare inset map –––
cv_counties <- map_data("county", region = "california") %>%
  filter(subregion %in% c("shasta","tehama","glenn","butte","colusa",
                          "sutter","yuba","yolo","sacramento","san joaquin",
                          "stanislaus","merced","madera","fresno","kings",
                          "tulare","kern"))
ca_state <- map_data("state", region = "california")

map_insert1=ggplot() +
  # a) Central Valley alluvial boundary base (your sf layer)
  geom_sf(data = cv_boundary,
          fill  = "#f0f0f0",
          color = "black",
          linewidth = 0.6) +
  # b) county subdivisions (your fortified polygon layer)
  geom_polygon(data = cv_counties,
               aes(x = long, y = lat, group = group),
               fill  = NA,
               color = "red",
               linewidth = 0.1)+
  theme_void() 

map_inset <- ggplot() +
  geom_polygon(data = ca_state, aes(long, lat, group = group),
               fill="grey90", color="black") +
  geom_polygon(data = cv_counties, aes(long, lat, group = group),
               fill="lightblue", color="blue") +
  coord_quickmap() +
  theme_void()


clean_theme <- theme_classic(base_size = 24) +
  theme(
    panel.border   = element_blank(),
    axis.line      = element_blank(),
    axis.ticks     = element_blank()
  )

#––– Panel A: Yearly with trend –––
pA <- ggplot(year_totals5, aes(year, count)) +
  geom_col(fill="black") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed", color="black") +
  ylim(0,100)+
  labs(title="(a) Yearly Dust Event (2005–2024)",
       x="Year", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

#––– Panel B: Monthly –––
pB <- ggplot(month_totals5, aes(month, count)) +
  geom_col(fill="black") +ylim(0,150)+
  labs(title="(b) Monthly Dust Event",
       x="Month", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

#––– Panel C: Duration + inset –––
pC_base <- ggplot(dur_counts5, aes(duration_group, count)) +
  geom_col(fill="black") +
  labs(title="(C) Event Count by Duration",
       x="Duration (Hours)", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

pC <- ggdraw(pC_base) +
  draw_plot(map_insert1, x=0.55, y=0.5, width=0.4, height=0.4)

#--------------------------------------------------------------------
pC_base <- ggplot(dur_counts5, aes(duration_group, count)) +
  geom_col(fill="black") +
  geom_text(aes(label = pct_lab), vjust = -0.35, size = 5) +
  coord_cartesian(ylim = c(0, ymax_c)) +
  labs(title="(d) Event Count by Duration",
       x="Duration (Hours)", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

pC <- ggdraw(pC_base) +
  draw_plot(map_insert1, x=0.55, y=0.5, width=0.4, height=0.4)
#---------------------------------------------------------------------

#––– Panel D: Diurnal –––
pD <- ggplot(diurnal_totals5, aes(hour, count)) +
  geom_col(fill="black", width=1) +
  scale_x_continuous(breaks=seq(0,24,2),
                     labels=sprintf("%02d", seq(0,24,2))) +
  coord_cartesian(xlim=c(0,23)) +ylim(0,100)+
  labs(title="(c) Diurnal Start-Time Distribution",
       x="Local Time (Hours)", y=NULL) +
  clean_theme +
  theme(plot.title = element_text(hjust = 0))

#––– 6) Arrange and add single Y label –––
png("plots9.png", width=1850, height=950)
grid.arrange(
  pA, pB,
  pC, pD,
  ncol = 2, nrow = 2,
  left = textGrob("Total Dust Events", rot = 90,
                  gp = gpar(fontsize = 26, fontface="bold"))
)

grid.arrange(
  pA, pB,
  pD, pC,
  ncol = 2, nrow = 2,
  left = textGrob("Total Dust Events", rot = 90,
                  gp = gpar(fontsize = 26, fontface="bold"))
)

grid.arrange(
  pA, pB,
  pD, pC_base,
  ncol = 2, nrow = 2,
  left = textGrob("Total Dust Events", rot = 90,
                  gp = gpar(fontsize = 26, fontface="bold"))
)
dev.off()


#---------------slope trend
# 


fit <- lm(count ~ year, data = year_totals5)
s   <- summary(fit)

slope    <- coef(fit)[["year"]]
se_slope <- s$coefficients["year","Std. Error"]
ci       <- confint(fit)["year", ]
r2       <- s$r.squared
mean_ev  <- mean(df$events)
sd_ev    <- sd(df$events)

cat(sprintf(
  "On average, this increasing trend is about %.2f (±%.2f) dust events per year.\n",
  slope, se_slope))
cat(sprintf("95%% CI: [%.2f, %.2f]; R² = %.2f; n = %d\n",
            ci[1], ci[2], r2, nrow(df)))
cat(sprintf("Overall mean ± SD: %.1f ± %.1f events/yr\n", mean_ev, sd_ev))


#-----------------plot duration fig s2--------------------------------------------------------
# ------------------------------------------------------
# libs
library(dplyr)
library(stringr)
library(lubridate)
library(tidyr)
library(forcats)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(sf)

# ------------------------------------------------------
# CONFIG
MAX_GAP_HOURS <- 2
DUR_LABELS <- c(as.character(1:9), ">10")
DUR_BREAKS <- c(-Inf, 1, 2, 3, 4, 5, 6, 7, 8, 9, Inf)  # ≤1, 2..9, ≥10

# ------------------------------------------------------
# 1) Build per-report PWC flags from dust_events_cv5 
base_cv <- dust_events_cv5 %>%
  mutate(
    txt_wx = paste(coalesce(as.character(wxcodes), ""),
                   coalesce(as.character(metar), "")),
    has_DU     = str_detect(txt_wx, regex("\\bDU\\b",     TRUE)),
    has_BLDU   = str_detect(txt_wx, regex("\\bBLDU\\b",   TRUE)),
    has_VCBLDU = str_detect(txt_wx, regex("\\bVCBLDU\\b", TRUE)),
    has_DS     = str_detect(txt_wx, regex("\\bDS\\b",     TRUE)),
    has_SS     = str_detect(txt_wx, regex("\\bSS\\b",     TRUE)),
    has_HZ     = str_detect(txt_wx, regex("\\bHZ\\b",     TRUE))
  ) %>%
  mutate(
    has_BLDU_plain = has_BLDU & !has_VCBLDU,
    n_pwc          = rowSums(cbind(has_DU, has_BLDU_plain, has_DS, has_SS, has_VCBLDU, has_HZ), na.rm = TRUE),
    ANY_PWC_report = n_pwc >= 1
  ) %>%
  filter(ANY_PWC_report) %>%
  arrange(station, local_time)

# ------------------------------------------------------
# 2) Segment events per station, compute durations & bins
dust_event_duration <- base_cv %>%
  group_by(station) %>%
  mutate(
    dt_hours  = as.numeric(difftime(local_time, lag(local_time), units = "hours")),
    new_event = is.na(dt_hours) | dt_hours > MAX_GAP_HOURS,
    event_id  = cumsum(new_event)
  ) %>%
  ungroup() %>%
  group_by(station, event_id) %>%
  summarise(
    start_time = min(local_time),
    end_time   = max(local_time),
    .groups    = "drop"
  ) %>%
  mutate(
    start_date   = as.Date(start_time),
    start_year   = year(start_time),
    start_month  = factor(month.abb[month(start_time)], levels = month.abb, ordered = TRUE),
    start_hour   = hour(start_time),
    start_period = if_else(start_hour < 12, "Morning", "Evening"),
    duration_hr  = as.numeric(difftime(end_time, start_time, units = "hours")),
    duration_group = cut(duration_hr, breaks = DUR_BREAKS, labels = DUR_LABELS,
                         right = TRUE, ordered_result = TRUE)
  ) %>%
  arrange(station, start_time)

dust_event_duration=events_dedup
# ------------------------------------------------------
# --- Station coordinates (order north → south) ---
coords_txt <- "
station,lat,lon
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
"
df_coords <- read.csv(text = coords_txt, stringsAsFactors = FALSE) %>%
  arrange(desc(lat))                     # north → south
station_ids <- intersect(df_coords$station, unique(dust_event_duration$station))
df_coords   <- df_coords %>% filter(station %in% station_ids)
station_ids <- df_coords$station

# --- Central Valley boundary (for the map inset) ---
library(sf); library(ggplot2); library(cowplot); library(gridExtra); library(forcats); library(tidyr); library(dplyr)

cv_boundary <- st_read("central_valley_alluvial_boundary_usgs.shp", quiet = TRUE) |>
  st_zm(drop = TRUE) |>
  st_transform(4326)

bb   <- sf::st_bbox(cv_boundary)
pad  <- 0.25
xlim <- c(bb["xmin"] - pad, bb["xmax"] + pad)
ylim <- c(bb["ymin"] - pad, bb["ymax"] + pad)

# one inset map per station
inset_maps <- setNames(vector("list", length(station_ids)), station_ids)
for (i in seq_len(nrow(df_coords))) {
  sc <- df_coords[i, ]
  inset_maps[[sc$station]] <-
    ggplot() +
    geom_sf(data = cv_boundary, fill = "#f0f0f0", color = "black", linewidth = 0.6) +
    geom_point(data = sc, aes(lon, lat), color = "red", size = 3, inherit.aes = FALSE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void()
}

# --- Build two-bin (≤1 hr, >1 hr) % bars per station ---
# duration_group 
plots_2bin <- list()

for (stn in station_ids) {
  stn_summary <- dust_event_duration %>%
    filter(station == stn) %>%
    mutate(
      bin = if_else(as.character(duration_bin) == "1", "≤ 1 hr", "> 1 hr"),
      bin = factor(bin, levels = c("≤ 1 hr", "> 1 hr"), ordered = TRUE)
    ) %>%
    count(bin, name = "count") %>%
    complete(bin, fill = list(count = 0)) %>%
    mutate(
      total  = sum(count, na.rm = TRUE),
      pct    = ifelse(total > 0, 100 * count / total, 0),
      pct_lab = sprintf("%.1f%%", pct)
    )
  
  p_main <- ggplot(stn_summary, aes(x = bin, y = pct)) +
    geom_col(fill = "black", width = 0.8) +
    geom_text(aes(label = pct_lab), vjust = -0.35, size = 16) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = stn, x = NULL, y = "Count(%)") +
    theme_minimal(base_size = 26) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 26, face = "bold"),
      axis.text.x = element_text(size = 26),
      axis.text.y = element_text(size = 26)
    )
  
  inset <- inset_maps[[stn]]
  if (is.null(inset)) inset <- ggplot() + theme_void()
  
  plots_2bin[[stn]] <- cowplot::ggdraw() +
    cowplot::draw_plot(p_main) +
    cowplot::draw_plot(inset, x = 0.68, y = 0.60, width = 0.4, height = 0.4)
}

# --- Arrange grid (north → south) and save ---
n <- length(plots_2bin)
ncol <- 4
nrow <- ceiling(n / ncol)

png("duration_plots_all_stations.png", width = 3000, height = 2000)
gridExtra::grid.arrange(grobs = plots_2bin, nrow = nrow, ncol = ncol)
dev.off()

