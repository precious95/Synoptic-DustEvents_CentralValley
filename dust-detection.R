#--------------THE DUST DETECTION CODE points us to where dust is []
#------------ note: we had to check manually to confirm those dates by cross-referencing with NWS report, satellite imagery, screening for MIST before any dust report
#----------- use this to find dust observations, but confirm the dates with extra sources, 

library(dplyr)
library(ggplot2)
library(lubridate)
library(riem)
library(gridExtra)
library(cowplot)
library(tidyr)

setwd('/Users/precious/downloads/ASOS-PAPER/stations/AWS')
asos_utc=read.csv('asos.csv')
# simplest: alphabetical by station
asos_utc_sorted <- asos_utc %>% arrange(station)

unique(asos_utc_sorted$station)

#------------------------------------calculate RH and add to the data------------------------------
library(dplyr)
library(lubridate)
library(stringr)

# --- helpers ---
safe_num <- function(x) suppressWarnings(as.numeric(x))

decode_M2 <- function(v){
  out <- rep(NA_real_, length(v))
  neg <- str_detect(v, "^M\\d{1,2}$")
  pos <- str_detect(v, "^\\d{1,2}$")
  out[neg] <- -safe_num(str_sub(v[neg], 2))
  out[pos] <-  safe_num(v[pos])
  out
}

# Parse T/Td: prefer "dd/dd" or "Mdd/Mdd"; fallback to RMK TsnTTTsdTdTdTd (tenths °C)
extract_T_Td <- function(metar_chr){
  n <- length(metar_chr)
  T_met  <- rep(NA_real_, n)
  Td_met <- rep(NA_real_, n)
  
  # 1) dd/dd with word boundaries (won’t match "1/2SM")
  m1 <- str_match(metar_chr, "(?<=\\s|^)(M?\\d{1,2})/(M?\\d{1,2})(?=\\s|$)")
  ok1 <- !is.na(m1[,1])
  if (any(ok1)) {
    T_met[ok1]  <- decode_M2(m1[ok1,2])
    Td_met[ok1] <- decode_M2(m1[ok1,3])
  }
  
  # 2) RMK T00000000 (tenths °C)
  m2 <- str_match(metar_chr, "\\bT(\\d{8})\\b")
  idx2 <- which(!is.na(m2[,2]))
  if (length(idx2)){
    v   <- m2[idx2,2]
    sT  <- as.integer(substr(v,1,1)); TTT <- as.integer(substr(v,2,4))
    sD  <- as.integer(substr(v,5,5)); DDD <- as.integer(substr(v,6,8))
    Tval <- (TTT/10) * ifelse(sT==1, -1, 1)
    Dval <- (DDD/10) * ifelse(sD==1, -1, 1)
    
    fill_idx <- idx2[is.na(T_met[idx2])]
    if (length(fill_idx)){
      sel <- is.na(T_met[idx2])
      T_met[fill_idx]  <- Tval[sel]
      Td_met[fill_idx] <- Dval[sel]
    }
  }
  tibble(T_met, Td_met)
}

# Magnus (Alduchov–Eskridge constants)
mag_rh <- function(Tc, Tdc){
  # RH = 100 * es(Td)/es(T)
  100 * exp((17.625*Tdc)/(243.04+Tdc) - (17.625*Tc)/(243.04+Tc))
}

# --- pipeline ---
parsed <- extract_T_Td(asos_utc_sorted$metar)

asos_data8 <- asos_utc_sorted %>%
  mutate(
    valid          = as.POSIXct(valid, tz = "America/Los_Angeles"),
    local_time     = with_tz(valid, tzone = "America/Los_Angeles"),
    date           = as.Date(local_time,tz = "America/Los_Angeles"),
    year           = year(local_time),
    month          = month(local_time, label = TRUE, abbr = TRUE),
    hour           = hour(local_time),
    wind_speed_mps = safe_num(sknt) * 0.514444,
    visibility_km  = safe_num(vsby) * 1.60934
  ) %>%
  bind_cols(parsed) %>%
  mutate(
    # prefer METAR-parsed T/Td; fallback to existing tmpc/dwpc; else NA
    T_use  = coalesce(T_met,  safe_num(tmpc)),
    Td_use = coalesce(Td_met, safe_num(dwpc)),
    RH_calc = if_else(is.na(T_use) | is.na(Td_use),
                      NA_real_,
                      round(mag_rh(T_use, Td_use), 2)),
    rh_source = case_when(
      !is.na(T_met) & !is.na(Td_met) ~ "metar_pair_or_Tgroup",
      is.na(T_met) & !is.na(tmpc) & !is.na(dwpc) ~ "tmpc_dwpc_fallback",
      TRUE ~ "NA"
    ),
    month = factor(month, levels = month.abb, ordered = TRUE)
  ) %>%
  select(-T_met, -Td_met, -T_use, -Td_use)



dust_events_CV <- asos_data8 %>%
  filter(
    # Dust codes in either field
    grepl("\\b(DU|BLDU|DS|SS|VCBLDU)\\b", wxcodes, ignore.case = TRUE) |
      grepl("\\b(DU|BLDU|DS|SS|VCBLDU)\\b", metar,   ignore.case = TRUE) |
      
      # HZ (no FU) with thresholds, from wxcodes
      (wind_speed_mps >= 6 & visibility_km <= 10 & RH_calc < 70 &
         grepl("\\bHZ\\b", wxcodes, ignore.case = TRUE) &
         !grepl("\\bFU\\b", wxcodes, ignore.case = TRUE)) |
      
      # HZ (no FU) with thresholds, from metar
      (wind_speed_mps >= 6 & visibility_km <= 10 & RH_calc < 70 &
         grepl("\\bHZ\\b", metar,   ignore.case = TRUE) &
         !grepl("\\bFU\\b", metar,   ignore.case = TRUE))
  )


dust_events_cv5 <- dplyr::filter(dust_events_CV, is.na(wind_speed_mps) | wind_speed_mps > 0)

#remove all SS+IEM_tags
library(dplyr)

dust_events_cv6 <- dust_events_cv5 %>%
  filter(!grepl("\\(SS\\)\\s*IEM_GHCNH", metar))

