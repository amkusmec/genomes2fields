library(tidyverse)
library(lubridate)


# Load the various data sets ----------------------------------------------
munge <- read_rds("data/weather/weather_munged.rds")
prism <- read_rds("data/weather/weather_prism.rds")
power <- read_rds("data/weather/weather_power.rds")
metadata <- read_csv("data/metadata/metadata_clean.csv")


# Get site elevation ------------------------------------------------------
sites <- with(prism, paste(Environment, Year, sep = "_")) %>% unique() %>% sort()
prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
metadata <- metadata %>%
  filter(paste(Environment, Year, sep = "_") %in% sites) %>%
  rename(x = Longitude, y = Latitude) %>%
  select(x, y, Year, Environment)
elev <- elevatr::get_elev_point(as.data.frame(metadata), prj = prj, src = "epqs")
metadata <- metadata %>%
  select(Year, Environment, x, y) %>%
  rename(Longitude = x, Latitude = y) %>%
  inner_join(., elev@data, by = c("Year", "Environment")) %>%
  select(-elev_units) %>%
  rename(Elevation = elevation)


# Summarize the G2F weather data ------------------------------------------
munge <- munge %>%
  group_by(Year, Environment, Month, Day) %>%
  summarise(TMIN = min(TEMP, na.rm = TRUE), 
            TMAX = max(TEMP, na.rm = TRUE), 
            PPT = sum(PPT, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(paste(Environment, Year, sep = "_") %in% sites)


# Impute data from PRISM --------------------------------------------------
prism <- prism %>%
  mutate(Key = paste(Environment, Year, Month, Day, sep = "_"))

munge <- munge %>%
  mutate(Key = paste(Environment, Year, Month, Day, sep = "_")) %>%
  filter(!(Key %in% prism$Key)) %>%
  bind_rows(., prism)


# Combine with data from POWER --------------------------------------------
munge <- munge %>%
  select(-Key) %>%
  inner_join(., power, by = c("Year", "Environment", "Month", "Day")) %>%
  inner_join(., metadata, by = c("Year", "Environment")) %>%
  rename(z = Elevation)


# Estimate missing radiation values ---------------------------------------
Gsc <- 0.0820 # Global solar constant [MJ/m^2/min]
munge <- munge %>%
  mutate(J = paste(Year, Month, Day, sep = "-") %>% ymd() %>% yday(), 
         phi = (pi/180)*Latitude, 
         dr = 1 + 0.033*cos(2*pi*J/365), 
         d = 0.409*sin(2*pi*J/365 - 1.39), 
         ws = acos(-tan(phi)*tan(d)), 
         Ra = (24*60/pi)*Gsc*dr*(ws*sin(phi)*sin(d) + cos(phi)*cos(d)*sin(ws)), 
         Rso = (0.75 + z*2e-5)*Ra, 
         N = 24*ws/pi, 
         Rs = (0.25 + 0.50*(SG_DAY_HOUR_AVG/N))*Ra) %>%
  mutate(ALLSKY_SFC_SW_DWN = if_else(is.na(ALLSKY_SFC_SW_DWN), Rs, ALLSKY_SFC_SW_DWN),
         CLRSKY_SFC_SW_DWN = if_else(is.na(CLRSKY_SFC_SW_DWN), Rso, CLRSKY_SFC_SW_DWN)) %>%
  select(Year:z)


# Estimate daily evapotranspiration ---------------------------------------
# Uses the Penman-Monteith equation
l <- 2.45
cp <- 1.013e-3
e <- 0.622
a <- 0.23
s <- 4.903e-9
munge <- munge %>%
  mutate(P = 101.3*((293 - 0.0065*z)/293)^5.26, 
         g = (cp*P)/(e*l), 
         Rns = (1 - a)*ALLSKY_SFC_SW_DWN, 
         ea = 0.6108*exp(17.27*T2MDEW/(T2MDEW + 237.3)), 
         Rnl = s*0.5*((T2M_MAX + 273.16)^4 + (T2M_MIN + 273.16)^4)*
           (0.34 - 0.14*sqrt(ea))*(1.35*(ALLSKY_SFC_SW_DWN/CLRSKY_SFC_SW_DWN) - 0.35),
         Rn = Rns - Rnl,
         es = 0.5*(0.6108*exp(17.27*T2M_MAX/(T2M_MAX + 237.3)) +
                     0.6108*exp(17.27*T2M_MIN/(T2M_MIN + 237.3))), 
         TAVG = 0.5*(T2M_MAX + T2M_MIN), 
         d = 4098*(0.6108*exp(17.27*TAVG/(TAVG + 237.3)))/(TAVG + 237.3)^2, 
         ET = (0.408*d*Rn + g*(900/(TAVG + 273))*WS2M*(es - ea))/(d + g*(1 + 0.34*WS2M))) %>%
  select(Year:ALLSKY_SFC_SW_DWN, SG_DAY_HOUR_AVG, ET) %>%
  rename(SR = ALLSKY_SFC_SW_DWN)


# Calculate developmental time in photothermal time -----------------------
# PTT = SG_DAY_HOUR_AVG*CHU
munge <- munge %>%
  mutate(Tm = if_else(TMIN < 4.4, 4.4, TMIN), 
         TM = if_else(TMAX < 10, 10, TMAX), 
         CHUm = 1.8*(Tm - 4.4), 
         CHUM = 3.33*(TM - 10) - 0.084*(TM - 10)^2, 
         CHU = 0.5*(CHUM + CHUm),
         PTT = SG_DAY_HOUR_AVG*CHU) %>%
  select(Year:SR, ET, PTT)


# Calculate net precipitation since planting ------------------------------
munge <- munge %>%
  mutate(NET = PPT - ET)

write_rds(munge, "data/weather/env_variables.rds")


# Some plots of environmental variables -----------------------------------
munge %>%
  mutate(Site = paste(Environment, Year, sep = "_")) %>%
  gather(Variable, Value, TMIN:NET) %>%
  ggplot(., aes(x = Site, y = Value)) + 
    theme_bw() + geom_boxplot() +
    facet_wrap(~ Variable, scales = "free_y", ncol = 1, strip.position = "r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/select/evar_all_boxplots.pdf", width = 14, height = 12, 
       units = "in", dpi = 300)
