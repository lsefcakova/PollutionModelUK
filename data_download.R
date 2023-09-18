library(openair)
library(tidyverse)

meta_data <- importMeta(source = c("aurn","local","kcl",'aqe','saqn','waqn','ni'), all = TRUE,year=2020)

selected_data <- meta_data %>%
  filter(variable == "NO")

selected_sites <- selected_data %>%
  dplyr::select(code) %>%
  mutate_all(.funs = tolower)

data <- importAURN(site = selected_sites$code, year = 2021)

filtered_data <- meta_data %>%
  filter(variable == "PM10") %>%
  right_join(data, "code") %>%
  dplyr::select(date, pm10, code, site.x, site_type, latitude, longitude, variable) %>%
  rename(site = site.x)

write_csv(filtered_data, "filtered_data.csv")

