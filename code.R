

# load data and libraries -------------------------------------------------
setwd("/Users/jemutai/Downloads/chicks")

pacman::p_load(readr,
               lubridate,
               dplyr,
               tidyverse)


khis <- readr::read_csv("data 2/CHIKV_data_KHIS.csv")
county_fever_data <- readr::read_csv("data 2/KHIS_data_fever_county.csv")
data_fever <- readr::read_csv("data 2/KHIS_data_fever.csv")


# data cleaning -----------------------------------------------------------

#remove columns with only NA
khis <- khis[, !(names(khis)%in% c("organisationunitdescription","perioddescription"))]
county_fever_data <- county_fever_data[, !(names(county_fever_data)%in% c("organisationunitdescription","perioddescription"))]

data_fever <- data_fever[, sapply(data_fever, function(x) sum(is.na(x) | x == "")) != 3395]


#clean dates 

khis <- khis |>
  mutate(period = parse_date_time(periodname, orders = "B Y"),
         period = format(period, "%Y-%m")) 

#visualize the epi curve for chikungunya 
# Rename column for easier use
khis <- khis %>%
  rename(chik_cases = `MOH 705A Rev 2020_ Chikungunya`)

epi_county <- khis %>%
  group_by(organisationunitname, period) %>%
  summarise(cases = sum(chik_cases, na.rm = TRUE)) %>%
  ungroup()

ggplot(epi_county, aes(x = period, y = organisationunitname, fill = cases)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Chikungunya Epi Curve by County (Heatmap)",
       x = "Period", y = "County",
       fill = "Cases") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Summarise total cases per county
county_cases <- khis %>%
  group_by(organisationunitname) %>%
  summarise(total_cases = sum(chik_cases, na.rm = TRUE)) %>%
  arrange(desc(total_cases))

# Horizontal bar chart
ggplot(county_cases, aes(x = total_cases, 
                         y = reorder(organisationunitname, total_cases))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Total Chikungunya Cases by County",
       x = "Total Cases",
       y = "County") +
  theme_minimal()
