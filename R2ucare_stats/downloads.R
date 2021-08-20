# packages
library(tidyverse)
theme_set(theme_light())
library(lubridate)
library(cranlogs)

# get downloads
Sys.setlocale(locale = "en_US.UTF-8")
R2ucare_downloads <- cran_downloads(packages = "R2ucare",
                                    from = "2017-04-13",
                                    to = Sys.Date()-1)

# dataviz
plot <- R2ucare_downloads %>%
  mutate(date = ymd(date)) %>%
  group_by(month = floor_date(date, unit = "month")) %>%
  mutate(stats = sum(count)) %>%
  ggplot() +
  aes(x = month, y = stats) +
  geom_line() +
  geom_smooth() + # add GAM for trend
  scale_x_date(date_labels = "%b %y", breaks = "month", expand=c(0,0)) +
  labs(x = NULL, y = "number of downloads per month", title = "R2ucare downloads") +
  theme(axis.text.x  = element_text(size = 9,
                                    angle = 45,
                                    colour = "black",
                                    vjust = 1,
                                    hjust = 1))
plot
ggsave(plot = plot, filename = "R2ucare_stats/R2ucare_downloads.png", width = 12, height = 8, dpi = 300)


