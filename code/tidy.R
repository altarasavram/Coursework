library (tidyverse)
gapminder_orig <- read.csv("https://raw.githubusercontent.com/swcarpentry/r-novice-gapminder/gh-pages/_episodes_rmd/data/gapminder-FiveYearData.csv")
gapminder <- gapminder_orig

gapminder %>% map_chr (class)
gapminder %>% map_dbl (n_distinct)

gapminder %>% map_df(~(data.frame (n_distinct = n_distinct (.x), 
                                   class = class(.x))),   
                                  .id = "variable")

# Do same incrementally: 1) let x be the first column of gapminder
.x <- gapminder %>% pluck (1)
data.frame (n_distinct = n_distinct(.x),
            class = class(.x))

# 2) generalize to all columns
gapminder %>% map_df(~(data.frame (n_distinct = n_distinct (.x), 
                                   class = class(.x))))

###  map2 (x, y, function)
continent_year <- gapminder %>% distinct(continent, year)
continents <- continent_year %>% pull (continent)
years      <- continent_year %>% pull (year)

# Again incrementally
.x <- continents[1]
.y <- years[1]

gapminder %>% filter (continent == .x, 
                      year == .y) %>%
  ggplot () +
  geom_point (aes(x=gdpPercap, y=lifeExp)) +
  labs (title = paste (.x, .y))
# 2) generalize

plot_list <- map2 (.x = continents,
                   .y = years,
                   .f=~(
                     gapminder %>% filter (continent == .x, 
                                           year == .y) %>%
                       ggplot () +
                       geom_point (aes(x=gdpPercap, y=lifeExp)) +
                       labs (title = paste (.x, .y))
              ))

