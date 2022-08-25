
# setup -------------------------------------------------------------------

# tidyverse includes dplyr, lubridate, stringr, and some other packages that are generally helpful. Allen normally loads them separately, but I'm lazy and the difference in computational work for the computer to load them all is pretty negligible.

library(tidyverse)


# read_csv() is a little faster than read.csv(), but it's basically the same thing.

full <- read_csv('data/fullDataset_2021-12-01.csv')


# functions ---------------------------------------------------------------

# start by looking at the data

View(full)

# filter() uses a logical statement to pull out specific rows of the data, so this statement selects only rows of the dataframe 'full' where the user ID is 66

filter(full, UserFKOfObserver == 66)

# the pipe (%>%) basically plugs whatever is in front of it into the first argument of the dataframe. All dplyr functions and most recently-written ones for working with dataframes take the dataframe itself as the first argument, so they're pipe-friendly. this is the same as the function above.

full %>% 
  filter(UserFKOfObserver == 66)

# dplyr functions don't make output, they only display it, so if you want to make the output from the function an object, make sure you assign it a name, for example.

user66 <- full %>% 
  filter(UserFKOfObserver == 66)

# most dplyr functions allow you to plug in multiple operations to a single one - good practice is to put each one on a new line. the as.Date functions switches character objects into date objects - if the format isn't yyyy-mm-dd, you need to specify the format.

full %>% 
  filter(
    Name == 'Prairie Ridge Ecostation',
    LocalDate %in% as.Date(c('2018-06-11', '2018-05-14', '2018-05-21')),
    PlantFK == 2757)

# you can find the help for any function like this

?filter

# to pick out specific columns from a dataframe, use select(). distinct() removes duplicate rows, so for example, if I wanted to see what years surveys happened at each site:

full %>% 
  select(Year, Name) %>% 
  distinct()

# pay attention to function chaining with dplyr. Like above, you can always connect the output from one function as the input for the next with a pipe.

# the mutate() function creates a new column in a dataframe. It can be based on existing columns, or not, and it's named with name =. I'm sending this one through to View() so you can see the new columns.

full %>% 
  filter(Name == 'Prairie Ridge Ecostation') %>% 
  mutate(
    # str_c() combines character values into a single string.
    newcolumn1 = str_c(SiteFK, Circle, Code, sep = '_'),
    # this is just  a sum of all the columns listed in each row
    newcolumn2 = Pupa + Hairy + Rolled + Tented,
    # rep repeats something as many times as times =, the '.' represents the existing dataframe as it goes into the function
    newcolumn3 = rep('whee', times = nrow(.))) %>% 
  View()

# summarize() is used to get summary statistics from particular groupings in combination with group_by().

full %>% 
  group_by(Year) %>% 
  summarize(
    # n() just returns the number in each group
    n_rows = n(),
    # total biomass each year - not meaningful but a good example
    total_biomass = sum(Biomass_mg, na.rm = T))


