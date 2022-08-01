library(dplyr)

# make sure data are ready ------------------------------------------------


i <- 0

while(nccovid::get_county_covid_demographics(demographic = 'age_group', region = 'Guilford') %>% pull(week_of) %>% max() <Sys.Date()-6){

  i <- i + 1

  if(i==1) eadeploy::job_done('NEW R estimate update paused; NCDHHS data not ready.',
                              to = c('jennifer.wenner@conehealth.com'))
  break()
  Sys.sleep(60*30)
}


# run script --------------------------------------------------------------


# script estimates R, writes to `output`
source(here::here('src', 'estimate-r.R'))
