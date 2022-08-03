# Purpose: Calculate the reproduction numbers for:
# * North Carolina
# * All NC Counties
# * Combined Cone Health Region

# libraries ---------------------------------------------------------------
library(EpiEstim)
library(data.table)
library(dplyr)



# configs -----------------------------------------------------------------



start_date <- as.Date('2021-12-18') #omicron became dominant


#what serial interval to use?  going with this publication for now:

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8832521/ mean = 3.5, std = 2.4
#Michael originally had mean = 5.8, std = 2

config_lit <- make_config(
  mean_si = 3.5,
  std_si = 2.4
)

#possible si for predictions:
# hist(rgamma(1000, 3, scale = 1))
si <- distcrete::distcrete("gamma", interval = 1L,
                           shape = 3,
                           scale = 1,
                           w = 0.5)



# Loop through and collect outputs

# all counties ------------------------------------------------------------


county_target <- nccovid::nc_county_fips$county
collector <-list()

cases_all <- nccovid::get_county_covid_demographics(demographic = 'age_group')
cases_all <- cases_all %>% filter(week_of >= start_date) %>%
  mutate(cases_weekly = cases, deaths_weekly = deaths,
         date = week_of, state = 'North Carolina') %>%
  group_by(state, county, date) %>%
  summarise(cases_weekly = sum(cases_weekly)) %>%
  mutate(date = as.Date(date))

#turn weekly time series into daily; assume equal cases rate during the week?
dates <- seq.Date(min(cases_all$date),
                  max(cases_all$date), '1 day')


combos <- expand.grid(state= 'North Carolina',
                      county = unique(cases_all$county),
                      date = dates)
cases_all <- left_join(combos, cases_all,
                       by = c('state', 'county', 'date')) %>%
  group_by(state, county) %>% arrange(date) %>%
  mutate(cases_daily = round(cases_weekly/7),
         cases_daily = zoo::na.locf(cases_daily, na.rm = F)) %>%
  setDT()

project_incidence <- T

for(i in county_target){

  cases_county <- cases_all[county==i]

  epi_start <- min(cases_county[cases_daily>1]$date)

  cases <- cases_county[date>=epi_start,list(cases_daily,date)]

  names(cases) <- c("I", "dates")

  epiestim_res_lit <- estimate_R(
    incid = cases,
    method = "parametric_si",
    config = config_lit
  ) %>% suppressMessages()
#suppressing the following:
    #   Default config will estimate R on weekly sliding windows.
    #   To change this change the t_start and t_end arguments.



  #simple R projection using arima:
  f <- forecast::auto.arima(epiestim_res_lit$R$`Median(R)`)
  f <- forecast::forecast(f,h=21,level = c(75,90,95)) # 3 weeks
  f <- as.data.frame(f) %>%
    transmute(t_start = NA, t_end = NA, `Mean(R)` = `Point Forecast`, `Std(R)`= NA,
              `Quantile.0.025(R)` = `Lo 95`, `Quantile.0.05(R)` = `Lo 90`,
              `Quantile.0.25(R)` = `Lo 75`, `Median(R)` = `Point Forecast`,
              `Quantile.0.75(R)` = `Hi 75`, `Quantile.0.95(R)` = `Hi 90`,
              `Quantile.0.975(R)` = `Hi 95`,
              dates = tail(epiestim_res_lit$dates,1)+1:21,
              type = 'forecast', variable = 'R'
              )

  if(project_incidence){

    mean_r <- tail(epiestim_res_lit$R$`Mean(R)`,1)
    sd_r <- tail(epiestim_res_lit$R$`Std(R)`,1)
    shapescale <- epitrix::gamma_mucv2shapescale(mu = mean_r, cv = sd_r/mean_r)
    plausible_r <- rgamma(1000, shape = shapescale$shape, scale = shapescale$scale)
  inc <- unlist(mapply(date = cases$dates, reps = cases$I,
                       FUN = function(date, reps) rep(date, reps)))
  class(inc) <- 'Date'
  inc <- incidence::incidence(inc)
  proj <- projections::project(x = inc,
                               R = plausible_r,
                               si = si,
                               n_days = 21, n_sim = 1000)
  #plot(proj)

  proj <- as.matrix(proj)
  proj <- proj %>% apply(MARGIN = 1,
                         FUN = function(x)list(mean = round(mean(x)),
                                               sd = round(sd(x)),
                                               q = quantile(x, c(q025 = 0.025,q05 = 0.05,q25 = 0.25,
                                                                 med = .5,q75 = .75,q95 = 0.95, q975 = 0.975))))
  proj <- lapply(proj, function(x) x %>% unlist() %>% rbind() %>% as_tibble())
  proj <- rbindlist(proj) %>%
    transmute(t_start = NA, t_end = NA, `Mean(R)` = mean, `Std(R)`= sd,
              `Quantile.0.025(R)` = `q.2.5%`, `Quantile.0.05(R)` = `q.5%`,
              `Quantile.0.25(R)` = `q.25%`, `Median(R)` = `q.50%`,
              `Quantile.0.75(R)` = `q.75%`, `Quantile.0.95(R)` = `q.95%`,
              `Quantile.0.975(R)` = `q.97.5%`,
              dates = tail(epiestim_res_lit$dates,1)+1:21,
              type = 'forecast', variable = 'reported_cases'
    ) %>%
    bind_rows(cases %>% transmute(dates, `Mean(R)` = I, `Median(R)` = I,
                                  type = 'actual', variable = 'reported_cases'))





  }

  tmp_out <- rbind(cbind(epiestim_res_lit$R, dates = tail(epiestim_res_lit$dates,-7),
                         type = 'estimate', variable = 'R'),
                   f)

  if(project_incidence) tmp_out <- rbind(tmp_out, proj)

  collector[[i]] <- tmp_out
  cat('\nFinished', i,'---------------------------------------------------------\n')
}

out <- rbindlist(collector, idcol = "county")


# one for nc --------------------------------------------------------------


cases_all_nc <- cases_all %>% group_by(state, date) %>%
  summarise(cases_daily = sum(cases_daily)) %>%
  ungroup() %>% setDT()

epi_start <- min(cases_all_nc[cases_daily>1]$date)

cases <- cases_all_nc[date>=epi_start,list(cases_daily,date)]

cases[,date:=as.Date(date)]

names(cases) <- c("I", "dates")

epiestim_res_lit <- estimate_R(
  incid = cases,
  method = "parametric_si",
  config = config_lit
)

#simple R projection using arima:
f <- forecast::auto.arima(epiestim_res_lit$R$`Median(R)`)
f <- forecast::forecast(f,h=21,level = c(75,90,95)) # 3 weeks
f <- as.data.frame(f) %>%
  transmute(t_start = NA, t_end = NA, `Mean(R)` = `Point Forecast`, `Std(R)`= NA,
            `Quantile.0.025(R)` = `Lo 95`, `Quantile.0.05(R)` = `Lo 90`,
            `Quantile.0.25(R)` = `Lo 75`, `Median(R)` = `Point Forecast`,
            `Quantile.0.75(R)` = `Hi 75`, `Quantile.0.95(R)` = `Hi 90`,
            `Quantile.0.975(R)` = `Hi 95`,
            dates = tail(epiestim_res_lit$dates,1)+1:21,
            type = 'forecast', variable = 'R'
  )

if(project_incidence){

  mean_r <- tail(epiestim_res_lit$R$`Mean(R)`,1)
  sd_r <- tail(epiestim_res_lit$R$`Std(R)`,1)
  shapescale <- epitrix::gamma_mucv2shapescale(mu = mean_r, cv = sd_r/mean_r)
  plausible_r <- rgamma(1000, shape = shapescale$shape, scale = shapescale$scale)
  inc <- unlist(mapply(date = cases$dates, reps = cases$I,
                       FUN = function(date, reps) rep(date, reps)))
  class(inc) <- 'Date'
  inc <- incidence::incidence(inc)
  proj <- projections::project(x = inc,
                               R = plausible_r,
                               si = si,
                               n_days = 21, n_sim = 1000)
  #plot(proj)

  proj <- as.matrix(proj)
  proj <- proj %>% apply(MARGIN = 1,
                         FUN = function(x)list(mean = round(mean(x)),
                                               sd = round(sd(x)),
                                               q = quantile(x, c(q025 = 0.025,q05 = 0.05,q25 = 0.25,
                                                                 med = .5,q75 = .75,q95 = 0.95, q975 = 0.975))))
  proj <- lapply(proj, function(x) x %>% unlist() %>% rbind() %>% as_tibble())
  proj <- rbindlist(proj) %>%
    transmute(t_start = NA, t_end = NA, `Mean(R)` = mean, `Std(R)`= sd,
              `Quantile.0.025(R)` = `q.2.5%`, `Quantile.0.05(R)` = `q.5%`,
              `Quantile.0.25(R)` = `q.25%`, `Median(R)` = `q.50%`,
              `Quantile.0.75(R)` = `q.75%`, `Quantile.0.95(R)` = `q.95%`,
              `Quantile.0.975(R)` = `q.97.5%`,
              dates = tail(epiestim_res_lit$dates,1)+1:21,
              type = 'forecast', variable = 'reported_cases'
    ) %>%
    bind_rows(cases %>% transmute(dates, `Mean(R)` = I, `Median(R)` = I,
                                  type = 'actual', variable = 'reported_cases'))





}

out_nc <- rbind(cbind(epiestim_res_lit$R, dates = tail(epiestim_res_lit$dates,-7),
                      type = 'estimate', variable = 'R'),
                f) %>% rbind(proj)

# cone service region -----------------------------------------------------

cases_all_cone <- cases_all %>% filter(county %in% nccovid::cone_region)

cases_all_cone <- cases_all_cone[,.(cases_reported = sum(cases_daily)), by = "date"]

epi_start <- min(cases_all_cone[cases_reported>1]$date)

cases <- cases_all_cone[date>=epi_start,list(cases_reported,date)]

cases[,date:=as.Date(date)]

names(cases) <- c("I", "dates")

epiestim_res_lit <- estimate_R(
  incid = cases,
  method = "parametric_si",
  config = config_lit
)


#simple R projection using arima:
f <- forecast::auto.arima(epiestim_res_lit$R$`Median(R)`)
f <- forecast::forecast(f,h=21,level = c(75,90,95)) # 3 weeks
f <- as.data.frame(f) %>%
  transmute(t_start = NA, t_end = NA, `Mean(R)` = `Point Forecast`, `Std(R)`= NA,
            `Quantile.0.025(R)` = `Lo 95`, `Quantile.0.05(R)` = `Lo 90`,
            `Quantile.0.25(R)` = `Lo 75`, `Median(R)` = `Point Forecast`,
            `Quantile.0.75(R)` = `Hi 75`, `Quantile.0.95(R)` = `Hi 90`,
            `Quantile.0.975(R)` = `Hi 95`,
            dates = tail(epiestim_res_lit$dates,1)+1:21,
            type = 'forecast', variable = 'R'
  )


if(project_incidence){

  mean_r <- tail(epiestim_res_lit$R$`Mean(R)`,1)
  sd_r <- tail(epiestim_res_lit$R$`Std(R)`,1)
  shapescale <- epitrix::gamma_mucv2shapescale(mu = mean_r, cv = sd_r/mean_r)
  plausible_r <- rgamma(1000, shape = shapescale$shape, scale = shapescale$scale)
  inc <- unlist(mapply(date = cases$dates, reps = cases$I,
                       FUN = function(date, reps) rep(date, reps)))
  class(inc) <- 'Date'
  inc <- incidence::incidence(inc)
  proj <- projections::project(x = inc,
                               R = plausible_r,
                               si = si,
                               n_days = 21, n_sim = 1000)
  #plot(proj)

  proj <- as.matrix(proj)
  proj <- proj %>% apply(MARGIN = 1,
                         FUN = function(x)list(mean = round(mean(x)),
                                               sd = round(sd(x)),
                                               q = quantile(x, c(q025 = 0.025,q05 = 0.05,q25 = 0.25,
                                                                 med = .5,q75 = .75,q95 = 0.95, q975 = 0.975))))
  proj <- lapply(proj, function(x) x %>% unlist() %>% rbind() %>% as_tibble())
  proj <- rbindlist(proj) %>%
    transmute(t_start = NA, t_end = NA, `Mean(R)` = mean, `Std(R)`= sd,
              `Quantile.0.025(R)` = `q.2.5%`, `Quantile.0.05(R)` = `q.5%`,
              `Quantile.0.25(R)` = `q.25%`, `Median(R)` = `q.50%`,
              `Quantile.0.75(R)` = `q.75%`, `Quantile.0.95(R)` = `q.95%`,
              `Quantile.0.975(R)` = `q.97.5%`,
              dates = tail(epiestim_res_lit$dates,1)+1:21,
              type = 'forecast', variable = 'reported_cases'
    ) %>%
    bind_rows(cases %>% transmute(dates, `Mean(R)` = I, `Median(R)` = I,
                                  type = 'actual', variable = 'reported_cases'))





}

out_cone <- rbind(cbind(epiestim_res_lit$R, dates = tail(epiestim_res_lit$dates,-7),
                        type = 'estimate', variable = 'R'),
                  f) %>% rbind(proj)


# put it all together -----------------------------------------------------


#keeping names consistent with previous R estimates for easy swap of EpiNow2 estimates with these EpiEstim estimates
final <- bind_rows(out,
                   out_cone %>% mutate(county = 'Cone Health'),
                   out_nc %>% mutate(county = 'North Carolina')
                   ) %>% mutate(last_update = Sys.Date(),
                                strat = '') %>%
  select(county, date = dates, last_update, variable, strat, type,
         median = `Median(R)`, mean = `Mean(R)`,
         sd = `Std(R)`, bottom = `Quantile.0.025(R)`,
         lower = `Quantile.0.05(R)`, central_lower = `Quantile.0.25(R)`,
         central_upper = `Quantile.0.75(R)`, upper = `Quantile.0.95(R)`,
         top = `Quantile.0.975(R)`)


# library(ggplot2)
# ggplot(final %>% filter(date > Sys.Date() - 90)) +
#   geom_line(aes(x = date, y = median, color = county, linetype = type))



# write data --------------------------------------------------------------


write.csv(final,
          here::here('output','latest_r_coviddata.csv'),
          row.names = F)

