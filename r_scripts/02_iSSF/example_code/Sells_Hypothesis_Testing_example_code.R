#### Sells Hypothesis Testing Example Code ####

# Author: Read Barbee

# Date:2023-07-07 

# Purpose:


#########################################################################
##
## 1. Fit iSSF model to each individual
##
##########################################################################

#Make tracks and resmple to 3 hours. Fit global model to each track
tracks <- vector(mode = "list", length = 1)
for (i in 1:length(unique(dat$bear_name))) {
  dat_ <- dat %>% dplyr::filter(sex.id == i) %>% dplyr::arrange(t)
  track <- dat_ %>% amt::make_track(x, y, t, id = bear_name,
                                    crs = sp::CRS("+proj=utm +zone=12 +datum=WGS84")) %>%
    amt::transform_coords(grizmove::raster.crs) %>%
    amt::track_resample(rate = hours(3), tolerance = minutes(45)) %>%
    amt::steps_by_burst() %>% 
    dplyr::mutate(bear_name = first(dat_$bear_name),
                  sex = first(dat_$sex), sex.id = i) %>% 
    dplyr::filter(sl_ >=100) %>% dplyr::filter(sl_ <= 15000)
  track <- sample_bear_steps(dat = track)
  mod <- run_global(track)
  dat.n <- track %>% tidyr::nest(data = -bear_name)
  dat.n <- dat.n %>% dplyr::mutate(mod = list(mod))
  tracks <- rbind(tracks, dat.n)
}

sample_bear_steps <- function (dat)
{
  dat <- dat %>% dplyr::arrange(burst_, t1_, .by_group = T) %>%
    amt::random_steps(10) %>% 
    amt::extract_covariates(cov_stack, where = "end") %>%
    dplyr::mutate(log_sl_ = log(sl_),    cos_ta_ = cos(ta_))
  return(dat)
}

run_global <- function (dat)
{
  mod <- dat %>% amt::fit_issf(case_ ~ ndvi_e + I(ndvi_e^2) +
                                 ruggedness_e + I(ruggedness_e^2) + d2forestedge_e +
                                 I(d2forestedge_e^2) + densforestedge_e + I(densforestedge_e^2) +
                                 densriparian_e + I(densriparian_e^2) + densbuildings_e +
                                 I(densbuildings_e^2) + d2core_e + I(d2core_e^2) + sl_ +
                                 log_sl_ + cos_ta_ + strata(step_id_), model = TRUE)
  return(mod)
}


#########################################################################
##
## 2. Calculate log-RSS and classify individual responses to each covariate
##
##########################################################################

# Functions for griz log-RSS:

set_df_s1 <- function ()
{
  s1 <- data.frame(sl_ = 100, log_sl_ = log(100), cos_ta_ = 1,
                   ruggedness_e = seq(from = -2, to = 2, length.out = 200),
                   d2forestedge_e = 0, densforestedge_e = 0, ndvi_e = 0,
                   densriparian_e = 0, densbuildings_e = 0, d2core_e = 0,
                   ruggedness_s = 0, d2forestedge_s = 0, densforestedge_s = 0,
                   ndvi_s = 0, densriparian_s = 0, densbuildings_s = 0,
                   d2core_s = 0, `(sl_ + log_sl_ + cos_ta_)` = 0)
  s1$ruggedness_e <- 0
  return(s1)
}

set_df_s2 <- function ()
{
  s2 <- data.frame(sl_ = 100, log_sl_ = log(100), cos_ta_ = 1,
                   ruggedness_e = 0, d2forestedge_e = 0, densforestedge_e = 0,
                   ndvi_e = 0, densriparian_e = 0, densbuildings_e = 0,
                   d2core_e = 0, ruggedness_s = 0, d2forestedge_s = 0,
                   densforestedge_s = 0, ndvi_s = 0, densriparian_s = 0,
                   densbuildings_s = 0, d2core_s = 0, `(sl_ + log_sl_ + cos_ta_)` = 0)
  return(s2)
}

classify_results <- function (tracks, curr_param)
{
  classifications <- data.frame(bear_name = NA, uncertain = NA,
                                positive = NA)
  classifications_ <- classifications
  for (i in 1:nrow(tracks)) {
    mod <- tracks$mod[[i]]
    data <- tracks$data[[i]]
    s1 <- set_df_s1()
    if (curr_param == "NDVI") {
      s1$ndvi_e <- seq(from = min(data$ndvi_e, na.rm = T),
                       to = max(data$ndvi_e, na.rm = T), length.out = 200)
    }
    if (curr_param == "Ruggedness") {
      s1$ruggedness_e <- seq(from = min(data$ruggedness_e,
                                        na.rm = T), to = max(data$ruggedness_e, na.rm = T),
                             length.out = 200)
    }
    if (curr_param == "D2forestedge") {
      s1$d2forestedge_e <- seq(from = min(data$d2forestedge_e,
                                          na.rm = T), to = max(data$d2forestedge_e, na.rm = T),
                               length.out = 200)
    }
    if (curr_param == "Densforestedge") {
      s1$densforestedge_e <- seq(from = min(data$densforestedge_e,
                                            na.rm = T), to = max(data$densforestedge_e,
                                                                 na.rm = T), length.out = 200)
    }
    if (curr_param == "Densriparian") {
      s1$densriparian_e <- seq(from = min(data$densriparian_e,
                                          na.rm = T), to = max(data$densriparian_e, na.rm = T),
                               length.out = 200)
    }
    if (curr_param == "Densbuildings") {
      s1$densbuildings_e <- seq(from = min(data$densbuildings_e,
                                           na.rm = T), to = max(data$densbuildings_e, na.rm = T),
                                length.out = 200)
    }
    if (curr_param == "D2core") {
      s1$d2core_e <- seq(from = min(data$d2core_e, na.rm = T),
                         to = max(data$d2core_e, na.rm = T), length.out = 200)
    }
    lr._ci_se <- amt::log_rss(mod, s1, s2, ci = "se", ci_level = 0.95)
    classifications_$bear_name <- tracks$bear_name[[i]]
    count_pos <- lr._ci_se$df %>% summarize(count_pos = sum(log_rss >
                                                              0))
    classifications_$positive <- ifelse(count_pos > 100,
                                        T, F)
    count_overlap <- lr._ci_se$df %>% rowwise %>% dplyr::summarize(count_overlap = sum(lwr <=
                                                                                         0 && upr >= 0)) %>% dplyr::filter(count_overlap ==
                                                                                                                             1)
    count_overlap <- nrow(count_overlap)
    classifications_$uncertain <- ifelse(count_overlap >
                                           100, T, F)
    classifications <- rbind(classifications, classifications_)
  }
  classifications <- classifications[-1, ]
  classifications <- classifications %>% dplyr::group_by(bear_name) %>%
    dplyr::arrange(bear_name)
  classifications$plot.order <- 1:nrow(classifications)
  return(classifications)
}






