################################################################################
#Function built by Josh to split long data into survey periods/intervals
#DO NOT EDIT
################################################################################


make_interval <- function(x, dt_col, time_int = "month", increment = 1) {
  stopifnot(tibble::is_tibble(x))
  tmp <- dplyr::pull(x, {{dt_col}})
  if (time_int == "hour") {
    stopifnot(inherits(tmp, c("POSIXct", "POSIXlt")))
  }
  stopifnot(all(!is.na(tmp)))
  stopifnot(time_int %in% c("hour", "day", "week"))
  if (dplyr::is_grouped_df(x)) {
    args_passed <- as.list(match.call())
    stop(
      "Data is grouped, try something like",
      "\n\t",
      paste0(
        "x %>% ",
        "dplyr::do(gbn_make_interval(",
        args_passed$dt_col, ", ",
        "'", args_passed$time_int, "', ",
        args_passed$increment, "))"
      )
    )
  }
  
  int <- lubridate::interval(min(tmp), tmp)
  
  x %>%
    dplyr::mutate(
      rnd_dt = lubridate::round_date(tmp, unit = paste(increment, time_int)),
      interval = switch(
        time_int,
        "week" = int %/% lubridate::weeks(increment),
        "day" = int %/% lubridate::days(increment),
        "hour" = int %/% lubridate::hours(increment)
      ) + 1
    )
}