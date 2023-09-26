library(crawl)
library(sfnetworks)
library(pathroutr)
library(dplyr)
library(sf)
library(ggplot2)
library(ggspatial)
library(dplyr)

data("akcoast")

harborSeal_sf <- harborSeal_sf %>% sf::st_transform(3338)

# l <- harborSeal_sf %>% 
#   dplyr::filter(!sf::st_is_empty(.)) %>% 
#   dplyr::summarise(do_union = FALSE) %>% 
#   sf::st_cast('LINESTRING')

land_region <- sf::st_buffer(harborSeal_sf, dist = 50000) %>% 
  sf::st_union() %>% 
  sf::st_convex_hull() %>% 
  sf::st_intersection(akcoast) %>% 
  st_collection_extract('POLYGON') %>% 
  st_sf()

vis_graph <- prt_visgraph(land_region)


fixPar = c(log(250), log(500), log(1500), rep(NA,5), 0)
displayPar( mov.model=~1, err.model=list(x=~Argos_loc_class-1),data=harborSeal_sf, 
            activity=~I(1-DryTime),fixPar=fixPar)

constr=list(
  lower=c(rep(log(1500),3), rep(-Inf,2)),
  upper=rep(Inf,5)
)

set.seed(123)
fit1 <- crwMLE(
  mov.model=~1, err.model=list(x=~Argos_loc_class-1), activity=~I(1-DryTime),
  data=harborSeal_sf,  Time.name="Time", 
  fixPar=fixPar, theta=c(rep(log(5000),3),log(3*3600), 0),
  constr=constr, method="L-BFGS-B",
  control=list(maxit=2000, trace=1, REPORT=1)
)


pred1 = crwPredict(fit1, predTime = '15 min')
pred1_sf <- pred1 %>% crw_as_sf("POINT","p")

segs_tbl <- get_barrier_segments(pred1_sf,land_region)
segs_tbl

segs_tbl <- segs_tbl %>% prt_shortpath(vis_graph)

track_pts_fix <- prt_reroute(pred1_sf, land_region, vis_graph)

track_pts_fix <- prt_update_points(track_pts_fix, pred1_sf)

track_line_fixed <- track_pts_fix %>% summarise(do_union = FALSE) %>% st_cast('LINESTRING')



