#### Fixed Rank Kriging Craigslist data for 2019 ------------------------------

library(dplyr)
library(ggplot2)
library(yaml)
library(DBI)
library(FRK)
library(sp)
library(sf)
library(patchwork)

#set wd to project base folder
setwd("H:/natrent-krige")

#store credentials at base dir of natrent-city-sub as YAML
cred <- read_yaml("./natrent0.yaml")

#create a database connection
natrent <- dbConnect(
  drv = RPostgres::Postgres(),
  dbname = "natrent",
  host = "natrent0.csde.washington.edu",
  user = names(cred),
  password = cred[[names(cred)]],
  bigint = "numeric"
)

#take advantage of the sim clusters CPUs
opts_FRK$set("parallel", 4L)

#read cbsa shapefile
cbsa <- read_sf("./input/ACS CBSA Polygon/US_cbsa_2017.shp") %>%
  st_set_crs(102008) %>% #first properly identify the CRS
  st_transform(4326) #then reproject

#### A. Query Neighborhood Estimates of CL Listing Activity --------------------

#query for Washington DC listings as test case
query <- "SELECT DISTINCT d.listing_date AS date, d.clean_beds AS beds, d.clean_sqft AS sqft, 
                          d.clean_rent AS rent, c.state, c.county,
                          ROUND(CAST(ST_X(ST_TRANSFORM(d.geometry, 4326)) as numeric), 3) as lng, 
                          ROUND(CAST(ST_Y(ST_TRANSFORM(d.geometry, 4326)) as numeric), 3) as lat
          FROM (
                SELECT a.gisjoin AS trt_id, a.geometry, b.cbsafp AS met_id, b.statefp AS state, 
                       b.countyfp AS county
                FROM tract17 a
                JOIN county17 b ON a.statefp = b.statefp AND a.countyfp = b.countyfp
                WHERE b.cbsafp = ?cbsa_code
          ) c
          LEFT JOIN clean d ON ST_Contains(c.geometry, d.geometry) 
          WHERE d.listing_date BETWEEN ?start AND ?end AND
                d.match_type NOT IN ('No Address Found', 'Google Maps Lat/Long') AND
                d.clean_beds = ?beds AND d.clean_rent IS NOT NULL
          ORDER BY d.listing_date DESC"


#set temporal cutoffs, bed size and cbsa (can be input args later)
start_date <- "2019-01-01"
cutoff_date <- "2019-12-31"
bed_size <- 1
cbsa_code <- "42660"

#work these into the query
query <- sqlInterpolate(natrent, query, 
                           start = start_date,
                           end = cutoff_date,
                           beds = bed_size,
                           cbsa_code = cbsa_code)

#submit the query
query <- dbGetQuery(natrent, query)

#reduce to unique lat/lng via aggregation
query <- query %>%
  group_by(lat, lng) %>%
  summarize(rent = median(rent),
            sqft = median(sqft),
            state = unique(state),
            county = unique(county))

#initialize as Spatial object
coordinates(query) <- ~ lng + lat

#set the projection to WGS84 lat/lng
raster::projection(query) <- raster::crs(4326)


#### B. Estimate Spatial Random Effects Model for listing data ----------------

#e.g., points for 1B colored by rent
#ggplot(sf::st_as_sf(query), 
#       aes(color = rent)) +
#  geom_sf() +
#  scale_color_viridis_c()

#run automated fixed rank kriging via FRK library
#frk_formula <- log(rent) ~ 1 + beds
#S <- FRK(f = frk_formula, 
#         data = list(query))

#predict values to obtain surface
#Pred <- predict(S, obs_fs = FALSE) %>%
#  st_as_sf()

#plot surface
#ggplot() +
#  geom_tile(data = Pred, aes(x = lng, y = lat, fill = exp(mu)),
#            color = "grey80", lwd = 0.01) +
#  geom_point(data = as.data.frame(query@coords), 
#                  aes(x = lng, y = lat),
#             pch = 21, alpha = .5) +
#  scale_fill_viridis_c() +
#  coord_fixed() +
#  theme_bw()

#fine grained estimation of SRE
set.seed(1)
GridBAUs <- auto_BAUs(manifold = plane(), # 2D plane
                      type = "hex", # hex grid
                      data = query) # data around which to create BAUs

GridBAUs$fs <- 1 # fine-scale variation at BAU level

G <- auto_basis(manifold = plane(), # 2D plane
                data = query, # queried data
                nres = 3, # number of resolutions
                type = "Gaussian", # type of basis function
                regular = 1) # place regularly in domain

#show_basis(G) + # illustrate basis functions
#  coord_fixed() + # fix aspect ratio
#  xlab("Easting (m)") + # x-label
#  ylab("Northing (m)") # y-label

f <- log(rent) ~ 1 # formula for SRE model

S <- SRE(f = f, # formula
         data = list(query), # list of datasets
         BAUs = GridBAUs, # BAUs
         basis = G, # basis functions
         est_error = TRUE, # estimation measurement error
         average_in_BAU = TRUE) # average data over BAUs

S <- SRE.fit(SRE_model = S, # SRE model
             n_EM = 10, # max. no. of EM iterations
             tol = 0.01, # tolerance at which EM is assumed to have converged
             print_lik=TRUE) # print log-likelihood at each iteration

GridBAUsPred <- predict(S, obs_fs = FALSE) %>%
  st_as_sf(coord = c("lng", "lat")) %>%
  st_set_crs(4326)

GridBAUsPred <- st_join(GridBAUsPred, cbsa, left = FALSE)


#### C. Plot the Kriging Grid -------------------------------------------------

theme_krige <- function(...){
  theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(1, "inch"),
          legend.key.height = unit(.1, "inch"),
          axis.text = element_blank()) +
    xlab("") + ylab("") +
    coord_fixed()
}

#plot predicted values and standard errors
g1 <- ggplot() +
  geom_hex(data = GridBAUsPred, 
           aes(x = lng, y = lat, fill = exp(mu)),
           stat = "identity", color = NA) +
  scale_fill_distiller(palette = "Spectral",
                       name = paste0("Predicted\n", bed_size, "B Rent"),
                       labels = scales::dollar) +
  geom_density_2d(data = as_tibble(query@coords), 
                  aes(x = lng, y = lat),
                  colour = "grey80") + 
  theme_krige()

#Similar to above but with s.e.
g2 <- ggplot() + 
  geom_hex(data= GridBAUsPred,
            aes(x = lng, y = lat, fill = sqrt(var)),
            stat = "identity") +
  scale_fill_distiller(palette = "BrBG",
                       name = "SE",
                       guide = guide_legend(title="SE")) +
  theme_krige()

g1 + g2 +
  ggsave(filename = paste0("./output/maps/", cbsa_code, "_kriged.pdf"),
         width = 14, height = 9, dpi = 300)

