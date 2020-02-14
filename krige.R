#### Fixed Rank Kriging Craigslist data for 2019 ------------------------------

#dependencies
library(dplyr)
library(yaml)
library(DBI)
library(FRK)
library(sp)
library(sf)
library(purrr)
library(stringr)

#set wd to project base folder
setwd("H:/natrent-krige")

#store credentials at base dir of natrent-city-sub as YAML
cred <- read_yaml("./natrent0.yaml")

#create a natrent database connection
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
  select(CSAFP, CBSAFP, NAME, geometry) %>%
  rename(METRO = NAME)

#read place shapefile
places <- read_sf("./input/ACS Place Polygon/US_place_2017.shp") %>%
  select(NAME, geometry) %>%
  rename(PLACE = NAME)


#### Function to run Fixed Rank Kriging on listings for a given metro ---------

kriger <- function(metro_code, bed_size = 1,
                   start_date = "2019-01-01", end_date = "2019-12-31"){
  
  ## Define SQL query and submit to natrent
  
  #basic query
  query <- "SELECT DISTINCT d.listing_date AS date, d.clean_beds AS beds, d.clean_sqft AS sqft, 
                            d.clean_rent AS rent, c.state, c.county, c.met_name,
                            ROUND(CAST(ST_X(ST_TRANSFORM(d.geometry, 4326)) as numeric), 3) as lng, 
                            ROUND(CAST(ST_Y(ST_TRANSFORM(d.geometry, 4326)) as numeric), 3) as lat
            FROM (
                  SELECT a.geometry, a.cbsafp AS met_id, a.name AS met_name, b.statefp AS state, 
                         b.countyfp AS county
                  FROM cbsa17 a
                  JOIN county17 b ON a.cbsafp = b.cbsafp
                  WHERE a.cbsafp = ?cbsa_code AND b.cbsafp = ?cbsa_code
            ) c
            LEFT JOIN clean d ON ST_Contains(c.geometry, d.geometry) 
            WHERE d.listing_date BETWEEN ?start AND ?end AND
                  d.match_type NOT IN ('No Address Found', 'Google Maps Lat/Long') AND
                  d.clean_beds = ?beds AND d.clean_rent IS NOT NULL
            ORDER BY d.listing_date DESC"
  
  #work fn args into the query
  query <- sqlInterpolate(natrent, query, 
                          start = start_date,
                          end = end_date,
                          beds = bed_size,
                          cbsa_code = metro_code)
  
  #submit the query
  query <- dbGetQuery(natrent, query)
  
  #reduce to unique lat/lng via aggregation
  query <- query %>%
    group_by(lat, lng) %>%
    summarize(rent = median(rent),
              sqft = median(sqft),
              state = max(state),
              county = max(county),
              metro = max(met_name))
  
  #assign a vector to capture the metro shorthand for output naming
  metro_nickname <- unique(str_split_fixed(query$metro, pattern = "-", n = 2)[1])
  
  #initialize as Spatial object
  coordinates(query) <- ~ lng + lat
  
  #set the projection to WGS84 lat/lng
  proj4string(query) <- CRS("+init=epsg:4326")
  
  ## Estimate Spatial Random Effects Model for listing data
  
  #fine grained estimation of SRE
  set.seed(1)
  
  #define Basic Areal Units
  GridBAUs <- auto_BAUs(manifold = plane(),
                        type = "hex", # hex grid
                        data = query) # data around which to create BAUs
  
  GridBAUs$fs <- 1 # fine-scale variation at BAU level
  
  #define basis for kriging
  G <- auto_basis(manifold = plane(), 
                  data = query, # queried data
                  nres = 3, # number of resolutions
                  type = "Gaussian", # type of basis function
                  regular = 1) # place regularly in domain
  
  #show_basis(G) + # illustrate basis functions
  #  coord_fixed() + # fix aspect ratio
  #  xlab("Easting (m)") + # x-label
  #  ylab("Northing (m)") # y-label
  
  #formula for SRE model
  f <- rent ~ 1 
  
  #run the SRE model routine
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
  
  #create predicted values for grid of hex cells
  GridBAUsPred <- predict(S, obs_fs = FALSE) %>%
    st_as_sf(coord = c("lng", "lat")) %>%
    st_set_crs(4326) %>%
    st_transform(st_crs(cbsa))
  
  #trim to CBSA bounds, rename some columns
  GridBAUsPred <- st_join(GridBAUsPred, cbsa, left = FALSE) %>%
    filter(CBSAFP == metro_code) %>%
    rename(predicted = mu,
           variance = var,
           std_dev = sd) 
  
  if(!dir.exists(file.path(paste0("./output/shp/", metro_nickname)))){
    dir.create(file.path(paste0("./output/shp/", metro_nickname)))
  }
  
  #write shapefile of hexes to disk (no )
  st_write(GridBAUsPred, 
           file.path(paste0("./output/shp/", metro_nickname, "/", 
                            bed_size, "B.shp")),
           update = TRUE)
  
  #create map of kriged hexes
  source("map.R", local = TRUE)
  
  #save model variogram/summary statistics (to do)
}


#test a few
kriger("42660")
kriger("42660", bed_size = 2)
kriger("37980")
kriger("37980", bed_size = 2)

#### Kriging for top 100 metros across bedroom sizes --------------------------

#query for n by cbsa
metro_sql <- "SELECT a.cbsafp, a.name, count(b.geometry) as n
              FROM cbsa17 a
              INNER JOIN clean b ON ST_CONTAINS(a.geometry, b.geometry)
              WHERE b.listing_date BETWEEN '2019-12-01' AND '2019-12-31' AND
                    a.memi = '1'
              GROUP BY a.cbsafp, a.name
              ORDER BY n DESC"

#fetch the query
metros <- dbGetQuery(natrent, metro_sql)

#filter the query result
metros <- metros %>%
  top_n(100, n)

#specify bedroom sizes for the kriging
bed_sizes <- 0:3

#make a list out of these values to pass to purr::cross_df()
args <- list(cbsafp = metros$cbsafp, bed_size = bed_sizes) %>%
  cross_df() %>%
  arrange(cbsafp)

#join the counts
args <- left_join(args, metros)

#arrange from most to least
args <- arrange(args, desc(n))

#map the columns to kriger for running the function on each combination
map2(args$cbsafp, args$bed_size, kriger)

