

setwd( R'(C:\Users\James.Thorson\Desktop\Work files\Meetings and Presentations\2025-05 -- Tunis)' )

# packages
library(tinyVAST)
library(sf)
library(fmesher)
library(googledrive)

# Read sampling data
# Data are on google drive here: https://drive.google.com/drive/folders/1EeaN-mJcWyjEyvQbzxk7bzlOELfGuuj5
data_dir = file.path(getwd(), "data")
if( !("merl_TATBTC.AMMED.rds" %in% list.files(data_dir)) ){
  Data_drive <- drive_get( "merl_TATBTC.AMMED.rds",
                           corpus = "allDrives" )
  drive_download( Data_drive$id,
                  path = file.path(data_dir,"merl_TATBTC.AMMED.rds") )
}
if( !("pape_TATBTC.AMMED.rds" %in% list.files(data_dir)) ){
  Data_drive <- drive_get( "pape_TATBTC.AMMED.rds",
                           corpus = "allDrives" )
  drive_download( Data_drive$id,
                  path = file.path(data_dir,"pape_TATBTC.AMMED.rds") )
}
if( !("grid_1516.rds" %in% list.files(data_dir)) ){
  Data_drive <- drive_get( "grid_1516.rds",
                           corpus = "allDrives" )
  drive_download( Data_drive$id,
                  path = file.path(data_dir,"grid_1516.rds") )
}

merl = readRDS( file.path("data", "merl_TATBTC.AMMED.rds") )
pape = readRDS( file.path("data","pape_TATBTC.AMMED.rds") )
grid_1516 = readRDS( file.path("data", "grid_1516.rds") )

# Hake
# for hake, we could cap the length bins at 60 cm, which would reduce the total to around 30.
# If that still feels like too many, we could consider increasing the bin width from 2 cm to 4 cm
# length_class in mm
table(merl$length_class)
# Doesn't have zeros
tapply( merl$nblon,
        INDEX = list( merl$year, merl$length_class ),
        FUN = length )

# Shrimp
# setting a maximum length of 3.8 cm and grouping into 2 mm bins
# length_class in mm
table(pape$length_class)
# Doesn't have zeros
tapply( pape$nblon,
        INDEX = list( pape$year, pape$length_class ),
        FUN = length )

#
plot( grid$X, grid$Y )
points(pape$X, pape$Y, col="blue")

library(ggplot2)
ggplot( pape ) +
  geom_point( aes( x=X, y=Y, col=month) ) +
  facet_wrap( ~ year )

#
merl$length_bin = cut(
  merl$length_class,
  breaks = c(0, seq( 60, 600, by = 40 ), Inf ),
  include.lowest = TRUE
)
table(merl$length_bin)

#
pape$length_bin = cut(
  pape$length_class,
  breaks = c(0, seq( 10, 38, by = 4 ), Inf ),
  include.lowest = TRUE
)
table(pape$length_bin)

#
tapply( pape$nblon,
        INDEX = list( pape$year, pape$length_bin ),
        FUN = sum )

#################
# QAQC and fix known errors
#################

# Show duplicates
which_rows = which( pape$code == "ITA_16_1999_11" & pape$length_class == 8 )
which_rows
pape[which_rows,]
# hist(table(code_length_class), breaks=seq(-0.5,10.5,by=1))

#
fix_duplicates <-
function( dat ){
  code_length_class = paste( dat$code, dat$length_class, sep="_" )
  keep_rows = match( unique(code_length_class), code_length_class )
  return( dat[keep_rows,] )
}

#
pape = fix_duplicates(pape)
merl = fix_duplicates(mer)

#
tapply( pape$nblon,
        INDEX = list( pape$year, pape$length_bin ),
        FUN = sum )

################
# Aggregate bins
################

combine_bins <-
function( dat ){
  total_catch = tapply( dat$nblon,
                        INDEX = list( code=dat$code, length_bin=dat$length_bin ),
                        FUN = sum )
  total_catch = ifelse( is.na(total_catch), 0, total_catch )

  dat2 = expand.grid( dimnames(total_catch) )
  dat2$nblon = as.vector(total_catch)

  #
  match_rows = match( dat2$code, dat$code )
  dat2 = cbind( dat2, dat[match_rows,c('X','Y','year','depth','tmp_bot','swept.area')] )
  if( any(is.na(dat2)) ) stop("Check")
  if( sum(dat2$nblon) != sum(dat$nblon) ) stop("Check")
  return(dat2)
}

pape = combine_bins( pape )
merl = combine_bins( merl )

#################
# Add UTM
#################

#
get_utm <-
function( dat, crs = 32633 ){
  dat_sf = st_as_sf( dat, coords = c("X","Y"), crs = st_crs(4326) )
  proj_sf = st_transform( dat_sf, crs = st_crs(crs) )
  coords = st_coordinates( proj_sf ) / 1000
  colnames(coords) = c("X_crs", "Y_crs")
  return(coords)
}

pape = cbind( pape, get_utm(pape) )
merl = cbind( merl, get_utm(merl) )

#################
# Rescale covariates
#################

# Rescale depth
pape$depth = pape$depth / 100
merl$depth = merl$depth / 100

# Rescale tmp_bot
pape$tmp_bot = (pape$tmp_bot - 10) / 10
merl$tmp_bot = (merl$tmp_bot - 10) / 10

# Get log-offset to match link
pape$log.swept.area = log( pape$swept.area )
merl$log.swept.area = log( merl$swept.area )


#################
# Rename length_bin levels
#################

rename_bins <-
function( length_bin ){
  length_bin = gsub( x = length_bin, pattern = ",", replace = "_", fixed=TRUE )
  length_bin = gsub( x = length_bin, pattern = "(", replace = "", fixed=TRUE )
  length_bin = gsub( x = length_bin, pattern = "]", replace = "", fixed=TRUE )
  length_bin = gsub( x = length_bin, pattern = "[", replace = "", fixed=TRUE )
  return( length_bin )
}

#
pape$length_bin = rename_bins( pape$length_bin )
merl$length_bin = rename_bins( merl$length_bin )

#################
# Exclude length_bin:year with all 0s
#################

remove_all_zeros <-
function( samples,
          dat = samples ){
  #
  prop_tc =
  tapply( samples$nblon,
          INDEX = list( year = samples$year, length_bin = samples$length_bin ),
          FUN = function(vec){round(mean(vec>0),3)} )
  prop_tc = ifelse( is.na(prop_tc), 0, prop_tc )
  if( any(colSums(prop_tc)==0) ) stop("Check_1")

  #
  prop_z = expand.grid(dimnames(prop_tc))
  prop_z$prop = as.vector(prop_tc)

  #
  which_intercepts = which(prop_z$prop == 0)
  if( length(which_intercepts) > 0 ){
    prop_z = prop_z[which_intercepts,,drop=FALSE]

    #
    which_rows = which( paste(dat$year,dat$length_bin) %in% paste(prop_z$year,prop_z$length_bin) )
    if( length(which_rows) > 0 ){
      dat = dat[-which_rows,]
    }
  }
  return(dat)
}

pape = remove_all_zeros( pape )
merl = remove_all_zeros( merl )

#################
# Final QAQC
#################

#



#################
# Fit model
#################

species = c("pape", "merl")[1]
covs = c("by_len", "const", "no")[1]

Date = Sys.Date()
  run_dir = file.path( getwd(), Date, paste0(species,"_covs=",covs) )
  dir.create(run_dir, recursive = TRUE)

if( species == "pape" ){
  # Data
  data_species = pape

  #
  space_term = "
    0_10 <-> 0_10, sd
    10_14 <-> 10_14, sd
    14_18 <-> 14_18, sd
    18_22 <-> 18_22, sd
    22_26 <-> 22_26, sd
    26_30 <-> 26_30, sd
    30_34 <-> 30_34, sd
    34_38 <-> 34_38, sd
    38_Inf <-> 38_Inf, sd
  "

  #
  spacetime_term = "
    0_10 <-> 0_10, 0, sd
    10_14 <-> 10_14, 0, sd
    14_18 <-> 14_18, 0, sd
    18_22 <-> 18_22, 0, sd
    22_26 <-> 22_26, 0, sd
    26_30 <-> 26_30, 0, sd
    30_34 <-> 30_34, 0, sd
    34_38 <-> 34_38, 0, sd
    38_Inf <-> 38_Inf, 0, sd

    0_10 -> 0_10, 1, ar
    10_14 -> 10_14, 1, ar
    14_18 -> 14_18, 1, ar
    18_22 -> 18_22, 1, ar
    22_26 -> 22_26, 1, ar
    26_30 -> 26_30, 1, ar
    30_34 -> 30_34, 1, ar
    34_38 -> 34_38, 1, ar
    38_Inf -> 38_Inf, 1, ar
  "

  # Define distribution for each category
  Family = list(
    "0_10" = poisson(),
    "10_14" = poisson(),
    "14_18" = poisson(),
    "18_22" = poisson(),
    "22_26" = poisson(),
    "26_30" = poisson(),
    "30_34" = poisson(),
    "34_38" = poisson(),
    "38_Inf" = poisson()
  )
}else{
  data_species = merl

  #
  space_term = "
    0_60 <-> 0_60, sd
    60_100 <-> 60_100, sd
    100_140 <-> 100_140, sd
    140_180 <-> 140_180, sd
    180_220 <-> 180_220, sd
    220_260 <-> 220_260, sd
    260_300 <-> 260_300, sd
    300_340 <-> 300_340, sd
    340_380 <-> 340_380, sd
    380_420 <-> 380_420, sd
    420_460 <-> 460_500, sd
    460_500 <-> 460_500, sd
    500_540 <-> 500_540, sd
    540_580 <-> 540_580, sd
    580_Inf <-> 580_Inf, sd
  "

  #
  spacetime_term = "
    0_60 <-> 0_60, 0, sd
    60_100 <-> 60_100, 0, sd
    100_140 <-> 100_140, 0, sd
    140_180 <-> 140_180, 0, sd
    180_220 <-> 180_220, 0, sd
    220_260 <-> 220_260, 0, sd
    260_300 <-> 260_300, 0, sd
    300_340 <-> 300_340, 0, sd
    340_380 <-> 340_380, 0, sd
    380_420 <-> 380_420, 0, sd
    420_460 <-> 460_500, 0, sd
    460_500 <-> 460_500, 0, sd
    500_540 <-> 500_540, 0, sd
    540_580 <-> 540_580, 0, sd
    580_Inf <-> 580_Inf, 0, sd

    0_60 -> 0_60, 1, ar
    60_100 -> 60_100, 1, ar
    100_140 -> 100_140, 1, ar
    140_180 -> 140_180, 1, ar
    180_220 -> 180_220, 1, ar
    220_260 -> 220_260, 1, ar
    260_300 -> 260_300, 1, ar
    300_340 -> 300_340, 1, ar
    340_380 -> 340_380, 1, ar
    380_420 -> 380_420, 1, ar
    420_460 -> 460_500, 1, ar
    460_500 -> 460_500, 1, ar
    500_540 -> 500_540, 1, ar
    540_580 -> 540_580, 1, ar
    580_Inf -> 580_Inf, 1, ar
  "

  # Define distribution for each category
  Family = list(
    "0_60" = poisson(),
    "60_100" = poisson(),
    "100_140" = poisson(),
    "140_180" = poisson(),
    "180_220" = poisson(),
    "220_260" = poisson(),
    "260_300" = poisson(),
    "300_340" = poisson(),
    "340_380" = poisson(),
    "380_420" = poisson(),
    "420_460" = poisson(),
    "460_500" = poisson(),
    "500_540" = poisson(),
    "540_580" = poisson(),
    "580_Inf" = poisson()
  )
}

# make SPDE mesh
mesh = fm_mesh_2d(
  loc = data_species[,c("X_crs","Y_crs")],
  cutoff = 20
)

# Add bin_num
data_species$bin_num = scale( as.numeric(as.factor(data_species$length_bin)), center=TRUE, scale = TRUE )

#
if( covs=="const" ){
  form = nblon ~ offset(log.swept.area) + 0 + interaction(year,length_bin) + s(tmp_bot) + s(depth)
}else if(covs == "no"){
  form = nblon ~ offset(log.swept.area) + 0 + interaction(year,length_bin)
}else if( covs == "by_len" ){
  #
  form = nblon ~ offset(log.swept.area) + 0 + interaction(year,length_bin) + s(tmp_bot) + s(depth) + ti(depth,bin_num)
  # Didn't converge
  form = nblon ~ offset(log.swept.area) + 0 + interaction(year,length_bin) + s(tmp_bot) + s(depth) + ti(tmp_bot,bin_num)
  # Full covariate:bin_num didn't work
  #form = nblon ~ offset(log.swept.area) + 0 + interaction(year,length_bin) + te(tmp_bot,bin_num) + te(depth,bin_num)
  # Marginal covariate + ti( covariate, bin_num) didn't converge
  #form = nblon ~ offset(log.swept.area) + 0 + interaction(year,length_bin) + s(tmp_bot) + s(depth) + ti(tmp_bot,bin_num) + ti(depth,bin_num)
}else{
  stop("Check")
}

# swept.area
fit = tinyVAST(
  # Specification
  data = data_species,
  formula = form,
  space_term = space_term,
  spacetime_term = spacetime_term,
  family = Family,

  # Indicators
  variable_column = "length_bin",
  time_column = "year",
  times = min(data_species$year):max(data_species$year),
  distribution_column = "length_bin",
  space_columns = c("X_crs", "Y_crs"),

  # Settings
  control = tinyVASTcontrol(
    trace = 1
  ),
  spatial_domain = mesh
)

#
saveRDS( fit,
      file = file.path(run_dir, "fit.RDS") )

##################
# Plots
##################
#library(visreg)
#visreg(fit, xvar="tmp_bot", what="p_g")

if( covs=="yes" ){
  # compute partial dependence plot
  library(pdp)

  Partial = partial( object = fit,
                     pred.var = "depth",
                     pred.fun = \(object,newdata) predict(object,newdata),
                     train = data_species,
                     approx = TRUE )
  png( file = file.path(run_dir, "depth.png"), width=4, height=4, res=200, units="in" )
    plotPartial( Partial )
  dev.off()

  Partial = partial( object = fit,
                     pred.var = "tmp_bot",
                     pred.fun = \(object,newdata) predict(object,newdata),
                     train = data_species,
                     approx = TRUE )
  png( file = file.path(run_dir, "tmp_bot.png"), width=4, height=4, res=200, units="in" )
    plotPartial( Partial )
  dev.off()
}

#
grid = as.data.frame(grid_1516)
grid$tmp_bot = (grid$tmp_bot - 10) / 10
grid$depth = grid$depth / 100

#
grid = subset( grid, year %in% unique(data_species$year) )
grid_sf = st_as_sf( x = grid, coords = c("X.utm","Y.utm"), crs = st_crs(32633) )
grid_sf = st_geometry(grid_sf)

#
for( bin_i in seq_along(unique(data_species$length_bin)) ){
  newdata = cbind( grid,
                   X_crs = grid$X.utm / 1000,
                   Y_crs = grid$Y.utm / 1000,
                   length_bin = unique(data_species$length_bin)[bin_i] )
  # Remove same zeros
  newdata = remove_all_zeros( data_species, newdata )
  newdata = cbind( newdata, log.swept.area = log(1) )

  #
  log_d = predict(fit, newdata = newdata, what = "p_g")
  plot_sf = st_sf( grid_sf, log_d = log_d, year = grid$year )

  #
  ggplot( plot_sf) +
    geom_sf( aes(col = log_d) ) +
    facet_wrap( vars(year) )
  ggsave( file = file.path(run_dir,paste0("logd_",unique(data_species$length_bin)[bin_i],".png")), width = 10, height = 10 )
}

# simulate new data conditional on fixed and random effects
y_ir = replicate( n = 100,
           expr = fit$obj$simulate()$y_i )
res = DHARMa::createDHARMa( simulatedResponse = y_ir,
                            observedResponse = data_species$nblon,
                            fittedPredictedResponse = fitted(fit) )
png( file.path(run_dir, "DHARMa.png"), width=5, height=5, res=200, units="in" )
  plot(res, form = data_species$length_bin)
dev.off()

# Get abundance
N_jz = expand.grid( length_bin = unique(data_species$length_bin), year=unique(data_species$year) )
N_jz = remove_all_zeros( data_species, N_jz )
N_jz = cbind( N_jz, "Biomass"=NA, "SE"=NA, "BiasCorr"=NA )
for( j in seq_len(nrow(N_jz)) ){
  message( "Integrating ", N_jz[j,'year'], " ", N_jz[j,'length_bin'], ": ", Sys.time() )
  if( is.na(N_jz[j,'Biomass']) ){
    grid0 = subset(grid, year==1999)
    newdata = cbind( grid0[,c('depth','tmp_bot')],
                     X_crs = grid0$X.utm / 1000,
                     Y_crs = grid0$Y.utm / 1000,
                     length_bin = N_jz[j,'length_bin'],
                     year = N_jz[j,'year'],
                     log.swept.area = log(1) )
    # Area-expansion
    index1 = integrate_output( fit,
                    area = rep(1,nrow(newdata)),
                    newdata = newdata,
                    apply.epsilon = TRUE,
                    bias.correct = FALSE,
                    intern = TRUE )
    N_jz[j,c('Biomass','SE','BiasCorr')] = index1[1:3] / 1e3
  }
}
N_ct = tapply( N_jz$BiasCorr,
               INDEX = list(length_bin = N_jz$length_bin, year = N_jz$year),
               FUN = sum )
N_ct = N_ct / outer( rep(1,nrow(N_ct)), colSums(N_ct,na.rm=TRUE) )
write.csv( N_ct, file = file.path(run_dir,"N_ct.csv") )

#
library(ggplot2)
long = cbind( expand.grid(dimnames(N_ct)), "p"=as.numeric(N_ct) )
ggplot( data=long  ) +
  facet_wrap( vars(length_bin), scales="free_y", ncol = 3, axis.labels = "margins" ) +
  geom_line( aes(x=year, y=p, group = length_bin) ) # + scale_y_log10()
ggsave( file = file.path(run_dir,"index_by_length.png"), width=7, height=7 )

library(ggplot2)
ggplot( data=long ) +
  facet_wrap( vars(year), scales="free_y", ncol = 3, axis.labels = "margins" ) +
  geom_line( aes(x=length_bin, y=p, group = year) ) # + scale_y_log10()
ggsave( file = file.path(run_dir,"comps_by_year.png"), width=7, height=7 )
