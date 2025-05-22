

setwd( R'(C:\Users\james\OneDrive\Desktop\Tunis trip)' )
data_dir = file.path(getwd(), "data")

merl = readRDS( file.path("data", "merl_TATBTC.AMMED.rds") )
pape = readRDS( file.path("data","pape_TATBTC.AMMED.rds") )
grid_1516 = readRDS( file.path("data", "grid_1516.rds") )

get_weights <- 
function( RDS_format ){
  match_row = match( unique(RDS_format$code), RDS_format$code )
  RDS_format[match_row,c("X", "Y", "year", "swept.area", "ptot", "tmp_bot", "SPECIE")]
}

DF = rbind(
  get_weights(merl), 
  get_weights(pape) 
)
