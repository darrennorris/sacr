# to do
# 1) checkNAs 637 ....
#1) Use merge in function "dataF5". Depends on sorting of multiple tables...
#2) check why number of females at start (Year 0) is not correct??
#3) add zero values for "propKM"


# run demography  ("forrandom5.R") for 
# each of the river lengths ("RiverLengthSummary.R")
#1) load file,
rlb <- "C:\\Users\\Darren\\Documents\\2018 Unifilis demography\\analysis\\riverl.RDS"
riverl <- readRDS(rlb)
popb <- "C:\\Users\\Darren\\Documents\\2018 Unifilis demography\\analysis\\dfpop.RDS"
dfpop <- readRDS(popb)


#2) apply function to each row, 
# result is data frame with projections for each row of "dfpop"
# 2.1) list with data and matrix for projection
prepop <- function(x){
  vpop <- unlist(x[ ,4:19])
  class(vpop)
  tracaja <- matrix(vpop, byrow = TRUE, ncol=4)
  dimnames(tracaja) <- list(c("a", "b", "c", "d"),
                            c(1,2,3,4))
  
  #numeric vector of individuals at different age stage
  library(reshape2)
  dft <- reshape2::melt(riverl, id.vars = c('BASIN_N', 'subbasn', 
                                         'sub_area_km2', 'accessible'), 
                 measure.vars = c('tot_km', 'tot_notPA', 'tot_PA',
                                  'tot_Ind', 'tot_SP', 'tot_use'), 
                 value.name = c('distKMa'))
  
  myq <- function(x){
    myquant <- seq(0.05,1,by=0.05)
    dfout <- data.frame(propKM = myquant, distKM = x$distKMa * myquant)
    dfout
  }
  library(plyr)
  dft2 <- plyr::ddply(dft, .(BASIN_N, subbasn, sub_area_km2, 
                             accessible, variable), myq)
  adultF.d <- 10 # adult female density per river km
  #distkm <- 100 # length of river
  dft2$adultF.n <- trunc(adultF.d * dft2$distKM)
  
  l1 <- list(rdata = dft2, tracajam = tracaja)
  
}
library(plyr)
l.gpop <- dlply(dfpop, .(species, type, increase), prepop)
#l.gpop$`Podocnemis unifilis.headstart.0`$rdata
#class(l.gpop$`Podocnemis unifilis.headstart.0`$tracajam)

# 2.2) project list
proj.rivl <- function(x){
  #tracaja <- l.gpop$`Podocnemis unifilis.headstart.0`$tracajam
  #x <- l.gpop$`Podocnemis unifilis.headstart.0`$rdata
  tracaja <- x$tracajam
  
  doproj <- function(x) {
    tracaja_n <-  x$adultF.n * c(11.1, 4, 2, 1) 
    
    # project PPM 
    library(popdemo)
    library(popbio)
    pr_tracaja <- project(tracaja, vector=tracaja_n, time=50)
    
    # data for plotting
    len <- length(pr_tracaja)
    Time.intervals <- 0:(len - 1)
    eggs <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[1])))
    eju <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[2])))
    lju <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[3])))
    ad.fe <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[4])))
    plambda = popbio::lambda(tracaja)
    
    # make dataframe 
    dfout <- data.frame(lambda = plambda,
                        Years = Time.intervals, Individuals = pr_tracaja,
                        ss_egghatchling = popbio::stable.stage(tracaja)[1],
                        ss_earlyjuven = popbio::stable.stage(tracaja)[2],
                        ss_latejuven = popbio::stable.stage(tracaja)[3],
                        ss_adultfemale = popbio::stable.stage(tracaja)[4],
                        egghatch = eggs,
                        early_juven = eju,
                        late_juven = lju,
                        adult_females = ad.fe
    )
    fem0 <- dfout[(dfout$Years == 0), 'adult_females']
    dft <- data.frame(dfout, fem_t0 = fem0)
    dft$adult_female_diff <- round(((dft$adult_females - dft$fem_t0) / dft$fem_t0), 3)
    dft$change50_flag <- as.integer(ifelse(abs(dft$adult_female_diff) > 0.499, 1, 0))
    dft$double_flag <- as.integer(ifelse(dft$adult_female_diff > 0.999, 1, 0))
    dft
  }
  
  dfin <- x$rdata
  #dfin <- l.gpop$`Podocnemis unifilis.headstart.0`$rdata
  dout <- ddply(dfin, .(BASIN_N, subbasn, sub_area_km2, accessible,
                        variable, propKM, distKM), doproj)
  
  dout
}

memory.limit(84000)
dfres <- ldply(l.gpop, proj.rivl)

saveRDS(dfres, "dfpopres.RDS")

# 3)
# Table 1: summaries and basic checks
riverl <- readRDS("riverl.RDS")
dfpop <- readRDS("dfpop.RDS")
dfpop.res <- readRDS("dfpopres.RDS")


# Business as usual 
# final year totals
# Accessible populations with nest collection (first year survival 0.1) and
# adult harvest (10%).
selBA <- which(dfpop.res$Years == max(dfpop.res$Years) & 
                 dfpop.res$accessible == "Yes" &
                 dfpop.res$variable == "tot_km" &
                 dfpop.res$type =="headstart, female-hunt 10%" &
                 dfpop.res$increase == "0.1") # 53 rows (total for each subbasin)
# Inaccessible at base rates.
selBNA <- which(dfpop.res$Years == max(dfpop.res$Years) & 
                 dfpop.res$accessible == "No" &
                dfpop.res$variable == "tot_km" &
                 dfpop.res$type =="headstart" &
                 dfpop.res$increase == "0.2") # 53 rows (total for each subbasin)

dft.BAU <-merge(dfpop.res[selBA, ], dfpop.res[selBNA, ], 
            by = c("BASIN_NAME","subbasin"))
# summarise totals
library(plyr)
dfBAU <-  ddply(dft.BAU, .(BASIN_NAME, subbasin), summarize,
       BAU_km_tot = distKM.x + distKM.y, 
       BAU_km_accessible = distKM.x,
       BAU_km_inaccessible = distKM.y,
       BAU_lambda_accessible = lambda.x, 
       BAU_lambda_inaccessible = lambda.y,
       BAU_female_current = (distKM.x + distKM.y) * 10,
       BAU_female = round(adult_females.x + adult_females.y, 3),
       BAU_female_diff = (round(adult_females.x + adult_females.y, 3) - 
         ((distKM.x + distKM.y) * 10)), 
       BAU_female_change = ((round(adult_females.x + adult_females.y, 3) - 
                              ((distKM.x + distKM.y) * 10)) / ((distKM.x + distKM.y) * 10)) * 100

       ) 

# Strict protection, 
#Accessible with PAs set to base
selSPA.pa <- which(dfpop.res$Years == max(dfpop.res$Years) & 
                 dfpop.res$accessible == "Yes" &
                 dfpop.res$variable == "tot_PA" &
                 dfpop.res$type =="headstart" &
                 dfpop.res$increase == "0.2") # 53 rows (total for each subbasin)
#Accessible not PAs set to nest collection (first year survival 0.1) and
# adult harvest (10%).
selSPA.npa <- which(dfpop.res$Years == max(dfpop.res$Years) & 
                  dfpop.res$accessible == "Yes" &
                  dfpop.res$variable == "tot_notPA" &
                  dfpop.res$type =="headstart, female-hunt 10%" &
                  dfpop.res$increase == "0.1") # 53 rows (total for each subbasin)
# Inaccessible at base rates.
selSPNA <- which(dfpop.res$Years == max(dfpop.res$Years) & 
                  dfpop.res$accessible == "No" &
                  dfpop.res$variable == "tot_km" &
                  dfpop.res$type =="headstart" &
                  dfpop.res$increase == "0.2") # 53 rows (total for each subbasin)

dft.SP <- merge(dfpop.res[selSPA.pa, ], dfpop.res[selSPA.npa, ],
            by = c("BASIN_NAME","subbasin"))
dft.SP <-  merge(dft.SP, dfpop.res[selSPNA, ], 
                 by = c("BASIN_NAME","subbasin"))

dfSP <-  ddply(dft.SP, .(BASIN_NAME, subbasin), summarize,
                SP_km_tot = distKM.x + distKM.y + distKM, 
                SP_km_accessible = distKM.x + distKM.y,
                SP_km_inaccessible = distKM,
                SP_lambda_accessible = lambda.y, 
                SP_lambda_inaccessible = lambda,
                SP_female_current = (distKM.x + distKM.y + distKM) * 10,
                SP_female = round(adult_females.x + adult_females.y + adult_females, 3),
                SP_female_diff = (round(adult_females.x + adult_females.y + adult_females, 3) - 
                                     ((distKM.x + distKM.y + distKM) * 10)), 
                SP_female_change = ((round(adult_females.x + adult_females.y + adult_females, 3) - 
                                        ((distKM.x + distKM.y + distKM) * 10)) / ((distKM.x + distKM.y + distKM) * 10)) * 100
                
) 

# Community management
#Accessible with PAs set set to nest collection (first year survival 0.1) and
# adult harvest (10%).
selCMA.pa <- which(dfpop.res$Years == max(dfpop.res$Years) & 
                     dfpop.res$accessible == "Yes" &
                     dfpop.res$variable == "tot_PA" &
                     dfpop.res$type =="headstart, female-hunt 10%" &
                     dfpop.res$increase == "0.1") # 53 rows (total for each subbasin)

#Accessible not PAs to headstart 0.5 and harvest 10%
selCMA.npa <- which(dfpop.res$Years == max(dfpop.res$Years) & 
                      dfpop.res$accessible == "Yes" &
                      dfpop.res$variable == "tot_notPA" &
                      dfpop.res$type =="headstart, female-hunt 10%" &
                      dfpop.res$increase == "0.5") # 53 rows (total for each subbasin)
# Inaccessible at base rates.
selCMNA <- which(dfpop.res$Years == max(dfpop.res$Years) & 
                   dfpop.res$accessible == "No" &
                   dfpop.res$variable == "tot_km" &
                   dfpop.res$type =="headstart" &
                   dfpop.res$increase == "0.2") # 53 rows (total for each subbasin)

dft.CM <- merge(dfpop.res[selCMA.pa, ], dfpop.res[selCMA.npa, ],
                by = c("BASIN_NAME","subbasin"))
dft.CM <-  merge(dft.CM, dfpop.res[selCMNA, ], 
                 by = c("BASIN_NAME","subbasin"))

dfCM <-  ddply(dft.CM, .(BASIN_NAME, subbasin), summarize,
               CM_km_tot = distKM.x + distKM.y + distKM, 
               CM_km_accessible = distKM.x + distKM.y,
               CM_km_inaccessible = distKM,
               CM_lambda_accessible = lambda.x, 
               CM_lambda_inaccessible = lambda,
               CM_female_current = (distKM.x + distKM.y + distKM) * 10,
               CM_female = round(adult_females.x + adult_females.y + adult_females, 3),
               CM_female_diff = (round(adult_females.x + adult_females.y + adult_females, 3) - 
                                   ((distKM.x + distKM.y + distKM) * 10)), 
               CM_female_change = ((round(adult_females.x + adult_females.y + adult_females, 3) - 
                                      ((distKM.x + distKM.y + distKM) * 10)) / ((distKM.x + distKM.y + distKM) * 10)) * 100
               
) 

dfsum <- merge(dfBAU, dfSP, by = c("BASIN_NAME","subbasin"))
dfsum <- merge(dfsum, dfCM, by = c("BASIN_NAME","subbasin"))  
amsub <- c("Abacaxis", "Amazon floodplain", "Putumayo", "Japurá - Caquetá", "Javari",
           "Juruá", "Madeira", "Marañón", "Curuá-una", "Guama", "Jari", "Jutai",
           "Madeirinha", "Manacapuru", "Nanay", "Pacajá", "Piorini", "Tefe", 
           "Uatumá", "Napo", "Negro", "Purus", "Tapajós", "Tocantins",
           "Trombetas", "Ucayali", "Xingu")
dfsum$subbasinT <- c(amsub, dfsum$subbasin[28:53])
dfsum$BAU_flag_50 <- ifelse(dfsum$BAU_female_change < -49.999, 1,0)
dfsum$BAU_flag_30 <- ifelse(dfsum$BAU_female_change < -29.999, 1,0)
dfsum$SP_flag_50 <- ifelse(dfsum$SP_female_change < -49.999, 1,0)
dfsum$SP_flag_30 <- ifelse(dfsum$SP_female_change < -29.999, 1,0)
dfsum$CM_flag_50 <- ifelse(dfsum$CM_female_change < -49.999, 1,0)
dfsum$CM_flag_30 <- ifelse(dfsum$CM_female_change < -29.999, 1,0)
            
write.csv2(dfsum, "dfpopsum.csv")
saveRDS(dfsum, "dfpopsum.RDS")

library(htmlTable)
outc <- c('BASIN_NAME', 'subbasin', 'subbasinT', 'BAU_female', 
          'BAU_female_change', 'SP_female', 'SP_female_change', 
          'CM_female', 'CM_female_change')
dft4 <- dfsum[, outc]
si <- with(dft4, order(BASIN_NAME, subbasinT))
rv <- c(0,1,0,1,0,1)
htmlTable(txtRound(dft4[si, ],  
                   digits = rv, excl.cols = 1:3 ))
library(plyr)
ddply(dfsum, .(BASIN_NAME), summarize,
      BAU_femc = round(mean(na.omit(BAU_female_change)),1), 
      SP_femc = round(mean(na.omit(SP_female_change)),1),
      CM_femc = round(mean(na.omit(CM_female_change)),1), 
      BAU_fem = round(sum(na.omit(BAU_female)),0), 
      SP_fem = round(sum(na.omit(SP_female)),0),
      CM_fem = round(sum(na.omit(CM_female)),0)
)

round(sum(na.omit(dfsum$BAU_female)),0) 
round(sum(na.omit(dfsum$SP_female)),0)
round(sum(na.omit(dfsum$CM_female)),0)

# overall loss in adult females
current <- round(sum(na.omit(dfsum$BAU_female_current)),0)
diff <- round(sum(na.omit(dfsum$BAU_female)),0) - current
(diff / current) * 100

selLoss <- which(dfsum$BAU_female_change < 0)
# 4) Figure 4
# use sf with gplot2 pretty maps
library(sf)
library(dplyr)
#library(lwgeom)
memory.limit(84000)
sf1 <- read_sf("shapes/amazon_orinoco/amazon_orinoco.shp")
sf1clean <- st_union(st_buffer(sf1,0), by_feature = TRUE) %>% 
  group_by(BASIN_NAME, subbasin) %>% 
  summarise(area = min(Area_km2)) %>% st_cast("MULTIPOLYGON")
# add subbasin as per article table 
amsub <- c("Abacaxis", "Amazon floodplain", "Putumayo", "Japurá - Caquetá", "Javari",
"Juruá", "Madeira", "Marañón", "Curuá-una", "Guama", "Jari", "Jutai",
"Madeirinha", "Manacapuru", "Nanay", "Pacajá", "Piorini", "Tefe", 
"Uatumá", "Napo", "Negro", "Purus", "Tapajós", "Tocantins",
"Trombetas", "Ucayali", "Xingu")
sf1clean$subbasinT <- c(amsub, sf1clean$subbasin[28:53])

sf1clean$BASIN_FLAG <- as.numeric(as.factor(sf1clean$BASIN_NAME))
sf1clean$SUBBASIN_FLAG <- as.numeric(as.factor(sf1clean$subbasin))
sf1clean2 <- merge(sf1clean, dfsum, by = c("BASIN_NAME","subbasin"))
sf13395 <- st_transform(sf1clean2, crs=3395)
sf13395$area_km2 <- units::set_units(st_area(sf13395), km^2)
sum(sf13395$area_km2) #8859761 km^2
sf1clean2$area_km2<- sf13395$area_km2

# overall loss in adult females
current <- round(sum(na.omit(sf1clean2$BAU_female_current)),0)
diff <- round(sum(na.omit(sf1clean2$BAU_female)),0) - current
(diff / current) * 100

selLoss <- which(sf1clean2$BAU_female_change < 0)
sum(sf1clean2$area_km2[selLoss]) 
sum(sf1clean2$area_km2[selLoss]) / sum(sf1clean2$area_km2)

selLoss30 <- which(sf1clean2$BAU_female_change < -29.99)
sum(sf1clean2$area_km2[selLoss30]) 
sum(sf1clean2$area_km2[selLoss30]) / sum(sf1clean2$area_km2)

selLoss50 <- which(sf1clean2$BAU_female_change < -49.99)
sum(sf1clean2$area_km2[selLoss50]) 
sum(sf1clean2$area_km2[selLoss50]) / sum(sf1clean2$area_km2)

library(ggplot2)
# use scale_fill_gradientn for finer conttrol of fill
fbau <- ggplot(sf1clean2) +
  geom_sf(aes(fill = BAU_female_change)) +
  scale_fill_gradient2("% change") +
  ggtitle("A) Business as usual")

fsp <- ggplot(sf1clean2) +
  geom_sf(aes(fill = SP_female_change)) +
  scale_fill_gradient2("% change") +
  ggtitle("B) Strict protection")

fcm <- ggplot(sf1clean2) +
  geom_sf(aes(fill = CM_female_change)) +
  scale_fill_gradient2("% change") +
  ggtitle("C) Community management")
source("multiplot.R")
windows(width=14, height=5)
multiplot(fbau, fsp, fcm, cols=3)

# map of BAU loss for IUCN
# It's recommended to use a named vector
cols <- c("0" = "darkgrey", "1" = "red", "NA" = "darkgrey")
f50BAU <- ggplot(sf1clean2) +
  geom_sf(aes(fill = factor(BAU_flag_50))) +
  scale_fill_manual("",
    values = cols,
    breaks = c("0", "1", "NA"),
    labels = c("No", "Loss (50%)", "NA")
  ) +
  ggtitle("D) Subbasins with 50% loss")

f50SP <- ggplot(sf1clean2) +
  geom_sf(aes(fill = factor(SP_flag_50))) +
  scale_fill_manual("",
                    values = cols,
                    breaks = c("0", "1", "NA"),
                    labels = c("No", "Loss (50%)", "NA")
  ) +
  ggtitle("E) Subbasins with 50% loss")

f50CM <- ggplot(sf1clean2) +
  geom_sf(aes(fill = factor(CM_flag_50))) +
  scale_fill_manual("",
                    values = cols,
                    breaks = c("0", "1", "NA"),
                    labels = c("No", "Loss (50%)", "NA")
  ) +
  ggtitle("F) Subbasins with 50% loss")

source("multiplot.R")
windows(width=14, height=5)
multiplot(f50BAU, f50SP, f50CM, cols=3)

f30BAU <- ggplot(sf1clean2) +
  geom_sf(aes(fill = factor(BAU_flag_30))) +
  scale_fill_manual("",
                    values = cols,
                    breaks = c("0", "1", "NA"),
                    labels = c("No", "Loss (30%)", "NA")
  ) +
  ggtitle("H) Subbasins with 30% loss")

f30SP <- ggplot(sf1clean2) +
  geom_sf(aes(fill = factor(SP_flag_30))) +
  scale_fill_manual("",
                    values = cols,
                    breaks = c("0", "1", "NA"),
                    labels = c("No", "Loss (30%)", "NA")
  ) +
  ggtitle("I) Subbasins with 30% loss")

f30CM <- ggplot(sf1clean2) +
  geom_sf(aes(fill = factor(CM_flag_30))) +
  scale_fill_manual("",
                    values = cols,
                    breaks = c("0", "1", "NA"),
                    labels = c("No", "Loss (30%)", "NA")
  ) +
  ggtitle("J) Subbasins with 30% loss")
source("multiplot.R")
windows(width=14, height=5)
multiplot(f30BAU, f30SP, f30CM, cols=3)


#Fig 5 minimum lengths for pop increase?
# headstart only and max years 
dfpop.res <- readRDS("dfpopres.RDS")
selH <- which(dfpop.res$type == "headstart" & 
                dfpop.res$Years == max(dfpop.res$Years))
dfpop.resF5 <- dfpop.res[selH, ]

# acessible  propKM between 0.05 and .95
# subbasin * levels propKm * levels variable
53 * 19 * 6 # 6042 rows
selAc <- which(dfpop.resF5$accessible == "Yes" & 
                 dfpop.resF5$propKM < 1)
# subbasin * levels increase * levels variable
with(dfpop.resF5[selAc, ], table(propKM)) # 3180 (53 * 10 * 6)

53 * 19 * 10 * 6 # 60420 rows
# sort by basin, subbasin, increase, propKM, variable
library(dplyr)
df5as <- dplyr::arrange(dfpop.resF5[selAc, ], 
                        increase,  propKM, BASIN_N, subbasn, variable)
# add inacessible
# inaccessible propKM = 1 and at base rates
selIn <- which(dfpop.resF5$accessible == "No" & 
                 dfpop.resF5$propKM == 1 &
                 dfpop.resF5$increase == "0.2")
# should be all 1
with(dfpop.resF5[selIn, ], table(propKM)) # 318 (53 * 6)
library(dplyr)
df5In <- dplyr::arrange(dfpop.resF5[selIn, ], 
                        BASIN_N, subbasn, increase, variable)
df5In[which(df5In$subbasn=="Abacaxis"), ]
df5as$subbasn_In <- rep(df5In$subbasn, (60420/318))
df5as$variable_In <- rep(df5In$variable, (60420/318))
df5as$In_prop <- rep(df5In$propKM, (60420/318))
df5as$distKM_In <- rep(df5In$distKM, (60420/318))
df5as$adult_females_baseIn <- rep(df5In$adult_females, (60420/318))
# test 
selB <- which(df5as$subbasn != df5as$subbasn_In)
# add missing proportion
# base
selAcBase <- which(dfpop.resF5$accessible == "Yes" & 
                 dfpop.resF5$propKM < 1 & 
                   dfpop.resF5$increase == "0.2" )
# BAU
selAcBaseBAU <- which(dfpop.res$type == "headstart, female-hunt 10%" & 
                        dfpop.res$Years == max(dfpop.res$Years) &
                        dfpop.res$accessible == "Yes" & 
                     dfpop.res$propKM < 1 & 
                     dfpop.res$increase == "0.1" )
table(dfpop.resF5$type)
# order right?
df5asBaseBAU <- dplyr::arrange(dfpop.res[selAcBaseBAU, ], 
                  desc(propKM), BASIN_N, subbasn, increase,   variable)
df5as$subbasn_baseAcBAU <- rep(df5asBaseBAU$subbasn, 10)
df5as$variable_baseAcBAU <- rep(df5asBaseBAU$variable, 10)
df5as$propKM_baseACBAU <- rep(df5asBaseBAU$propKM, 10)
df5as$distKM_baseACBAU <- rep(df5asBaseBAU$distKM, 10)
df5as$adult_females_baseACBAU <- rep(df5asBaseBAU$adult_females, 10)
df5as$acc_propBAU <- df5as$propKM_baseACBAU + df5as$propKM
table(df5as$acc_propBAU) # should be all 1

# test 
selB <- which(df5as$subbasn != df5as$subbasn_baseAC)

# integers for comparison
df5as$adult_females_total <- trunc(df5as$adult_females_baseIn + df5as$adult_females_baseAC + df5as$adult_females)
df5as$adult_females_current <- trunc((df5as$distKM + df5as$distKM_In + df5as$distKM_baseAC) * 10)
df5as$prop_change <- (df5as$adult_females_total - df5as$adult_females_current) / df5as$adult_females_current
df5as$prop_change_clean <- ifelse(df5as$prop_change > 1, 1, df5as$prop_change)

# only accessible
df5as$adult_females_total_acc <- trunc(df5as$adult_females_baseAC + df5as$adult_females)
df5as$adult_females_current_acc <- trunc((df5as$distKM + df5as$distKM_baseAC) * 10)
df5as$prop_change_acc <- (df5as$adult_females_total_acc - df5as$adult_females_current_acc) / df5as$adult_females_current_acc
df5as$prop_change_acc_clean <- ifelse(df5as$prop_change_acc > 1, 1, df5as$prop_change_acc)

# only accessible BAU
df5as$adult_females_total_accBAU <- trunc(df5as$adult_females_baseACBAU + df5as$adult_females)
df5as$adult_females_current_accBAU <- trunc((df5as$distKM + df5as$distKM_baseACBAU) * 10)
df5as$prop_change_accBAU <- (df5as$adult_females_total_accBAU - df5as$adult_females_current_acc) / df5as$adult_females_current_acc
df5as$prop_change_accBAU_clean <- ifelse(df5as$prop_change_accBAU > 1, 1, df5as$prop_change_accBAU)

summary(df5as$prop_change)
summary(df5as$prop_change_clean)
summary(df5as$prop_change_acc)
summary(df5as$prop_change_acc_clean)
summary(df5as$prop_change_accBAU_clean)



# now function for all hunting levels
library(plyr)
library(dplyr)
gc()
# make tables
dfpop.res <- readRDS("dfpopres.RDS")
selH50 <- which(dfpop.res$Years == max(dfpop.res$Years))
dfpop.res50 <- dfpop.res[selH50, ]

# add inacessible
# inaccessible propKM = 1 and at base rates
selIn <- which(dfpop.res50$type == "headstart" & 
                 dfpop.res50$accessible == "No" & 
                 dfpop.res50$propKM == 1 &
                 dfpop.res50$increase == "0.2")

df5In <- dplyr::arrange(dfpop.res50[selIn, ], 
                        BASIN_N, subbasn, increase, variable)

# add missing proportion
# base
selAcBase <- which(dfpop.res50$type == "headstart" & 
                     dfpop.res50$accessible == "Yes" & 
                     dfpop.res50$propKM < 1 & 
                     dfpop.res50$increase == "0.2" )
df5asBase <- dplyr::arrange(dfpop.res50[selAcBase, ], 
                               desc(propKM), BASIN_N, subbasn, increase,   variable)
df5asBase$propKM_join <- 1 - df5asBase$propKM
# BAU
selAcBaseBAU <- which(dfpop.res$type == "headstart, female-hunt 10%" & 
                        dfpop.res$Years == max(dfpop.res$Years) &
                        dfpop.res$accessible == "Yes" & 
                        dfpop.res$propKM < 1 & 
                        dfpop.res$increase == "0.1" )
# order right?
df5asBaseBAU <- dplyr::arrange(dfpop.res[selAcBaseBAU, ], 
                               desc(propKM), BASIN_N, subbasn, increase,   variable)
df5asBaseBAU$propKM_join <- 1 - df5asBaseBAU$propKM

dataF5 <- function(x) {
# acessible  propKM between 0.05 and .95
selAc <- which(x$accessible == "Yes" & x$propKM < 1)
dft <- x[selAc, ]
# sort by basin, subbasin, increase, propKM, variable
library(dplyr)
df5as <- dplyr::arrange(dft, 
                        increase, propKM, BASIN_N, subbasn, variable)


#Add missing levels
df5as$subbasn_In <- rep(df5In$subbasn, (60420/318))
df5as$variable_In <- rep(df5In$variable, (60420/318))
df5as$In_prop <- rep(df5In$propKM, (60420/318))
df5as$distKM_In <- rep(df5In$distKM, (60420/318))
df5as$adult_females_baseIn <- rep(df5In$adult_females, (60420/318))

df5as$subbasn_baseAC <- rep(df5asBase$subbasn, 10)
df5as$variable_baseAC <- rep(df5asBase$variable, 10)
df5as$propKM_baseAC <- rep(df5asBase$propKM, 10)
df5as$distKM_baseAC <- rep(df5asBase$distKM, 10)
df5as$adult_females_baseAC <- rep(df5asBase$adult_females, 10)
df5as$acc_propbaseAC <- df5as$propKM_baseAC + df5as$propKM

df5as$subbasn_baseACBAU <- rep(df5asBaseBAU$subbasn, 10)
df5as$variable_baseACBAU <- rep(df5asBaseBAU$variable, 10)
df5as$propKM_baseACBAU <- rep(df5asBaseBAU$propKM, 10)
df5as$distKM_baseACBAU <- rep(df5asBaseBAU$distKM, 10)
df5as$adult_females_baseACBAU <- rep(df5asBaseBAU$adult_females, 10)
df5as$acc_propBAU <- df5as$propKM_baseACBAU + df5as$propKM

# integers for comparison
df5as$adult_females_total <- trunc(df5as$adult_females_baseIn + 
                                     df5as$adult_females_baseAC + df5as$adult_females)
df5as$adult_females_current <- trunc((df5as$distKM + 
                                        df5as$distKM_In + df5as$distKM_baseAC) * 10)
df5as$prop_change <- ((df5as$adult_females_total - df5as$adult_females_current) / 
                              df5as$adult_females_current)
df5as$prop_change_clean <- ifelse(df5as$prop_change > 1, 1, df5as$prop_change)

# only accessible
df5as$adult_females_total_acc <- trunc(df5as$adult_females_baseAC + df5as$adult_females)
df5as$adult_females_current_acc <- trunc((df5as$distKM + df5as$distKM_baseAC) * 10)
df5as$prop_change_acc <- (df5as$adult_females_total_acc - df5as$adult_females_current_acc) / df5as$adult_females_current_acc
df5as$prop_change_acc_clean <- ifelse(df5as$prop_change_acc > 1, 1, df5as$prop_change_acc)

# only accessible BAU
df5as$adult_females_total_accBAU <- trunc(df5as$adult_females_baseACBAU + df5as$adult_females)
df5as$adult_females_current_accBAU <- trunc((df5as$distKM + df5as$distKM_baseACBAU) * 10)
df5as$prop_change_accBAU <- (df5as$adult_females_total_accBAU - df5as$adult_females_current_acc) / df5as$adult_females_current_acc
df5as$prop_change_accBAU_clean <- ifelse(df5as$prop_change_accBAU > 1, 1, df5as$prop_change_accBAU)

df5as
}

library(plyr)
 df5hs <- plyr::ddply(dfpop.res50, .(species, type), .fun = dataF5)
 table(df5hs$type) # 60420
names(df5hs)
levels(df5hs$variable)
levels(df5hs$variable) <- c("All", "Not protected", "Protected", 
                        "Indigenous", "Strict", "Use")
levels(df5hs$type)
levels(df5hs$type) <- c("No hunt", "Hunt 2.5%", "Hunt 10%",
                            "Hunt 25%", "Hunt 50%")

write.csv2(dft1, "basins.csv")
write.csv2(df5hs, "popchangeF5.csv")
# checking
#examples where pop increased with current at 0?
selE <- which(df5hs$adult_females_current == 0 & 
                df5hs$adult_females_total > 0) #none
# NAs when there are currently females?
# If none then Nas may represent 0/0 = nothing to project and no change
selNA <- which(is.na(df5hs$prop_change)==TRUE & 
                 df5hs$adult_females_current > 0) #none
# examples with large values
class(df5hs)
selL <- which(df5hs$prop_change > 200)
sam <- sample(selL, 10)
df5hs[sam, ]

memory.limit(84000)
# Fig 2 plot = How much?
library(ggplot2)
table(df5hs$variable)
mycol <- c("#FF00FF", "#CC33CC", "#FF00CC", 
           "#FFFF33", "#FF9933", "#CC6600", "#993300", 
           "#00FF00", "#339900", "#336600")
selKM <- which(df5hs$variable %in% c("All", "Not protected", "Protected"))
selKM1 <- which(df5hs$variable %in% c("All"))
# 173 mm (7inch) width, http://onlinelibrary.wiley.com/journal/10.1111/(ISSN)1755-263X/homepage/ForAuthors.html

saveRDS(df5hs[selKM1, ], "Fig2dat.RDS") 
getf <- "C:\\Users\\Darren\\Documents\\2018 Unifilis demography\\analysis\\Fig2dat.RDS"
f2d <- readRDS(getf)
library(ggplot2)
pdf("fig5.pdf", width= 7, height = 3.5, useDingbats = FALSE)
ggplot(f2d, aes(propKM, prop_change_clean, color = increase)) +
  geom_hline(yintercept = 0) +
  geom_jitter(width = 0.1, height = 0.1, alpha=0.1) +
  stat_smooth(se=FALSE) +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  #facet_grid(variable~type, labeller = label_wrap_gen(width=19)) +
  facet_wrap(~type, nrow=1, labeller = label_wrap_gen(width=19)) +
  ylab("Relative population change\n(50 year projection)") + 
  xlab("Scenario cover (Proportion of subbasin river length)") + 
  scale_color_manual(name="Hatchling\nGraduation", values = mycol)

dev.off()

png("fig5.png", type="cairo", 
    width= 7, height = 3.5, units="in", pointsize=12, res=300)
ggplot(df5hs[selKM1, ], aes(propKM, prop_change_clean, color = increase)) +
  geom_hline(yintercept = 0) +
  geom_jitter(width = 0.1, height = 0.1, alpha=0.1) +
  stat_smooth(se=FALSE) +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  #facet_grid(variable~type, labeller = label_wrap_gen(width=19)) +
  facet_wrap(~type, nrow=1, labeller = label_wrap_gen(width=19)) +
  ylab("Relative population change\n(50 year projection)") + 
  xlab("Scenario cover (Proportion of subbasin river length)") + 
  scale_color_manual(name="Hatchling\nGraduation", values = mycol)
dev.off()


## Now how long ?
memory.limit(84000)
dfpop.res <- readRDS("dfpopres.RDS")
table(dfpop.res$lambda) # 0.465 to  1.1635
names(dfpop.res)
sel50All <- which(dfpop.res$change50_flag==1 & 
                    dfpop.res$accessible=="Yes",
                    dfpop.res$variable=="tot_km")
mycoln <- c("type" , "increase" , "BASIN_N" ,"subbasn", 
            "accessible", "propKM",  "distKM", "lambda" )           
      
dfpop.resY <- plyr::ddply(dfpop.res[sel50All, ], (mycoln), 
                          summarise, 
                          Years_50 = min(Years))
dfpop.resY$Years_50c <- ifelse(as.numeric(dfpop.resY$lambda) < 1, 
                               dfpop.resY$Years_50 * -1, 
                               dfpop.resY$Years_50)
levels(dfpop.resY$type)
levels(dfpop.resY$type) <- c("No hunt", "Hunt 2.5%", "Hunt 10%",
                        "Hunt 25%", "Hunt 50%")
selP1 <- which(dfpop.resY$propKM==1)
library(ggplot2)
ggplot(dfpop.resY[selP1, ], aes(x = increase, y = Years_50c, color=lambda)) +
  geom_jitter(alpha=0.3) +
  scale_color_gradientn("lambda", 
                       colours = c("darkred","tomato1", 
                                   "lightblue","darkblue"), 
                       values = c(0, 0.76,0.77, 1)) +
  scale_y_continuous("Years to change", limits = c(-50, 50), 
                     labels = c("50", "25", "0", "25", "50")) +
  scale_x_discrete("Hatchling graduation", breaks = c(0, 0.3, 0.6, 0.9)) +
  facet_wrap(~type, nrow=1, labeller = label_wrap_gen(width=19)) 
 
# 15/12/2018 accessibility analysis for sensitivity
#load("~/2018 Unifilis demography/analysis/testl2.RData")
load("~/2018 Unifilis demography/analysis/acess_sens.RData")
#5) Get results from 3 scenarios
library(cmartr)
source("C:/Users/Darren/Documents/2018 Unifilis demography/analysis/cmartr/R/PopScenAcess.R")
lscen <- PopScenAcess(dflup)
table(lscen$dft.BAU$namekey)  
library(ggplot2)
sel0 <- which(lscen$dft.BAU$dist_km.x > 0)
bau_clean <- lscen$dft.BAU[sel0, ]

#which(bau_clean$variable=="tot_notPA" & bau_clean$prop_km.x >0.75 & 
#        bau_clean$adult_female_diff > 0.5) # 2182 2183
#bau_clean[c(2182, 2183), ]
bau_clean$flag_use <- 1
bau_clean[which((bau_clean$prop_km.x > bau_clean$prop_km.y) & 
                  (bau_clean$dist_km.x < bau_clean$dist_km.y)), 'flag_use'] <- NA
bau_clean[which((bau_clean$prop_km.x < bau_clean$prop_km.y) & 
                  (bau_clean$dist_km.x > bau_clean$dist_km.y)), 'flag_use'] <- NA
bau_clean <- bau_clean[which(bau_clean$flag_use==1), ]
bau_clean <- droplevels(bau_clean)
levels(bau_clean$variable) <- c("Not protected", "Protected")

# see why still have increase in some with > 50% aceessible
ggplot(bau_clean, aes(prop_km.x, adult_female_diff, color = variable)) +
  geom_point(alpha=0.3) + stat_smooth() + 
  ylab("Relative population change\n(50 year projection)") + 
  xlab("Accessible (Proportion of catchment river length)") +
  coord_cartesian(ylim = c(-1, 2)) +
  theme(legend.title = element_blank(), 
        legend.justification = c(1,1),
        legend.position = c(0.9,0.9))

#SP
sel0sp <- which(lscen$dft.SP$dist_km.x > 0)
SP_clean <- lscen$dft.SP[sel0sp, ]

#which(SP_clean$variable=="tot_notPA" & SP_clean$prop_km.x >0.75 & 
#        SP_clean$adult_female_diff > 0.5) # 2182 2183
#SP_clean[c(2182, 2183), ]
SP_clean$flag_use <- 1
SP_clean[which((SP_clean$prop_km.x > SP_clean$prop_km.y) & 
                  (SP_clean$dist_km.x < SP_clean$dist_km.y)), 'flag_use'] <- NA
SP_clean[which((SP_clean$prop_km.x < SP_clean$prop_km.y) & 
                  (SP_clean$dist_km.x > SP_clean$dist_km.y)), 'flag_use'] <- NA
SP_clean <- SP_clean[which(SP_clean$flag_use==1), ]
SP_clean <- droplevels(SP_clean)
levels(SP_clean$variable) <- c("Not protected", "Protected")

ggplot(SP_clean, aes(prop_km.x, adult_female_diff, color = variable)) +
  geom_point(alpha=0.3) + stat_smooth() + 
  ylab("Relative population change\n(50 year projection)") + 
  xlab("Accessible (Proportion of catchment river length)") +
  coord_cartesian(ylim = c(-1, 2)) +
  theme(legend.title = element_blank(), 
        legend.justification = c(0.1,0),
        legend.position = c(0.1,0.1))

# community
sel0CM <- which(lscen$dft.CM$dist_km.x > 0)
CM_clean <- lscen$dft.CM[sel0CM, ]
CM_clean$flag_use <- 1
CM_clean[which((CM_clean$prop_km.x > CM_clean$prop_km.y) & 
                 (CM_clean$dist_km.x < CM_clean$dist_km.y)), 'flag_use'] <- NA
CM_clean[which((CM_clean$prop_km.x < CM_clean$prop_km.y) & 
                 (CM_clean$dist_km.x > CM_clean$dist_km.y)), 'flag_use'] <- NA
CM_clean <- CM_clean[which(CM_clean$flag_use==1), ]
CM_clean <- droplevels(CM_clean)
levels(CM_clean$variable) <- c("Not protected", "Protected")
summary(CM_clean$adult_female_diff) # great increase
CM_clean$afd <- CM_clean$adult_female_diff
CM_clean$afd <- ifelse(CM_clean$variable == "Not protected", CM_clean$adult_female_diff -1.4, 
                       CM_clean$adult_female_diff)  
CM_clean$afd <- ifelse(CM_clean$variable == "Not protected" & CM_clean$adult_female_diff > 2, 
                       CM_clean$adult_female_diff * 0.3, 
                       CM_clean$adult_female_diff)  

#CM_clean[which(CM_clean$variable == "Not protected" & CM_clean$adult_female_diff > 2), 'afd'] <- 2

ggplot(CM_clean, aes(prop_km.x, afd, color = variable)) +
  geom_point(alpha=0.3) + stat_smooth() + 
  ylab("Relative population change\n(50 year projection)") + 
  xlab("Accessible (Proportion of catchment river length)") +
  coord_cartesian(ylim = c(-1, 2)) +
  theme(legend.title = element_blank(), 
        legend.justification = c(0.1,0),
        legend.position = c(0.1,0.1))

# 26/12/2018
# lambda values for each 
source("C:/Users/Darren/Documents/2018 Unifilis demography/analysis/cmartr/R/PopScenLambda.R")
library(plyr)
#21:11-21:28
dfl <- ddply(dflup, .(akey), .fun = PopScenLambda)
min(dfl$lambda_min); max(dfl$lambda_max)


#27/12/2018
# general point based estimates
# create look up tables with projections
dfpop <- cmartr::PopParam(species = "Podocnemis unifilis", make_rds = FALSE)
PopProjF <- function(x, fed = 10){
vpop <- unlist(x[ ,4:19])
tracaja <- matrix(vpop, byrow = TRUE, ncol=4)
dimnames(tracaja) <- list(c("a", "b", "c", "d"),
                          c(1,2,3,4))
adultF.d <- fed # adult female density per river km
dist_km <- 1
adultF.n <- trunc(adultF.d * dist_km)
tracaja_n <-  adultF.n * c(11.1, 4, 2, 1) 

# project PPM 
pr_tracaja <- popdemo::project(tracaja, vector=tracaja_n, time=50)
# data for plotting
len <- length(pr_tracaja)
Time.intervals <- 0:(len - 1)
eggs <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[1])))
eju <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[2])))
lju <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[3])))
ad.fe <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[4])))
plambda = popbio::lambda(tracaja)

# make dataframe 
dfout <- data.frame(lambda = plambda,
                    ayear = Time.intervals, 
                    individuals = as.integer(trunc(pr_tracaja)),
                    ss_egghatchling = round(as.numeric(popbio::stable.stage(tracaja)[1]),3),
                    ss_earlyjuven = round(as.numeric(popbio::stable.stage(tracaja)[2]),3),
                    ss_latejuven = round(as.numeric(popbio::stable.stage(tracaja)[3]),3),
                    ss_adultfemale = round(as.numeric(popbio::stable.stage(tracaja)[4]),3),
                    egghatch = eggs,
                    early_juven = eju,
                    late_juven = lju,
                    adult_females = ad.fe
)
fem0 <- adultF.n
dft <- data.frame(species = x$species, hunt = x$type, 
                  increase = x$increase, dfout, fem_t0 = fem0)
dft$adult_female_diff <- round(((dft$adult_females - dft$fem_t0) / dft$fem_t0), 3)
dft$change50_flag <- as.integer(ifelse(abs(dft$adult_female_diff) > 0.499, 1, 0))
dft$double_flag <- as.integer(ifelse(dft$adult_female_diff > 0.999, 1, 0))
dft
}
library(plyr)
dff <- ddply(dfpop, .(akey), .fun = PopProjF)
# limit final to 10 * original
sel50 <- which(dff$ayear == max(dff$ayear))
dffinal <- dff[sel50, ]
dffinal$adult_females_clean <- ifelse(dffinal$adult_females > (10*dffinal$fem_t0), 
                                      10*dffinal$fem_t0, dffinal$adult_females)

dffinal$adult_female_diff_clean <- ifelse(dffinal$adult_females > (10*dffinal$fem_t0), 
                                      9, dffinal$adult_female_diff)
# now acessibility
# from resTab.R
rp <- system.file("shape/shapes_rivers3395", package="cmartr")
myfiles <- file.info(list.files(rp, pattern = "\\.shp$", full.names =TRUE))

# 215866
library(sf)
sf.rivp <- rbind(read_sf(rownames(myfiles)[1]), 
                 read_sf(rownames(myfiles)[2]),
                 read_sf(rownames(myfiles)[3]),
                 read_sf(rownames(myfiles)[4]),
                 read_sf(rownames(myfiles)[5]),
                 read_sf(rownames(myfiles)[6]),
                 read_sf(rownames(myfiles)[7]),
                 read_sf(rownames(myfiles)[8]),
                 read_sf(rownames(myfiles)[9]),
                 read_sf(rownames(myfiles)[10]),
                 read_sf(rownames(myfiles)[11]),
                 read_sf(rownames(myfiles)[12]),
                 read_sf(rownames(myfiles)[13]),
                 read_sf(rownames(myfiles)[14]),
                 read_sf(rownames(myfiles)[15]),
                 read_sf(rownames(myfiles)[16]),
                 read_sf(rownames(myfiles)[17]),
                 read_sf(rownames(myfiles)[18]),
                 read_sf(rownames(myfiles)[19]),
                 read_sf(rownames(myfiles)[20]),
                 read_sf(rownames(myfiles)[21]),
                 read_sf(rownames(myfiles)[22]),
                 read_sf(rownames(myfiles)[23]),
                 read_sf(rownames(myfiles)[24]),
                 read_sf(rownames(myfiles)[25]),
                 read_sf(rownames(myfiles)[26]),
                 read_sf(rownames(myfiles)[27]),
                 read_sf(rownames(myfiles)[28]),
                 read_sf(rownames(myfiles)[29]),
                 read_sf(rownames(myfiles)[30]),
                 read_sf(rownames(myfiles)[31]),
                 read_sf(rownames(myfiles)[32]),
                 read_sf(rownames(myfiles)[33]),
                 read_sf(rownames(myfiles)[34]),
                 read_sf(rownames(myfiles)[35]),
                 read_sf(rownames(myfiles)[36]),
                 read_sf(rownames(myfiles)[37]),
                 read_sf(rownames(myfiles)[38]),
                 read_sf(rownames(myfiles)[39]),
                 read_sf(rownames(myfiles)[40]),
                 read_sf(rownames(myfiles)[41]),
                 read_sf(rownames(myfiles)[42]),
                 read_sf(rownames(myfiles)[43]),
                 read_sf(rownames(myfiles)[44]),
                 read_sf(rownames(myfiles)[45]),
                 read_sf(rownames(myfiles)[46]),
                 read_sf(rownames(myfiles)[47]),
                 read_sf(rownames(myfiles)[48]),
                 read_sf(rownames(myfiles)[49]),
                 read_sf(rownames(myfiles)[50]),
                 read_sf(rownames(myfiles)[51]),
                 read_sf(rownames(myfiles)[52])
)

# writing as a single layer gpkg does not work
sf.rivp[3030, ]
selAp <- which(sf.rivp$SUBBASI %in% c(2, 31, 4))
fp <- "C:/Users/Darren/Documents/2018 Unifilis demography/analysis/cmartr/inst/other"
st_write(sf.rivp[selAp, ], dsn = file.path(fp, "rcall.shp"), delete_dsn = TRUE)
spr <- as(sf.rivp, "Spatial")
outgpg <- "C:/Users/Darren/Documents/2018 Unifilis demography/analysis/cmartr/inst/other/rcall.gpkg"
library(rgdal)
writeOGR(spr, outgpg, driver="GPKG", layer = 'riverpoints', delete_dsn = TRUE)

table(sf.rivp$urban)
ddply(sf.rivp, .(BASIN_N, subbasn), summarise, 
      bid = min(na.omit(SUBBASI)) 
)

# write as multilayer geopackage. 28.8 MB
mgt <- function(x){
  fp <- "C:/Users/Darren/Documents/2018 Unifilis demography/analysis/cmartr/inst/other"
   sf1 <- read_sf(rownames(x))
  layid <- min(sf1$SUBBASI)
  layn <- paste("layer", layid, sep="")
  st_write(sf1, dsn = file.path(fp, "riversec.gpkg"), layer = layn) 
}
plyr::a_ply(myfiles, .margins = 1, .fun = mgt)

library(sf)
inshp <- "C:/Users/Darren/Documents/2018 Unifilis demography/analysis/shapes/apriver.shp"
s2 <- read_sf(inshp)
s2$fa_t0 <- 10
plot(s2["All"])

# intersect with Amapa to get ottobacias and municipios
# ottobacias cover more so intersect with Amapa municipios, then points
inshp2 <- "C:/Users/Darren/Documents/2018 Unifilis demography/analysis/ap_municipios.shp"
s3 <- read_sf(inshp2)
inshp3 <- "C:/Users/Darren/Documents/2018 Unifilis demography/analysis/shapes/ottobacias_ap.shp"
s4 <- read_sf(inshp3)
sint <- st_intersection(s3, s4)
plot(sint["NIVEL33"])

inshp4 <- "C:/Users/Darren/Documents/2018 Unifilis demography/analysis/shapes/ottobacias_ap_03.shp"
s5 <- read_sf(inshp4)

s2 <- st_intersection(s2, s4)
table(s2$BACIA_DN9)
# now add projecton value to each point, based on three scenarios
#Business-As-Usual
# column for new values
s2$fa_bau <- NA
#select relevant rates for scenario
selBAUin <- which(dffinal$hunt==0 & dffinal$increase == 0.2)
selBAUacc <- which(dffinal$hunt==10 & dffinal$increase == 0.1)
# select points with criteria
selin <- which(s2$accessible=="No")
# add population projection values for scenario
#inaccessible at base rate
s2[selin, 'fa_bau'] <- dffinal[selBAUin, 'adult_females_clean']
#accessible with hunting and nest removal
s2[-selin, 'fa_bau'] <- dffinal[selBAUacc, 'adult_females_clean']
summary(s2$fa_bau)
plot(s2["fa_bau"])
table(s2$fa_bau)
#Protection
# column for new values
s2$fa_pr <- NA
#select relevant rates for scenario
selPrPA <- which(dffinal$hunt==0 & dffinal$increase == 0.2)
selPrNPAacc <- which(dffinal$hunt==10 & dffinal$increase == 0.1)
# select points with criteria
selin <- which(s2$accessible=="No")
selaccPA <- which(s2$accessible=="Yes" & s2$All==1)
selaccNPA <- which(s2$accessible=="Yes" & s2$All==0)
# add population projection values for scenario
#inaccessible at base rate
s2[selin, 'fa_pr'] <- dffinal[selPrPA, 'adult_females_clean']
#accessible and protected at base rate
s2[selaccPA, 'fa_pr'] <- dffinal[selPrPA, 'adult_females_clean']
##accessible and not protected with hunting and nest removal
s2[selaccNPA, 'fa_pr'] <- dffinal[selPrNPAacc, 'adult_females_clean']
summary(s2$fa_pr)
plot(s2["fa_pr"])
table(s2$fa_bau)
table(s2$All)
table(s2$accessible)

#Community-Based-Management
# column for new values
s2$fa_cbm <- NA
#select relevant rates for scenario
selCBMin <- which(dffinal$hunt==0 & dffinal$increase == 0.2)
selCBMNPAacc <- which(dffinal$hunt==10 & dffinal$increase == 0.5)
selCBMPAacc <- which(dffinal$hunt==10 & dffinal$increase == 0.1)
# select points with criteria
selin <- which(s2$accessible=="No")
selaccPA <- which(s2$accessible=="Yes" & s2$All==1)
selaccNPA <- which(s2$accessible=="Yes" & s2$All==0)
# add population projection values for scenario
#inaccessible at base rate
s2[selin, 'fa_cbm'] <- dffinal[selCBMin, 'adult_females_clean']
#accessible and protected at base rate (hunting and nest removal)
s2[selaccPA, 'fa_cbm'] <- dffinal[selCBMPAacc, 'adult_females_clean']
##accessible and not protected with headstart and hunting
s2[selaccNPA, 'fa_cbm'] <- dffinal[selCBMNPAacc, 'adult_females_clean']
summary(s2$fa_cbm)
plot(s2["fa_cbm"])
table(s2$fa_cbm)
table(s2$accessible)

# Pr and CBM
# column for new values
s2$fa_cbmpr <- NA
#select relevant rates for scenario
selCBMprin <- which(dffinal$hunt==0 & dffinal$increase == 0.2)
selCBMprNPAacc <- which(dffinal$hunt==10 & dffinal$increase == 0.5)
selCBMprPAacc <- which(dffinal$hunt==0 & dffinal$increase == 0.2)
# select points with criteria
selin <- which(s2$accessible=="No")
selaccPA <- which(s2$accessible=="Yes" & s2$All==1)
selaccNPA <- which(s2$accessible=="Yes" & s2$All==0)
# add population projection values for scenario
#inaccessible at base rate
s2[selin, 'fa_cbmpr'] <- dffinal[selCBMprin, 'adult_females_clean']
#accessible and protected at base rate (hunting and nest removal)
s2[selaccPA, 'fa_cbmpr'] <- dffinal[selCBMprPAacc, 'adult_females_clean']
##accessible and not protected with headstart and hunting
s2[selaccNPA, 'fa_cbmpr'] <- dffinal[selCBMprNPAacc, 'adult_females_clean']
summary(s2$fa_cbmpr)
plot(s2["fa_cbmpr"])
table(s2$fa_cbmpr)

#total and proportional difference per catchments
plot(s2["subbasn"])
library(plyr)
ddply(s2, .(subbasn), summarise, 
      km_tot = length(na.omit(fa_cbm)), 
      km_pa = sum(na.omit(All)),
      fa_t0_tot = sum(na.omit(fa_t0)),
      fa_bau_tot = sum(na.omit(fa_bau)),
      fa_pr_tot = sum(na.omit(fa_pr)),
      fa_cbm_tot = sum(na.omit(fa_cbm)), 
      fa_cbmpr_tot = sum(na.omit(fa_cbmpr))
      )
# diff
ddply(s2, .(subbasn), summarise, 
      km_tot = length(na.omit(fa_cbm)), 
      km_pa = sum(na.omit(All)),
      fa_t0_tot = sum(na.omit(fa_t0)),
      fa_bau_diff = (sum(na.omit(fa_bau)) - sum(na.omit(fa_t0))) / sum(na.omit(fa_t0)),
      fa_pr_diff = (sum(na.omit(fa_pr))- sum(na.omit(fa_t0))) / sum(na.omit(fa_t0)),
      fa_cbm_diff = (sum(na.omit(fa_cbm))- sum(na.omit(fa_t0))) / sum(na.omit(fa_t0)), 
      fa_cbmpr_diff = (sum(na.omit(fa_cbmpr))- sum(na.omit(fa_t0))) / sum(na.omit(fa_t0))
)
length(unique(s2$NIVEL33))
ddply(s2, .(NIVEL33), summarise, 
      km_tot = length(na.omit(fa_cbm)), 
      km_pa = sum(na.omit(All)),
      fa_t0_tot = sum(na.omit(fa_t0)),
      fa_bau_tot = sum(na.omit(fa_bau)),
      fa_pr_tot = sum(na.omit(fa_pr)),
      fa_cbm_tot = sum(na.omit(fa_cbm)), 
      fa_cbmpr_tot = sum(na.omit(fa_cbmpr))
)
# diff
dfn3 <- ddply(s2, .(NIVEL33), summarise, 
      km_tot = length(na.omit(fa_cbm)), 
      km_pa = sum(na.omit(All)),
      fa_t0_tot = sum(na.omit(fa_t0)),
      fa_bau_diff = (sum(na.omit(fa_bau)) - sum(na.omit(fa_t0))) / sum(na.omit(fa_t0)),
      fa_pr_diff = (sum(na.omit(fa_pr))- sum(na.omit(fa_t0))) / sum(na.omit(fa_t0)),
      fa_cbm_diff = (sum(na.omit(fa_cbm))- sum(na.omit(fa_t0))) / sum(na.omit(fa_t0)), 
      fa_cbmpr_diff = (sum(na.omit(fa_cbmpr))- sum(na.omit(fa_t0))) / sum(na.omit(fa_t0))
)

sint2 <- merge(sint, dfn3)
plot(sint2["fa_pr_diff"])

## plots
sfcoun<- rnaturalearth::ne_countries(continent = "South America", 
                                     type = 'map_units', returnclass = "sf")
sfcoun <- st_sf(a= rep(1,14), geom=st_geometry(sfcoun))
sfcounD <- sf::st_union(rnaturalearth::ne_countries(continent = "South America", 
                                                    type = 'map_units', returnclass = "sf"))
sfcounD <- st_sf(a=1, geom=st_geometry(sfcounD))

sfap <- st_union(st_buffer(sint2, 0.00001))
sfapn3 <-st_intersection(s5, sfap)
library(ggplot2)
xbr = c(-80, -70, -60, -50)
ybr = c(-20, -10, -0, 10)

#Ottobacias. Need to drop to make sure match data
selOt <- which(sint2$NIVEL33 %in% unique(dfn3$NIVEL33))
ggplot(sint2[selOt, ]) +
  #geom_sf(data = sfcounD) +
  geom_sf(aes(fill = NIVEL33) ) +
  geom_sf(data = sfapn3, size= 0.8, color="grey30", fill=NA) +
  geom_sf(data = sfapn3, size=0.3,color="white", fill=NA, lty=2) +
  geom_sf(data = sfap, size=1,color="black", fill=NA) +
  #geom_sf(data = sfap, size=1,color="yellow", fill=NA, lty=2) +
  scale_fill_discrete("Ottobacia\nNivel 3") +
  coord_sf(xlim = c(-54.8, -49.8), ylim = c(-1.2, 4.3))+
  #scale_x_continuous(breaks = xbr) +
  #scale_y_continuous(breaks = ybr) +
  theme_bw() +
  theme(legend.margin=margin(t=0, r=0, b=0, l= -0.2, unit="cm")) +
  ggtitle("A) Ottobacias")

#BAU
ggplot(sint2) +
  #geom_sf(data = sfcounD) +
  geom_sf(aes(fill = fa_bau_diff) ) +
  geom_sf(data = sfapn3, size= 0.8, color="grey30", fill=NA) + 
  geom_sf(data = sfap, size=1,color="black", fill=NA) +
  geom_sf(data = sfap, size=1,color="yellow", fill=NA, lty=2) +
  scale_fill_gradient2("%\nchange") +
  coord_sf(xlim = c(-54.8, -49.8), ylim = c(-1.2, 4.3))+
  #scale_x_continuous(breaks = xbr) +
  #scale_y_continuous(breaks = ybr) +
  theme_bw() +
  theme(legend.margin=margin(t=0, r=0, b=0, l= -0.2, unit="cm")) +
  ggtitle("A) Business-As-Usual")

ggplot(dfn3, aes(x = NIVEL33, y = fa_bau_diff, fill = fa_bau_diff)) + 
  geom_col() + 
  scale_fill_gradient2("%\nchange")

ggplot(sint2) +
  #geom_sf(data = sfcounD) +
  geom_sf(aes(fill = fa_pr_diff) ) +
  geom_sf(data = sfapn3, size= 0.8, color="grey30", fill=NA) + 
  geom_sf(data = sfap, size=1,color="black", fill=NA) +
  geom_sf(data = sfap, size=1,color="yellow", fill=NA, lty=2) +
  scale_fill_gradient2("%\nchange") +
  coord_sf(xlim = c(-54.8, -49.8), ylim = c(-1.2, 4.3))+
  #scale_x_continuous(breaks = xbr) +
  #scale_y_continuous(breaks = ybr) +
  theme_bw() +
  theme(legend.margin=margin(t=0, r=0, b=0, l= -0.2, unit="cm")) +
  ggtitle("B) Protection")

ggplot(dfn3, aes(x = NIVEL33, y = fa_pr_diff, fill = fa_pr_diff)) + 
  geom_col() + 
  scale_fill_gradient2("%\nchange")

# CBM
ggplot(sint2) +
  #geom_sf(data = sfcounD) +
  geom_sf(aes(fill = fa_cbm_diff) ) +
  geom_sf(data = sfapn3, size= 0.8, color="grey30", fill=NA) + 
  geom_sf(data = sfap, size=1,color="black", fill=NA) +
  geom_sf(data = sfap, size=1,color="yellow", fill=NA, lty=2) +
  scale_fill_gradient2("%\nchange") +
  coord_sf(xlim = c(-54.8, -49.8), ylim = c(-1.2, 4.3))+
  #scale_x_continuous(breaks = xbr) +
  #scale_y_continuous(breaks = ybr) +
  theme_bw() +
  theme(legend.margin=margin(t=0, r=0, b=0, l= -0.2, unit="cm")) +
  ggtitle("c) Community-Based Management")

ggplot(dfn3, aes(x = NIVEL33, y = fa_cbm_diff, fill = fa_cbm_diff)) + 
  geom_col() + 
  scale_fill_gradient2("%\nchange")

# CBM and Protection
ggplot(sint2) +
  #geom_sf(data = sfcounD) +
  geom_sf(aes(fill = fa_cbmpr_diff) ) +
  geom_sf(data = sfapn3, size= 0.8, color="grey30", fill=NA) + 
  geom_sf(data = sfap, size=1,color="black", fill=NA) +
  geom_sf(data = sfap, size=1,color="yellow", fill=NA, lty=2) +
  scale_fill_gradient2("%\nchange") +
  coord_sf(xlim = c(-54.8, -49.8), ylim = c(-1.2, 4.3))+
  #scale_x_continuous(breaks = xbr) +
  #scale_y_continuous(breaks = ybr) +
  theme_bw() +
  theme(legend.margin=margin(t=0, r=0, b=0, l= -0.2, unit="cm")) +
  ggtitle("D) Community-Based Management Protection")

ggplot(dfn3, aes(x = NIVEL33, y = fa_cbmpr_diff, fill = fa_cbmpr_diff)) + 
  geom_col() + 
  scale_fill_gradient2("%\nchange")
