#' Title
#'
#' @param x Data.frame with parameters from cmartr::PopParam
#' @param fed Density of females per km
#'
#' @description Projects populations through different scenarios across 
#' species range.
#' 
#' @return Data.frame with projections for each parameter provided in x.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' dfpop <- cmartr::PopParam(species = "Podocnemis unifilis", make_rds = FALSE)
#' library(plyr)
#' dff <- ddply(dfpop, .(akey), .fun = PopProjF)
#' }
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