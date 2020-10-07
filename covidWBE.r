
#### covidWBE.r
#### Nick Crosbie, May 2020
#### source("./bin/covidWBE.r")

library(EnvStats)

# WASTEWATER CATCHMENT (calculate dry weather flow in millilitres per day)
catchmentPop <- 10^4 # number of people in catchment
dryWeatherFlowPerPerson <- 2.12799e-4 # dry weather flow contributed by each person, in ** megalitres ** per day (average derived from all participating ColoSSoS WWTP catchments)
dryWeatherFlow_MLD <- catchmentPop * dryWeatherFlowPerPerson # catchment dry weather flow, in ** megalitres ** per day
dryWeatherFlow_mLD <- dryWeatherFlow_MLD * 1e9 # catchment dry weather flow, in ** millilitres ** per day

# SARS-CoV-2 GENOMES SHED PER DAY TO CATCHMENT BY DEFECATION (number of genomes)
modeShed <- 10^3 # mode number of genomes shed in faecies by an infected person (genomes per millilitre of faecies)
minShed <- 10^2 # min number of genomes shed in faecies by an infected person (genomes per millilitre of faecies)
maxShed <- 10^7 # max number of genomes shed in faecies by an infected person (genomes per millilitre of faecies)

defVol <- 200 # volume of faecies defacated per day by an individual, in millilitres (no distinction made between volume defecated by infected and non-infected individuals)
numberShedding <- 1 # number of people shedding in catchment
defVolInfCat <- defVol * numberShedding #  combined volume of faecies defacated by shedders per day, in millilitres

modeShedCat <- modeShed * defVolInfCat # mode number of genomes shed per day to catchment by infected population
minShedCat <- minShed * defVolInfCat # minimum number of genomes shed per day to catchment by infected population
maxShedCat <- maxShed  * defVolInfCat # maximum number of genomes shed per day to catchment by infected population

# TRIANGULAR DISTRIBUTION OF SARS-CoV-2 GENOMES SHED PER DAY TO CATCHMENT BY DEFECATION (number of genomes)
triagShedCat <- simulateVector(n = 10000, distribution = "tri", param.list = list(min = minShedCat, mode = modeShedCat, max = maxShedCat),
sample.method = "LHS", seed = 47)

pdfPlot(distribution = "tri", param.list = list(min = minShedCat, mode = modeShedCat, max = maxShedCat), n.points = 10000)

# SARS-CoV-2 CONCENTRATION AT POINT OF WASTEWATER SUBSAMPLING (genomes per mL)
catchmentDilution <- defVolInfCat / dryWeatherFlow_mLD
wastewaterConc <- triagShedCat * catchmentDilution # genomes per mL of catchment wastewater, ** assuming fully mixed catchment **

# plot(ecdf(wastewaterConc))

# SARS-CoV-2 CONCENTRATION OF FINAL CONCENTRATED PREPARATION
concFactor <- 4  # concentration factor from point of wastewater sampling to final concentrated preparation
finalPrepConc <- wastewaterConc * concFactor # virus concentration prior to subsampling for PCR (genomes per millilitre)

#plot(ecdf(finalPrepConc))

# SARS-CoV-2 GENOMES PER PCR REACTION (number of genomes)
finalConcVol <- 2 # volume of sample prior to subsampling for PCR reaction, in millilitres
genomesFinalPrep <- finalPrepConc * finalConcVol # number of genomes contained in the final prep (step immediately prior to PCR)

templateVol <- 0.005 # volume in which the final preparation is resuspended prior to pipetting to PCR reaction
volRxn <- 0.05 # volume of one PCR reaction, in millilitres (i.e. 50 microlitres)
genomesRxn <- genomesFinalPrep * (templateVol/volRxn) # genomes per reaction
LOD <- 10 # PCR limit of detection as genomes per reaction

n <- 10000 # number of wastewater samples analysed for SARS-CoV-2
detects <- sum(genomesRxn > LOD) # number of SARS-CoV-2 detects
trueDetectProportion <- 0.90 # proportion of detects that are true positives
percentRealDetects <- ((detects * trueDetectProportion) / n) * 100 # percentage of samples where SARS-CoV-2 is detected 

print(percentRealDetects)

# plot(ecdf(genomesRxn))


