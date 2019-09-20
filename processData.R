library(PharmacoGx)
library(readr)
library(tximport)
library(rhdf5)
library(gdata)
library(readxl)
library(openxlsx)
library(CoreGx)
library(data.table)


prefix <- "/pfs/"

options(stringsAsFactors = FALSE)
gCSI_GR_AOC <- read.csv(file.path(prefix,"gCSI2018_data/gCSI_GRvalues_v1.2.tsv"), sep="\t")

## Using data.table for efficency and consistency with other psets. If you 
## arent familiar with syntax (it implements inplace modifications), please 
## read up on it.
gCSI_GR_AOC <- data.table(gCSI_GR_AOC)

gCSI_GR_AOC[,expid := paste(CellLineName, DrugName, ExperimentNumber,sep="_")]
max_conc <- max(gCSI_GR_AOC[,.N, expid][,N])
uids <- unique(gCSI_GR_AOC$expid)

## Create matrices for the data
doses_final <- matrix(NA_real_,nrow=length(uids), ncol = max_conc, dimnames= list(uids, paste0("dose", seq_len(max_conc))))
viability_final <- matrix(NA_real_,nrow=length(uids), ncol = max_conc, dimnames= list(uids, paste0("dose", seq_len(max_conc))))

setorder(gCSI_GR_AOC, expid, log10Concentration)
gCSI_GR_AOC_list <- split(gCSI_GR_AOC, by="expid")

for (exp in names(gCSI_GR_AOC_list)) {
  xx <- gCSI_GR_AOC_list[[exp]]
  concentrations=10^(xx$log10Concentration) #remove log10
  concentrations= concentrations / 1.0e-6 #convert molar to uM
  viability = xx$relative_cell_count * 100 
  
  doses_final[exp,1:length(concentrations)] <- concentrations
  viability_final[exp,1:length(viability)] <- viability

}


raw.sensitivity <- array(c(as.matrix(doses_final), as.matrix(viability_final)),
                             c(length(uids), max_conc, 2),
                             dimnames=list(rownames(doses_final),
                                           colnames(doses_final),
                                           c("Dose", "Viability")))
sensitivityInfo_2018 <- as.data.frame(gCSI_GR_AOC[,c("expid", "DrugName","CellLineName","ExperimentNumber","TrtDuration","doublingtime")])
sensitivityInfo_2018 <- unique(sensitivityInfo_2018)
rownames(sensitivityInfo_2018) <- sensitivityInfo_2018$expid

## Read in published data.

gCSI_GR_AOC_Pub <- read.csv("2018/raw_data/gCSI_GRmetrics_v1.2.tsv", sep="\t")



uids2 <- unique(sprintf("%s_%s_%s",gCSI_GR_AOC_Pub$CellLineName,
                       gCSI_GR_AOC_Pub$DrugName,
                       gCSI_GR_AOC_Pub$ExperimentNumber))

rownames(gCSI_GR_AOC_Pub) <- uids2

sensitivity.info <- sensitivityInfo_2018
sensitivity.published <- gCSI_GR_AOC_Pub

save(raw.sensitivity, sensitivity.info, sensitivity.published, file="/pfs/out/sens.data.RData")

raw.sensitivity.x <- parallel::splitIndices(nrow(raw.sensitivity), floor(nrow(raw.sensitivity)/1000))




dir.create("/pfs/out/slices/")

for(i in seq_along(raw.sensitivity.x)){

  slce <- raw.sensitivity[raw.sensitivity.x[[i]],,]
  saveRDS(slce, file=paste0("/pfs/out/slices/gcsi2017_raw_sens_", i, ".rds"))

}

