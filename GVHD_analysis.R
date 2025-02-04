library(survival)
library(mstate)
library(data.table)
library(LTRCtrees)
library(partykit)
library(rpart)
library(openxlsx)

set.seed(1)

path_to_data <- "/path/to/data"
pathtofigures <- "/path/to/figures"

dataSet_org <- read.xlsx(paste0(path_to_data,"Data.xlsx"), sheet = 'Data')

dataSet_org$Osdays = dataSet_org$DFSdays/365
dataSet_org$Start =  dataSet_org$Start/365
dataSet_org$Stop =  dataSet_org$Stop/365

covariate_org = c("Diag_ALL", "Diag_CML", "Diag_MDS", "Disease_risk", "Ptsex","Dosex", 
                 "FtoM", 'trans_year', "PB1vsother0", "CB1vsother0", 
                "PtCMVAb", "DoCMVAb", "PSgood0poor1", "R0U1", 
                "HLA8_8dis", "TBI", "TCD", "MAC0RIC1", "Fkbased", "MMF", "CMVearlyTx",
                "Diag_ALL_CMV", "Diag_CML_CMV", "Diag_MDS_CMV", "Disease_risk_CMV", "Ptsex_CMV","Dosex_CMV", 
                 "FtoM_CMV", 'trans_year_CMV', "PB1vsother0_CMV", "CB1vsother0_CMV",
                "PtCMVAb_CMV", "DoCMVAb_CMV", "PSgood0poor1_CMV", "R0U1_CMV", 
                "HLA8_8dis_CMV", "TBI_CMV", "TCD_CMV", "MAC0RIC1_CMV", "Fkbased_CMV", "MMF_CMV", "G2_4GVHD")

covariate_Factor = covariate_org

UsingColumn = c(covariate_Factor, "Ptage", "Doage", "Ptage_CMV", "Doage_CMV", "Start", "Stop", "Event", "ID2", "DFSalive0RL1NRM2")
dataSet_org2 = dataSet_org[, UsingColumn]
dataSet = na.omit(dataSet_org2)


ix = c()
Num_covariate = length(covariate_Factor)
for (i in 1:Num_covariate){
  ix = c(ix, match(covariate_Factor[i], colnames(dataSet)))
}

dataSet[ix] = lapply(dataSet[ix], as.factor)


covariate_N = c("Ptage", "Doage", "Ptage_CMV", "Doage_CMV")

ix = c()
Num_covariate = length(covariate_N)
for (i in 1:Num_covariate){
  ix = c(ix, match(covariate_N[i], colnames(dataSet)))
}

dataSet[ix] = lapply(dataSet[ix], as.numeric)
covariate = c(covariate_Factor, covariate_N)

str(dataSet)
head(dataSet)
summary(dataSet)

dataSet <- as.data.table(dataSet)

dataSet_w <- crprep(Tstart = "Start", Tstop = "Stop", status = "Event", data = dataSet, 
                   trans = 2, cens = 0, id = "ID2", 
                   keep = covariate)

dataSet_w <- as.data.table(dataSet_w)

dataSet_w[, event := ifelse(status == 2, 1, 0)]
weight <- dataSet_w$weight.cens * dataSet_w$weight.trunc


# Cox regression
model_OS_w <- coxph(as.formula(paste0("Surv(Tstart, Tstop, status == \"2\") ~ ", paste(covariate, collapse = " + "))),
                                                     data = dataSet_w,
                                                     weight = weight.cens * weight.trunc)

model_OS_w_fit <- survfit(model_OS_w)
summary(model_OS_w_fit)
summary(model_OS_w)


# LTRCtrees
covariate_tree <- c("Diag_ALL", "Diag_CML", "Diag_MDS", "Disease_risk", "Ptsex","Dosex", "Ptage", "Doage", 
                    "FtoM", 'trans_year', "PB1vsother0", "CB1vsother0", "PtCMVAb", "DoCMVAb", "PSgood0poor1", "R0U1", 
                    "HLA8_8dis", "TBI", "TCD", "MAC0RIC1", "Fkbased", "MMF", "CMVearlyTx", "G2_4GVHD")

cp = 9e-04
pathtofigures_N = paste0(pathtofigures, "cp=", as.character(cp))
dir.create(pathtofigures_N)
pathtofigures_N2 = paste0(pathtofigures_N, "/")


minsplit=20
rpart_control <- list(minsplit = 20, minbucket = round(minsplit/3), cp = cp, 
                                    maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10,
                                    surrogatestyle = 0, maxdepth = 30)

DF = as.data.frame(do.call(rbind, rpart_control))
write.table(DF, paste0(pathtofigures_N2, "Parm.txt"), quote=F, sep="\t", row.names=T, col.names=F)

mlab <- function(id, nobs) id = paste0(id, " (n = ", round(nobs), ")")


for (i in 1:100){
# model 1
set.seed(i)
model_tree_1a <- LTRCART(as.formula(paste0("Surv(Tstart, Tstop, event) ~ ", paste(covariate_tree, collapse = " + "))),
                                                     data = dataSet_w,
                                                     weights = weight,
                                                     control = rpart_control)
model_tree_1a_party <- as.party(model_tree_1a)
model_tree_1a_party$fitted[["(response)"]] <- Surv(dataSet_w$Tstart, dataSet_w$Tstop, dataSet_w$event)


#########
output_file <- "output_Prob.xlsx"
wb <- createWorkbook()
output_file2 <- "output_point.xlsx"
wb2 <- createWorkbook()
wb2_sheet_name <- "DataSheet"  
addWorksheet(wb2, sheetName = wb2_sheet_name)

node_table = table(predict(model_tree_1a_party, type = "node")) 
branch_ids <- as.numeric(names(node_table))
branch_size <- length(branch_ids)

data_row <- t(branch_ids)
writeData(wb2, sheet = wb2_sheet_name, data_row, startRow = 1, startCol = 1)


for(bi in 1:branch_size){
branch_id <- branch_ids[bi]
branch_id_string <- as.character(branch_id)

branch_data <- dataSet_w[model_tree_1a$where == branch_id, ]

surv_fit  <- survfit(Surv(branch_data$Tstart, branch_data$Tstop, branch_data$event) ~ 1, data = branch_data, weights = weight.cens * weight.trunc) 

result_data <- data.frame(
  surv = surv_fit$time,
  time = surv_fit$surv
)

sheet_name <- paste("branch_id_", branch_id, sep = "")

addWorksheet(wb, sheetName = sheet_name)

writeData(wb, sheet = sheet_name, result_data)

surv_prob_years <- numeric(11)
for(time_point in 1:11){
  surv_prob_years[time_point] <- surv_fit$surv[which(surv_fit$time >= time_point)][1]
}

writeData(wb2, sheet = wb2_sheet_name, surv_prob_years, startRow = 3, startCol = bi)
}

saveWorkbook(wb, output_file, overwrite = TRUE)
saveWorkbook(wb2, output_file2, overwrite = TRUE)


png(paste0(pathtofigures_N2, "model1a_All_", as.character(i), ".png"), width = 80, height = 42, units = "cm", res = 500)
plot(model_tree_1a_party, 
     tp_args = list(mainlab = mlab, gp = gpar(fontsize=17), at=c(0,1,11)), 
     ip_args = list(gp = gpar(fontsize=20)),
     gp = gpar(fontsize=17),
     tnex = 1.5)
dev.off()

}
