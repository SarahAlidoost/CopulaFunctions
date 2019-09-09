# a script to test copulas on the data: ELF_human
library(RSQLite, quietly = TRUE)
library(fcros, quietly = TRUE)
library(devtools, quietly = TRUE)
library(copula, quietly = TRUE)
library(VineCopula, quietly = TRUE)

dbfile = "/home/sarah/MaxQuantData/ELF_human.sqlite"
dts = 'VH10'

query<-paste0("SELECT
    A.grp_id grp_id,
    IFNULL(GROUP_CONCAT(DISTINCT gene), '-') genes,
   ", dts, "_L0_M0_H1_norm_ratio_HL AS ratio_H1L0, -- norm. ratio ON/OFF  (treat1)
   ", dts, "_L1_M1_H0_norm_ratio_LH AS ratio_L1H0, -- norm. ratio ON/OFF  (treat2)
   ", dts, "_L0_M0_H1_norm_ratio_HM AS ratio_H1M0, -- norm. ratio ON/OFF  (treat3)
   ", dts, "_L1_M1_H0_norm_ratio_MH AS ratio_M1H0, -- norm. ratio ON/OFF  (treat4)
   ", dts, "_L0_M0_H1_norm_ratio_LM AS ratio_L0M0, -- norm. ratio OFF/OFF (ctrl1)
   ", dts, "_L1_M1_H0_norm_ratio_LM AS ratio_L1M1  -- norm. ratio ON/ON   (ctrl2)
FROM
   VVV_PGROUP_QUANT A, PROT2GRP B, V_PROTEIN C
WHERE
   A.grp_id = B.grp_id
   AND B.prot_acc = C.acc
   AND ((", dts, "_L1_M1_H0_norm_ratio_HL > 1
   AND   ", dts, "_L1_M1_H0_norm_ratio_HM > 1
   AND   ", dts, "_L0_M0_H1_norm_ratio_LH > 1
   AND   ", dts, "_L0_M0_H1_norm_ratio_MH > 1)
   OR   (", dts, "_L1_M1_H0_norm_ratio_LH > 1
   AND   ", dts, "_L1_M1_H0_norm_ratio_MH > 1
   AND   ", dts, "_L0_M0_H1_norm_ratio_HL > 1
   AND   ", dts, "_L0_M0_H1_norm_ratio_HM > 1))
GROUP BY A.grp_id")

# fetch data from db
drv<-dbDriver('SQLite')
conn<-tryCatch(dbConnect(drv, dbname = dbfile),
               error = function(e) {
                 cat(paste0("Error: Database file '", dbfile, "' not found.\n"))
                 q(save = 'no', status = 1)
               })
res<-tryCatch(dbSendQuery(conn, query), error = function(e) {
  cat('Error: Selected data set not found in dbfile.\n')
  q(save = 'no', status = 1)})
tbl<-fetch(res, n = -1)

# columns excluded from further processing
skipcols<-c('grp_id', 'genes')

# dynamically generate column names with log2ratio_* prefix
rcols<-names(tbl)[!(names(tbl) %in% skipcols)] # select columns with SILAC ratios

# dynamically generate column names with log2ratio_* prefix
rcols<-names(tbl)[!(names(tbl) %in% skipcols)] # select columns with SILAC ratios
log2cols<-paste0('log2', rcols)

# log2-transform SILAC ratios, append new columns to the table
tbl[log2cols]<-log2(tbl[rcols])

# columns with log2 ratios corresponding to control conditions
ctrl<-c('log2ratio_L0M0', 'log2ratio_L1M1')

# select columns with log2 ratios corresponding to treatment conditions
treat<-log2cols[!(log2cols %in% ctrl)]

# ranking data: empirical probabilities
tblRn<-pobs(tbl, ties.method = "average") 

# selecting two variables: ratio_L1H0 and ratio_H1M0
u1<-c(2)
u2<-c(3)
plot(tblRn[,treat[u1:u2]])

# correlations and bivariate copula
tblTau<-cor(tblRn[,rcols[u1:u2]], use ="everything", method ="kendall") # kendall correlation
tblrho<-cor(tblRn[,rcols[u1:u2]], use ="everything", method ="pearson")
tblCop<-BiCopSelect(tblRn[,rcols[u1]], tblRn[,rcols[u2]], familyset = c(1:5), selectioncrit = "AIC", indeptest = T, rotations = F)

