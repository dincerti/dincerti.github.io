
# This creates the dataset loadclean_meps.csv
setwd("C:/Users/Devin/Dropbox/Projects/dincerti.github.io")

# set up
for(con in dbListConnections(MySQL())) dbDisconnect(con)
library("RMySQL")
library("data.table")
con <- dbConnect(MySQL(), dbname="meps", host="localhost",
                 default.file = "C:/Program Files/MySQL/MySQL Server 5.5/my.ini")

# full-year consolidated data files
fyc <- data.table(dbGetQuery(conn = con, "SELECT * FROM fyc WHERE panel >= 10"))
pe <- c("angidx", "asthdx", "chddx","diabdx", "emphdx", "hibpdx", "ohrtdx",
        "midx", "strkdx", "choldx")
for (col in pe) fyc[, (col) := ifelse(get(col) == 1, 1, ifelse(get(col) > 0, 0, NA))]
fyc[, dsa1c := ifelse(dsa1c < 0, NA, ifelse(dsa1c == 96, 0, dsa1c))]
ds <- c("dsdiet", "dseypr", "dskidn", "dsinsu")
for (col in ds) fyc[, (col) := ifelse(get(col) == 1, 1, ifelse(get(col) > 0, 0, NA))]

# medical conditions
mc <- data.table(dbGetQuery(conn = con, "SELECT * FROM mc WHERE panel >= 10"))
mc <- mc[as.numeric(cccodex) > 0]
mc <- mc[, .N, by = .(dupersid, panel, year, cccodex)]
mc$N <- 1
mc$ccs <- paste0("ccs", mc$cccodex)
mc <- dcast.data.table(mc, year + dupersid + panel ~ ccs, value.var = "N")
mc[is.na(mc)] <- 0

# meps data set
meps <- merge(fyc, mc, by = c("year", "dupersid", "panel"), all.x = TRUE)
ccsnames <- colnames(meps)[grep("ccs", colnames(meps))]
for (col in ccsnames) meps[is.na(get(col)), (col) := 0]
meps <- meps[ccs049 == 1 | ccs050 == 1]
meps <- meps[inscov == 1 | inscov == 2]
meps <- meps[order(dupersid, panel, year)]
meps[, l_totexp :=c(NA, totexp[-.N]), by = dupersid]
meps[, age := ifelse(age < 0, NA, age)]
meps[, age2 := age^2]
meps[, srh := as.factor(ifelse(rthlth53 < 0, NA, rthlth53))] # health status
meps <- cbind(meps, model.matrix(~ srh, model.frame(~ srh, meps, na.action = NULL))[, -1])
meps[, totexp_ecdf := ecdf(totexp)(totexp), by = year]
meps[, topexp := ifelse(totexp_ecdf > .9, 1, 0)]
meps[, l_topexp :=c(NA, topexp[-.N]), by = dupersid]
meps[, log_dsa1c:= log(dsa1c + 1)]

# save
write.csv(meps, "data/mepsdiab.csv")