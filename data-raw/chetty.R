# Process Chetty et al. (2014) data for the spatialunitroot package
library(readxl)

# Read raw data
labels <- readxl::read_excel("inst/extdata/Chetty_Data_Labels.xlsx",
                              sheet = "Sheet1", col_names = TRUE)
d <- readxl::read_excel("inst/extdata/Chetty_Data_1.xlsx",
                         sheet = "Sheet1", col_names = TRUE)

# Rename coordinates to s_* for Stata/SPUR compatibility
names(d)[names(d) == "Lat"] <- "s_1"
names(d)[names(d) == "Lon"] <- "s_2"

# Drop non-contiguous states
d <- d[!d$State %in% c("HI", "AK"), ]

# Select key variables from the example.do file
keep_vars <- c("State", "CZ", "CZName", "s_1", "s_2",
               "AM", "FracBlack", "RacSeg", "SegPov25",
               "FracCom15", "HIPC", "Gini", "IncSh1",
               "TSR", "TSPerc", "HSDrop", "SCInd",
               "FracRel", "CrimeR", "FracSM", "FracDiv",
               "FracMar", "LocTR", "ColPC", "ColTui",
               "ColGrad", "ManShare", "ChImp", "TLFPR",
               "MigIRate", "MigORate", "FracFor",
               "SegAfl75", "IncSeg", "Frac2575",
               "GiniBot99", "LocGexp", "SchExp", "LFPR")

# Keep only columns that exist
keep_vars <- intersect(keep_vars, names(d))
chetty <- d[, keep_vars]

# Convert to lowercase names for R convention
names(chetty) <- tolower(names(chetty))

# Remove rows with missing coordinates
chetty <- chetty[!is.na(chetty$s_1) & !is.na(chetty$s_2), ]

# Save
usethis::use_data(chetty, overwrite = TRUE)

message("Chetty data processed: ", nrow(chetty), " obs, ", ncol(chetty), " vars")
