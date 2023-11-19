############################################
### Code for: Deer Sex Ratios Revisited  ###
### By: Jane Dentinger                   ###
### Last Updated: 1/10/2023              ###
############################################

# Clear the workspace
rm(list = ls(all.names = TRUE))
gc()
# Clear all loaded packages
lapply(paste("package:", names(sessionInfo()$otherPkgs), sep=""), 
       detach, 
       character.only = TRUE, 
       unload = TRUE)

set.seed(1)

# ##################
# # Clean the Data #
# ##################
# # Import the capture data
# deer = read.csv("Deer Fetal Sex Ratio Database.csv")
# 
# # Cleaning libraries
# library(dplyr)
# library(tidyverse)
# library(hablar)
# 
# # Clean dataset
# deer = deer %>%
#   dplyr::na_if(".") %>%
#   # Remove NAs
#   filter(!is.na(KFI) &
#         !is.na(D_WEIGHT) &
#         !is.na(AGE)) %>%
#   # Remove unwanted columns
#     dplyr::select(-c(DEER_NO, SEX, KIDNEYS, KID_FAT, CL_NO,
#                      COMMENTS, DATEYEAR, JULIAN, D, M, Y,
#                      BREEDDAY, BREEDJUL, JUL, JBREEDAY, INTERVAL)) %>%
#   # Convert to correct file types
#     hablar::convert(num(AGE, D_WEIGHT, KFI, FETUS_NO, FET_AGE, FETUS_T, FETUS_M, FETUS_F, YEAR),
#             chr(CODE, SITE, FET_SEX),
#             dte(COLLDATE, CON_DATE, .args=list(format = "%m/%d/%Y")),
#             fct(SOILCODE)) %>%
#     dplyr::rename(fetus_sex = FET_SEX,
#            fetus_total = FETUS_T,
#            fetus_age = FET_AGE,
#            dress_wt = D_WEIGHT,
#            concept_date = CON_DATE,
#            collect_date = COLLDATE,
#            fetus_num = FETUS_NO) %>%
#     rename_with(toupper)
# 
# # Clean fetus sex groups
# deer = deer %>% mutate(FETUS_SEX = case_when(
#                 (FETUS_TOTAL == 1) ~
#                         dplyr::recode(FETUS_SEX,
#                         "1 F" = "F",
#                         "1 M" = "M",
#                         "1 M" = "M",
#                         "1F" = "F",
#                         "1M" = "M",
#                         "f" = "F",
#                         "m" = "M"),
#                (FETUS_TOTAL==2) ~
#                         dplyr::recode(FETUS_SEX,
#                         "M/M" = "MM",
#                         "M/F" = "MF",
#                         "1M 1F" = "MF",
#                         "2F" = "FF",
#                         "F/F" = "FF",
#                         "2M" = "MM",
#                         "F/M" = "MF",
#                         "mf" = "MF",
#                         "mm" = "MM",
#                         "ff" = "FF",
#                         "2 M" = "MM",
#                         "2 F" = "FF",
#                         "1 F 1 M" = "MF",
#                         "M" = "MM",
#                         "1M1 F" = "MF",
#                         "1F 1M" = "MF",
#                         "1M,1F" = "MF",
#                         "F" = "FF",
#                         "M,F" = "MF",
#                         "F,F" = "FF",
#                         "M,M" = "MM",
#                         "1M, 1F" = "MF",
#                         "2 MALES" = "MM",
#                         "1m1f" = "MF",
#                         "F,M" = "MF",
#                         "1 M 1 F" = "MF",
#                         "1F,1M" = "MF",
#                         "FM" = "MF",
#                         "m" = "MM"),
#                  (FETUS_TOTAL==3) ~
#                         dplyr::recode(FETUS_SEX,
#                         "M/M/M" = "MMM",
#                         "3M" = "MMM",
#                         "2F 1M" = "MFF",
#                         "mmf" = "MFF",
#                         "fff" = "FFF",
#                         "mff" = "MFF",
#                         "3F" = "FFF",
#                         "1M 2F" = "MFF",
#                         "2M 1F" = "MMF",
#                         "FFF" = "FFF",
#                         "ffm" = "MFF",
#                         "mmm" = "MMM",
#                         "FMM" = "MMF",
#                         "F.F.M" = "MFF",
#                         "M.M.M" = "MMM",
#                         "M,M,M" = "MMM",
#                         "1M, 2F" = "MFF",
#                         "2M, 1F" = "MMF",
#                         "1F, 2M" = "MMF",
#                         "2F, 1M" = "MFF",
#                         "F" = "FFF",
#                         "3 M" = "MMM",
#                         "M,F,F" = "MFF",
#                         "F,M,F" = "MFF",
#                         "M" = "MMM",
#                         "FFM" = "MFF")))
# 
# # Filter by year
# deer = deer %>% filter(YEAR > "1990")
# 
# # Calculate proportion male
# deer$PROP_MALE = deer$FETUS_M/deer$FETUS_TOTAL
# 
# # Import soil shape file
# library(sf)
# soil_type_layer = st_read(dsn = "soils/mlra_a_ms.shp")
# # Extract coordinate reference system information
# crs1 = st_crs(soil_type_layer)
# 
# # Import site locations (lat long)
# latlong = read_csv("Deer_Site_LatLong.csv")
# 
# # Clean site locations file
# latlong = latlong %>%
#   # Remove blank last two columns
#   dplyr::select(Site, `Number of Captures`, X, Y, SOILCODE) %>%
#   # Rename column
#   dplyr::rename(Num_Captures = `Number of Captures`) %>%
#   # Convert to uppercase
#   rename_with(toupper) %>%
#   # Remove rows with NAs
#   na.omit()
# 
# # Convert to spatial object
# latlong = st_as_sf(latlong, coords = c("X","Y"), crs = crs1)
# 
# # Intersect the points to the polygon and extract some variables
# soil_df = st_intersection(latlong, soil_type_layer)
# 
# # Clean merged dataset
# soil_df = soil_df %>%
#   # Remove unwanted columns
#   dplyr::select(SITE, SOILCODE, MLRA_ID, MLRA_NAME) %>%
#   # Rename columns
#   dplyr::rename(REGION_ID = MLRA_ID,
#          FULL_NAME = MLRA_NAME) %>%
#   # Convert to data frame
#   as.data.frame() %>%
#   # Remove geometry column - will not cooperate
#   subset(select = 1:4) %>%
#   # Create abbreviation column
#   mutate(ABBREV = dplyr::recode(FULL_NAME,
#                          "Southern Mississippi River Alluvium" = "SMRV",
#                          "Gulf Coastal Plain" = "GCP",
#                          "Southern Mississippi Valley Loess" = "SVL",
#                          "Alabama and Mississippi Blackland Prairie" = "BP",
#                          "Gulf Coast Marsh" = "GCM",
#                          "Eastern Gulf Coast Flatwoods" = "GCF"))
# # Subset site
# site = soil_df %>%
#   dplyr::select(SITE, REGION_ID, FULL_NAME, ABBREV)
# 
# # Join the regional information to the site
# deer = left_join(deer, site, by = c("SITE"))
# 
# # remove regional NAs
# deer = deer[!is.na(deer$ABBREV),]
# 
# # Export the cleaned data
# write_csv(deer, "Deer Fetal Sex Ratio Database Cleaned.csv")
# 
# ################################################################################
# ################################################################################

##########
# Models #
##########
library(tidyverse)
library(ggplot2)
deer = read_csv("Deer Fetal Sex Ratio Database Cleaned.csv")

# Shift legend function
library(gtable)
library(cowplot)
library(grid)

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

################################################################
# Calculate Relative Values for Age, Region and Dressed Weight #
################################################################

#######
# KFI #
#######

# Raw Unscaled KFI by Region
KFI_region = ggplot(data = deer, aes(x = KFI, fill = FULL_NAME))+
  geom_histogram(color = "black")+
  xlab("Kidney Fat Index")+
  ylab("Count")+
  ggtitle("Kidney Fat Index by Region")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, "cm"))+
  facet_wrap(~ABBREV)+
  geom_vline(xintercept = mean(deer$KFI))+
  labs(fill = "Soil Type")

# Put legend as grob
grid.draw(shift_legend(KFI_region))

# Scale KFI by regional mean
deer = deer %>% 
  group_by(ABBREV) %>% 
  mutate(KFI_rel =  KFI - mean(KFI)) %>% 
  ungroup()

# plot scaled KFI
hist_KFI_rel =
  ggplot(data = deer, aes(x = KFI_rel))+
    geom_histogram(color = "black")+
    xlab("Relative Kidney Fat Index")+
    ylab("Count")+
    ggtitle("Relative Kidney Fat Index")+
    theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, "cm"))
hist_KFI_rel

##################
# Dressed Weight #
##################

# Raw Unscaled KFI by Region
DW_region = ggplot(data = deer, aes(x = DRESS_WT, fill = FULL_NAME))+
  geom_histogram(color = "black")+
  xlab("Dressed Weight (kg)")+
  ylab("Count")+
  ggtitle("Dressed Weight by Region")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, "cm"))+
  facet_wrap(~ABBREV)+
  geom_vline(xintercept = mean(deer$DRESS_WT))+
  labs(fill = "Soil Type")

# Put legend as grob
grid.draw(shift_legend(DW_region))

# Scale Dressed Weight by regional mean
deer = deer %>% 
  group_by(ABBREV) %>% 
  mutate(DRESS_WT_rel = DRESS_WT - mean(DRESS_WT)) %>% 
  ungroup()

# plot scaled dressed weight
hist_DW_rel =
  ggplot(data = deer, aes(x = DRESS_WT_rel))+
  geom_histogram(color = "black")+
  xlab("Relative Dressed Weight")+
  ylab("Count")+
  ggtitle("Relative Dressed Weight")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, unit = "cm"))
hist_DW_rel

#######
# AGE #
#######

AGE_region = ggplot(data = deer, aes(x = AGE, fill = FULL_NAME))+
  geom_histogram(color = "black")+
  xlab("Age (years)")+
  ylab("Count")+
  ggtitle("Age by Region")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, "cm"))+
  facet_wrap(~ABBREV)+
  geom_vline(xintercept = mean(deer$AGE))+
  labs(fill = "Soil Type")

# Put legend as grob
grid.draw(shift_legend(AGE_region))

# Scale Age by regional mean
deer = deer %>% 
  group_by(ABBREV) %>% 
  mutate(AGE_rel = AGE - mean(AGE)) %>% 
  ungroup()

# plot scaled age
hist_AGE_rel =
  ggplot(data = deer, aes(x = AGE_rel))+
  geom_histogram(color = "black")+
  xlab("Relative Age")+
  ylab("Count")+
  ggtitle("Relative Maternal Age")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, "cm"))
hist_AGE_rel

###################
# CONCEPTION DATE #
###################
# Convert conception date to Julian days
deer = deer %>% 
  mutate(CONCEPT_DATE_JUL = format(CONCEPT_DATE, "%j") %>% 
           as.numeric())

# View Conception Days
hist_CONCEPT_DATE_JUL =
  ggplot(data = deer, mapping = aes(x = CONCEPT_DATE_JUL))+
  geom_histogram(col = "black")+
  xlab("Conception Day (Julian Days)")+
  ylab("Count")+
  ggtitle("Conception Days (Raw)")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, unit = "cm"))
hist_CONCEPT_DATE_JUL

# Check
# Start 01-01 = 1
min(deer$CONCEPT_DATE_JUL) #1
deer$CONCEPT_DATE[deer$CONCEPT_DATE_JUL == min(deer$CONCEPT_DATE_JUL)] #01-01
# End 12-31 = 366
max(deer$CONCEPT_DATE_JUL) #366
deer$CONCEPT_DATE[deer$CONCEPT_DATE_JUL == max(deer$CONCEPT_DATE_JUL)] #12-31
# Upper min 11-02 = 306 
min(deer$CONCEPT_DATE_JUL[deer$CONCEPT_DATE_JUL > 300]) #306
deer$CONCEPT_DATE[deer$CONCEPT_DATE_JUL == min(deer$CONCEPT_DATE_JUL[deer$CONCEPT_DATE_JUL > 300])] #11-02
# Lower max 03-01 = 60
max(deer$CONCEPT_DATE_JUL[deer$CONCEPT_DATE_JUL < 300]) #60 
deer$CONCEPT_DATE[deer$CONCEPT_DATE_JUL == max(deer$CONCEPT_DATE_JUL[deer$CONCEPT_DATE_JUL < 300])] #03-01

# Test JDate Conversion
df = data.frame(original = deer$CONCEPT_DATE,
                jdate = as.numeric(format(deer$CONCEPT_DATE, "%j")))

# If on larger end of scale subtract max and add max of low end + 1
df$jdate_adj = ifelse(df$jdate>300, df$jdate-366+61, df$jdate+61)

# Check
df[df$jdate_adj == min(df$jdate_adj),] # 1 (jday_sc), 306 (jday), 11-02
df[df$jdate_adj == max(df$jdate_adj),] # 121 (jday_sc), 60 (jday), 03-01
df[df$jdate == 1,] # 62 (jday_sc), 1 (jday), 01-01
df[df$jdate == 366,] #61 (jday_sc), 366 (jday) 12-31

# Time is circular, put all on "linear plane" so can do math
# Add +61 to all, if above 300 (max of lower is 60) subtract max 366 so that 
# 12-31 (366) and 1-1 (1) are next to each other (61, 62)
deer$CONCEPT_DATE_JUL = ifelse(deer$CONCEPT_DATE_JUL > 300, 
                                   deer$CONCEPT_DATE_JUL - 366 + 61,
                                   deer$CONCEPT_DATE_JUL + 61)

# View Conception Days linearized
hist_CONCEPT_DATE_JUL =
  ggplot(data = deer, mapping = aes(x = CONCEPT_DATE_JUL))+
  geom_histogram(col = "black")+
  xlab("Conception Day (Adjusted Julian Days)")+
  ylab("Count")+
  ggtitle("Conception Days (Adjusted Julian Days)")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, unit = "cm"))
hist_CONCEPT_DATE_JUL
# JDAY 1 (11-02), JDAY 121 (03-01)

# MEAN LINERIZED JDATE
mean_jdate_region = deer %>% 
  group_by(ABBREV) %>% 
  summarize(Lin_JDay_mean = round(mean(CONCEPT_DATE_JUL))) %>% 
  ungroup()

# Convert back by -61 and if negative, + 366
mean_jdate_region$JDay_mean = ifelse(mean_jdate_region$Lin_JDay_mean - 61 <= 0, 
                                     mean_jdate_region$Lin_JDay_mean - 61 + 366,
                                     mean_jdate_region$Lin_JDay_mean - 61)
# Convert to Date (approximate assuming 1 = 1-1 and 366 = 12-31)
mean_jdate_region$Date = c("12-31", "1-8", "1-6", "12-22", "12-24")

# Calculate mean breeding date per region
deer = deer %>% 
  group_by(ABBREV) %>% 
  mutate(JDAY_SC = CONCEPT_DATE_JUL - mean(CONCEPT_DATE_JUL)) %>% 
  ungroup()
# Sample size per region
tapply(deer$CONCEPT_DATE_JUL, deer$ABBREV, length)

# View conception days scaled by region
hist_JDAY_SC =
  ggplot(data = deer, mapping = aes(x = JDAY_SC))+
  geom_histogram(col = "black")+
  xlab("Relative Conception Day")+
  ylab("Count")+
  ggtitle("Relative Conception Days")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = 0.7, unit = "cm"))
hist_JDAY_SC

#############################################################
# Question 1: How do sex-ratios change across MS over time? #
#############################################################
# Aggregate by time (year)
yearly = deer %>% 
  dplyr::group_by(YEAR) %>% 
  dplyr::summarize(
         FETUS_M = sum(FETUS_M),
         FETUS_F = sum(FETUS_F),
         FETUS_T = sum(FETUS_TOTAL))

# Calculate the sex ratio
yearly$SR = yearly$FETUS_M/yearly$FETUS_F

# Plot SR by Year
yearly_SR = ggplot(data = yearly, mapping = aes(x = YEAR, y = SR))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = mean(yearly$SR), linetype = 3, linewidth = 1.2)+
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1.2, col = 2)+
  xlim(1990,2001)+
  xlab("Year")+
  ylab("Sex Ratio (M/F)")+
  ggtitle("Mississippi Sex Ratios by Year")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = 0.7, unit = "cm"))
yearly_SR
# Does seem to be consistently greater than 1 but appears random

# Linear model for sex ratio over time
summary(lm(yearly$SR ~ yearly$YEAR))
# 0.02507, 0.142 p-value, 0.01556 se
# R2 = 0.1377
par(mfrow=c(2,2))
plot(lm(yearly$SR ~ yearly$YEAR))

# Aggregate by region
regionally = deer %>% 
  dplyr::group_by(ABBREV) %>% 
  dplyr::summarize(
    FETUS_M = sum(FETUS_M),
    FETUS_F = sum(FETUS_F),
    FETUS_TOTAL = sum(FETUS_TOTAL))

# Calculate the sex ratio
regionally$SR = regionally$FETUS_M/regionally$FETUS_F

#Plot SR by Region
regional_SR = ggplot(data = regionally, mapping = aes(x = ABBREV, y = SR))+
  geom_point()+
  geom_hline(yintercept = mean(regionally$SR), linetype = 3, linewidth = 1.2)+
  geom_hline(yintercept = 1, linetype = 3, linewidth = 1.2, col = 2)+
  xlab("Region")+
  ylab("Sex Ratio (M/F)")+
  ggtitle("Mississippi Sex Ratios by Region")+
  theme(plot.margin = margin(t = 1, r = 1.5, b = 1, l = .7, unit = "cm"))
regional_SR

library(Rmisc)
# Confidence interval for sex ratio
# By year
CI(yearly$SR) #1.3061
# Sample size per year
tapply(deer$YEAR, deer$YEAR, length)
# By region
CI(regionally$SR) #1.6321, lots of variation in sex ratio by region
# Sample size per region
tapply(deer$ABBREV, deer$ABBREV, length)
# Global sex ratio - regardless of year/region
sum(deer$FETUS_M)/sum(deer$FETUS_F) #1.34

################################################################
# Question 2: What is the relative contribution of body size,  #
# age and kidney fat to the probability of having a male fawn? #
################################################################

# Stratify by number of fawns
# No way to compare how 1 M compares to 2 F etc...

#############
# Singlets ##
#############
# Subset by number of offspring
singlets = deer %>% filter(FETUS_TOTAL == 1)

# Are variables correlated?
singlets %>% 
  dplyr::select(KFI_rel, AGE_rel, DRESS_WT_rel) %>% 
  cor() # Do not seem correlated

# Singlet Model
singlets_mod = glm(data = singlets, formula = FETUS_M ~ KFI_rel + AGE_rel + DRESS_WT_rel, family = binomial(link="logit"))
# View results
summary(singlets_mod)
# KFI_rel, -0.002949, negative and not precise
# AGE_rel, 0.007983, positive and not precise
# DRESS_WT_rel, 0.009624, positive and not precise
# View explanatory power of variables
library(fmsb)
NagelkerkeR2(singlets_mod) #0.008387994
# None of these variables are likely good predictors of whether a female will have a male

###############
## Twinsies ###
###############
# Subset by number of offspring
twinsies = deer %>% filter(FETUS_TOTAL == 2)

# Are variables correlated?
twinsies %>% 
  dplyr::select(KFI_rel, AGE_rel, DRESS_WT_rel) %>% 
  cor() # Do not seem correlated

# Ordinal logistic regression
library(MASS)
twinsies_mod = polr(as.factor(FETUS_M) ~ KFI_rel + AGE_rel + DRESS_WT_rel, data = twinsies, Hess=TRUE)
# View results
summary(twinsies_mod)
# KFI_rel, -0.001187, negative and not precise
# AGE_rel, -0.027420, negative and not precise
# DRESS_WT_rel, 0.003148, positive and not precise
# View confidence intervals
confint(twinsies_mod)
library("DescTools")
PseudoR2(twinsies_mod, "Nagelkerke") #0.001509573

##############
## Triplets ##
##############
# Subset by number of offpsring
triplets = deer %>% filter(FETUS_TOTAL == 3)

# Are variables correlated?
triplets %>% 
  dplyr::select(KFI_rel, AGE_rel, DRESS_WT_rel) %>% 
  cor() # Do not seem correlated

# Ordinal logistic regression
triplets_mod = polr(as.factor(FETUS_M) ~ KFI_rel + AGE_rel + DRESS_WT_rel, data=triplets, Hess=TRUE)
# View results
summary(triplets_mod)
# KFI_rel, 0.012363, positive and precise, although very close to 0
# AGE_rel, -0.249044, negative and not precise, although close to 0
# DRESS_WT_rel, 0.001269, positive and not precise, although close to 0
# View confidence intervals
confint(triplets_mod)
PseudoR2(triplets_mod, "Nagelkerke") #0.1217719 not so bad

#######################################################
# Question 3: Does body condition, age or kidney fat  #
# predict number of fawns?                            #
#######################################################
# Define as factor
deer$FETUS_TOTAL_FCT = as.factor(deer$FETUS_TOTAL)
# Ordinal logistic regression 
offspring_mod = polr(FETUS_TOTAL_FCT ~ KFI_rel + AGE_rel + DRESS_WT_rel, data=deer, Hess=TRUE)
# View results
summary(offspring_mod)
# KFI_rel, -0.0004183, negative and not precise, very close to 0
# AGE_rel, 0.1532373, positive and precise
# DRESS_WT_rel, 0.0689343, positive and precise
# View confidence intervals
confint(offspring_mod)
# Calculate R2
PseudoR2(offspring_mod, "Nagelkerke") #0.1946 not bad

########################################################
# Question 4: Does the probability of having a male    #
# change with when in the season the female conceived? #
########################################################
# expect more males towards end of breeding season
# alt: expect more males towards mean conception

#Define as factor
singlets$MALE = as.factor(singlets$FETUS_M)# later effects plot will not play nice with factor defined within the model
# model
singlets_concept = glm(data = singlets, formula = MALE ~ JDAY_SC, family = binomial(link="logit"))
# results
summary(singlets_concept)
# View explanatory power of variables
NagelkerkeR2(singlets_concept) #very very low

# define as factor
twinsies$MALE = as.factor(twinsies$FETUS_M) 
# ordinal logistic regression
twinsies_concept = polr(data = twinsies, formula = MALE ~ JDAY_SC, Hess=TRUE)
# results
summary(twinsies_concept)
# confidence intervals
confint(twinsies_concept)
# very small, negative but significant
# R2
PseudoR2(twinsies_concept, "Nagelkerke") #0.01990363

# define as factor
triplets$MALE = as.factor(triplets$FETUS_M)
# ordinal logistic regression
triplets_concept = polr(data = triplets, formula = MALE ~ JDAY_SC, Hess=TRUE)
# results
summary(triplets_concept)
# confidence intervals
confint(triplets_concept)

# not signif and close to 0
# R2
PseudoR2(triplets_concept, "Nagelkerke") #0.0255794

#################
# Figure Making #
#################

#########################
### Variable Boxplots ###
#########################
# Don't do with facet wrap/facet grid b/c of levels in Sex
library(gridExtra)
library(ggpubr)

# Age Boxplots
g1 = ggplot(data = singlets, mapping = aes(x = FETUS_SEX, y = AGE_rel))+
  geom_boxplot()+
  ggtitle("Single Fawns")+
  ylab("Relative Doe Age\n")+
  xlab("")+
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 14, face = "bold"))+
  ylim(c(-4, 12))
g1

g2 = ggplot(data = twinsies, mapping = aes(x = FETUS_SEX, y = AGE_rel))+
  geom_boxplot()+
  ggtitle("Twins")+
  xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 14, face = "bold"))+
  ylim(c(-4, 12))
g2

g3 = ggplot(data = triplets, mapping = aes(x = FETUS_SEX, y = AGE_rel))+
  geom_boxplot()+
  ggtitle("Triplets")+
  xlab("")+ylab("")+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 14, face = "bold"))+
  ylim(c(-4, 12))
g3

# Dressed Weight Boxplots
g4 = ggplot(data = singlets, mapping = aes(x = FETUS_SEX, y = DRESS_WT_rel))+
  geom_boxplot()+
  xlab("")+
  ylab("Relative Dressed Weight\n")+
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_blank())+
  ylim(c(-40, 60))
g4

g5 = ggplot(data = twinsies, mapping = aes(x = FETUS_SEX, y = DRESS_WT_rel))+
  geom_boxplot()+
  xlab("")+
  ylab("")+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_blank())+
  ylim(c(-40, 60))
g5

g6 = ggplot(data = triplets, mapping = aes(x = FETUS_SEX, y = DRESS_WT_rel))+
  geom_boxplot()+
  xlab("")+
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 12))+
  ylim(c(-40, 60))
g6

# Third is kidney fat
# Has x-axes
# No titles
g7 = ggplot(data = singlets, mapping = aes(x = FETUS_SEX, y = KFI_rel))+
  geom_boxplot()+
  xlab("\nSex")+
  ylab("Relative Kidney Fat \n")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))+
  ylim(c(-110, 430))
g7

g8 = ggplot(data = twinsies, mapping = aes(x = FETUS_SEX, y = KFI_rel))+
  geom_boxplot()+
  xlab("\nSex of Twins")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 12))+
  ylim(c(-110, 430))
g8

g9 = ggplot(data = triplets, mapping = aes(x = FETUS_SEX, y = KFI_rel))+
  geom_boxplot()+
  xlab("\nSex of Triplets")+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 12))+
  ylim(c(-110, 430))
g9

grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9, nrow = 3, ncol = 3)

# Histograms for scaled variables
grid.arrange(hist_AGE_rel, hist_DW_rel, hist_KFI_rel, hist_JDAY_SC, nrow = 2, ncol = 2,
             top=text_grob("\n\nHistograms of Scaled Variables",
                           size = 18, face = "bold"))

###############
## Offspring ##
###############
library(effects)
library(car)
library(splines)
library(MASS)

# Multi-variate plot (Age and Dressed Weight)
plot(effects::Effect(focal.predictors = c("AGE_rel", "DRESS_WT_rel"), offspring_mod),
     main = "Fetal Number by Relative Age and Dressed Weight",
     xlab = "Relative Age", ylab = "Probability of Fetal Number (1, 2 or 3)")

# Multi-variate plot (KFI and Julian Day)
plot(effects::Effect(focal.predictors = c("KFI_rel"), offspring_mod),
     main = "Fetal Number by Relative KFI",
     xlab = "Relative KFI", ylab = "Probability of Fetal Number (1, 2 or 3)")


################
## Julian Day ##
################

## Global ##
sum(deer$FETUS_M)/sum(deer$FETUS_F) #1.34
# Global probability for males
p_male = sum(deer$FETUS_M)/sum(deer$FETUS_TOTAL)
p_male #0.5718

## Singlets
# Global sex ratio for singlets (1.5026:1)
sum(singlets$FETUS_M)/sum(singlets$FETUS_F)

# Global probability of male for singlets ()
# 1 = p + q
prob_1m_1 = sum(singlets$FETUS_M)/sum(singlets$FETUS_TOTAL)
prob_1m_1 #0.6004
prob_1f_1 = 1 - prob_1m_1
# Check
prob_1m_1+prob_1f_1

## Twins ##
# Global sex ratio for twins (1.3207:1)
sum(twinsies$FETUS_M)/sum(twinsies$FETUS_F)

# Global probability of male for twinsies
# 1 = p^2 + 2pq + q^2
prob_2m_2 = sum(twinsies$FETUS_M[twinsies$FETUS_M==2])/sum(twinsies$FETUS_TOTAL) #0.3432
prob_1m_2 = 2 * sqrt(prob_2m_2)*(1-sqrt(prob_2m_2)) #0.4853
prob_0m_2 = (1-sqrt(prob_2m_2))^2 #0.1715
# check
prob_0m_2 + prob_1m_2 + prob_2m_2
# implied probability for single male given twins
prob_1m_d2 = sqrt(prob_2m_2) #0.5858

## Triplets ##
# Global sex ratio for triplets (1.13)
sum(triplets$FETUS_M)/sum(triplets$FETUS_F)
# Global sex prob for triplets (0.5294:1)
sum(triplets$FETUS_M)/sum(triplets$FETUS_TOTAL)

# Global probability of male for triplets
# 1 = p^3 + 3p^2q + 3pq^2 + q^3
prob_3m_3 = sum(triplets$FETUS_M[triplets$FETUS_M==3])/sum(triplets$FETUS_TOTAL) #0.1961
prob_2m_3 = 3 * prob_3m_3^(2/3) * (1 - prob_3m_3^(1/3)) #0.4243
prob_1m_3 = 3 * prob_3m_3^(1/3) * (1 - prob_3m_3^(1/3))^2 #0.3360
prob_0m_3 = sum(triplets$FETUS_F[triplets$FETUS_F==3])/sum(triplets$FETUS_TOTAL) #0.1176
# check
prob_3m_3 + prob_2m_3 + prob_1m_3 + prob_0m_3
# implied probability for a single male given triplets
prob_1m_d3 = (prob_3m_3)^(1/3) #0.581

# Sex Ratios
# Global (1.34), Singlets (1.50), Twins (1.32), Triplets (1.13)
# Sex ratios appear to become less male biased with increasing fetal number

# Probability of Male
# Global (0.57), Singlets (0.60), Twins (0.59), Triplets (0.58)
# Probability of having a male seems independent of fetal number

# singlets
p1 = plot(effects::Effect(focal.predictors = "JDAY_SC", singlets_concept), 
          xlab = "Relative Conception Day", 
          ylab = "Probability of Male Fetus",
          main = "Singlets: Sex and Relative Conception Day")

p2 = plot(effects::Effect(focal.predictors = "JDAY_SC", twinsies_concept),
       xlab = "Relative Conception Day", 
       ylab = "Probability of Male Fetuses (0, 1, 2)",
       main = "Twins: Male Fetuses and Relative Conception Day")

p3 = plot(effects::Effect(focal.predictors = "JDAY_SC", triplets_concept),
       xlab = "Relative Conception Day", 
       ylab = "Probability of Male Fetuses (0, 1, 2, 3)",
       main = "Triplets: Male Fetuses and Relative Conception Day")
grid.arrange(p1, p2, p3, nrow = 1)

# What is the baseline probability of a male fetus?
prob_m = mean(as.numeric(deer$FETUS_M)/as.numeric(deer$FETUS_TOTAL))
prob_m #0.576

