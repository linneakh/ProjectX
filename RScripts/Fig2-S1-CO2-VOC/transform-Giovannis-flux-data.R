##Linnea Honeker
##3/20/23
##linneah@arizona.edu
##Objective: To write a script that will transform VOC data tables from Giovanni into useable data frames in R for downstream analysis

####load libraries####
library(dplyr)
library(tidyr)
library(stringr)



########C13/C12 data- corrected version from Giovanni###############
#Load 13C data
data.13C12C.P <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetate-pre-drought-13C-12C.csv") 
data.13C12C.D <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetate-drought-13C-12C.csv") 
data.acetone.13C12C.P <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetone-pre-drought-13C-12C.csv") 
data.acetone.13C12C.D <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetone-drought-13C-12C.csv") 
data.diacetyl.13C12C.P <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/diacetyl-pre-drought-13C-12C.csv") 
data.diacetyl.13C12C.D <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/diacetyl-drought-13C-12C.csv") 

#######ACETATE############
#process pre-drought dataset
# gather dataframe 
data.13.P.long <- data.13C12C.P %>%
  pivot_longer(cols = contains("P"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))


# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows

df.p.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.13.P.long)) {
  if (isTRUE(data.13.P.long[i,1] == data.13.P.long[(i+1),1]) == TRUE) 
  {df.p.0[i,] <- spread(data.13.P.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.p <- df.p.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.p.f <- df.p %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Pre-Drought")

### process drought dataset
# gather dataframe 
data.13.D.long <- data.13C12C.D %>%
  pivot_longer(cols = contains("S"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))

data.13.D.long$Label
# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows
df.d.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.13.D.long)) {
  if (isTRUE(data.13.D.long[i,1] == data.13.D.long[(i+1),1]) == TRUE) 
  {df.d.0[i,] <- spread(data.13.D.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.d <- df.d.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.d.f <- df.d %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  mutate_at("Label", str_replace, ".W", "_W") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Drought")

#combine pre-drought and drought tables together
df.acetate <- rbind(df.p.f, df.d.f)

#filter to betwee 0 and 48 hours, and only keep labeled fluxes
acetate_flux <- df.acetate %>%
  mutate(Time.h = Time * 24) %>%
  dplyr::select(-Time) %>%
  filter(Time.h < 48 &
           Time.h > 0) %>%
  filter(Pyr != "W")

acetate_flux$Time.h
acetate_flux$Pyr
acetate_flux$Condition
acetate_flux$Label

write.csv(df.acetate, "./Output/CO2-VOCs/acetate_flux-13C-12C-all.csv")
write.csv(acetate_flux, "./Output/CO2-VOCs/acetate_flux-13C-12C-0-48.csv")

######ACETONE####
# process pre-drought dataset
# gather dataframe 
data.acetone.13.P.long <- data.acetone.13C12C.P %>%
  pivot_longer(cols = contains("P"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))



# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows

df.acetone.p.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.acetone.13.P.long)) {
  if (isTRUE(data.acetone.13.P.long[i,1] == data.acetone.13.P.long[(i+1),1]) == TRUE) 
  {df.acetone.p.0[i,] <- spread(data.acetone.13.P.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.acetone.p <- df.acetone.p.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.acetone.p.f <- df.acetone.p %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Pre-Drought")

### process drought dataset
# gather dataframe 
data.acetone.13.D.long <- data.acetone.13C12C.D %>%
  pivot_longer(cols = contains("S"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))

data.acetone.13.D.long$Label
# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows
df.acetone.d.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.acetone.13.D.long)) {
  if (isTRUE(data.acetone.13.D.long[i,1] == data.acetone.13.D.long[(i+1),1]) == TRUE) 
  {df.acetone.d.0[i,] <- spread(data.acetone.13.D.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.acetone.d <- df.acetone.d.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.acetone.d.f <- df.acetone.d %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  mutate_at("Label", str_replace, ".W", "_W") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Drought")

#combine pre-drought and drought tables together
df.acetone <- rbind(df.acetone.p.f, df.acetone.d.f)

#filter to betwee 0 and 48 hours, and only keep labeled fluxes
acetone_flux <- df.acetone %>%
  mutate(Time.h = Time * 24) %>%
  dplyr::select(-Time) %>%
  filter(Time.h < 48 &
           Time.h > 0) %>%
  filter(Pyr != "W")

acetone_flux$Time.h
acetone_flux$Pyr
acetone_flux$Condition
acetone_flux$Label

write.csv(df.acetone, "./Output/CO2-VOCs/acetone_flux-13C-12C-all.csv")
write.csv(acetone_flux, "./Output/CO2-VOCs/acetone_flux-13C-12C-0-48.csv")

######DIACETYL####
# process pre-drought dataset
# gather dataframe 
data.diacetyl.13.P.long <- data.diacetyl.13C12C.P %>%
  pivot_longer(cols = contains("P"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))



# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows

df.diacetyl.p.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.diacetyl.13.P.long)) {
  if (isTRUE(data.diacetyl.13.P.long[i,1] == data.diacetyl.13.P.long[(i+1),1]) == TRUE) 
  {df.diacetyl.p.0[i,] <- spread(data.diacetyl.13.P.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.diacetyl.p <- df.diacetyl.p.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.diacetyl.p.f <- df.diacetyl.p %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Pre-Drought")

### process drought dataset
# gather dataframe 
data.diacetyl.13.D.long <- data.diacetyl.13C12C.D %>%
  pivot_longer(cols = contains("S"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))

data.diacetyl.13.D.long$Label
# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows
df.diacetyl.d.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.diacetyl.13.D.long)) {
  if (isTRUE(data.diacetyl.13.D.long[i,1] == data.diacetyl.13.D.long[(i+1),1]) == TRUE) 
  {df.diacetyl.d.0[i,] <- spread(data.diacetyl.13.D.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.diacetyl.d <- df.diacetyl.d.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.diacetyl.d.f <- df.diacetyl.d %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  mutate_at("Label", str_replace, ".W", "_W") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Drought")

#combine pre-drought and drought tables together
df.diacetyl <- rbind(df.diacetyl.p.f, df.diacetyl.d.f)

#filter to betwee 0 and 48 hours, and only keep labeled fluxes
diacetyl_flux <- df.diacetyl %>%
  mutate(Time.h = Time * 24) %>%
  dplyr::select(-Time) %>%
  filter(Time.h < 48 &
           Time.h > 0) %>%
  filter(Pyr != "W")

diacetyl_flux$Time.h
diacetyl_flux$Pyr
diacetyl_flux$Condition
diacetyl_flux$Label

write.csv(df.diacetyl, "./Output/CO2-VOCs/diacetyl_flux-13C-12C-all.csv")
write.csv(diacetyl_flux, "./Output/CO2-VOCs/diacetyl_flux-13C-12C-0-48.csv")

########C13 data- corrected version from Giovanni###############
#Load 13C data
data.13C.P <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetate-pre-drought-13C.csv") 
data.13C.D <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetate-drought-13C.csv") 
data.acetone.13C.P <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetone-pre-drought-13C.csv") 
data.acetone.13C.D <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetone-drought-13C.csv") 
data.diacetyl.13C.P <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/acetone-pre-drought-13C.csv") 
data.diacetyl.13C.D <- read.csv("./Data/CO2-VOCs/voc_raw/corrected_raw_voc_revision/diacetyl-drought-13C.csv") 

#####ACETATE######
# process pre-drought dataset
# gather dataframe 
data.13.P.long <- data.13C.P %>%
  pivot_longer(cols = contains("P"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))


# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows

df.p.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.13.P.long)) {
  if (isTRUE(data.13.P.long[i,1] == data.13.P.long[(i+1),1]) == TRUE) 
  {df.p.0[i,] <- spread(data.13.P.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.p <- df.p.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.p.f <- df.p %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Pre-Drought")

### process drought dataset
# gather dataframe 
data.13.D.long <- data.13C.D %>%
  pivot_longer(cols = contains("S"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))

data.13.D.long$Label
# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows
df.d.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.13.D.long)) {
  if (isTRUE(data.13.D.long[i,1] == data.13.D.long[(i+1),1]) == TRUE) 
  {df.d.0[i,] <- spread(data.13.D.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.d <- df.d.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.d.f <- df.d %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  mutate_at("Label", str_replace, ".W", "_W") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Drought")

#combine pre-drought and drought tables together
df.acetate <- rbind(df.p.f, df.d.f)

#filter to betwee 0 and 48 hours, and only keep labeled fluxes
acetate_flux <- df.acetate %>%
  mutate(Time.h = Time * 24) %>%
  dplyr::select(-Time) %>%
  filter(Time.h < 48 &
           Time.h > 0) %>%
  filter(Pyr != "W")

acetate_flux$Time.h
acetate_flux$Pyr
acetate_flux$Condition
acetate_flux$Label

write.csv(df.acetate, "./Output/CO2-VOCs/acetate_flux-13C-all.csv")
write.csv(acetate_flux, "./Output/CO2-VOCs/acetate_flux-13C-0-48.csv")

######ACETONE######
# process pre-drought dataset
# gather dataframe 
data.acetone.13.P.long <- data.acetone.13C.P %>%
  pivot_longer(cols = contains("P"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))



# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows

df.acetone.p.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.acetone.13.P.long)) {
  if (isTRUE(data.acetone.13.P.long[i,1] == data.acetone.13.P.long[(i+1),1]) == TRUE) 
  {df.acetone.p.0[i,] <- spread(data.acetone.13.P.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.acetone.p <- df.acetone.p.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.acetone.p.f <- df.acetone.p %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Pre-Drought")

### process drought dataset
# gather dataframe 
data.acetone.13.D.long <- data.acetone.13C.D %>%
  pivot_longer(cols = contains("S"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))

data.acetone.13.D.long$Label
# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows
df.acetone.d.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.acetone.13.D.long)) {
  if (isTRUE(data.acetone.13.D.long[i,1] == data.acetone.13.D.long[(i+1),1]) == TRUE) 
  {df.acetone.d.0[i,] <- spread(data.acetone.13.D.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.acetone.d <- df.acetone.d.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.acetone.d.f <- df.acetone.d %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  mutate_at("Label", str_replace, ".W", "_W") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Drought")

#combine pre-drought and drought tables together
df.acetone <- rbind(df.acetone.p.f, df.acetone.d.f)

#filter to betwee 0 and 48 hours, and only keep labeled fluxes
acetone_flux <- df.acetone %>%
  mutate(Time.h = Time * 24) %>%
  dplyr::select(-Time) %>%
  filter(Time.h < 48 &
           Time.h > 0) %>%
  filter(Pyr != "W")

acetone_flux$Time.h
acetone_flux$Pyr
acetone_flux$Condition
acetone_flux$Label

write.csv(df.acetone, "./Output/CO2-VOCs/acetone_flux-13C-all.csv")
write.csv(acetone_flux, "./Output/CO2-VOCs/acetone_flux-13C-0-48.csv")

######DIACETYL######
# process pre-drought dataset
# gather dataframe 
data.diacetyl.13.P.long <- data.diacetyl.13C.P %>%
  pivot_longer(cols = contains("P"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))



# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows

df.diacetyl.p.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.diacetyl.13.P.long)) {
  if (isTRUE(data.diacetyl.13.P.long[i,1] == data.diacetyl.13.P.long[(i+1),1]) == TRUE) 
  {df.diacetyl.p.0[i,] <- spread(data.diacetyl.13.P.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.diacetyl.p <- df.diacetyl.p.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.diacetyl.p.f <- df.diacetyl.p %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Pre-Drought")

### process drought dataset
# gather dataframe 
data.diacetyl.13.D.long <- data.diacetyl.13C.D %>%
  pivot_longer(cols = contains("S"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))

data.diacetyl.13.D.long$Label
# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows
df.diacetyl.d.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.diacetyl.13.D.long)) {
  if (isTRUE(data.diacetyl.13.D.long[i,1] == data.diacetyl.13.D.long[(i+1),1]) == TRUE) 
  {df.diacetyl.d.0[i,] <- spread(data.diacetyl.13.D.long[i:(i+1),], Time, Time_Flux)[1,]
  }
}

# erase rows with NA and name columns
df.diacetyl.d <- df.diacetyl.d.0 %>%
  dplyr::rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.diacetyl.d.f <- df.diacetyl.d %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  mutate_at("Label", str_replace, ".W", "_W") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Drought")

#combine pre-drought and drought tables together
df.diacetyl <- rbind(df.diacetyl.p.f, df.diacetyl.d.f)

#filter to betwee 0 and 48 hours, and only keep labeled fluxes
diacetyl_flux <- df.diacetyl %>%
  mutate(Time.h = Time * 24) %>%
  dplyr::select(-Time) %>%
  filter(Time.h < 48 &
           Time.h > 0) %>%
  filter(Pyr != "W")

diacetyl_flux$Time.h
diacetyl_flux$Pyr
diacetyl_flux$Condition
diacetyl_flux$Label

write.csv(df.diacetyl, "./Output/CO2-VOCs/diacetyl_flux-13C-all.csv")
write.csv(diacetyl_flux, "./Output/CO2-VOCs/diacetyl_flux-13C-0-48.csv")

