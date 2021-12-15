# Required libraries
library(tidyverse)
library(corrplot)

# Read in clinical data
clues_time1 <- read.csv("clues_time1all.csv")
clues_time3 <- read.csv("CLUES_Time3all.csv") # Zach said we probably won't touch this file for the project, but I should still clean it up

# Read in the sample metadata from the RNA-seq
rnaseq_samples <- read.csv("samples.csv", header = F)

# Extract just the patient IDs from the RNA-seq data
tmp1 <- as.data.frame(str_split_fixed(rnaseq_samples$V1, "_", 2))
rnaseq_samples_mod <- as.data.frame(str_split_fixed(tmp1$V1, "-", 2))
rnaseq_samples_mod$V2 <- as.integer(rnaseq_samples_mod$V2)
rm(tmp1)

# Join RNA-seq sample IDs with clinical data
clues_time1_rnaseq_join <- left_join(clues_time1, rnaseq_samples_mod, by = c("subjectid" = "V2")) # Number of rows increased because many patients had RNAseq of multiple cell types, but not all of them

clues_time3_rnaseq_join <- left_join(clues_time3, rnaseq_samples_mod, by = c("subjectid" = "V2")) # Same thing as the above

# Filter out individuals for which we don't have RNA-seq data
clues_time1_rnaseq_join_filter <- clues_time1_rnaseq_join %>% drop_na(V1) # 480/4 = 120 samples we're supposed to have! Yay!
clues_time3_rnaseq_join_filter <- clues_time3_rnaseq_join %>% drop_na(V1)

# Count patients per level of RAPA activity
rapa_light_patients <- clues_time1_rnaseq_join_filter %>% group_by(subjectid) %>% summarise(count = sum(rapalight == 1))
rapa_rarely_patients <- clues_time1_rnaseq_join_filter %>% group_by(subjectid) %>% summarise(count = sum(rapararely == 1))
rapa_moderate_patients <- clues_time1_rnaseq_join_filter %>% group_by(subjectid) %>% summarise(count = sum(rapamod == 1))

rapa_light_patients <- rapa_light_patients %>% drop_na(count)
rapa_rarely_patients <- rapa_rarely_patients %>% drop_na(count)
rapa_moderate_patients <- rapa_moderate_patients %>% drop_na(count)

rapa_light_patients <- rapa_light_patients %>% subset(.$count != 0)
rapa_rarely_patients <- rapa_rarely_patients %>% subset(.$count != 0)
rapa_moderate_patients <- rapa_moderate_patients %>% subset(.$count != 0)


# Only one row per patient, easier for data wrangling. Thanks Jackie!
ppl_w_rnaseq <- subset(clues_time1, subjectid %in% clues_time1_rnaseq_join_filter$subjectid)
sum(ppl_w_rnaseq$rapararely, na.rm = T) 

ppl_w_rnaseq_time3 <- subset(clues_time3, subjectid %in% clues_time3_rnaseq_join_filter$subjectid)
sum(ppl_w_rnaseq_time3$rapararely)

ppl_w_rnaseq <- ppl_w_rnaseq %>% 
  mutate(exercisers = case_when(rapamod == 1 | rapalight == 1 | rapamodlt30 == 1 | rapaviglt20 == 1 | rapamodge30 == 1 | rapavigge20 ~ 1, # This is evaluated in order
                                rapamod == 0 | rapalight == 0 | rapamodlt30 == 0 | rapaviglt20 == 0 | rapamodge30 == 0 | rapavigge20 ~ 0,
                                TRUE ~ 3))

ppl_w_rnaseq$exercisers[ppl_w_rnaseq$exercisers == 3] <- NA

table(ppl_w_rnaseq$rapararely, ppl_w_rnaseq$exercisers)
table(ppl_w_rnaseq$rapararely, ppl_w_rnaseq$exercisers, useNA = "always")


# Plot distributions of some of the features
ggplot(ppl_w_rnaseq, aes(x = age)) +
  geom_histogram(aes(y = ..density..), 
                  binwidth = 1 ,colour = "blue", fill = "white") +
  geom_density(alpha = .2, fill="#FF6655") +
  geom_vline(aes(xintercept = mean(age, na.rm = T)),
             colour = "red", linetype ="longdash", size = .8) +
  theme_minimal()

ggplot(ppl_w_rnaseq, aes(x = sledaiscore)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 1 ,colour = "blue", fill = "white") +
  geom_density(alpha = .2, fill="#FF6655") +
  geom_vline(aes(xintercept = mean(sledaiscore, na.rm = T)),
             colour = "red", linetype ="longdash", size = .8) +
  theme_minimal()

ggplot(ppl_w_rnaseq, aes(x = acrcsum)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 1 ,colour = "blue", fill = "white") +
  geom_density(alpha = .2, fill="#FF6655") +
  geom_vline(aes(xintercept = mean(acrcsum, na.rm = T)),
             colour = "red", linetype ="longdash", size = .8) +
  theme_minimal()

ggplot(ppl_w_rnaseq, aes(x = agedx)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 1 ,colour = "blue", fill = "white") +
  geom_density(alpha = .2, fill="#FF6655") +
  geom_vline(aes(xintercept = mean(agedx, na.rm = T)),
             colour = "red", linetype ="longdash", size = .8) +
  theme_minimal()


# Logistic regression models to determine if certain features are associated with physical activity outcomes
glm1 <- glm(rapararely ~ age + sledaiscore , data = ppl_w_rnaseq, family = "binomial")
summary(glm1)

glm2 <- glm(exercisers ~ age + sledaiscore , data = ppl_w_rnaseq, family = "binomial")
summary(glm2)

glm3 <- glm(exercisers ~ age + sledaiscore + acrcsum , data = ppl_w_rnaseq, family = "binomial")
summary(glm3)

glm4 <- glm(rapararely ~ sf36mcs + sf36menthlth , data = ppl_w_rnaseq, family = "binomial")
summary(glm4)


# Correlation plots
ppl_w_rnaseq_subset <- ppl_w_rnaseq %>% select(subjectid, rapararely, rapalight, rapamod, age, sledaiscore, acrcsum, exercisers)

corr_mat <-cor(ppl_w_rnaseq_subset, method = "spearman")

corrplot(corr_mat, method="number")



# Distribution of the RAPA activity categorizations
ftable(ppl_w_rnaseq$rapararely == 1, ppl_w_rnaseq$rapalight == 1) # 11 patients in rarely, but 9 of them are also in light or other categories
ftable(ppl_w_rnaseq$rapalight == 1, ppl_w_rnaseq$rapavigge20 == 0) # 56 patients in light, but not vigorous
ftable(ppl_w_rnaseq$rapalight == 1, ppl_w_rnaseq$rapamod == 0) # 50 patients in light, but not moderate
ftable(ppl_w_rnaseq$rapararely == 1, ppl_w_rnaseq$rapamod == 0)
ftable(ppl_w_rnaseq$rapararely == 1, ppl_w_rnaseq$rapavigge20 == 0)

length(which(ppl_w_rnaseq$rapararely == 1 | ppl_w_rnaseq$rapalight == 1) & which(ppl_w_rnaseq$rapamod == 0 & ppl_w_rnaseq$rapavigge20 == 0)) # "Under-exercisers"
length(which(ppl_w_rnaseq$rapamod == 1 | ppl_w_rnaseq$rapavigge20 == 1) & which(ppl_w_rnaseq$rapararely == 0 & ppl_w_rnaseq$rapalight == 0))


