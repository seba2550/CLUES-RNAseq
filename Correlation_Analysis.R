# Required libraries
library(tidyverse)
library(nnet)
library(table1)
library(reshape2)

# Read in clinical data
clues_time1 <- read.csv("clues_time1all.csv")
clues_time3 <- read.csv("CLUES_Time3all.csv") 

# Read in the sample metadata from the RNA-seq
rnaseq_samples <- read.csv("samples.csv", header = F)

# Extract just the patient IDs from the RNA-seq data
tmp1 <- as.data.frame(str_split_fixed(rnaseq_samples$V1, "_", 2))
rnaseq_samples_mod <- as.data.frame(str_split_fixed(tmp1$V1, "-", 2))
rnaseq_samples_mod$V2 <- as.integer(rnaseq_samples_mod$V2)
rm(tmp1)
rm(rnaseq_samples)

samples_only <- rnaseq_samples_mod %>% distinct(V2)
write.table(samples_only, file = "RNAseq_samples.txt", quote = F, col.names = F, row.names = F) # Write the samples out to a txt file so that I can load it into AWS for RNA-seq sample processing

# Get the patients for which we have both RNAseq and clinical data (one row per patient)
ppl_w_rnaseq <- subset(clues_time1, subjectid %in% rnaseq_samples_mod$V2)
ppl_w_rnaseq_3 <- subset(clues_time3, subjectid %in% rnaseq_samples_mod$V2)

# Classify patients according to their activity level (answers to RAPA questionnaire)
ppl_w_rnaseq <- ppl_w_rnaseq %>% 
  mutate(exercising_classification = case_when((rapararely == 1 | rapalight == 1) & (rapamod == 0 & rapamodlt30 == 0 & rapaviglt20 == 0 & rapamodge30 == 0 & rapavigge20 == 0) ~ 0, # This is evaluated in order
                                               (rapararely == 0 | rapalight == 0) & (rapamod == 1 | rapamodlt30 == 1 | rapaviglt20 == 1 | rapamodge30 == 1 | rapavigge20) ~ 1))



ppl_w_rnaseq_subset <- ppl_w_rnaseq %>%
  select(female, subjectid, bmi, sledaiscore, phqscore, smokenow, age, steroralnow, steroral12mo, crfstroke, crfmi, crfhtn, crfdiabetes1, crfdiabetes2, crfmalig, crflupusneph, crfesrddial)

ppl_w_rnaseq_3_subset <- ppl_w_rnaseq_3 %>%
  select(subjectid, physactlevel, physactdays)

ppl_w_rnaseq_subset_followups <- left_join(ppl_w_rnaseq_subset, ppl_w_rnaseq_3_subset, by = "subjectid")



# ALternate way of classifying patients (based off the physactivitylevel column)
ppl_w_rnaseq_subset_followups <- ppl_w_rnaseq_subset_followups %>%
  mutate(physact_classification = case_when(physactlevel == 1 ~ 0,
                                            physactlevel == 2 | physactlevel == 3 ~ 1))

write.csv(ppl_w_rnaseq_subset_followups, "clues_time_rnaseq_physact.csv")

# Logistic regression to find clinical features correlated with physical activity classification
#glm1 <- glm(exercising_classification ~ 
            #bmi + sledaiscore + phqscore + smokenow + age + steroralnow + steroral12mo
            #+ crfstroke + crfmi + crfhtn + crfdiabetes1 + crfdiabetes2 + crfmalig + crflupusneph + crfesrddial,  #some of these variables are binary, might remove them for other models
            #data = ppl_w_rnaseq, family = "binomial")
#summary(glm1)

glm_all <- glm(physact_classification ~ 
            bmi + sledaiscore + phqscore + smokenow + age + steroralnow + steroral12mo
            + crfstroke + crfmi + crfhtn + crfdiabetes1 + crfdiabetes2 + crfmalig + crflupusneph + crfesrddial,
            data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_all)

glm_bmi_age <- glm(physact_classification ~
            age + bmi,  #Only two continuous variables, the easiest ones
            data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_bmi_age)

glm_continuous_vars <- glm(physact_classification ~ # All of the continuous variables
              age + bmi + sledaiscore + phqscore,
            data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_continuous_vars)

glm_bmi <- glm(physact_classification ~  #Single variable model
              bmi ,
            data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_bmi)

glm_steroral12mo <- glm(physact_classification ~  #Single variable model
              steroral12mo ,
            data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_steroral12mo)

# A few more models that Marina said I should do
glm_sledaiscore <- glm(physact_classification ~
                         sledaiscore,
                       data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_sledaiscore)

glm_phqscore <- glm(physact_classification ~
                      phqscore,
                    data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_phqscore)

glm_physactdays <- glm(physact_classification ~
                         physactdays,
                       data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_physactdays) # The mother of all sanity checks: physical activity days is significantly associated with physical activity levels

glm_sex <- glm(physact_classification ~
                 female,
               data = ppl_w_rnaseq_subset_followups, family = "binomial")
summary(glm_sex)

# Jackie tests
###ppl_w_rnaseq <- ppl_w_rnaseq %>% 
  mutate(exercisers = case_when(rapamod == 1 | rapalight == 1 | rapamodlt30 == 1 | rapaviglt20 == 1 | rapamodge30 == 1 | rapavigge20 ~ 1, # This is evaluated in order
                                rapamod == 0 | rapalight == 0 | rapamodlt30 == 0 | rapaviglt20 == 0 | rapamodge30 == 0 | rapavigge20 ~ 0,
                                TRUE ~ 3))
###ppl_w_rnaseq$exercisers[ppl_w_rnaseq$exercisers == 3] <- NA

###glmn <- glm(exercisers ~ 
              ###bmi + sledaiscore + phqscore + smokenow + age + steroralnow + steroral12mo
            ###+ crfstroke + crfmi + crfhtn + crfdiabetes1 + crfdiabetes2 + crfmalig + crflupusneph + crfesrddial,
            ###data = ppl_w_rnaseq, family = "binomial")
###summary(glmn)

###glmn2 <- glm(exercisers ~
               ###age + bmi + sledaiscore + phqscore,
             ###data = ppl_w_rnaseq, family = "binomial")
###summary(glmn2)

# Jackie suggested a table1 examination. This tells us ,for each column, the distribution for whether they are physically active or not
table1(~ bmi
       + age
       + sledaiscore 
       + phqscore
       + smokenow
       + steroralnow
       + steroral12mo
       + crfstroke
       + crfmi
       + crfhtn
       + crfdiabetes1
       + crfdiabetes2
       + crfmalig
       + crflupusneph
       + crfesrddial
       | physact_classification, data=ppl_w_rnaseq_subset_followups)



# Plot some of these features to compare their distributions
ppl_w_rnaseq_subset_followups <- ppl_w_rnaseq_subset_followups %>% drop_na(physact_classification)

# Age
ggplot(ppl_w_rnaseq_subset_followups, aes(x=as.factor(physact_classification), y=age, color = as.factor(physact_classification))) + 
  geom_violin(trim = F) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  xlab("Physical Activity Classification") +
  ylab("Age") +
  labs(color='Physical Activity Classification') + 
  stat_summary(fun=median, geom="point", shape=23, size=4, color = "red") 

# BMI
ggplot(ppl_w_rnaseq_subset_followups, aes(x=as.factor(physact_classification), y=bmi, color = as.factor(physact_classification))) + 
  geom_violin(trim = F) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  xlab("Physical Activity Classification") +
  ylab("BMI") +
  labs(color='Physical Activity Classification') + 
  stat_summary(fun=median, geom="point", shape=23, size=4, color = "red")

# SLEDAI score
ggplot(ppl_w_rnaseq_subset_followups, aes(x=as.factor(physact_classification), y=sledaiscore, color = as.factor(physact_classification))) + 
  geom_violin(trim = F) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  xlab("Physical Activity Classification") +
  ylab("SLEDAI Score") +
  labs(color='Physical Activity Classification') + 
  stat_summary(fun=median, geom="point", shape=23, size=2, color = "red") 

# PHQ Score
ggplot(ppl_w_rnaseq_subset_followups, aes(x=as.factor(physact_classification), y=phqscore, color = as.factor(physact_classification))) + 
  geom_violin(trim = F) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  xlab("Physical Activity Classification") +
  ylab("PHQ Score") +
  labs(color='Physical Activity Classification') + 
  stat_summary(fun=median, geom="point", shape=23, size=2, color = "red")

# Physical Activity Days
ggplot(ppl_w_rnaseq_subset_followups, aes(x=as.factor(physact_classification), y=physactdays, color = as.factor(physact_classification))) + 
  geom_violin(trim = F) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  xlab("Physical Activity Classification") +
  ylab("Days of Physical Activity (Time 3)") +
  labs(color='Physical Activity Classification') + 
  stat_summary(fun=median, geom="point", shape=23, size=2, color = "red")


# Test for enrichment of certain clinical features
chisq.test(ppl_w_rnaseq_subset_followups$physact_classification, ppl_w_rnaseq_subset_followups$steroral12mo) # Steroid use in the last 12 months
chisq.test(ppl_w_rnaseq_subset_followups$physact_classification, ppl_w_rnaseq_subset_followups$crfdiabetes1) # Warning on this and T2D. Probably because of missingness in data
chisq.test(ppl_w_rnaseq_subset_followups$physact_classification, ppl_w_rnaseq_subset_followups$crfdiabetes2)
chisq.test(ppl_w_rnaseq_subset_followups$physact_classification, ppl_w_rnaseq_subset_followups$crfstroke)
chisq.test(ppl_w_rnaseq_subset_followups$physact_classification, ppl_w_rnaseq_subset_followups$crfmalig)
chisq.test(ppl_w_rnaseq_subset_followups$physact_classification, ppl_w_rnaseq_subset_followups$smokenow)

# Which individuals have missing physical activity data?
missing_phys_act <- ppl_w_rnaseq_subset_followups[which(is.na(ppl_w_rnaseq_subset_followups$physact_classification)),]

melted_missing_phys_act <- melt(missing_phys_act, id = "subjectid")

melted_missing_phys_act_ages <- melted_missing_phys_act[which(melted_missing_phys_act$variable == "age"),]

melted_missing_phys_act_bmi <- melted_missing_phys_act[which(melted_missing_phys_act$variable == "bmi"),]

ggplot(melted_missing_phys_act_ages, aes(x = subjectid, y = value), fill = variable) +
  geom_col(width = 2) +
  ylab("Ages of Patients with Missing Physical Activity Data")

ggplot(melted_missing_phys_act_bmi, aes(x = subjectid, y = value), fill = variable) +
  geom_col(width = 2) +
  ylab("BMIs of Patients with Missing Physical Activity Data")


