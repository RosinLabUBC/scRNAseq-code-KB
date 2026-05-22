
library(ARTool)
library(emmeans)


# Dlx2 puncta frequency data frame ----------------------------------------
# male Dlx2 data
m_con_data <- c(
  12.77258567,	14.75409836,	2.605863192,
  29.90654206,	34.42622951,	7.166123779,
  16.19937695,	14.3442623,	6.51465798,
  14.01869159,	12.70491803,	9.771986971,
  9.657320872,	9.836065574,	9.771986971,
  6.542056075,	5.327868852,	14.33224756,
  4.049844237,	4.098360656,	13.02931596,
  3.738317757,	0.819672131,	8.469055375,
  0.934579439,	1.639344262,	10.09771987,
  1.246105919,	1.229508197,	4.234527687,
  0.934579439,	0.409836066,	5.863192182,
  0,	0,	2.280130293,
  0,	0.409836066,	1.954397394,
  0,	0,	1.954397394,
  0,	0,	0.977198697,
  0,	0,	0.651465798,
  0,	0,	0,
  0,	0,	0.325732899,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0)

m_cold_data <- c(
  8.813559322,	13.46153846,	2.762430939,
  29.15254237,	33.07692308,	17.12707182,
  21.3559322,	14.61538462,	16.57458564,
  15.25423729,	11.92307692,	16.29834254,
  11.86440678,	11.15384615,	11.60220994,
  8.13559322,	8.076923077,	9.944751381,
  2.033898305,	5,	7.182320442,
  1.694915254,	1.538461538,	7.73480663,
  0,	0.384615385,	5.248618785,
  0.677966102,	0.384615385,	2.209944751,
  0.677966102,	0.384615385,	1.657458564,
  0.338983051,	0,	0.828729282,
  0,	0,	0.276243094,
  0,	0,	0.276243094,
  0,	0,	0.276243094,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0)

# combine data frames
all_m_scores <- c(m_con_data, m_cold_data)

bin_names <- c("0", "1-5", "6-10", "11-15", "16-20", "21-25", "26-30", 
              "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", 
              "61-65", "66-70", "71-75", "76-80", "81-85", "86-90", 
              "91-95", "96-100")

bin_labels <- rep(bin_names, each = 3, times = 2)

treatment_labels <- rep(c("Control", "Cold"), each = 63)

df_male_dlx2 <- data.frame(
  Response = all_m_scores,
  Treatment = factor(treatment_labels),
  Bin = factor(bin_labels, levels = bin_names))

saveRDS(df_male_dlx2, '/Users/kodeebao/Documents/ART ANOVA/df_male_dlx2.rds')


# female Dlx2 data + n=4 
f_con_data <- c(
  82.38095238,	40.85106383,	16.796875,	40.78431373,
  11.42857143,	40.42553191,	44.140625,	38.43137255,
  3.80952381,	10.63829787,	16.015625,	13.33333333,
  1.904761905,	3.829787234,	10.546875,	3.529411765,
  0.476190476,	2.127659574,	5.078125,	1.568627451,
  0,	1.70212766,	3.515625,	0.784313725,
  0,	0,	1.953125,	0.392156863,
  0,	0.425531915,	1.5625,	0.784313725,
  0,	0,	0.390625,	0.392156863,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0)

f_cold_data <- c(
  8.4375,	38.11881188,	4.273504274,	4.658385093,
  32.8125,	37.12871287,	12.82051282,	15.2173913,
  12.1875,	8.415841584,	15.0997151,	10.86956522,
  11.5625,	4.95049505,	12.82051282,	14.9068323,
  11.25,	2.475247525,	13.39031339,	11.49068323,
  9.375,	3.96039604,	11.96581197,	8.695652174,
  5,	2.475247525,	9.686609687,	8.695652174,
  4.375,	0.495049505,	5.698005698,	8.695652174,
  2.5,	0.495049505,	3.703703704,	4.658385093,
  1.25,	0.99009901,	2.849002849,	4.658385093,
  0.9375,	0,	1.994301994,	2.173913043,
  0,	0.495049505,	2.564102564,	1.863354037,
  0.3125,	0,	1.13960114,	0.310559006,
  0,	0,	0.56980057,	1.552795031,
  0,	0,	0.56980057,	0.931677019,
  0,	0,	0.284900285,	0,
  0,	0,	0.284900285,	0.310559006,
  0,	0,	0.284900285,	0.310559006,
  0,	0,	0,	0,
  0,	0,	0,	0,
  0,	0,	0,	0)

# combine data frames
all_f_scores <- c(f_con_data, f_cold_data)

bin_names <- c("0", "1-5", "6-10", "11-15", "16-20", "21-25", "26-30", 
               "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", 
               "61-65", "66-70", "71-75", "76-80", "81-85", "86-90", 
               "91-95", "96-100")

bin_labels <- rep(bin_names, each = 4, times = 2)

treatment_labels <- rep(c("Control", "Cold"), each = 84)

df_female_dlx2 <- data.frame(
  Response = all_f_scores,
  Treatment = factor(treatment_labels),
  Bin = factor(bin_labels, levels = bin_names))

saveRDS(df_female_dlx2, '/Users/kodeebao/Documents/ART ANOVA/df_female_dlx2.rds')


# loading in data frames --------------------------------------------------


df_male_dlx2 <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_male_dlx2.rds")
df_female_dlx2 <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_female_dlx2.rds")


# ART on male dlx2 ------------------------------------------------------
# fit ART model
male_dlx2 <- art(Response ~ Treatment * Bin, data = df_male_dlx2)

# verify ART appropriateness
summary(male_dlx2)

# run the ART ANOVA
anova_results_male_dlx2 <- anova(male_dlx2)

write.csv(as.data.frame(anova_results_male_dlx2),
          "/Users/kodeebao/Documents/ART ANOVA/results dlx2/male_dlx2_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
male_dlx2_lm_treatment <- artlm(male_dlx2, "Treatment")
male_dlx2_lm_bin <- artlm(male_dlx2, "Bin")

# treatment differences within each Bin (Tukey adjustment)
treatment_posthoc_male_dlx2 <- emmeans(male_dlx2_lm_treatment,
                                          pairwise ~ Treatment | Bin,
                                          adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_male_dlx2$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results dlx2/male_dlx2_posthoc_treatment(tukey).csv",
          row.names = FALSE)

# puncta per cell bin differences within each Treatment (Tukey adjustment)
bin_posthoc_male_dlx2 <- emmeans(male_dlx2_lm_bin,
                                     pairwise ~ Bin | Treatment,
                                     adjust = "tukey")

write.csv(as.data.frame(bin_posthoc_male_dlx2$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results dlx2/male_dlx2_posthoc_bin(tukey).csv",
          row.names = FALSE)


# ART on female dlx2 ------------------------------------------------------
# fit ART model
female_dlx2 <- art(Response ~ Treatment * Bin, data = df_female_dlx2)

# verify ART appropriateness
summary(female_dlx2)

# run the ART ANOVA
anova_results_female_dlx2 <- anova(female_dlx2)

write.csv(as.data.frame(anova_results_female_dlx2),
          "/Users/kodeebao/Documents/ART ANOVA/results dlx2/female_dlx2_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
female_dlx2_lm_treatment <- artlm(female_dlx2, "Treatment")
female_dlx2_lm_bin <- artlm(female_dlx2, "Bin")

# treatment differences within each Bin (Tukey adjustment)
treatment_posthoc_female_dlx2 <- emmeans(female_dlx2_lm_treatment,
                                       pairwise ~ Treatment | Bin,
                                       adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_female_dlx2$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results dlx2/female_dlx2_posthoc_treatment(tukey).csv",
          row.names = FALSE)

# puncta per cell bin differences within each Treatment (Tukey adjustment)
bin_posthoc_female_dlx2 <- emmeans(female_dlx2_lm_bin,
                                 pairwise ~ Bin | Treatment,
                                 adjust = "tukey")

write.csv(as.data.frame(bin_posthoc_female_dlx2$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results dlx2/female_dlx2_posthoc_bin(tukey).csv",
          row.names = FALSE)


