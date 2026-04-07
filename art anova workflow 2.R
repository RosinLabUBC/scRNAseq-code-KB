
library(ARTool)
library(emmeans)


# primary sphere data frame -----------------------------------------------
# male primary spheres
m_con_45 <- c(17.5, 18.5, 32, 24.5, 40.5, 38.5, 43.5, 25.5)
m_con_100 <- c(14, 17, 23, 18.5, 24, 29, 17, 16)
m_con_200 <- c(2.5, 6.5, 12, 4.5, 11, 10, 13.5, 6.5)

m_cold_45 <- c(37.5, 29, 25.5, 13, 32.5, 35, 24.5, 46.5, 25, 25.5, 26, 25,	34,	37, 43.5)
m_cold_100 <- c(16.5,	12,	11,	19,	24.5,	24,	15,	20.5,		12.5,	19,	16,	13,	16.5,	12.5,	16)
m_cold_200 <- c(4,	7.5,	30,	12.5,	20,	11.5,	9,	1.5,	7,	10.5,	5.5,	10,	10,	3,	11)

# female primary spheres
f_con_45 <- c(17,	20,	24.5,	20.5,	26,	25.5,	23.5,	13.5,	35.5,	39,	27.5,	35,	25)
f_con_100 <- c(7.5,	22.5,	16.5,	14,	15.5,	8,	18.5,	9,	33,	21.5,	12,	16.5,	18)
f_con_200 <- c(3.5,	7,	8.5,	3,	6,	1,	7,	0,	4.5,	6.5,	8.5,	11.5,	5.5)

f_cold_45 <- c(15.5,	23.5,	33.5,	29.5,	25.5,	41,	32,	36.5,	26)
f_cold_100 <- c(12.5,	12.5,	16.5,	16.5,	21,	22.5,	15.5,	20.5,	14.5)
f_cold_200 <- c(10,	5,	13.5,	15.5,	7,	5,	4,	9.5,	6.5)


# combine male numbers
all_male_scores <- c(m_con_45, m_con_100, m_con_200, m_cold_45, m_cold_100, m_cold_200)

all_female_scores <- c(f_con_45, f_con_100, f_con_200, f_cold_45, f_cold_100, f_cold_200)


# create treatment labels
male_treatment_labels <- c(rep("Control", 24), rep("Cold", 45))

female_treatment_labels <- c(rep("Control", 39), rep("Cold", 27))


# create treatment labels
male_size_labels <- c(rep("45-100", 8), rep("100-200", 8), rep(">200", 8),
                      rep("45-100", 15), rep("100-200", 15), rep(">200", 15))

female_size_labels <- c(rep("45-100", 13), rep("100-200", 13), rep(">200", 13),
                      rep("45-100", 9), rep("100-200", 9), rep(">200", 9))


# make data frame
df_male_primary_spheres <- data.frame(
  Response = all_male_scores,
  Treatment = factor(male_treatment_labels),
  Size = factor(male_size_labels, levels = c("45-100", "100-200", ">200")))

df_female_primary_spheres <- data.frame(
  Response = all_female_scores,
  Treatment = factor(female_treatment_labels),
  Size = factor(female_size_labels, levels = c("45-100", "100-200", ">200")))


# save data frame
saveRDS(df_male_primary_spheres, '/Users/kodeebao/Documents/ART ANOVA/df_male_primary_sphere.rds')

saveRDS(df_female_primary_spheres, '/Users/kodeebao/Documents/ART ANOVA/df_female_primary_sphere.rds')


# secondary sphere data frame -----------------------------------------------
# male secondary spheres
m_con_45 <- c(2.5,	4,	12,	3.5,	16.5,	3.5,	9)
m_con_100 <- c(1,	3,	9,	6,	17.5,	4.5,	6.5)
m_con_200 <- c(4,	23,	26.5,	15,	18,	16,	11.5)

m_cold_45 <- c(3,	3,	7.5,	2,	4,	9,	14.5,	19.5,	4,	4)
m_cold_100 <- c(6,	5.5,	15,	13.5,	4,	10.5,	13,	24.5,	8,	9)
m_cold_200 <- c(12,	24.5,	20.5,	22.5,	27.5,	29,	20,	16,	16,	23)

# female secondary spheres
f_con_45 <- c(0.5,	3,	2.5,	3,	4,	10.5,	5,	3.5,	8,	4.5,	13.5,	14)
f_con_100 <- c(1,	2.5,	2,	3,	11,	16,	8.5,	7,	7.5,	6,	15,	8)
f_con_200 <- c(11.5,	4,	6.5,	1,	12.5,	29.5,	14,	17.5,	14,	17,	11.5,	17.5)

f_cold_45 <- c(4.5,	12.5,	7.5,	3.5,	3.5,	1.5,	6,	8,	9,	4)
f_cold_100 <- c(2.5,	6,	3.5,	2,	11,	2.5,	12,	3,	10,	7.5)
f_cold_200 <- c(14,	20,	13.5,	7.5,	18.5,	19,	29,	1.5,	11,	17.5)


# combine male numbers
all_male_scores <- c(m_con_45, m_con_100, m_con_200, m_cold_45, m_cold_100, m_cold_200)

all_female_scores <- c(f_con_45, f_con_100, f_con_200, f_cold_45, f_cold_100, f_cold_200)


# create treatment labels
male_treatment_labels <- c(rep("Control", 21), rep("Cold", 30))

female_treatment_labels <- c(rep("Control", 36), rep("Cold", 30))


# create treatment labels
male_size_labels <- c(rep("45-100", 7), rep("100-200", 7), rep(">200", 7),
                      rep("45-100", 10), rep("100-200", 10), rep(">200", 10))

female_size_labels <- c(rep("45-100", 12), rep("100-200", 12), rep(">200", 12),
                        rep("45-100", 10), rep("100-200", 10), rep(">200", 10))


# make data frame
df_male_secondary_spheres <- data.frame(
  Response = all_male_scores,
  Treatment = factor(male_treatment_labels),
  Size = factor(male_size_labels, levels = c("45-100", "100-200", ">200")))

df_female_secondary_spheres <- data.frame(
  Response = all_female_scores,
  Treatment = factor(female_treatment_labels),
  Size = factor(female_size_labels, levels = c("45-100", "100-200", ">200")))


# save data frame
saveRDS(df_male_secondary_spheres, '/Users/kodeebao/Documents/ART ANOVA/df_male_secondary_sphere.rds')

saveRDS(df_female_secondary_spheres, '/Users/kodeebao/Documents/ART ANOVA/df_female_secondary_sphere.rds')



# Dlx5 puncta frequency data frame ----------------------------------------
# male Dlx5 data
m_con_data <- c(
  0.71428571,	1.03092784,	0.79365079,
  15,	15.1202749,	9.12698413,
  21.0714286,	22.6804124,	16.2698413,
  15,	21.3058419,	19.047619,
  18.2142857,	17.5257732,	18.2539683,
  8.21428571,	7.21649485,	14.2857143,
  7.14285714,	6.52920962,	8.33333333,
  5.35714286,	4.46735395,	7.14285714,
  5,	1.71821306,	3.17460317,
  2.14285714,	1.03092784,	1.98412698,
  1.42857143,	0.34364261,	1.19047619,
  0,	0.34364261,	0.3968254,
  0.35714286,	0,	0,
  0,	0,	0,
  0.35714286,	0,	0,
  0,	0.34364261,	0,
  0,	0.34364261,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0)

m_cold_data <- c(
  2.30263158,	1.17647059,	1.42857143,
  13.1578947,	19.6078431,	10.7142857,
  12.1710526,	27.4509804,	17.8571429,
  13.1578947,	21.1764706,	22.5,
  16.4473684,	13.7254902,	15.7142857,
  13.1578947,	9.01960784,	15.7142857,
  10.1973684,	3.52941176,	8.21428571,
  7.56578947,	2.74509804,	4.64285714,
  4.27631579,	1.56862745,	2.85714286,
  3.94736842,	0,	0,
  1.97368421,	0,	0.35714286,
  0.65789474,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0.65789474,	0,	0,
  0.32894737,	0,	0,
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

df_male_dlx5 <- data.frame(
  Response = all_m_scores,
  Treatment = factor(treatment_labels),
  Bin = factor(bin_labels, levels = bin_names))

saveRDS(df_male_dlx5, '/Users/kodeebao/Documents/ART ANOVA/df_male_dlx5.rds')


# female Dlx5 data
f_con_data <- c(
  1.24481328,	9.42408377,	1.43369176,
  45.6431535,	49.7382199,	16.4874552,
  24.8962656,	21.9895288,	28.6738351,
  19.9170124,	8.90052356,	21.8637993,
  5.39419087,	4.71204188,	17.2043011,
  1.65975104,	2.61780105,	6.09318996,
  0.82987552,	1.04712042,	4.65949821,
  0,	0.52356021,	2.15053763,
  0,	0.52356021,	0.71684588,
  0,	0,	0.71684588,
  0.41493776,	0.52356021,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0)

f_cold_data <- c(
  0.55096419,	12.5,	0,
  7.98898072,	39.5833333,	10.0271003,
  18.4573003,	18.4027778,	18.4281843,
  17.630854,	10.0694444,	19.2411924,
  15.7024793,	7.29166667,	19.5121951,
  12.3966942,	5.90277778,	10.298103,
  7.71349862,	1.73611111,	8.67208672,
  7.16253444,	2.43055556,	5.14905149,
  4.4077135,	1.04166667,	3.25203252,
  4.1322314,	0.69444444,	2.7100271,
  0.82644628,	0.34722222,	1.08401084,
  1.65289256,	0,	0.27100271,
  0.82644628,	0,	0,
  0.27548209,	0,	0.54200542,
  0,	0,	0.54200542,
  0.27548209,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0.27100271,
  0,	0,	0)

# combine data frames
all_f_scores <- c(f_con_data, f_cold_data)

bin_names <- c("0", "1-5", "6-10", "11-15", "16-20", "21-25", "26-30", 
               "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", 
               "61-65", "66-70", "71-75", "76-80", "81-85", "86-90", 
               "91-95", "96-100")

bin_labels <- rep(bin_names, each = 3, times = 2)

treatment_labels <- rep(c("Control", "Cold"), each = 63)

df_female_dlx5 <- data.frame(
  Response = all_f_scores,
  Treatment = factor(treatment_labels),
  Bin = factor(bin_labels, levels = bin_names))

saveRDS(df_female_dlx5, '/Users/kodeebao/Documents/ART ANOVA/df_female_dlx5.rds')


# Epha4 puncta frequency data frame ----------------------------------------
# male Epha4 data
m_con_data <- c(
  21.6216216,	13.2930514,	4.08163265,
  52.027027,	42.9003021,	33.4693878,
  17.9054054,	28.7009063,	20.4081633,
  5.06756757,	10.2719033,	13.4693878,
  2.7027027,	2.11480363,	13.0612245,
  0.67567568,	2.71903323,	8.16326531,
  0,	0,	4.48979592,
  0,	0,	2.04081633,
  0,	0,	0.81632653,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0)

m_cold_data <- c(
  6.43274854,	10.4615385,	3.85674931,
  33.0409357,	40.9230769,	33.6088154,
  36.2573099,	26.7692308,	28.9256198,
  15.7894737,	10.7692308,	17.0798898,
  4.67836257,	6.15384615,	8.815427,
  1.4619883,	3.69230769,	4.4077135,
  1.16959064,	0.92307692,	1.65289256,
  0.58479532,	0.30769231,	1.37741047,
  0.29239766,	0,	0,
  0.29239766,	0,	0,
  0,	0,	0.27548209,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
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

df_male_epha4 <- data.frame(
  Response = all_m_scores,
  Treatment = factor(treatment_labels),
  Bin = factor(bin_labels, levels = bin_names))

saveRDS(df_male_epha4, '/Users/kodeebao/Documents/ART ANOVA/df_male_epha4.rds')


# female Epha4 data
f_con_data <- c(
  49.7005988,	36.7424242,	19.7368421,
  45.8083832,	50.3787879,	50.6578947,
  4.19161677,	9.84848485,	22.6973684,
  0.2994012,	1.13636364,	4.93421053,
  0,	1.13636364,	1.31578947,
  0,	0.37878788,	0.65789474,
  0,	0.37878788,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0)

f_cold_data <- c(
  18.5676393,	18.6440678,	6.94444444,
  58.6206897,	51.5254237,	35,
  17.5066313,	22.0338983,	27.7777778,
  4.77453581,	4.40677966,	15.5555556,
  0.26525199,	3.05084746,	7.5,
  0.26525199,	0.33898305,	4.16666667,
  0,	0,	2.5,
  0,	0,	0.27777778,
  0,	0,	0.27777778,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0,
  0,	0,	0)

# combine data frames
all_f_scores <- c(f_con_data, f_cold_data)

bin_names <- c("0", "1-5", "6-10", "11-15", "16-20", "21-25", "26-30", 
               "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", 
               "61-65", "66-70", "71-75", "76-80", "81-85", "86-90", 
               "91-95", "96-100")

bin_labels <- rep(bin_names, each = 3, times = 2)

treatment_labels <- rep(c("Control", "Cold"), each = 63)

df_female_epha4 <- data.frame(
  Response = all_f_scores,
  Treatment = factor(treatment_labels),
  Bin = factor(bin_labels, levels = bin_names))

saveRDS(df_female_epha4, '/Users/kodeebao/Documents/ART ANOVA/df_female_epha4.rds')


# loading in data frames --------------------------------------------------

df_male_primary_sphere <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_male_primary_sphere.rds")
df_female_primary_sphere <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_female_primary_sphere.rds")

df_male_secondary_sphere <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_male_secondary_sphere.rds")
df_female_secondary_sphere <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_female_secondary_sphere.rds")


df_male_dlx5 <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_male_dlx5.rds")
df_female_dlx5 <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_female_dlx5.rds")

df_male_epha4 <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_male_epha4.rds")
df_female_epha4 <- readRDS("/Users/kodeebao/Documents/ART ANOVA/df_female_epha4.rds")


# ART on male primary sphere size ------------------------------------------------------
# fit ART model
male_primary <- art(Response ~ Treatment * Size, data = df_male_primary_sphere)

# verify ART appropriateness
summary(male_primary)

# run the ART ANOVA
anova_results_male_primary <- anova(male_primary)

write.csv(as.data.frame(anova_results_male_primary),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_primary_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
male_primary_lm_treatment <- artlm(male_primary, "Treatment")
male_primary_lm_size <- artlm(male_primary, "Size")
male_primary_lm_interaction <- artlm(male_primary, "Treatment:Size")

# treatment differences within each Size (Tukey adjustment)
treatment_posthoc_male_primary <- emmeans(male_primary_lm_treatment,
                             pairwise ~ Treatment | Size,
                             adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_male_primary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_primary_posthoc_treatment.csv",
          row.names = FALSE)

# size differences within each Treatment (Tukey adjustment)
size_posthoc_male_primary <- emmeans(male_primary_lm_size,
                        pairwise ~ Size | Treatment,
                        adjust = "tukey")

write.csv(as.data.frame(size_posthoc_male_primary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_primary_posthoc_size.csv",
          row.names = FALSE)

# interaction post hoc (Tukey adjustment)
interaction_posthoc_male_primary <- emmeans(male_primary_lm_interaction,
                                          pairwise ~ Treatment*Size,
                                          adjust = "tukey")

write.csv(as.data.frame(interaction_posthoc_male_primary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_primary_posthoc_interaction.csv",
          row.names = FALSE)


# ART on female primary sphere size ------------------------------------------------------
# fit ART model
female_primary <- art(Response ~ Treatment * Size, data = df_female_primary_sphere)

# verify ART appropriateness
summary(female_primary)

# run the ART ANOVA
anova_results_female_primary <- anova(female_primary)

write.csv(as.data.frame(anova_results_female_primary),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_primary_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
female_primary_lm_treatment <- artlm(female_primary, "Treatment")
female_primary_lm_size <- artlm(female_primary, "Size")
female_primary_lm_interaction <- artlm(female_primary, "Treatment:Size")

# treatment differences within each Size (Tukey adjustment)
treatment_posthoc_female_primary <- emmeans(female_primary_lm_treatment,
                                          pairwise ~ Treatment | Size,
                                          adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_female_primary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_primary_posthoc_treatment.csv",
          row.names = FALSE)

# size differences within each Treatment (Tukey adjustment)
size_posthoc_female_primary <- emmeans(female_primary_lm_size,
                                     pairwise ~ Size | Treatment,
                                     adjust = "tukey")

write.csv(as.data.frame(size_posthoc_female_primary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_primary_posthoc_size.csv",
          row.names = FALSE)

# interaction post hoc (Tukey adjustment)
interaction_posthoc_female_primary <- emmeans(female_primary_lm_interaction,
                                            pairwise ~ Treatment*Size,
                                            adjust = "tukey")

write.csv(as.data.frame(interaction_posthoc_female_primary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_primary_posthoc_interaction.csv",
          row.names = FALSE)


# ART on male secondary sphere size ------------------------------------------------------
# fit ART model
male_secondary <- art(Response ~ Treatment * Size, data = df_male_secondary_sphere)

# verify ART appropriateness
summary(male_secondary)

# run the ART ANOVA
anova_results_male_secondary <- anova(male_secondary)

write.csv(as.data.frame(anova_results_male_secondary),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_secondary_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
male_secondary_lm_treatment <- artlm(male_secondary, "Treatment")
male_secondary_lm_size <- artlm(male_secondary, "Size")
male_secondary_lm_interaction <- artlm(male_secondary, "Treatment:Size")

# treatment differences within each Size (Tukey adjustment)
treatment_posthoc_male_secondary <- emmeans(male_secondary_lm_treatment,
                                          pairwise ~ Treatment | Size,
                                          adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_male_secondary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_secondary_posthoc_treatment.csv",
          row.names = FALSE)

# size differences within each Treatment (Tukey adjustment)
size_posthoc_male_secondary <- emmeans(male_secondary_lm_size,
                                     pairwise ~ Size | Treatment,
                                     adjust = "tukey")

write.csv(as.data.frame(size_posthoc_male_secondary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_secondary_posthoc_size.csv",
          row.names = FALSE)

# interaction post hoc (Tukey adjustment)
interaction_posthoc_male_secondary <- emmeans(male_secondary_lm_interaction,
                                            pairwise ~ Treatment*Size,
                                            adjust = "tukey")

write.csv(as.data.frame(interaction_posthoc_male_secondary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_secondary_posthoc_interaction.csv",
          row.names = FALSE)


# ART on female secondary sphere size ------------------------------------------------------
# fit ART model
female_secondary <- art(Response ~ Treatment * Size, data = df_female_secondary_sphere)

# verify ART appropriateness
summary(female_secondary)

# run the ART ANOVA
anova_results_female_secondary <- anova(female_secondary)

write.csv(as.data.frame(anova_results_female_secondary),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_secondary_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
female_secondary_lm_treatment <- artlm(female_secondary, "Treatment")
female_secondary_lm_size <- artlm(female_secondary, "Size")
female_secondary_lm_interaction <- artlm(female_secondary, "Treatment:Size")

# treatment differences within each Size (Tukey adjustment)
treatment_posthoc_female_secondary <- emmeans(female_secondary_lm_treatment,
                                            pairwise ~ Treatment | Size,
                                            adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_female_secondary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_secondary_posthoc_treatment.csv",
          row.names = FALSE)

# size differences within each Treatment (Tukey adjustment)
size_posthoc_female_secondary <- emmeans(female_secondary_lm_size,
                                       pairwise ~ Size | Treatment,
                                       adjust = "tukey")

write.csv(as.data.frame(size_posthoc_female_secondary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_secondary_posthoc_size.csv",
          row.names = FALSE)

# interaction post hoc (Tukey adjustment)
interaction_posthoc_female_secondary <- emmeans(female_secondary_lm_interaction,
                                              pairwise ~ Treatment*Size,
                                              adjust = "tukey")

write.csv(as.data.frame(interaction_posthoc_female_secondary$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_secondary_posthoc_interaction.csv",
          row.names = FALSE)



# ART on male dlx5 ------------------------------------------------------
# fit ART model
male_dlx5 <- art(Response ~ Treatment * Bin, data = df_male_dlx5)

# verify ART appropriateness
summary(male_dlx5)

# run the ART ANOVA
anova_results_male_dlx5 <- anova(male_dlx5)

write.csv(as.data.frame(anova_results_male_dlx5),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_dlx5_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
male_dlx5_lm_treatment <- artlm(male_dlx5, "Treatment")
male_dlx5_lm_bin <- artlm(male_dlx5, "Bin")

# treatment differences within each Bin (Tukey adjustment)
treatment_posthoc_male_dlx5 <- emmeans(male_dlx5_lm_treatment,
                                          pairwise ~ Treatment | Bin,
                                          adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_male_dlx5$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_dlx5_posthoc_treatment(tukey).csv",
          row.names = FALSE)

# puncta per cell bin differences within each Treatment (Tukey adjustment)
bin_posthoc_male_dlx5 <- emmeans(male_dlx5_lm_bin,
                                     pairwise ~ Bin | Treatment,
                                     adjust = "tukey")

write.csv(as.data.frame(bin_posthoc_male_dlx5$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_dlx5_posthoc_bin(tukey).csv",
          row.names = FALSE)


# ART on female dlx5 ------------------------------------------------------
# fit ART model
female_dlx5 <- art(Response ~ Treatment * Bin, data = df_female_dlx5)

# verify ART appropriateness
summary(female_dlx5)

# run the ART ANOVA
anova_results_female_dlx5 <- anova(female_dlx5)

write.csv(as.data.frame(anova_results_female_dlx5),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_dlx5_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
female_dlx5_lm_treatment <- artlm(female_dlx5, "Treatment")
female_dlx5_lm_bin <- artlm(female_dlx5, "Bin")

# treatment differences within each Bin (Tukey adjustment)
treatment_posthoc_female_dlx5 <- emmeans(female_dlx5_lm_treatment,
                                       pairwise ~ Treatment | Bin,
                                       adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_female_dlx5$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_dlx5_posthoc_treatment(tukey).csv",
          row.names = FALSE)

# puncta per cell bin differences within each Treatment (Tukey adjustment)
bin_posthoc_female_dlx5 <- emmeans(female_dlx5_lm_bin,
                                 pairwise ~ Bin | Treatment,
                                 adjust = "tukey")

write.csv(as.data.frame(bin_posthoc_female_dlx5$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_dlx5_posthoc_bin(tukey).csv",
          row.names = FALSE)


# ART on male epha4 ------------------------------------------------------
# fit ART model
male_epha4 <- art(Response ~ Treatment * Bin, data = df_male_epha4)

# verify ART appropriateness
summary(male_epha4)

# run the ART ANOVA
anova_results_male_epha4 <- anova(male_epha4)

write.csv(as.data.frame(anova_results_male_epha4),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_epha4_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
male_epha4_lm_treatment <- artlm(male_epha4, "Treatment")
male_epha4_lm_bin <- artlm(male_epha4, "Bin")

# treatment differences within each Bin (Tukey adjustment)
treatment_posthoc_male_epha4 <- emmeans(male_epha4_lm_treatment,
                                       pairwise ~ Treatment | Bin,
                                       adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_male_epha4$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_epha4_posthoc_treatment(tukey).csv",
          row.names = FALSE)

# puncta per cell bin differences within each Treatment (Tukey adjustment)
bin_posthoc_male_epha4 <- emmeans(male_epha4_lm_bin,
                                 pairwise ~ Bin | Treatment,
                                 adjust = "tukey")

write.csv(as.data.frame(bin_posthoc_male_epha4$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/male_epha4_posthoc_bin(tukey).csv",
          row.names = FALSE)


# ART on female epha4 ------------------------------------------------------
# fit ART model
female_epha4 <- art(Response ~ Treatment * Bin, data = df_female_epha4)

# verify ART appropriateness
summary(female_epha4)

# run the ART ANOVA
anova_results_female_epha4 <- anova(female_epha4)

write.csv(as.data.frame(anova_results_female_epha4),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_epha4_anova.csv",
          row.names = TRUE)

# post hoc 
# convert ART object to lm for emmeans
female_epha4_lm_treatment <- artlm(female_epha4, "Treatment")
female_epha4_lm_bin <- artlm(female_epha4, "Bin")

# treatment differences within each Bin (Tukey adjustment)
treatment_posthoc_female_epha4 <- emmeans(female_epha4_lm_treatment,
                                         pairwise ~ Treatment | Bin,
                                         adjust = "tukey")

write.csv(as.data.frame(treatment_posthoc_female_epha4$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_epha4_posthoc_treatment(tukey).csv",
          row.names = FALSE)

# puncta per cell bin differences within each Treatment (Tukey adjustment)
bin_posthoc_female_epha4 <- emmeans(female_epha4_lm_bin,
                                   pairwise ~ Bin | Treatment,
                                   adjust = "tukey")

write.csv(as.data.frame(bin_posthoc_female_epha4$contrasts),
          "/Users/kodeebao/Documents/ART ANOVA/results 2/female_epha4_posthoc_bin(tukey).csv",
          row.names = FALSE)


