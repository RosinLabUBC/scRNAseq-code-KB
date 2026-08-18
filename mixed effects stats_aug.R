
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)


# primary spheres ----------------------------------------------------------

primary_sphere_numbers <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    34, 42,
    # cold M (2 embryos)
    58, 48.5,
    # control F (1 embryo)
    28,
    # cold F (2 embryos)
    38, 41,
    
    # --- LITTER 2 ---
    # con M (1)
    86, 
    # cold M (6)
    66.5, 44.5, 77, 70.5, 48.5, 68.5, 
    # con F (7)
    49.5, 49.5, 37.5, 47.5, 34.5, 49, 22.5,
    # cold F (4)
    63.5, 61.5, 53.5, 68.5,
    
    # --- LITTER 3 ---
    # con M (5)
    47.5, 75.5, 77.5, 74, 48,
    # cold M (7)
    44.5, 55, 47.5, 48, 60.5, 52.5, 70.5, 
    # con F (5)
    73, 67, 48, 63, 48.5, 
    # cold F (3)
    51.5, 66.5, 47),
  
  sex = c(
    rep("M", 4), rep("F", 3),
    rep("M", 7), rep("F", 11),
    rep("M", 12), rep("F", 8)),
  
  condition = c(
    rep("control", 2), rep("cold", 2), rep("control", 1), rep("cold", 2), 
    rep("control", 1), rep("cold", 6), rep("control", 7), rep("cold", 4), 
    rep("control", 5), rep("cold", 7), rep("control", 5), rep("cold", 3)),
  
  litter = c(
    rep("L1", 7),
    rep("L2", 18),
    rep("L3", 20)))

# convert to factors
primary_sphere_numbers$sex <- factor(primary_sphere_numbers$sex)
primary_sphere_numbers$condition <- factor(primary_sphere_numbers$condition)
primary_sphere_numbers$litter <- factor(primary_sphere_numbers$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = primary_sphere_numbers)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_primary_sphere_numb.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_primary_sphere_numb.csv")



# prim spheres + size -----------------------------------------------------

sphere_values = c(
    # --- LITTER 1 ---
    # Embryo E1 (Control M): counts across size bins
    17.5, 14, 2.5,   # 45-100, 100-200, >200
    # Embryo E2 (Control M)
    18.5, 17, 6.5,
    
    # Embryo E3 (Cold M)
    37.5, 16.5, 4,   # 45-100, 100-200, >200
    # Embryo E4 (Cold M)
    29, 12, 7.5,
    
    # Embryo E5 (Con F)
    17, 7.5, 3.5,   # 45-100, 100-200, >200
    
    # Embryo E6 (Cold F)
    15.5, 12.5, 10,   # 45-100, 100-200, >200
    # Embryo E7 (Cold F)
    23.5, 12.5, 5,   
    
    # --- LITTER 2 ---
    # Embryo E1 (Control M): counts across size bins
    32, 23, 12,   # 45-100, 100-200, >200
    
    # Embryo E2 (Cold M)
    25.5, 11, 30,   # 45-100, 100-200, >200
    # Embryo E3 (Cold M)
    13, 19, 12.5,
    # Embryo E4 (Cold M)
    32.5, 24.5, 20,
    # Embryo E5 (Cold M)
    35, 24, 11.5,
    # Embryo E6 (Cold M)
    24.5, 15, 9,
    # Embryo E7 (Cold M)
    46.5, 20.5, 1.5,
    
    # Embryo E8 (Con F)
    20, 22.5, 7,   # 45-100, 100-200, >200
    # Embryo E9 (Con F)
    24.5, 16.5, 8.5,   
    # Embryo E10 (Con F)
    20.5, 14, 3,   
    # Embryo E11 (Con F)
    26, 15.5, 6,
    # Embryo E12 (Con F)
    25.5, 8, 1,   
    # Embryo E13 (Con F)
    23.5, 18.5, 7,   
    # Embryo E14 (Con F)
    13.5, 9, 0, 
    
    # Embryo E15 (Cold F)
    33.5, 16.5, 13.5,   # 45-100, 100-200, >200
    # Embryo E16 (Cold F)
    29.5, 16.5, 15.5,
    # Embryo E17 (Cold F)
    25.5, 21, 7,  
    # Embryo E18 (Cold F)
    41, 22.5, 5,  
    
    # --- LITTER 3 ---
    # Embryo E1 (Control M): counts across size bins
    24.5, 18.5, 4.5,   # 45-100, 100-200, >200
    # Embryo E2 (Control M)
    40.5, 24, 11,
    # Embryo E3 (Control M)
    38.5, 29, 10,
    # Embryo E4 (Control M)
    43.5, 17, 13.5,
    # Embryo E5 (Control M)
    25.5, 16, 6.5,
    
    # Embryo E6 (Cold M)
    25, 12.5, 7,   # 45-100, 100-200, >200
    # Embryo E7 (Cold M)
    25.5, 19, 10.5,
    # Embryo E8 (Cold M)
    26, 16, 5.5,
    # Embryo E9 (Cold M)
    25, 13, 10,   
    # Embryo E10 (Cold M)
    34, 16.5, 10,   
    # Embryo E11 (Cold M)
    37, 12.5, 3,
    # Embryo E12 (Cold M)
    43.5, 16, 11,   
    
    # Embryo E13 (Con F)
    35.5, 33, 4.5,   # 45-100, 100-200, >200   
    # Embryo E14 (Con F)
    39, 21.5, 6.5, 
    # Embryo E15 (Con F)
    27.5, 12, 8.5,
    # Embryo E16 (Con F)
    35, 16.5, 11.5, 
    # Embryo E17 (Con F)
    25, 18, 5.5,
    
    # Embryo E18 (Cold F)
    32, 15.5, 4,   # 45-100, 100-200, >200
    # Embryo E19 (Cold F)
    36.5, 20.5, 9.5,  
    # Embryo E20 (Cold F)
    26, 14.5, 6.5)
  
sex = c(
    rep("M", 4), rep("F", 3),
    rep("M", 7), rep("F", 11),
    rep("M", 12), rep("F", 8))
  
condition = c(
    rep("control", 2), rep("cold", 2), rep("control", 1), rep("cold", 2), 
    rep("control", 1), rep("cold", 6), rep("control", 7), rep("cold", 4), 
    rep("control", 5), rep("cold", 7), rep("control", 5), rep("cold", 3))
  
litter = c(
    rep("L1", 7),
    rep("L2", 18),
    rep("L3", 20))
  
# Automatically generate sequential embryo IDs (L1_E1, L1_E2...)
embryo_id <- paste0(litter, "_E", c(1:7, 1:18, 1:20))
  
# 3. Build the final long-format data frame
primary_sphere_numbers_sized <- data.frame(
  value = sphere_values,
  # Repeat each embryo's metadata 3 times (once for each bin)
  sex_sized = factor(rep(sex, each = 3)),
  condition_sized = factor(rep(condition, each = 3)),
  litter_sized = factor(rep(litter, each = 3)),
  embryo_id_sized = factor(rep(embryo_id, each = 3)),
  
  # Repeat the 3 size bins 45 times (once per embryo)
  bin_size = factor(
    rep(c("45-100", "100-200", ">200"), times = 45), 
    levels = c("45-100", "100-200", ">200")))


# 4. Run mixed-effects model
model <- lmer(value ~ sex_sized * condition_sized * bin_size + (1 | litter_sized/embryo_id_sized), data = primary_sphere_numbers_sized)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex_sized * condition_sized * bin_size)
posthoc_results <- bind_rows(
  "Condition (Control vs Cold)" = as.data.frame(pairs(emm, by = c("bin_size", "sex_sized"), adjust = "bonferroni")),
  "Sex (Male vs Female)" = as.data.frame(pairs(emm, by = c("bin_size", "condition_sized"), adjust = "bonferroni")),
  "Bin (Between Bins)" = as.data.frame(pairs(emm, by = c("condition_sized", "sex_sized"), adjust = "bonferroni")),
  
  "Bin and Condition Interaction" = as.data.frame(pairs(emm, by = "sex_sized", adjust = "bonferroni")), .id = "Comparison_Type")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_primary_sphere_numbers_sized.csv")

write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_primary_sphere_numbers_sized.csv")


# secondary spheres ----------------------------------------------------------

secondary_sphere_numbers <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    7.5, 30,
    # cold M (2 embryos)
    21, 33,
    # control F (4 embryo)
    13, 9.5, 11, 7,
    # cold F (3 embryos)
    21, 38.5, 24.5,
    
    # --- LITTER 2 ---
    # con M (1)
    47.5, 
    # cold M (4)
    43, 38, 35.5, 48.5, 
    # con F (4)
    27.5, 56, 27.5, 28, 
    # cold F (4)
    13, 33, 23, 47,
    
    # --- LITTER 3 ---
    # con M (4)
    24.5, 52, 24, 27, 
    # cold M (4)
    47.5, 60, 28, 36,  
    # con F (4)
    29.5, 27.5, 40, 39.5,  
    # cold F (3)
    12.5, 30, 29),
  
  sex = c(
    rep("M", 4), rep("F", 7),
    rep("M", 5), rep("F", 8),
    rep("M", 8), rep("F", 7)),
  
  condition = c(
    rep("control", 2), rep("cold", 2), rep("control", 4), rep("cold", 3), 
    rep("control", 1), rep("cold", 4), rep("control", 4), rep("cold", 4), 
    rep("control", 4), rep("cold", 4), rep("control", 4), rep("cold", 3)),
  
  litter = c(
    rep("L1", 11),
    rep("L2", 13),
    rep("L3", 15)))

# convert to factors
secondary_sphere_numbers$sex <- factor(secondary_sphere_numbers$sex)
secondary_sphere_numbers$condition <- factor(secondary_sphere_numbers$condition)
secondary_sphere_numbers$litter <- factor(secondary_sphere_numbers$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = secondary_sphere_numbers)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_secondary_sphere_numb.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_secondary_sphere_numb.csv")



# secondary spheres + size ------------------------------------------------

sphere_values = c(
  # --- LITTER 1 ---
  # Embryo E1 (Control M): counts across size bins
  2.5, 1, 4,   # 45-100, 100-200, >200
  # Embryo E2 (Control M)
  4, 3, 23,
  
  # Embryo E3 (Cold M)
  3, 6, 12,   # 45-100, 100-200, >200
  # Embryo E4 (Cold M)
  3, 5.5, 24.5,
  
  # Embryo E5 (Con F)
  0.5, 1, 11.5,   # 45-100, 100-200, >200
  # Embryo E6 (Con F)
  3, 2.5, 4,
  # Embryo E7 (Con F)
  2.5, 2, 6.5,
  # Embryo E8 (Con F)
  3, 3, 1,
  
  # Embryo E9 (Cold F)
  4.5, 2.5, 14,   # 45-100, 100-200, >200
  # Embryo E10 (Cold F)
  12.5, 6, 20,
  # Embryo E11 (Cold F)
  7.5, 3.5, 13.5,  
  
  # --- LITTER 2 ---
  # Embryo E1 (Control M): counts across size bins
  12, 9, 26.5,   # 45-100, 100-200, >200
  
  # Embryo E2 (Cold M)
  7.5, 15, 20.5,   # 45-100, 100-200, >200
  # Embryo E3 (Cold M)
  2, 13.5, 22.5,
  # Embryo E4 (Cold M)
  4, 4, 27.5,
  # Embryo E5 (Cold M)
  9, 10.5, 29,
  
  # Embryo E6 (Con F)
  4, 11, 12.5,   # 45-100, 100-200, >200
  # Embryo E7 (Con F)
  10.5, 16, 29.5,   
  # Embryo E8 (Con F)
  5, 8.5, 14,   
  # Embryo E9 (Con F)
  3.5, 7, 17.5,
  
  # Embryo E10 (Cold F)
  3.5, 2, 7.5,   # 45-100, 100-200, >200
  # Embryo E11 (Cold F)
  3.5, 11, 18.5,
  # Embryo E12 (Cold F)
  1.5, 2.5, 19,  
  # Embryo E13 (Cold F)
  6, 12, 29,  
  
  # --- LITTER 3 ---
  # Embryo E1 (Control M): counts across size bins
  3.5, 6, 15,   # 45-100, 100-200, >200
  # Embryo E2 (Control M)
  16.5, 17.5, 18,
  # Embryo E3 (Control M)
  3.5, 4.5, 16,
  # Embryo E4 (Control M)
  9, 6.5, 11.5,

  # Embryo E5 (Cold M)
  14.5, 13, 20,   # 45-100, 100-200, >200
  # Embryo E6 (Cold M)
  19.5, 24.5, 16,
  # Embryo E7 (Cold M)
  4, 8, 16,
  # Embryo E8 (Cold M)
  4, 9, 23,

  # Embryo E9 (Con F)
  8, 7.5, 14,   # 45-100, 100-200, >200   
  # Embryo E10 (Con F)
  4.5, 6, 17, 
  # Embryo E11 (Con F)
  13.5, 15, 11.5,
  # Embryo E12 (Con F)
  14, 8, 17.5, 

  # Embryo E13 (Cold F)
  8, 3, 1.5,   # 45-100, 100-200, >200
  # Embryo E14 (Cold F)
  9, 10, 11,  
  # Embryo E15 (Cold F)
  4, 7.5, 17.5)

sex = c(
  rep("M", 4), rep("F", 7),
  rep("M", 5), rep("F", 8),
  rep("M", 8), rep("F", 7))

condition = c(
  rep("control", 2), rep("cold", 2), rep("control", 4), rep("cold", 3), 
  rep("control", 1), rep("cold", 4), rep("control", 4), rep("cold", 4), 
  rep("control", 4), rep("cold", 4), rep("control", 4), rep("cold", 3))

litter = c(
  rep("L1", 11),
  rep("L2", 13),
  rep("L3", 15))

# Automatically generate sequential embryo IDs (L1_E1, L1_E2...)
embryo_id <- paste0(litter, "_E", c(1:11, 1:13, 1:15))

# 3. Build the final long-format data frame
secondary_sphere_numbers_sized <- data.frame(
  value = sphere_values,
  # Repeat each embryo's metadata 3 times (once for each bin)
  sex_sized = factor(rep(sex, each = 3)),
  condition_sized = factor(rep(condition, each = 3)),
  litter_sized = factor(rep(litter, each = 3)),
  embryo_id_sized = factor(rep(embryo_id, each = 3)),
  
  # Repeat the 3 size bins 45 times (once per embryo)
  bin_size = factor(
    rep(c("45-100", "100-200", ">200"), times = 39), 
    levels = c("45-100", "100-200", ">200")))


# 4. Run mixed-effects model
model <- lmer(value ~ sex_sized * condition_sized * bin_size + (1 | litter_sized/embryo_id_sized), data = secondary_sphere_numbers_sized)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex_sized * condition_sized * bin_size)
posthoc_results <- bind_rows(
  "Condition (Control vs Cold)" = as.data.frame(pairs(emm, by = c("bin_size", "sex_sized"), adjust = "bonferroni")),
  "Sex (Male vs Female)" = as.data.frame(pairs(emm, by = c("bin_size", "condition_sized"), adjust = "bonferroni")),
  "Bin (Between Bins)" = as.data.frame(pairs(emm, by = c("condition_sized", "sex_sized"), adjust = "bonferroni")),
  
  "Bin and Condition Interaction" = as.data.frame(pairs(emm, by = "sex_sized", adjust = "bonferroni")), .id = "Comparison_Type")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_secondary_sphere_numbers_sized.csv")

write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_secondary_sphere_numbers_sized.csv")

# TUJ1 neuron numbers -----------------------------------------------------

tuj1_neuron_numbers <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    2.86, 2.78,
    # cold M
    8.24, 9.77,
    # control F
    0.613, 0.94,
    # cold F
    5.52, 1.91,
    
    # --- LITTER 2 ---
    2.07, 4.49,
    9.65, 5.14,
    1.98, 8.08,
    3.51, 5.8,
    
    # --- LITTER 3 ---
    3.49, 2.59,
    9.87, 21.43,
    7.8, 2.63,
    19.42, 7.59),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
tuj1_neuron_numbers$sex <- factor(tuj1_neuron_numbers$sex)
tuj1_neuron_numbers$condition <- factor(tuj1_neuron_numbers$condition)
tuj1_neuron_numbers$litter <- factor(tuj1_neuron_numbers$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = tuj1_neuron_numbers)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_tuj1_neuron_num.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_tuj1_neuron_num.csv")

# OLIG2 OPC numbers -----------------------------------------------------

olig2_opc_numbers <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    14.71, 19.45,
    # cold M
    20.46, 31.37,
    # control F
    27.8, 21.79,
    # cold F
    22.13, 22.56,
    
    # --- LITTER 2 ---
    29.63, 31.25,
    41.58, 15.45,
    29.94, 39.86,
    27.09, 21.61,
    
    # --- LITTER 3 ---
    32.47, 26.19,
    30.25, 17.39,
    13.17, 21.76,
    24.05, 24),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
olig2_opc_numbers$sex <- factor(olig2_opc_numbers$sex)
olig2_opc_numbers$condition <- factor(olig2_opc_numbers$condition)
olig2_opc_numbers$litter <- factor(olig2_opc_numbers$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = olig2_opc_numbers)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_olig2_opc_num.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_olig2_opc_num.csv")

# OLIG2/PDGFRa oligo numbers -----------------------------------------------------

olig2_pdgfra_oligo_numbers <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    6.47, 4.99,
    # cold M
    11.55, 17.84,
    # control F
    12.64, 7.16,
    # cold F
    7.76, 8.42,
    
    # --- LITTER 2 ---
    9.88, 10.63,
    12.63, 4.55,
    10.19, 18.91,
    16.75, 13.92,
    
    # --- LITTER 3 ---
    17.01, 10.12,
    10.49, 7.97,
    7.8, 11.18,
    5.06, 8.57),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
olig2_pdgfra_oligo_numbers$sex <- factor(olig2_pdgfra_oligo_numbers$sex)
olig2_pdgfra_oligo_numbers$condition <- factor(olig2_pdgfra_oligo_numbers$condition)
olig2_pdgfra_oligo_numbers$litter <- factor(olig2_pdgfra_oligo_numbers$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = olig2_pdgfra_oligo_numbers)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_olig2_pdgfra_oligo_num.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_olig2_pdgfra_oligo_num.csv")

# GFAP astrocytes numbers -----------------------------------------------------

gfap_astrocyte_numbers <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    74.83, 63.38,
    # cold M
    53.22, 44.09,
    # control F
    48.28, 55.05,
    # cold F
    67.99, 58.76,
    
    # --- LITTER 2 ---
    58.52, 37.11,
    55.5, 72.44,
    59.54, 31.74,
    60.29, 65.51,
    
    # --- LITTER 3 ---
    58.9, 48.99,
    54.55, 61.59,
    62.2, 37.5,
    41.67, 43.56),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
gfap_astrocyte_numbers$sex <- factor(gfap_astrocyte_numbers$sex)
gfap_astrocyte_numbers$condition <- factor(gfap_astrocyte_numbers$condition)
gfap_astrocyte_numbers$litter <- factor(gfap_astrocyte_numbers$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = gfap_astrocyte_numbers)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_gfap_astrocyte_num.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_gfap_astrocyte_num.csv")


# total number dendrites ---------------------------------------------
# input data
total_num_dendrites <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    7.6, 3.8,
    # cold M
    7.6, 2.6,
    # control F
    2.2, 3.8,
    # cold F
    3, 4.8,
    
    # --- LITTER 2 ---
    3.4, 7.6,
    5.8, 7.2,
    1.6, 2.8,
    3.4, 6.8,
    
    # --- LITTER 3 ---
    4, 4,
    4, 5.4,
    4, 2.8,
    4.4, 6),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
total_num_dendrites$sex <- factor(total_num_dendrites$sex)
total_num_dendrites$condition <- factor(total_num_dendrites$condition)
total_num_dendrites$litter <- factor(total_num_dendrites$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = total_num_dendrites)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_total_num_dendrites.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_total_num_dendrites.csv")


# cumulative dendritic length ---------------------------------------------
# input data
cumulative_dendritic_length <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    151.7046, 97.2642,
    # cold M
    202.9116, 54.7154,
    # control F
    60.5548, 64.1272,
    # cold F
    93.5394, 112.9394,
    
    # --- LITTER 2 ---
    123.5388, 163.6448,
    123.0862, 185.1732,
    60.3858, 94.752,
    103.554, 153.783,
    
    # --- LITTER 3 ---
    164.4986, 132.5422,
    113.9104, 95.0378,
    68.1888, 91.7796,
    97.6992, 125.5952),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
cumulative_dendritic_length$sex <- factor(cumulative_dendritic_length$sex)
cumulative_dendritic_length$condition <- factor(cumulative_dendritic_length$condition)
cumulative_dendritic_length$litter <- factor(cumulative_dendritic_length$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = cumulative_dendritic_length)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_total_dendritic_length.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_total_dendritic_length.csv")


# primary dendrite num ---------------------------------------------
# input data
primary_dendrite_num <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    2.4, 1.8,
    # cold M
    3.2, 2,
    # control F
    2, 1.6,
    # cold F
    2, 2.2,
    
    # --- LITTER 2 ---
    2, 3.2,
    3.6, 3,
    2, 1.6,
    1.8, 3.2,
    
    # --- LITTER 3 ---
    2, 2.4,
    1.6, 3.4,
    2.4, 2,
    2, 2.6),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
primary_dendrite_num$sex <- factor(primary_dendrite_num$sex)
primary_dendrite_num$condition <- factor(primary_dendrite_num$condition)
primary_dendrite_num$litter <- factor(primary_dendrite_num$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = primary_dendrite_num)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_primary_dendrite_num.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_primary_dendrite_num.csv")


# primary dendrite length ---------------------------------------------
# input data
primary_dendrite_length <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    37.919, 35.64,
    # cold M
    45.475, 23.94,
    # control F
    24.194, 15.715,
    # cold F
    42.966, 41.543,
    
    # --- LITTER 2 ---
    34.388, 35.735,
    33.441, 42.951,
    34.73, 45.743,
    32.701, 27.431,
    
    # --- LITTER 3 ---
    56.125, 33.767,
    43.659, 24.048,
    19.309, 32.351,
    35.495, 33.075),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
primary_dendrite_length$sex <- factor(primary_dendrite_length$sex)
primary_dendrite_length$condition <- factor(primary_dendrite_length$condition)
primary_dendrite_length$litter <- factor(primary_dendrite_length$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = primary_dendrite_length)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_primary_dendrite_length.csv")



posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_primary_dendrite_length.csv")


# secondary dendrite num ---------------------------------------------
# input data
secondary_dendrite_num <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    4.8, 1.8,
    # cold M
    4, 0.6,
    # control F
    0.2, 2,
    # cold F
    1, 2.6,
    
    # --- LITTER 2 ---
    1.4, 4.2,
    2.2, 4,
    0, 1.2,
    1.6, 3,
    
    # --- LITTER 3 ---
    1.8, 1.6,
    2.4, 1.8,
    1.4, 0.8,
    2.4, 3.2),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
secondary_dendrite_num$sex <- factor(secondary_dendrite_num$sex)
secondary_dendrite_num$condition <- factor(secondary_dendrite_num$condition)
secondary_dendrite_num$litter <- factor(secondary_dendrite_num$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = secondary_dendrite_num)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_secondary_dendrite_num.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_secondary_dendrite_num.csv")


# secondary dendrite length ---------------------------------------------
# input data
secondary_dendrite_length <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    15.176, 11.356,
    # cold M
    16.945, 4.361,
    # control F
    1.101, 10.487,
    # cold F
    15.757, 11.979,
    
    # --- LITTER 2 ---
    11.329, 15.178,
    7.017, 12.213,
    0, 6.038,
    13.658, 10.748,
    
    # --- LITTER 3 ---
    19.175, 12.806,
    13.136, 7.645,
    7.239, 13.338,
    7.385, 12.202),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
secondary_dendrite_length$sex <- factor(secondary_dendrite_length$sex)
secondary_dendrite_length$condition <- factor(secondary_dendrite_length$condition)
secondary_dendrite_length$litter <- factor(secondary_dendrite_length$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = secondary_dendrite_length)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_secondary_dendrite_length.csv")



posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_secondary_dendrite_length.csv")


# tertiary dendrite num ---------------------------------------------
# input data
tertiary_dendrite_num <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    0.4, 0.2,
    # cold M
    0.4, 0,
    # control F
    0, 0.2,
    # cold F
    0, 0,
    
    # --- LITTER 2 ---
    0, 0.2,
    0, 0.2,
    0, 0,
    0, 0.6,
    
    # --- LITTER 3 ---
    0.2, 0,
    0, 0.2,
    0.2, 0,
    0, 0.2),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
tertiary_dendrite_num$sex <- factor(tertiary_dendrite_num$sex)
tertiary_dendrite_num$condition <- factor(tertiary_dendrite_num$condition)
tertiary_dendrite_num$litter <- factor(tertiary_dendrite_num$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = tertiary_dendrite_num)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_tertiary_dendrite_num.csv")

posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_tertiary_dendrite_num.csv")


# tertiary dendrite length ---------------------------------------------
# input data
tertiary_dendrite_length <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    3.14, 2.233,
    # cold M
    4.054, 0,
    # control F
    0, 1.245,
    # cold F
    0, 0,
    
    # --- LITTER 2 ---
    0, 0.87,
    0, 2.216,
    0, 0,
    0, 5.666,
    
    # --- LITTER 3 ---
    1.447, 0,
    0, 0.437,
    0.81, 0,
    0, 4.369),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
tertiary_dendrite_length$sex <- factor(tertiary_dendrite_length$sex)
tertiary_dendrite_length$condition <- factor(tertiary_dendrite_length$condition)
tertiary_dendrite_length$litter <- factor(tertiary_dendrite_length$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = tertiary_dendrite_length)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_tertiary_dendrite_length.csv")



posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_tertiary_dendrite_length.csv")


# average axon length ---------------------------------------------

ave_axon_length <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    140.378, 190.895,
    # cold M
    413.589, 262.697,
    # control F
    141.628, 196.893,
    # cold F
    520.54, 217.863,
    
    # --- LITTER 2 ---
    146.35, 144.133,
    400.404, 327.615,
    151.705, 208.395,
    306.634, 498.815,
    
    # --- LITTER 3 ---
    284.124, 189.922,
    516.038, 401.581,
    325.054, 158.698,
    312.418, 285.311),
  
  sex = c(
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4),
    rep("M", 4), rep("F", 4)),
  
  condition = c(
    rep(c("control","control","cold","cold"), 3),
    rep(c("control","control","cold","cold"), 3)),
  
  litter = c(
    rep("L1", 8),
    rep("L2", 8),
    rep("L3", 8)))

# convert to factors
ave_axon_length$sex <- factor(ave_axon_length$sex)
ave_axon_length$condition <- factor(ave_axon_length$condition)
ave_axon_length$litter <- factor(ave_axon_length$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = ave_axon_length)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_ave_axon_length.csv")



posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_ave_axon_length.csv")



# dlx2 frequency ----------------------------------------------------------

puncta_values = c(
  # --- LITTER 1 ---
  # Litter 1, embryo 1 (Control M): counts across size bins (0 - 100)
  12.7725857,
  29.9065421,
  16.1993769,
  14.0186916,
  9.65732087,
  6.54205607,
  4.04984424,
  3.73831776,
  0.93457944,
  1.24610592,
  0.93457944,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  # Litter 1, embryo 2 (Control M)
  2.60586319,
  7.16612378,
  6.51465798,
  9.77198697,
  9.77198697,
  14.3322476,
  13.029316,
  8.46905537,
  10.0977199,
  4.23452769,
  5.86319218,
  2.28013029,
  1.95439739,
  1.95439739,
  0.9771987,
  0.6514658,
  0,
  0.3257329,
  0,
  0,
  0,
  
  # Litter 1, embryo 1 (Cold M)
  8.81355932,
  29.1525424,
  21.3559322,
  15.2542373,
  11.8644068,
  8.13559322,
  2.03389831,
  1.69491525,
  0,
  0.6779661,
  0.6779661,
  0.33898305,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,  
  
  # Litter 1, embryo 2 (Cold M)
  2.76243094,
  17.1270718,
  16.5745856,
  16.2983425,
  11.6022099,
  9.94475138,
  7.18232044,
  7.73480663,
  5.24861878,
  2.20994475,
  1.65745856,
  0.82872928,
  0.27624309,
  0.27624309,
  0.27624309,
  0,
  0,
  0,
  0,
  0,
  0,
  
  # Litter 1, embryo 1 (Control F)
  82.3809524,
  11.4285714,
  3.80952381,
  1.9047619,
  0.47619048,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  
  # Litter 1, embryo 1 (Cold F)
  8.4375,
  32.8125,
  12.1875,
  11.5625,
  11.25,
  9.375,
  5,
  4.375,
  2.5,
  1.25,
  0.9375,
  0,
  0.3125,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  # Litter 1, embryo 2 (Cold F)
  4.65838509,
  15.2173913,
  10.8695652,
  14.9068323,
  11.4906832,
  8.69565217,
  8.69565217,
  8.69565217,
  4.65838509,
  4.65838509,
  2.17391304,
  1.86335404,
  0.31055901,
  1.55279503,
  0.93167702,
  0,
  0.31055901,
  0.31055901,
  0,
  0,
  0,
  
  # --- LITTER 2 ---
  # Litter 2, embryo 1 (Control M): counts across size bins (0 - 100)
  14.7540984,
  34.4262295,
  14.3442623,
  12.704918,
  9.83606557,
  5.32786885,
  4.09836066,
  0.81967213,
  1.63934426,
  1.2295082,
  0.40983607,
  0,
  0.40983607,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  # Litter 2, embryo 2 (Control M)
  8.2278481,
  27.8481013,
  24.3670886,
  12.9746835,
  10.443038,
  7.27848101,
  1.89873418,
  3.48101266,
  1.26582278,
  1.26582278,
  0.3164557,
  0.3164557,
  0.3164557,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  
  # Litter 2, embryo 1 (Cold M)
  13.4615385,
  33.0769231,
  14.6153846,
  11.9230769,
  11.1538462,
  8.07692308,
  5,
  1.53846154,
  0.38461538,
  0.38461538,
  0.38461538,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  # Litter 2, embryo 2 (Cold M)
  15.1162791,
  41.8604651,
  22.6744186,
  9.30232558,
  4.06976744,
  2.61627907,
  2.61627907,
  0.29069767,
  1.1627907,
  0,
  0,
  0.29069767,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,

  # Litter 2, embryo 1 (Control F)
  40.7843137,
  38.4313725,
  13.3333333,
  3.52941176,
  1.56862745,
  0.78431373,
  0.39215686,
  0.78431373,
  0.39215686,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  # Litter 2, embryo 2 (Control F)
  16.796875,
  44.140625,
  16.015625,
  10.546875,
  5.078125,
  3.515625,
  1.953125,
  1.5625,
  0.390625,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  
  # Litter 2, embryo 1 (Cold F)
  38.1188119,
  37.1287129,
  8.41584158,
  4.95049505,
  2.47524752,
  3.96039604,
  2.47524752,
  0.4950495,
  0.4950495,
  0.99009901,
  0,
  0.4950495,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  
  # --- LITTER 3 ---
  # Litter 3, embryo 1 (Control F)
  40.85106383,
  40.42553191,
  10.63829787,
  3.829787234,
  2.127659574,
  1.70212766,
  0,
  0.425531915,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  
  # Litter 3, embryo 1 (Cold F)
  4.27350427,
  12.8205128,
  15.0997151,
  12.8205128,
  13.3903134,
  11.965812,
  9.68660969,
  5.6980057,
  3.7037037,
  2.84900285,
  1.99430199,
  2.56410256,
  1.13960114,
  0.56980057,
  0.56980057,
  0.28490028,
  0.28490028,
  0.28490028,
  0,
  0,
  0)

sex = c(
  rep("M", 4), rep("F", 3),
  rep("M", 4), rep("F", 3),
  rep("M", 0), rep("F", 2))

condition = c(
  rep("control", 2), rep("cold", 2), rep("control", 1), rep("cold", 2), 
  rep("control", 2), rep("cold", 2), rep("control", 2), rep("cold", 1), 
  rep("control", 0), rep("cold", 0), rep("control", 1), rep("cold", 1))

litter = c(
  rep("L1", 7),
  rep("L2", 7),
  rep("L3", 2))

# Automatically generate sequential embryo IDs (L1_E1, L1_E2...)
embryo_id <- paste0(litter, "_E", c(1:7, 1:7, 1:2))

# 3. Build the final long-format data frame
dlx2_frequency <- data.frame(
  value = puncta_values,
  # Repeat each embryo's metadata for each bin
  sex_sized = factor(rep(sex, each = 21)),
  condition_sized = factor(rep(condition, each = 21)),
  litter_sized = factor(rep(litter, each = 21)),
  embryo_id_sized = factor(rep(embryo_id, each = 21)),
  
  # Repeat the 21 size bins 16 times (once per embryo)
  bin_size = factor(
    rep(c(
      "≤0",
      "1-5",
      "6-10",
      "11-15",
      "16-20",
      "21-25",
      "26-30",
      "31-35",
      "36-40",
      "41-45",
      "46-50",
      "51-55",
      "56-60",
      "61-65",
      "66-70",
      "71-75",
      "76-80",
      "81-85",
      "86-90",
      "91-95",
      "96-100"), times = 16), 
    levels = c(
      "≤0",
      "1-5",
      "6-10",
      "11-15",
      "16-20",
      "21-25",
      "26-30",
      "31-35",
      "36-40",
      "41-45",
      "46-50",
      "51-55",
      "56-60",
      "61-65",
      "66-70",
      "71-75",
      "76-80",
      "81-85",
      "86-90",
      "91-95",
      "96-100")))


# 4. Run mixed-effects model
model <- lmer(value ~ sex_sized * condition_sized * bin_size + (1 | litter_sized/embryo_id_sized), data = dlx2_frequency)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex_sized * condition_sized * bin_size)
# Bind all three types of comparisons into one large dataframe
posthoc_results <- bind_rows(
  
  # 1. Compares Control vs Cold (within each bin and sex)
  "Condition (Control vs Cold)" = as.data.frame(pairs(emm, by = c("bin_size", "sex_sized"), adjust = "bonferroni")),
  
  # 2. Compares Male vs Female (within each bin and condition)
  "Sex (Male vs Female)" = as.data.frame(pairs(emm, by = c("bin_size", "condition_sized"), adjust = "bonferroni")),
  
  # 3. Compares Bin vs Bin (e.g., 0 vs 1-5, within each condition and sex)
  "Bin (Between Bins)" = as.data.frame(pairs(emm, by = c("condition_sized", "sex_sized"), adjust = "bonferroni")), 
  
  # 4. Requested by JR
  "Bin and Sex Interaction" = as.data.frame(pairs(emm, by = "condition_sized", adjust = "bonferroni")),
  "Bin and Condition Interaction" = as.data.frame(pairs(emm, by = "sex_sized", adjust = "bonferroni")),
  "All vs All" = as.data.frame(pairs(emm, adjust = "bonferroni")),
  .id = "Comparison_Type")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_dlx2_frequency.csv")

write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_dlx2_frequency.csv")

# Dlx2 puncta per cell ---------------------------------------------
Dlx2_puncta <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    16.5109034, 50.1628664,
    # cold M
    12.5423729, 32.320442,
    # control F
    0, 
    # cold F
    22.5, 35.4037267,
    
    # --- LITTER 2 ---
    13.1147541, 15.18987342,
    15.3846154, 6.686046512,
    2.35294118, 7.421875,
    8.41584158, 
    
    # --- LITTER 3 ---
    #con F
    2.127659574, 
    #cold F
    33.9031339),
  
  sex = c(
    rep("M", 4), rep("F", 3),
    rep("M", 4), rep("F", 3),
    rep("M", 0), rep("F", 2)),
  
  condition = c(
    rep("control", 2), rep("cold", 2), rep("control", 1), rep("cold", 2), 
    rep("control", 2), rep("cold", 2), rep("control", 2), rep("cold", 1), 
    rep("control", 0), rep("cold", 0), rep("control", 1), rep("cold", 1)),
  
  litter = c(
    rep("L1", 7),
    rep("L2", 7),
    rep("L3", 2)))

# convert to factors
Dlx2_puncta$sex <- factor(Dlx2_puncta$sex)
Dlx2_puncta$condition <- factor(Dlx2_puncta$condition)
Dlx2_puncta$litter <- factor(Dlx2_puncta$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = Dlx2_puncta)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_Dlx2_puncta_per_cell.csv")


posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_Dlx2_puncta_per_cell.csv")

# dlx5 frequency ----------------------------------------------------------

puncta_values = c(
  # --- LITTER 1 ---
  # Litter 1, embryo 1 (Control M): counts across size bins (0 - 100)
  0.71428571,
  15,
  21.0714286,
  15,
  18.2142857,
  8.21428571,
  7.14285714,
  5.35714286,
  5,
  2.14285714,
  1.42857143,
  0,
  0.35714286,
  0,
  0.35714286,
  0,
  0,
  0,
  0,
  0,
  0,
  # Litter 1, embryo 2 (Control M)
  0.79365079,
  9.12698413,
  16.2698413,
  19.047619,
  18.2539683,
  14.2857143,
  8.33333333,
  7.14285714,
  3.17460317,
  1.98412698,
  1.19047619,
  0.3968254,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,

  # Litter 1, embryo 1 (Cold M)
  2.30263158,
  13.1578947,
  12.1710526,
  13.1578947,
  16.4473684,
  13.1578947,
  10.1973684,
  7.56578947,
  4.27631579,
  3.94736842,
  1.97368421,
  0.65789474,
  0,
  0,
  0.65789474,
  0.32894737,
  0,
  0,
  0,
  0,
  0,
  # Litter 1, embryo 2 (Cold M)
  1.42857143,
  10.7142857,
  17.8571429,
  22.5,
  15.7142857,
  15.7142857,
  8.21428571,
  4.64285714,
  2.85714286,
  0,
  0.35714286,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,

  # Litter 1, embryo 1 (Control F)
  1.24481328,
  45.6431535,
  24.8962656,
  19.9170124,
  5.39419087,
  1.65975104,
  0.82987552,
  0,
  0,
  0,
  0.41493776,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,

  # Litter 1, embryo 1 (Cold F)
  0.55096419,
  7.98898072,
  18.4573003,
  17.630854,
  15.7024793,
  12.3966942,
  7.71349862,
  7.16253444,
  4.4077135,
  4.1322314,
  0.82644628,
  1.65289256,
  0.82644628,
  0.27548209,
  0,
  0.27548209,
  0,
  0,
  0,
  0,
  0,
  # Litter 1, embryo 2 (Cold F)
  2.05278592,
  25.5131965,
  22.8739003,
  15.542522,
  11.143695,
  8.2111437,
  4.98533724,
  2.93255132,
  2.93255132,
  1.75953079,
  0.8797654,
  1.17302053,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,

  # --- LITTER 2 ---
  # Litter 2, embryo 1 (Control M): counts across size bins (0 - 100)
  1.03092784,
  15.1202749,
  22.6804124,
  21.3058419,
  17.5257732,
  7.21649485,
  6.52920962,
  4.46735395,
  1.71821306,
  1.03092784,
  0.34364261,
  0.34364261,
  0,
  0,
  0,
  0.34364261,
  0.34364261,
  0,
  0,
  0,
  0,
  # Litter 2, embryo 1 (Cold M)
  1.17647059,
  19.6078431,
  27.4509804,
  21.1764706,
  13.7254902,
  9.01960784,
  3.52941176,
  2.74509804,
  1.56862745,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,

  # Litter 2, embryo 1 (Control F)
  1.43369176,
  16.4874552,
  28.6738351,
  21.8637993,
  17.2043011,
  6.09318996,
  4.65949821,
  2.15053763,
  0.71684588,
  0.71684588,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  # Litter 2, embryo 2 (Control F)
  2.23880597,
  24.2537313,
  27.6119403,
  23.880597,
  9.32835821,
  6.34328358,
  4.10447761,
  0.37313433,
  0.74626866,
  0,
  0.37313433,
  0.37313433,
  0.37313433,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  # Litter 2, embryo 1 (Cold F)
  12.5,
  39.5833333,
  18.4027778,
  10.0694444,
  7.29166667,
  5.90277778,
  1.73611111,
  2.43055556,
  1.04166667,
  0.69444444,
  0.34722222,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  
  # --- LITTER 3 ---
  # Litter 3, embryo 1 (Control F)
  9.42408377,
  49.7382199,
  21.9895288,
  8.90052356,
  4.71204188,
  2.61780105,
  1.04712042,
  0.52356021,
  0.52356021,
  0,
  0.52356021,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,

  # Litter 3, embryo 1 (Cold F)
  0,
  10.0271003,
  18.4281843,
  19.2411924,
  19.5121951,
  10.298103,
  8.67208672,
  5.14905149,
  3.25203252,
  2.7100271,
  1.08401084,
  0.27100271,
  0,
  0.54200542,
  0.54200542,
  0,
  0,
  0,
  0,
  0.27100271,
  0)

sex = c(
  rep("M", 4), rep("F", 3),
  rep("M", 2), rep("F", 3),
  rep("M", 0), rep("F", 2))

condition = c(
  rep("control", 2), rep("cold", 2), rep("control", 1), rep("cold", 2), 
  rep("control", 1), rep("cold", 1), rep("control", 2), rep("cold", 1), 
  rep("control", 0), rep("cold", 0), rep("control", 1), rep("cold", 1))

litter = c(
  rep("L1", 7),
  rep("L2", 5),
  rep("L3", 2))

# Automatically generate sequential embryo IDs (L1_E1, L1_E2...)
embryo_id <- paste0(litter, "_E", c(1:7, 1:5, 1:2))

# 3. Build the final long-format data frame
dlx5_frequency <- data.frame(
  value = puncta_values,
  # Repeat each embryo's metadata for each bin
  sex_sized = factor(rep(sex, each = 21)),
  condition_sized = factor(rep(condition, each = 21)),
  litter_sized = factor(rep(litter, each = 21)),
  embryo_id_sized = factor(rep(embryo_id, each = 21)),
  
  # Repeat the 21 size bins 16 times (once per embryo)
  bin_size = factor(
    rep(c(
      "≤0",
      "1-5",
      "6-10",
      "11-15",
      "16-20",
      "21-25",
      "26-30",
      "31-35",
      "36-40",
      "41-45",
      "46-50",
      "51-55",
      "56-60",
      "61-65",
      "66-70",
      "71-75",
      "76-80",
      "81-85",
      "86-90",
      "91-95",
      "96-100"), times = 14), 
    levels = c(
      "≤0",
      "1-5",
      "6-10",
      "11-15",
      "16-20",
      "21-25",
      "26-30",
      "31-35",
      "36-40",
      "41-45",
      "46-50",
      "51-55",
      "56-60",
      "61-65",
      "66-70",
      "71-75",
      "76-80",
      "81-85",
      "86-90",
      "91-95",
      "96-100")))


# 4. Run mixed-effects model
model <- lmer(value ~ sex_sized * condition_sized * bin_size + (1 | litter_sized/embryo_id_sized), data = dlx5_frequency)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex_sized * condition_sized * bin_size)
# Bind all three types of comparisons into one large dataframe
posthoc_results <- bind_rows(
  
  # 1. Compares Control vs Cold (within each bin and sex)
  "Condition (Control vs Cold)" = as.data.frame(pairs(emm, by = c("bin_size", "sex_sized"), adjust = "bonferroni")),
  
  # 2. Compares Male vs Female (within each bin and condition)
  "Sex (Male vs Female)" = as.data.frame(pairs(emm, by = c("bin_size", "condition_sized"), adjust = "bonferroni")),
  
  # 3. Compares Bin vs Bin (e.g., 0 vs 1-5, within each condition and sex)
  "Bin (Between Bins)" = as.data.frame(pairs(emm, by = c("condition_sized", "sex_sized"), adjust = "bonferroni")), 
  
  # 4. Requested by JR
  "Bin and Sex Interaction" = as.data.frame(pairs(emm, by = "condition_sized", adjust = "bonferroni")),
 .id = "Comparison_Type")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_dlx5_frequency.csv")

write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_dlx5_frequency.csv")

# Dlx5 puncta per cell ---------------------------------------------
Dlx5_puncta <- data.frame(
  value = c(
    # --- LITTER 1 ---
    # control M (2 embryos)
    27.8571429, 34.9206349,
    # cold M
    39.1447368, 31.4285714,
    # control F
    2.48962656, 
    # cold F
    35.8126722, 20.8211144,
    
    # --- LITTER 2 ---
    20.9621993, 
    16.8627451, 
    14.3369176, 11.5671642,
    11.8055556, 
    
    # --- LITTER 3 ---
    #con F
    4.71204188, 
    #cold F
    30.0813008),
  
  sex = c(
    rep("M", 4), rep("F", 3),
    rep("M", 2), rep("F", 3),
    rep("M", 0), rep("F", 2)),
  
  condition = c(
    rep("control", 2), rep("cold", 2), rep("control", 1), rep("cold", 2), 
    rep("control", 1), rep("cold", 1), rep("control", 2), rep("cold", 1), 
    rep("control", 0), rep("cold", 0), rep("control", 1), rep("cold", 1)),
  
  litter = c(
    rep("L1", 7),
    rep("L2", 5),
    rep("L3", 2)))

# convert to factors
Dlx5_puncta$sex <- factor(Dlx5_puncta$sex)
Dlx5_puncta$condition <- factor(Dlx5_puncta$condition)
Dlx5_puncta$litter <- factor(Dlx5_puncta$litter)

# 4. Run mixed-effects model
model <- lmer(value ~ sex * condition + (1 | litter), data = Dlx5_puncta)

# 5. View results
summary(model)
anova(model)

# Post-hoc comparisons
emm <- emmeans(model, ~ sex * condition)
posthoc <- pairs(emm, adjust = "bonferroni")

# results
anova_results <- as.data.frame(anova(model))
write.csv(anova_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/anova_results_Dlx5_puncta_per_cell.csv")


posthoc_results <- as.data.frame(posthoc)
write.csv(posthoc_results,
          "/Users/kodeebao/Desktop/ROSIN LAB/manuscript v.2/aug revisions/results/posthoc_results_Dlx5_puncta_per_cell.csv")
