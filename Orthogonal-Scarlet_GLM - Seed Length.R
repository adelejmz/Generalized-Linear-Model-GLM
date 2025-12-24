# ----------------------------------------------------------------------------------
# 0. SETUP: Load Packages, Read, and Clean Data
# ----------------------------------------------------------------------------------

# Load required libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("emmeans", quietly = TRUE)) install.packages("emmeans")

library(ggplot2)
library(emmeans)
library(dplyr)

# Load the dataset
scarlet_data <- read.csv("C:/Users/adele.jamalzei/Desktop/Dissertation/HighBiomass Project/CIMMYT/Raw Data/combined_cimmyt_data_without Kelse.csv")

# Check column names
colnames(scarlet_data)

# Inspect the data
table(scarlet_data$Planting, useNA = "ifany")
table(scarlet_data$Year, useNA = "ifany")
table(scarlet_data$Planting, scarlet_data$Year)

# ----------------------------------------------------------------------------------
# 1. Function to Detect Outliers (Tukey's Rule: 1.5 * IQR) for Seed Length
# ----------------------------------------------------------------------------------
detect_outliers <- function(data, year_value, planting_value, qtl_status) {
  
  # 🔹 Filter dataset for specific Year, Planting & QTL Status
  df_subset <- data %>% filter(Year == year_value, Planting == planting_value, QTL.Status == qtl_status)
  
  # 🔹 Compute Q1, Q3 & IQR for each High Biomass Level
  outliers_df <- df_subset %>%
    group_by(High.bio.Levels) %>%
    summarise(
      Q1 = quantile(Seed.length, 0.25, na.rm = TRUE),
      Q3 = quantile(Seed.length, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      Lower_Bound = Q1 - 1.5 * IQR,
      Upper_Bound = Q3 + 1.5 * IQR
    ) %>%
    left_join(df_subset, by = "High.bio.Levels") %>%
    filter(Seed.length < Lower_Bound | Seed.length > Upper_Bound) %>%
    select(Year, Planting, QTL.Status, High.bio.Levels, Seed.length)
  
  return(outliers_df)
}

# ----------------------------------------------------------------------------------
# 1-2. Apply the Function to All Year × Planting × QTL Combinations
# ----------------------------------------------------------------------------------

# 🔹 List of unique Year & Planting conditions
year_planting_combinations <- unique(scarlet_data %>% select(Year, Planting))

# 🔹 Initialize empty dataframe for storing all outliers
all_outliers <- data.frame()

# 🔹 Loop through each Year × Planting condition
for (i in 1:nrow(year_planting_combinations)) {
  year_val <- year_planting_combinations$Year[i]
  planting_val <- year_planting_combinations$Planting[i]
  
  # 🔹 Find outliers for WITH QTL
  outliers_with_qtl <- detect_outliers(scarlet_data, year_val, planting_val, "WITH QTL")
  
  # 🔹 Find outliers for WITHOUT QTL
  outliers_without_qtl <- detect_outliers(scarlet_data, year_val, planting_val, "WITHOUT QTL")
  
  # 🔹 Combine results
  all_outliers <- bind_rows(all_outliers, outliers_with_qtl, outliers_without_qtl)
}
print(all_outliers)

# Clean and standardize data
scarlet_data$Planting <- trimws(tolower(scarlet_data$Planting))
scarlet_data$Year <- trimws(tolower(scarlet_data$Year))

# Convert to factors and define levels
scarlet_data$Planting <- factor(scarlet_data$Planting, levels = c("normal", "late"))
scarlet_data$QTL.Status <- factor(scarlet_data$QTL.Status, levels = c("WITH QTL", "WITHOUT QTL"))
scarlet_data$Year <- factor(scarlet_data$Year, levels = c("2021-2022", "2022-2023", "2023-2024"))
scarlet_data$High.bio.Levels <- factor(scarlet_data$High.bio.Levels, levels = c("500", "501", "502", "503", "504"))

# ----------------------------------------------------------------------------------
# 2. Define orthogonal contrasts
# ----------------------------------------------------------------------------------

contrasts(scarlet_data$Planting) <- cbind(Normal_vs_Late = c(1, -1))

contrasts(scarlet_data$Year) <- cbind(
  Year_2021_2022_vs_2022_2023 = c(1, -1, 0),
  Year_2022_2023_vs_2023_2024 = c(0, 1, -1)
)

contrasts(scarlet_data$QTL.Status) <- cbind(WITH_vs_WITHOUT = c(1, -1))

contrasts(scarlet_data$High.bio.Levels) <- cbind(
  `500 vs 501` = c(1, -1, 0, 0, 0),
  `(500_501) vs 502` = c(1, 1, -2, 0, 0),
  `(500_501_502) vs 503` = c(1, 1, 1, -3, 0),
  `(500_501_502_503) vs 504` = c(1, 1, 1, 1, -4)
)


# ----------------------------------------------------------------------------------
# 3. Fit the GLM Model for Seed Length (Gamma)
# ----------------------------------------------------------------------------------

# Fit GLM with Gamma distribution (for continuous positive data)
model_glm_SeedLength <- glm(
  Seed.length ~ Planting * QTL.Status * Year * High.bio.Levels, 
  data = scarlet_data, 
  family = Gamma(link = "log")
)

# ----------------------------------------------------------------------------------
# 3-1. Model Evaluation and Summary
# ----------------------------------------------------------------------------------

# Get the summary of GLM
summary_glm_SeedLength <- summary(model_glm_SeedLength)
print(summary_glm_SeedLength)

# Akaike Information Criterion (AIC) & Bayesian Information Criterion (BIC)
AIC_value <- AIC(model_glm_SeedLength)
BIC_value <- BIC(model_glm_SeedLength)
cat("AIC:", AIC_value, "\nBIC:", BIC_value, "\n")

# ----------------------------------------------------------------------------------
# 3-2. Save the Significant Results (Corrected)
# ----------------------------------------------------------------------------------

# Extract coefficients table from GLM summary
coefficients_table_glm_SeedLength <- as.data.frame(summary_glm_SeedLength$coefficients)

# Rename columns for clarity
colnames(coefficients_table_glm_SeedLength) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")

# Add a significance column based on p-values
coefficients_table_glm_SeedLength$Significance <- cut(
  round(coefficients_table_glm_SeedLength$`Pr(>|z|)`, 3), 
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", ".")
)

# Filter only significant results (p < 0.05)
significant_coefficients_glm_SeedLength <- subset(coefficients_table_glm_SeedLength, `Pr(>|z|)` < 0.05)

# Define the file path
save_path <- "C:/Users/adele.jamalzei/Desktop/Dissertation/HighBiomass Project/Data from CIMMYT/Significant_GLM_SeedLength_Below_0.05.csv"

# Ensure the directory exists
if (!dir.exists(dirname(save_path))) {
  dir.create(dirname(save_path), recursive = TRUE)
}

# Save the significant coefficients table to a CSV file
write.csv(significant_coefficients_glm_SeedLength, file = save_path, row.names = TRUE)

cat("Significant GLM regression summary for Seed Length (p < 0.05) saved successfully to:", save_path)

# ----------------------------------------------------------------------------------
# 3-3. Model Diagnostics (Checking Assumptions)
# ----------------------------------------------------------------------------------

par(mfrow = c(2, 2))  # Arrange plots
plot(model_glm_SeedLength)  # Diagnostic plots



# ----------------------------------------------------------------------------------
# 4. Visualization of GLM Result (Seed Length)
# ----------------------------------------------------------------------------------
library(ggplot2)
library(emmeans)
library(dplyr)

# ----------------------------------------------------------------------------------
# 4-1. Extract EMMs for Seed Length
# ----------------------------------------------------------------------------------

# ✅ Extract estimated marginal means (EMMs) for Seed Length
emm_results <- emmeans(model_glm_SeedLength, ~ QTL.Status * High.bio.Levels * Year * Planting, type = "response")

# ✅ Convert EMMs to a dataframe
emm_df <- as.data.frame(emm_results)

# ✅ Remove non-estimable (NA) values
emm_df <- emm_df %>% filter(!is.na(response))

# ----------------------------------------------------------------------------------
# 4-2. Define Conditions and Plot Function
# ----------------------------------------------------------------------------------

# ✅ Define correct Year × Planting conditions
valid_conditions <- data.frame(
  Year = c("2021-2022", "2022-2023", "2022-2023", "2023-2024"),
  Planting = c("late", "normal", "late", "normal")
)

# ✅ Function to plot each Year × Planting condition separately, allowing different y-axis limits
plot_qtl_effect <- function(year_value, planting_value) {
  
  # 🔹 Filter for the specific Year & Planting condition
  df_subset <- emm_df %>% filter(Year == year_value, Planting == planting_value)
  
  # 🔹 Skip plotting if no data exists for this condition
  if (nrow(df_subset) == 0) {
    message("Skipping: No data for ", year_value, " - ", planting_value)
    return(NULL)
  }
  
  # 🔹 Compute y-axis range for each plot separately
  y_min_local <- min(df_subset$response - df_subset$SE, na.rm = TRUE)  
  y_max_local <- max(df_subset$response + df_subset$SE, na.rm = TRUE) * 1.05  # Add a small buffer
  
  # 🔹 Generate the plot
  ggplot(df_subset, aes(x = factor(QTL.Status), y = response, color = factor(High.bio.Levels))) +
    
    # ✅ Points with error bars
    geom_point(position = position_dodge(width = 0.4), size = 4, alpha = 0.8) +  
    geom_errorbar(aes(ymin = response - SE, ymax = response + SE), 
                  position = position_dodge(width = 0.4), width = 0.2, size = 1.1) +
    
    # ✅ Labels & theme
    labs(title = paste("Predicted Seed Length -", year_value, "-", planting_value),
         x = "QTL Presence",
         y = "Predicted Seed Length",
         color = "High Biomass Levels") +
    
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    
    # ✅ Apply local y-axis range for each plot
    ylim(y_min_local, y_max_local)
}

# ----------------------------------------------------------------------------------
# 4-3. Generate and Print Plots for Each Condition
# ----------------------------------------------------------------------------------

for (i in 1:nrow(valid_conditions)) {
  year_val <- valid_conditions$Year[i]
  planting_val <- valid_conditions$Planting[i]
  
  plot <- plot_qtl_effect(year_val, planting_val)
  if (!is.null(plot)) print(plot)
}



# --------------------------------------------------------------------------------------
# 5. Visualization of GLM Result: Each Factor Separately (Seed Length)
# ----------------------------------------------------------------------------------
library(ggplot2)
library(emmeans)
library(dplyr)

# ✅ Convert estimated marginal means (EMMs) from GLM to a dataframe
emm_results <- emmeans(model_glm_SeedLength, ~ High.bio.Levels * Year * Planting, type = "response")
emm_df <- as.data.frame(emm_results) %>% filter(!is.na(response))  # Remove NAs

# ✅ Define "Year × Planting" column
emm_df <- emm_df %>%
  mutate(Year_Planting = factor(paste(Year, Planting, sep = " - "),
                                levels = c("2021-2022 - late", "2022-2023 - normal",
                                           "2022-2023 - late", "2023-2024 - normal")))

# ✅ Generate a combined plot for **High Biomass Levels** effect on Predicted Seed Length
ggplot(emm_df, aes(x = High.bio.Levels, y = response, color = High.bio.Levels)) +
  
  # 🔹 Points with error bars
  geom_point(position = position_dodge(width = 0.4), size = 4, alpha = 0.8) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE), 
                position = position_dodge(width = 0.4), width = 0.2, size = 1.1) +
  
  # 🔹 Facet by Year × Planting
  facet_grid(~ Year_Planting, scales = "free_x", space = "free_x") +
  
  # 🔹 Labels & Theme
  labs(
    title = "Effect of High Biomass Levels on Predicted Seed Length Across Years",
    x = "High Biomass Levels",
    y = "Predicted Seed Length",
    color = "High Biomass Levels"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_blank(),  # 🔥 Remove the 500, 501 labels from X-axis
    axis.ticks.x = element_blank(),  # 🔥 Remove X-axis ticks
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  
  # 🔹 Custom colors for High Biomass Levels
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD"))  # Adjust colors


#----------------------------------------------------------------------------------
library(ggplot2)
library(emmeans)
library(dplyr)

# ✅ Convert estimated marginal means (EMMs) from GLM to a dataframe
emm_results <- emmeans(model_glm_SeedLength, ~ QTL.Status * Year * Planting, type = "response")
emm_df <- as.data.frame(emm_results) %>% filter(!is.na(response))  # Remove NAs

# ✅ Define "Year × Planting" column
emm_df <- emm_df %>%
  mutate(Year_Planting = factor(paste(Year, Planting, sep = " - "),
                                levels = c("2021-2022 - late", "2022-2023 - normal",
                                           "2022-2023 - late", "2023-2024 - normal")))

# ✅ Generate a combined plot for **QTL Status Effect on Predicted Seed Length**
ggplot(emm_df, aes(x = QTL.Status, y = response, color = QTL.Status)) +
  
  # 🔹 Points with error bars
  geom_point(position = position_dodge(width = 0.4), size = 4, alpha = 0.8) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE), 
                position = position_dodge(width = 0.4), width = 0.2, size = 1.1) +
  
  # 🔹 Facet by Year × Planting
  facet_grid(~ Year_Planting, scales = "free_x", space = "free_x") +
  
  # 🔹 Labels & Theme
  labs(
    title = "Effect of QTL Status on Predicted Seed Length Across Years",
    x = "QTL Status",
    y = "Predicted Seed Length",
    color = "QTL Status"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_blank(),  # 🔥 Remove "WITH QTL" & "WITHOUT QTL" from x-axis
    axis.ticks.x = element_blank(), # 🔥 Remove x-axis ticks
    axis.title.y = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  
  # 🔹 Custom colors for QTL Status
  scale_color_manual(values = c("WITH QTL" = "#D62728", "WITHOUT QTL" = "#2CA02C"))  # Red & Green

# ----------------------------------------------------------------------------------
# 6. Visualizations: Violin Plots for Seed Length - Raw Data
# ----------------------------------------------------------------------------------
# 1. Load Required Libraries
library(ggplot2)
library(dplyr)

# 6-1. Function to Create Violin Plots for Seed Length

plot_qtl_biomass_per_year <- function(data, year_value, planting_value) {
  ggplot(data %>% filter(Year == year_value & Planting == planting_value), 
         aes(x = QTL.Status, y = Seed.length, fill = High.bio.Levels)) +
    
    # Adjust violin width and separation
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8, 
                color = "black", lwd = 1, 
                position = position_dodge(width = 0.9)) +  
    
    # Align boxplots properly inside violins
    geom_boxplot(width = 0.15, outlier.shape = NA, 
                 color = "black", lwd = 1, 
                 position = position_dodge(width = 0.9)) +  
    
    theme_classic(base_size = 16) +  
    labs(
      title = paste("Seed Length by QTL & High Biomass Levels -", year_value, "-", planting_value),
      x = "QTL Status",
      y = "Seed Length",
      fill = "High Biomass Level"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  
      axis.title.x = element_text(size = 16, face = "bold"),  
      axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),  
      axis.title.y = element_text(size = 16, face = "bold"),  
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold")
    ) +
    scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD"))  # Custom colors for biomass
}

# 6-2. Generate and Print Plots for Each Year × Planting

p_2021_late   <- plot_qtl_biomass_per_year(scarlet_data, "2021-2022", "late")
p_2022_normal <- plot_qtl_biomass_per_year(scarlet_data, "2022-2023", "normal")
p_2022_late   <- plot_qtl_biomass_per_year(scarlet_data, "2022-2023", "late")
p_2023_normal <- plot_qtl_biomass_per_year(scarlet_data, "2023-2024", "normal")

# ✅ Print the separate plots
print(p_2021_late)   # 2021-2022 Late
print(p_2022_normal) # 2022-2023 Normal
print(p_2022_late)   # 2022-2023 Late
print(p_2023_normal) # 2023-2024 Normal




#------------------------------------------------------------------------------------
# 7 - Plot on Each Factor Separately (Seed Length)
#------------------------------------------------------------------------------------

# Load Required Libraries
library(ggplot2)
library(dplyr)

# Define Conditions

# ✅ Define correct Year × Planting conditions
valid_conditions <- data.frame(
  Year = c("2021-2022", "2022-2023", "2022-2023", "2023-2024"),
  Planting = c("late", "normal", "late", "normal")
)

# ✅ Create a "Year × Planting" column with better formatting
scarlet_data <- scarlet_data %>%
  mutate(Year_Planting = factor(paste(Year, Planting, sep = " - "), 
                                levels = paste(valid_conditions$Year, valid_conditions$Planting, sep = " - ")))

# Create Violin Plot with Proper Faceting

ggplot(scarlet_data, aes(x = QTL.Status, y = Seed.length, fill = QTL.Status)) +
  
  # ✅ Violin plot
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black", lwd = 1) +
  
  # ✅ Boxplot inside violin plot
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", lwd = 1) +
  
  # ✅ Facet by Year × Planting with properly formatted labels
  facet_grid(~ Year_Planting, scales = "free_x", space = "free_x") + 
  
  # ✅ Customize labels and appearance
  theme_classic(base_size = 16) +  
  labs(
    title = "Effect of QTL Status on Seed Length Across Years",
    x = "QTL Status",  # 🔥 Adds the x-axis title
    y = "Seed Length",
    fill = "QTL Status"
  ) +
  
  # ✅ Improve text and layout
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  
    axis.title.x = element_text(size = 16, face = "bold"),  # 🔥 Ensures "QTL Status" appears as x-axis title
    axis.text.x = element_blank(),   # 🔥 Removes x-axis text (WITH/WITHOUT QTL) to avoid redundancy
    axis.ticks.x = element_blank(),  # 🔥 Removes x-axis ticks
    axis.title.y = element_text(size = 16, face = "bold"),  
    strip.text = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5, margin = margin(4, 4, 4, 4)), # 🔥 Ensures text fits in boxes
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  
  # ✅ Custom colors for QTL Status
  scale_fill_manual(values = c("WITH QTL" = "#D62728", "WITHOUT QTL" = "#2CA02C"))  # Red for WITH QTL, Green for WITHOUT QTL


# ----------------------------------------------------------------------------------
# Load Required Libraries
# ----------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)

# ----------------------------------------------------------------------------------
# Define Conditions
# ----------------------------------------------------------------------------------

# ✅ Define correct Year × Planting conditions
valid_conditions <- data.frame(
  Year = c("2021-2022", "2022-2023", "2022-2023", "2023-2024"),
  Planting = c("late", "normal", "late", "normal")
)

# ✅ Create a "Year × Planting" column with better formatting
scarlet_data <- scarlet_data %>%
  mutate(Year_Planting = factor(paste(Year, Planting, sep = " - "), 
                                levels = paste(valid_conditions$Year, valid_conditions$Planting, sep = " - ")))

# ----------------------------------------------------------------------------------
# Create Violin Plot for High Biomass Levels (Seed Length)
# ----------------------------------------------------------------------------------

ggplot(scarlet_data, aes(x = High.bio.Levels, y = Seed.length, fill = High.bio.Levels)) +
  
  # ✅ Violin plot
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black", lwd = 1) +
  
  # ✅ Boxplot inside violin plot
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black", lwd = 1) +
  
  # ✅ Facet by Year × Planting with properly formatted labels
  facet_grid(~ Year_Planting, scales = "free_x", space = "free_x") + 
  
  # ✅ Customize labels and appearance
  theme_classic(base_size = 16) +  
  labs(
    title = "Effect of High Biomass Levels on Seed Length Across Years",
    x = "High Biomass Level",  # 🔥 Adds the x-axis title
    y = "Seed Length",
    fill = "High Biomass Level"
  ) +
  
  # ✅ Improve text and layout
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  
    axis.title.x = element_text(size = 16, face = "bold"),  
    axis.text.x = element_blank(),   # 🔥 Removes x-axis labels (500, 501, etc.)
    axis.ticks.x = element_blank(),  # 🔥 Removes x-axis ticks
    axis.title.y = element_text(size = 16, face = "bold"),  
    strip.text = element_text(size = 14, face = "bold", hjust = 0.5, vjust = 0.5, margin = margin(4, 4, 4, 4)), # 🔥 Ensures text fits in boxes
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1.5),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  
  # ✅ Custom colors for High Biomass Levels
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD"))  # Custom colors
