library(readr)   
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)   
library(survival)
setwd("/fs/ess/PCON0022/Weidong/Short/XY/os")

library(readxl)

# Folder path
folder_path <- "/fs/ess/PCON0022/Weidong/Short/XY/os/Data"

# Get all .xlsx files
files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)

# Read all into a named list
data_list <- lapply(files, function(f) {
  read_excel(f, na = c("", "NA", "0"))
})

names(data_list) <- gsub("_", " ", tools::file_path_sans_ext(basename(files)))

# Example: show loaded
for (nm in names(data_list)) {
  cat("Loaded:", nm, "→ dim", dim(data_list[[nm]]), "\n")
}

i = 3
low_p = 0.5
df = data_list[[i]]

title = names(data_list)[i]

label = c("With Liraglutide", "Without Liraglutide")
pvalue_list = c("= 0.669","= 0.307","< 0.001",
                "< 0.001","= 0.152","= 0.022"
                ,"= 0.046","= 0.145","= 0.351",
                "= 0.581","< 0.001","< 0.001",
                "= 0.596","= 0.583","= 0.253",
                "= 0.473","= 0.620","= 0.012",
                "= 0.003","= 0.054")
pvalue_current = pvalue_list[i]


df_long <- df %>% 
  pivot_longer(cols = starts_with("Cohort"),
               names_to   = "Cohort",
               values_to  = "Survival") %>%   
  filter(!is.na(Survival))                    


my_cols <- c(Cohort1 = "#FFB6C1",   
             Cohort2 = "#B0E0E6")   


p <- ggplot(df_long,
            aes(x = Months/12, y = Survival/100, colour = Cohort)) +   # 👈 Months → Years
  geom_step(size = 1.1, direction = "hv") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(low_p, 1),
    breaks = seq(low_p, 1, 0.1)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(df_long$Months)/12, 1),   # 👈 tick every year
    labels = function(x) paste0(x)         # 👈 label as years
  ) +
  labs(x = "Years", y = "Overall Survival", title = title, colour = "Cohort") +
  theme_minimal(base_size = 14) +
  theme(
    axis.line        = element_line(colour = "black", size = 0.8),
    axis.ticks       = element_line(colour = "black", size = 0.6),
    axis.ticks.length= unit(0.25, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, margin = margin(b = 15))   # 👈 center the title
  ) + labs(caption = title) +
  theme(
    plot.caption = element_text(hjust = 0.5, size = 14, face = "bold")  # centered under plot
  ) +
  annotate("text", 
           x      = 0.2,   
           y      = low_p+0.05,
           label  = paste0("Log–rank P ", pvalue_current),
           hjust  = 0,
           size   = 5) +
  scale_color_manual(values = my_cols,
                     name   = "Cohort",
                     labels = label)
p
ggsave(paste0("figure/",title,".pdf"), plot = p, width = 7.5, height = 4, units = "in")
