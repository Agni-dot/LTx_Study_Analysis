###################### mIF Analysis Script #############################
#                                                                      #
# In the same folder, this script expects to find:                     #
# a. WilcoxPlotFunctions 24.03.R functions script file                 #
# b. The output folder                                                 #
# c. The ALL-RAW.xlsx data file with:                                  #
#     One sheet for each panel, named after the panel                  #
#          1. T-cell panel (Panel 1)  = "TC"                           #
#          2. Mono/ Macro Panel (Panel 2) = "MM"                       #
#     Each sheet must have:                                            #
#       one column for class "Interface" and/or "Diagnosis"            #
#       one column "aCohort" with the cohort code (here "AC", "ES")    #
#       any number of columns starting with "C" corresponding          #
#       to the different mIF markers (CD16, CD68 etc)                  #
#       any other columns are ignored                                  #
#       (no other column should start with "C"!)                       #
#       any rows with blank first cell are ignored                     #
#                                                                      #
# For different cohorts, panels, diagnosis and interface, adjust the   #
# content of the vectors in lines 49-146                               #
#                                                                      #
########################################################################

# Required packages:
library(readxl)
library(ggpubr)
library(data.table)
library(rstatix)
library(rstudioapi)
library(dplyr)
library(openxlsx)

############################ SETTINGS #################################
# Adjust as required                                                  #
plotRes <- 300
plotWidth <- 3330 / plotRes
plotHeight <- 2400 / plotRes
subFolder <- "./output/" # it must exist before execution!!!
#                                                                     #
#######################################################################

# Set current file folder as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Caution: The WilcoxPlotFunctions 24.03.R has to be ALWAYS in the same folder
source("WilcoxPlotFunctions 24.03.R")

# Load Excel File
CELL.DATA <- "./ALL-RAW.xlsx"

# We create vectors with the panel and cohort group and the class levels to analyze
vector1 = list(
  panel = "TC",
  cohort = "AC",
  class = "Diagnosis",
  levels = c("No rejection", "Indeterminate", "Rejection"),
  colors = c("#1984c5", "#e2e2e2", "#c23728")
)
vector2 = list(
  panel = "MM",
  cohort = "AC",
  class = "Diagnosis",
  levels = c("No rejection", "Indeterminate", "Rejection"),
  colors = c("#1984c5", "#e2e2e2", "#c23728")
)
vector3 = list(
  panel = "TC",
  cohort = "ES",
  class = "Diagnosis",
  levels = c("No rejection", "Rejection"),
  colors = c("#1984c5", "#c23728")
)
vector4 = list(
  panel = "MM",
  cohort = "ES",
  class = "Diagnosis",
  levels = c("No rejection", "Rejection"),
  colors = c("#1984c5", "#c23728")
)
vector5 = list(
  panel = "TC",
  cohort = "AC",
  class = "Interface",
  levels = c("no", "mild / minimal", "moderate / severe"),
  colors = c("#ffb400", "#9080ff", "#524e7e")
)
vector6 = list(
  panel = "TC",
  cohort = "AC",
  class = "Interface",
  levels = c("no", "yes"),
  colors = c("#ffb400", "#524e7e")
)
vector7 = list(
  panel = "MM",
  cohort = "AC",
  class = "Interface",
  levels = c("no", "mild / minimal", "moderate / severe"),
  colors = c("#ffb400", "#9080ff", "#524e7e")
)
vector8 = list(
  panel = "MM",
  cohort = "AC",
  class = "Interface",
  levels = c("no", "yes"),
  colors = c("#ffb400", "#524e7e")
)
vector9 = list(
  panel = "TC",
  cohort = "ES",
  class = "Interface",
  levels = c("no", "mild / minimal", "moderate / severe"),
  colors = c("#ffb400", "#9080ff", "#524e7e")
)
vector10 = list(
  panel = "TC",
  cohort = "ES",
  class = "Interface",
  levels = c("no", "yes"),
  colors = c("#ffb400", "#524e7e")
)
vector11 = list(
  panel = "MM",
  cohort = "ES",
  class = "Interface",
  levels = c("no", "mild / minimal", "moderate / severe"),
  colors = c("#ffb400", "#9080ff", "#524e7e")
)
vector12 = list(
  panel = "MM",
  cohort = "ES",
  class = "Interface",
  levels = c("no", "yes"),
  colors = c("#ffb400", "#524e7e")
)

# Vectors are added to a list...
iteration <- list(vector1,
                  vector2,
                  vector3,
                  vector4,
                  vector5,
                  vector6,
                  vector7,
                  vector8,
                  vector9,
                  vector10,
                  vector11,
                  vector12)

# ...and we iterate the list to repeat the same steps in every vector (=group and class levels):
for (vector in iteration) {
  panel = vector$panel
  selCohort = vector$cohort
  class = vector$class
  levels = vector$levels
  custom_colors = vector$colors
  
  ClassifierCD <- read_excel(CELL.DATA, sheet = panel) # sheet name = panel MM or TC
  # count per mm2 and filter Cohort and column with Diagnosis or Interface
  ClassifierCD <- prepare_excel(ClassifierCD, selCohort, class, levels)
  
  
  # filter each ROI and create class levels for markers
  markers.whole <- count_markersB(ClassifierCD, "Whole", class, levels)
  markers.portal <- count_markersB(ClassifierCD, "Portal Tract", class, levels)
  markers.cvein <- count_markersB(ClassifierCD, "Central Vein", class, levels)
  
  # Wilcoxon
  stat.whole <- do_wilcoxon3(markers.whole, class, levels)
  stat.portal <- do_wilcoxon3(markers.portal, class, levels)
  stat.cvein <- do_wilcoxon3(markers.cvein, class, levels)
  
  # Plot
  fig.whole <- do_plot_stars(stat.whole,
                             markers.whole,
                             "Whole Slide",
                             custom_colors,
                             class)
  fig.portal <- do_plot_stars(stat.portal,
                              markers.portal,
                              "Portal Tract",
                              custom_colors,
                              class)
  fig.cvein <- do_plot_stars(stat.cvein,
                             markers.cvein,
                             "Central Vein",
                             custom_colors,
                             class)
  
  # Save Plots
  fileNameWithPath = paste0(subFolder,
                            selCohort,
                            "-",
                            panel,
                            "-",
                            class,
                            "-",
                            length(levels))
  ggsave(
    paste0(fileNameWithPath, "-whole.tiff"),
    plot = fig.whole,
    device = "tiff",
    width = plotWidth,
    height = plotHeight,
    dpi = plotRes
  )
  ggsave(
    paste0(fileNameWithPath, "-portalT.tiff"),
    plot = fig.portal,
    device = "tiff",
    width = plotWidth,
    height = plotHeight,
    dpi = plotRes
  )
  ggsave(
    paste0(fileNameWithPath, "-centralV.tiff"),
    plot = fig.cvein,
    device = "tiff",
    width = plotWidth,
    height = plotHeight,
    dpi = plotRes
  )
  
  print(paste(selCohort, panel, class, sep = "-"))
  print(stat.whole)
  print(stat.portal)
  print(stat.cvein)
  
  # Save XL with Wilcoxon results - Do keep order of stats...
  saveWS(paste0(fileNameWithPath, ".xlsx"),
         stat.whole,
         stat.portal,
         stat.cvein)
  
  # Merge plots
  print(
    ggarrange(
      fig.whole,
      fig.portal,
      fig.cvein,
      nrow = 3,
      labels = c("A", "B", "C"),
      common.legend = TRUE,
      legend = "top"
    )
  )
  
} # next vector iteration

print(warnings())