# functions to be used by "withFunctions-Data_Analysis.R" script


#' prepare_excel
#'
#' Preproccessing of the dataframe from Excel file - selects specified Cohort Rows
#' and Class Column and counts Markers per mm2, according to Area mm2
#'
#' @param ClassifierCD 
#' @param selCohort 
#' @param class 
#'
#' @return
#'
prepare_excel <- function(ClassifierCD, selCohort, class="Diagnosis", levels)
{
  # skip rows with blank first cell to remove last lines with totals etc:
  ClassifierCD <- ClassifierCD[!(is.na(ClassifierCD$Image) |
                                   ClassifierCD$Image == ""), ]
  
  columns_to_remove <- c("CD8+CD3-", "CD8+", "CD163+CD68-", "C4d+", "CK7+",
                         "CD56+", "CD14+CD16-", "CD16+CD14-", "CD16+CD56-", "CD56+CD16-")
  existing_columns <- intersect(columns_to_remove, colnames(ClassifierCD))
  
  ClassifierCD <- ClassifierCD[, setdiff(colnames(ClassifierCD), existing_columns)]
  
  # select columns with cd...
  ClassifierCD.cd <- select(ClassifierCD, starts_with('C'))
  
  # calculate cd... per mm2
  ClassifierCD.cd <- (1000 * 1000 * ClassifierCD.cd / ClassifierCD$Area_µm2)
  
  # rebuild the dataframe with the necessary columns
  ClassifierCD <- cbind(ClassifierCD[, c("Image", "aCohort", "ROI", class)], ClassifierCD.cd)
  
  ClassifierCD <- filter(ClassifierCD, aCohort == selCohort)
  ClassifierCD <- subset(ClassifierCD, select = -c(aCohort))
  
  
  if (length(levels) == 2 & class == "Interface")
    ClassifierCD <- ClassifierCD %>%
    mutate(Interface = case_when(
      Interface %in% c("yes, minimal", "yes, many") ~ "yes",
      # Map both to "yes"
      TRUE ~ Interface  # Keep all other values as they are
    ))
  
  if (length(levels) == 2 & class=="Diagnosis")
    ClassifierCD <- filter(ClassifierCD, Diagnosis != "Borderline")
  
  
  if (length(levels) == 3 & class == "Diagnosis")
    ClassifierCD <- ClassifierCD %>%
    mutate(Diagnosis = case_when(
      Diagnosis %in% c("Borderline") ~ "Indeterminate",
      # Map both to "yes"
      TRUE ~ Diagnosis  # Keep all other values as they are
    ))
  
  if (length(levels) == 3 & class == "Interface")
    ClassifierCD <- ClassifierCD %>%
    mutate(
      Interface = case_when(
        Interface %in% c("yes, minimal") ~ "mild / minimal",
        Interface %in% c("yes, many") ~ "moderate / severe",
        TRUE ~ Interface # Keep original value if it doesn't match any condition
    ))
  
  
  return(ClassifierCD)
}



#' count_markersB 
#' 
#' New version - can handle any Class with 2-3 levels
#'
#' @param data 
#' @param selROI 
#' @param class 
#' @param levels 
#'
#' @return dataframe with Markers for the selected ROI and Class
#'
count_markersB <- function(data, selROI,
                           class = "Diagnosis",
                           levels = c("Rejection", "Borderline", "No rejection")) {
  
  data <- filter(data, ROI == selROI )
  data <- subset(data, select = -c(ROI))
  
  markers <- as.data.frame(melt(
    setDT(data),
    id.vars = c("Image", class),
    variable.name = "Marker",
    value.name = "Value"
  ))
  markers[[class]] <- factor(markers[[class]],
                              levels = levels)
  return(markers)
}


#' convert_p_to_stars
#'
#' Function to convert p-values to stars
#'
#' @param p 
#'
#' @return a string with 1-3 stars
#'
convert_p_to_stars <- function(p) {
  if (p < 0.001)
    "***"
  else if (p < 0.01)
    "**"
  else if (p < 0.05)
    "*"
  else
    ""
}

#' do_plot_stars
#'
#' Creates the boxplot with stars - can adjust size and colour
#' Needs annotation positioning info stored at stat dataframe.
#'
#' @param stat 
#' @param markers 
#' @param title 
#' @param custom_colors 
#' @param class 
#'
#' @return plot fig
#'
#' @examples fig.cvein <- do_plot_stars(stat.cvein, markers.cvein, "Central Vein", custom_colors, class)
do_plot_stars <- function(stat, markers, title, custom_colors, class="Diagnosis") {
  plot_outliers=19 # set to NA to remove them
  
  max_value <- max(markers$Value, na.rm = TRUE)
  max_annotation <- max(stat$y.position, na.rm = TRUE) * 1.15
  
  # Update stat for labels
  stat <- stat %>%
    mutate(
      label = sapply(p, convert_p_to_stars),
      empty_label = "",
      x = as.numeric(as.factor(Marker)),
      # Assuming 'Marker' is a factor
      y = y.position
    )
  
  fig <- ggboxplot(
    markers,
    x = "Marker",
    y = "Value",
    fill = class,
    palette = custom_colors,
    outlier.shape = plot_outliers 
  ) +
    stat_pvalue_manual(
      stat %>% filter(label != ""), # Filter out empty labels
      label = "empty_label" ,
      # Not using the label here for p-value
      tip.length = 0.03,
      bracket.nudge.y = -3
    ) +
    geom_text(
      data = stat,
      aes((xmin+xmax)/2, y.position, label = label),
      color = "red",
      size = 5,
      vjust = +0.1
    ) +
    theme(axis.text.x = element_text(angle = 90),
          #axis.text.y = element_text(size = 20),
          plot.title = element_text(face = "bold"),
          #legend.text = element_text(size = 18), # Adjust this value to increase the legend text size
          #legend.key.size = unit(3, 'lines'), # Adjust this value to change the size of the legend keys
          legend.title = element_text(face = "bold"), # Increase legend title size and make it bold
          axis.title.y = element_text(face = "bold"), # Enlarge and bold the y-axis label
          axis.title.x = element_text(face = "bold") # Enlarge and bold the x-axis title
    ) +
    labs(y = "Number of cells per mm^2",
         title = paste(title, panel, selCohort, sep = " - ")) +
    ylim(0, max(max_value, max_annotation))
  
  return(fig)
}


#' do_plot - version that plots p values instead of stars
#' 
#' Creates the boxplot with black stars or p-values (default)
#' Needs annotation positioning info stored at stat dataframe.
#'
#' @param stat 
#' @param markers 
#' @param title 
#' @param custom_colors 
#' @param plot_stars 
#'
#' @return plot figure
#'
#' @examples
do_plot <- function(stat,
                    markers,
                    title,
                    custom_colors,
                    plot_stars = FALSE) {
  plot_outliers=19 # set to NA to remove them
  
  if (plot_stars) {
    stat <- stat %>%
      mutate(label = sapply(p, convert_p_to_stars))
    sel_label = "label"
  } else{
    sel_label = "p"
  }
  
  max_value <- max(markers$Value, na.rm = TRUE)
  max_annotation <- max(stat$y.position, na.rm = TRUE) * 1.15
  
  fig <- ggboxplot(
    markers,
    x = "Marker",
    y = "Value",
    fill = "Diagnosis",
    palette = custom_colors,
    outlier.shape = plot_outliers  
  ) +
    stat_pvalue_manual(
      stat,
      label = sel_label, #Plot p-value or stars
      tip.length = 0.01,
      bracket.nudge.y = -2
    ) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(y = "Number of cells per mm^2",
         title = paste(title, panel, selCohort, sep = " - ")) +
    ylim(0, max(max_value, max_annotation)) # extend y-axis to accommodate p-values
  
  return(fig)
}


#' do_wilcoxon3
#'
#' Counts the p-values for each marker 
#' Counts dynamically the annotation positioning for the plot
#' Returns dataframe with statistical important p-values - may return an empty dataframe
#' Handles up to 3 levels per class
#' If there's insufficient non-zero variance in any Marker, it does not proceed with wilcox.
#'
#' @param markers 
#' @param class 
#' @param levels 
#'
#' @return statistics dataframe with markers with p-values <0.05 and annotation positioning locations
#'
#' @examples stat.cvein <- do_wilcoxon3(markers.cvein, class, levels)
do_wilcoxon3 <- function(markers, class, levels) {
  
  stat <- markers %>%
    group_by(Marker) %>%
    # Filter out groups where all values are zero or there's insufficient non-zero variance
    filter(any(Value != 0) && n_distinct(Value[Value != 0]) > 2) %>%
    # Now apply the Wilcoxon test
    wilcox_test(as.formula(paste("Value", class, sep = "~")))
#  %>% filter(p < 0.05)
  
  y_range <- max(markers$Value, na.rm = TRUE)
  
  # Calculate positions for each level within each Marker
  stat <- stat %>%
    #  group_by(Marker) %>%
    mutate(
      xmin = as.numeric(as.factor(Marker))  + (  as.numeric(factor(group1, levels = levels)) -2)*
        (0.07 * length(levels) +0.06),
      # Adjust based on layout
      xmax = as.numeric(as.factor(Marker))  + ( as.numeric(factor(group2, levels = levels))-length(levels)+1) *
        (0.07 * length(levels) +0.06)
    ) %>%
    ungroup()
  
  # Calculate positions dynamically
  position_data <- markers %>%
    group_by(Marker) %>%
    summarise(
      Q1 = quantile(Value, 0.25, na.rm = TRUE),
      Q3 = quantile(Value, 0.75, na.rm = TRUE),
      IQR = IQR(Value, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      upper_whisker = 0.9 * Q3 + 0 * IQR,
      y.position = 1.1 * y_range  # Add % of the y-axis range
    )
  
  # Set fixed annotation height, according to max
  position_data$y.position=max(position_data$y.position)
  
  # Merge position data back to stat.whole
  stat <- merge(stat, position_data, by = "Marker")
  
  # Adjust y.position based on the number of p-values per Marker
  # Using cumsum to correctly handle increments within each group
  stat <- stat %>%
    group_by(Marker) %>%
    mutate(increment = (cumsum(p < 0.05) - 1) * y_range/8,
           # Increment for each subsequent p-value
           y.position = y.position + increment) %>%
    ungroup()
  return(stat)
}


#' Saves 3 dataframes as sheets in a Excel Workbook
#' 
#' @param fileWithPath 
#' @param df1 
#' @param df2
#' @param df3
#'
#' @examples saveWS(paste0(fileNameWithPath, ".xlsx"), stat.whole, stat.portal, stat.cvein)
#' 
saveWS<-function(fileWithPath, df1, df2, df3) {
  # Create a new Excel workbook
  wb <- createWorkbook()
  
  # Add data frames to different sheets
  addWorksheet(wb, "whole")  # Sheet name: Sheet1
  writeData(wb, "whole", df1)  # Write df1 to Sheet1
  
  addWorksheet(wb, "portalT")  # Sheet name: Sheet2
  writeData(wb, "portalT", df2)  # Write df2 to Sheet2
  
  addWorksheet(wb, "centralV")  # Sheet name: Sheet3
  writeData(wb, "centralV", df3)  # Write df3 to Sheet3
  
  # Save the workbook as an Excel file
  saveWorkbook(wb, fileWithPath, overwrite = TRUE)
  
}
