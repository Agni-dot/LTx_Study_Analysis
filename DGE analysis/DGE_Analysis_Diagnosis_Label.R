
########### NANOSTRING NORMALIZATION AND ANALYSIS PIPELINE #############
#                                                                      #
# Place this R script in the analysis folder.                          #
# This folder should also contain:                                     #
#   a) a "Nanostring_Data" subfolder containing all RCC files          #
#   b) an annotation file named "Clinical_Data.xlsx"                   #
#      reporting sample names as rownames                              #
#      and other additional columns for the study group and the        #
#      batch label (Here we used "Diagnosis")                          #
#   c) a file called "BHOT_entrez_mapping.csv" including the           #
#      BHOT panel entrezIDs list                                       #
#   d) a Results folder                                                #
#                                                                      #
#   Note: NanoString Technologies has since been acquired              #
#   by Bruker Corporation                                              #
#                                                                      #
########################################################################


######## LOAD OR INSTALL NECESSARY PACKAGES ########  


library(MASS)
  library(readxl)
  library(dplyr)
  library(tidyverse)
  library(nanostringr)
  library(RUVSeq)
  library(NanoNormIter)
  library(EDASeq)
  library(DESeq2)
  library(plotly)
  library(ggfortify)
  library(cowplot)
  library(ggrepel)
  library(gage)
  

  NANOSTRING.DATA.DIR <- "./Nanostring_Data"
  RESULTS.FOLDER <- "./Results-Diagnosis"
  CLIN.DATA <- "./Clinical_Data.xlsx"
  RUVSEQ.K <- 2 # Determine this first by manual iteration and visualization of results (integer between 1-5)
  RUVSEQ.EMPIRICAL <- T # Determine this first by manual iteration and visualization of results
  LFC.THRESHOLD <- 0.8
  PVAL.THRESHOLD <- 0.05
  BHOT.ENTREZ.IDS <- "./BHOT_entrez_mapping.csv"
  
  #Indeterminate for rejection group is sometimes handled as "borderline" in this script
  
  ######## DATA IMPORT AND PREPROCESSING ########
  
  # Reading RCC files to retrieve p_data and f_data
  rcc_files <- read_rcc(NANOSTRING_DATA.DIR)
  
  raw <- rcc_files$raw
  rownames(raw) <- raw$Name
  
  p_data <- NanoStringQC(rcc_files$raw, rcc_files$exp) %>% # Import Patient Data
    column_to_rownames(var = "File.Name")
  
  
  
  f_data <- rcc_files$raw %>% # Import probeset information
    dplyr::select(c("Name", "Code.Class", "Accession")) %>%
    column_to_rownames(var = "Name")
  colnames(f_data)[1] <- "Class"
  
  # Read previously normalized nanostring data and clinical data
  clin_data <- read_excel(CLIN.DATA, sheet="Sheet2")
  clin_data$File <- gsub(".RCC", "", clin_data$File)
  rownames(clin_data) <- clin_data$File
  
  # Filter nanostring data on samples present in clin_data
  raw <- as.data.frame(raw[,colnames(raw) %in% row.names(clin_data) | colnames(raw) %in% c("Code.Class", "Name", "Accession")]) # Filter samples that are present in the dataset
  
  # Change order of p_data and clin_data to match order in raw
  p_data <- p_data[colnames(raw[, -c(1:3)]),]
  clin_data <- clin_data[colnames(raw[, -c(1:3)]),]
  clin_data$Rejection[clin_data$Rejection=="Borderline"]<-"Indeterminate"
  p_data$Diagnosis <- clin_data$Rejection
  
  #### Flag samples with low values close to the limit of detection ####
  neg_raw <- raw[f_data$Class == "Negative", -c(1:3)] # Creates Negative Controls data frame
  lod <- colMeans(neg_raw) - apply(neg_raw, 2, sd) # Calculate LOD
  p_data$num_endogenous_blod <- colSums(raw[f_data$Class == "Endogenous", -c(1:3)] < lod) # Count endogenous genes per sample below LOD
  p_data$num_hk_blod <- colSums(raw[f_data$Class == "Housekeeping", -c(1:3)] < lod) # Count housekeeping genes per sample below LOD
  
  # Remove unwanted samples by selecting number of housekeeping genes below lod to dump the sample
  raw <- raw[, ! (colnames(raw) %in% names(p_data$num_hk_blod[p_data$num_hk_blod > 7]))] # Remove from raw counts
  p_data <- p_data[!(rownames(p_data) %in% names(p_data$num_hk_blod[p_data$num_hk_blod > 7])), ] # Remove from p_data
  clin_data <- clin_data[!(row.names(clin_data) %in% names(p_data$num_hk_blod[p_data$num_hk_blod > 7])), ] # Remove from annotations
  
  #### Verify if HK genes are associated with phenotypes or biological conditions ####
  
  # Create list of HK genes and identify samples with no expression of HK genes
  c_idx <- rownames(f_data[f_data$Class == "Housekeeping", ])
  p_data$HK_Gene_Miss <- colSums(raw[c_idx, -c(1:3)] == 0)
  
  hk_raw <- raw[c_idx, -c(1:3)] # HK genes raw matrix
  
  # Calculate association between HK and the investigated variable
  pval <- vector(length = nrow(hk_raw)) # Create p-values vector for associations
  for (i in seq_len(nrow(hk_raw))){
    reg <- glm.nb(as.numeric(hk_raw[i, ]) ~ as.factor(p_data$Diagnosis))
    pval[i] <- coef(summary(reg))[2, 4]
  }
  
  # Add p-values to the hk table
  hk_raw <- add_column(hk_raw,
                       pval,
                       .before = 1)
  
  sum(pval <= .05) # if > 0 some HK needs to be removed from normalization
  
  # Create df of the HK to remove from Normalization
  biased_hk <- hk_raw[hk_raw$pval < 0.05, ]
  
  # Eventually label associated HK genes as endogenous
  c_idx <- c_idx[!(c_idx %in% row.names(biased_hk))] # Remove from c_idx biased HK genes
  if(nrow(biased_hk) > 0){
    f_data[row.names(f_data) %in% row.names(biased_hk) & f_data$Class == "Housekeeping",]$Class <- "Endogenous"
  }
  
  # Start of RUVSeq
  norm_dat <- RUV_total(raw[,-c(1:3)], p_data, f_data, k = RUVSEQ.K) # First normalization
  
  if(RUVSEQ.EMPIRICAL == T){
    # Steps 1-3 Should be followed only at the first iteration. Successively, work only on 4 to select control genes.
    ## Step 1: Create DESeqDataSet object without including RUV W_1 correction variable
    dds <- DESeqDataSetFromMatrix(countData = counts(norm_dat$set),
                                  colData = pData(norm_dat$set),
                                  design = ~ Diagnosis)
    ## Step 2: Run DESeq without continuous variable
    dds <- DESeq(dds)
    ## Step 3: Identify and replace outliers.
    # Quote it once you run it
    dds <- replaceOutliers(dds)
    ## Step 4: Identify empirical control genes. Tweak the lfc and padj thresholds to taste
    # Determine all unique pairwise comparisons
    comparisons <- sapply( as.data.frame(combn(unique(p_data$Diagnosis),2)), function(x) as.character(x), simplify = FALSE)
    # Initialize an empty list to store results for each comparison
    results.list <- list()
    # Loop through pairwise comparisons
    for (comp in comparisons) {
      # Extract results for the current comparison
      result <- as.data.frame(results(dds,
                                      contrast = c("Diagnosis", comp[1], comp[2]),
                                      altHypothesis = "lessAbs",
                                      lfcThreshold = LFC.THRESHOLD))
      print(head(result))
      nondegs <- rownames(result[result$padj < PVAL.THRESHOLD & !(is.na(result$padj)),])
      # Store results in the list
      results.list[[paste(comp[1], ".vs.", comp[2])]] <- nondegs
    }
    # Get the genes that are nonsignificant in all three comparisons
    empirical <- Reduce(intersect, results.list)
    # Changes the gene class of the newly selected HK genes to HK
    c_idx <- unique(c(c_idx, empirical))
    f_data[row.names(f_data) %in% empirical & f_data$Class == "Endogenous",]$Class <- "Housekeeping"
    print(f_data[f_data$Class == "Housekeeping",])
    print(nrow(f_data[f_data$Class == "Housekeeping",]))
    ## Step 5: Re-normalize
    norm_dat <- RUV_total(raw[,-c(1:3)], p_data, f_data, k = RUVSEQ.K)
  }
  
  # Log transform nanostring data
  log_dat <- log(norm_dat$set@assayData$normalizedCounts+1)
  
  # Create color coding vector for batch variable
  col_batch <- vector()
  col_batch[p_data$cartridgeID == unique(p_data$cartridgeID)[1]] <- "red" # Insert Factor name/s
  col_batch[p_data$cartridgeID == unique(p_data$cartridgeID)[2]] <- "blue" # Insert Factor name/s
  col_batch[p_data$cartridgeID == unique(p_data$cartridgeID)[3]] <- "green" # Insert Factor name/s
  
  # Create color coding vector for biological variable
  col_diagnosis <- vector()
  col_diagnosis[p_data$Diagnosis == unique(p_data$Diagnosis)[1]] <- "#c23728" # Insert Factor name/s
  col_diagnosis[p_data$Diagnosis == unique(p_data$Diagnosis)[2]] <- "#1984c5" # Insert Factor name/s
  col_diagnosis[p_data$Diagnosis == unique(p_data$Diagnosis)[3]] <- "#333333" # Insert Factor name/s
  #PINK
  #col_diagnosis[p_data$Diagnosis == unique(p_data$Diagnosis)[1]] <- "#c80064" # Insert Factor name/s
  #col_diagnosis[p_data$Diagnosis == unique(p_data$Diagnosis)[2]] <- "#54bebe" # Insert Factor name/s
  #col_diagnosis[p_data$Diagnosis == unique(p_data$Diagnosis)[3]] <- "#333333" # Insert Factor name/s
  
  raw_seq <- newSeqExpressionSet(as.matrix(raw[,-c(1:3)]), phenoData=p_data, featureData=f_data)
  plotRLE(raw_seq,
          col = col_batch,
          style = "full",
          main = "",
          xlab = "Sample",
          xaxt = "n")
  plotRLE(log_dat,
          col = col_batch,
          style = "full",
          main = "",
          xlab = "Sample",
          xaxt = "n")
  plotRLE(raw_seq,
          col = col_diagnosis,
          style = "full",
          main = "",
          xlab = "Sample",
          xaxt = "n")
  plotRLE(log_dat,
          col = col_diagnosis,
          style = "full",
          main = "",
          xlab = "Sample",
          xaxt = "n")
  
  # PCA of raw counts colored by batch
  pca_raw <- prcomp(t(raw[,-c(1:3)]), scale. = T)
  autoplot(pca_raw,
           data = p_data,
           colour = "cartridgeID",
           label = F) + theme_bw()
  
  # PCA of log counts colored by batch
  pca_log <- prcomp(t(log_dat), scale. = T)
  autoplot(pca_log,
           data = p_data,
           colour = "cartridgeID",
           label = F) + theme_bw()
  
  # PCA of raw counts colored by diagnosis
  
  p_data$Diagnosis <- as.factor(p_data$Diagnosis)
  cols <- c('No' = "#1984c5", 'Indeterminate' = "#333333", 'Yes' = "#c23728")
  #cols <- c('Borderline' = "#333333", 'No' = "#54bebe", 'Yes' = "#c80064")
  
  autoplot(pca_raw,
           data = p_data,
           colour = "Diagnosis",
           label = F) + theme_bw() + scale_colour_manual(values = cols)
  
  
  # PCA of log counts colored by diagnosis
  autoplot(pca_log,
           data = p_data,
           colour = "Diagnosis",
           label = F) + theme_bw() + scale_colour_manual(values = cols)
  
  
  
  
  ######## DIFFERENTIAL GENE EXPRESSION ANALYSIS #########
  
  clin_data$Rejection[clin_data$Rejection=="Borderline"]<-"Indeterminate"
  
  #### PCA plot #########################################
  PlotPCA <- function(data, clin_data){
    pca <- prcomp(data, scale.= T) # Perform PCA. Try to change "scale." to see the effect on the plot
    df_pca <- cbind(pca$x[,1:2], clin_data$Rejection) %>% as.data.frame()
    
    df_pca$PC1 <- as.numeric(df_pca$PC1)
    df_pca$PC2 <- as.numeric(df_pca$PC2)
    #  df_pca$V3 <- as.factor(df_pca$V3)
    # Explicitly setting the factor levels to control the order in the legend
    df_pca$V3 <- factor(df_pca$V3, levels = c("Yes", "Indeterminate", "No"))
    
    var_explained_pc1 <- summary(pca)$importance[2,][[1]]
    var_explained_pc2 <- summary(pca)$importance[2,][[2]]
    p1 <- ggplot(df_pca, aes(PC1, PC2, colour = V3)) +
      geom_point(size = 3) +
      stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))),
                   data = df_pca[df_pca$V3 == "Yes" |  df_pca$V3 == "Indeterminate" | df_pca$V3 == "No" ,], size = 1) +
      labs(x=paste0("PC1 (", var_explained_pc1*100, "%)"), y=paste0("PC2 (", var_explained_pc2*100, "%)"), color="Rejection") +
      scale_colour_manual(values = cols)+
      theme(
        legend.text = element_text(size = 14), # Enlarge legend text
        legend.title = element_text(size = 16), # Enlarge legend title
        legend.key.size = unit(2, "lines"), # Increase the size of the legend keys
        axis.title.x = element_text(size = 16, face = "bold"), # Enlarge and bold the x-axis title
        axis.title.y = element_text(size = 16, face = "bold")  # Enlarge and bold the y-axis title
      )
    
    # Add density curves to y and x axis
    xdens <-
      axis_canvas(p1, axis = "x") +
      geom_density(data = df_pca, aes(x = PC1, fill = V3, colour = V3), alpha = 0.3) +
      scale_colour_manual(values = cols, aesthetics = c("colour", "fill"))
    ydens <-
      axis_canvas(p1, axis = "y", coord_flip = TRUE) +
      geom_density(data = df_pca, aes(x = PC2, fill = V3, colour = V3), alpha = 0.3) +
      scale_colour_manual(values = cols, aesthetics = c("colour", "fill")) +
      coord_flip()
    print(p1 %>%
            insert_xaxis_grob(xdens, grid::unit(1, "in"), position = "top") %>%
            insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
            ggdraw())
  }
  
  # Create new DESeqDataSet object for downstream analysis
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(norm_dat$set@assayData$normalizedCounts)),
                                colData = p_data,
                                design = ~ Diagnosis)
  
  # Run DESeq on normalized data with full design
  dds <- DESeq(dds)
  
  head(dds)
  class(dds)

  #Save normalized data for GEO submission
  count_data = dds@assays@data[[1]]
  count_data = cbind('Gene' = rownames(count_data), count_data) |> as.data.frame()
  help(write.csv)
  
  csv_path <- paste0(RESULTS.FOLDER, "GEO_Data", ".csv")
  tryCatch({
    write.csv(count_data, file = csv_path, row.names = FALSE)
    print(paste("CSV saved for", "GEO", "at:", csv_path))
  }, error = function(e) {
    print(paste("Error saving CSV for", "GEO", e$message))
  })
  
  
  
  
  
  ####################################### TASK 1 ########################################################
  #              Save significant genes from each comparison in separate lists                           #
  #                   Make PCA plots for sig genes and top 15 sig genes                                  # 
 
  # Initialize lists for storing results
  resDEGs <- list()
  sigDEGs <- list()
  topDEGs <- list()
  
  # Task 1a & 1b: Save significant genes and their log fold changes + PCA plots for significant genes
  for (i in 1:length(comparisons)) {
    comp <- comparisons[[i]]
    print(paste("Processing comparison:", comp[1], "vs", comp[2]))
    
    # Extract DEG results for the current comparison
    resDEG <- as.data.frame(results(dds, contrast = c("Diagnosis", comp[1], comp[2])))
    resDEGs[[i]] <- resDEG
    
    # Identify significant genes (based on padj and LFC thresholds)
    # sig_genes <- resDEG[resDEG$padj < PVAL.THRESHOLD & abs(resDEG$log2FoldChange) > LFC.THRESHOLD & !is.na(resDEG$padj), ]
    sig_genes <- resDEG[resDEG$padj < PVAL.THRESHOLD & !is.na(resDEG$padj), ]
    sigDEGs[[i]] <- rownames(sig_genes)
    
    if (nrow(sig_genes) > 0) {
      # Save significant genes along with their log2 fold changes and adjusted p-values
      sig_genes_with_logfc <- sig_genes[, c("log2FoldChange", "padj")]
      sig_genes_with_logfc$Gene <- rownames(sig_genes_with_logfc)
      sig_genes_with_logfc <- sig_genes_with_logfc[, c("Gene", "log2FoldChange", "padj")]
      
      # Debugging: Check the data before saving
      print(paste("Number of significant genes in comparison", comp[1], "vs", comp[2], ":", nrow(sig_genes_with_logfc)))
      print(head(sig_genes_with_logfc))  # Print a few rows for confirmation
      
      # Save the results to a CSV file
      csv_path <- paste0(RESULTS.FOLDER, "/Significant_Genes_", comp[1], "_vs_", comp[2], ".csv")
      tryCatch({
        write.csv(sig_genes_with_logfc, file = csv_path, row.names = FALSE)
        print(paste("CSV saved for", comp[1], "vs", comp[2], "at:", csv_path))
      }, error = function(e) {
        print(paste("Error saving CSV for", comp[1], "vs", comp[2], ":", e$message))
      })
      
      # Save PCA plot for significant genes
      exp_sig_genes <- t(log_dat[rownames(log_dat) %in% rownames(sig_genes), ])
      pca_plot_path <- paste0(RESULTS.FOLDER, "/PCA_sig_genes_", comp[1], "_vs_", comp[2], ".jpeg")
      jpeg(pca_plot_path, units = "cm", width = 28, height = 28, res = 300)
      PlotPCA(exp_sig_genes, clin_data)
      dev.off()
      print(paste("PCA plot saved for", comp[1], "vs", comp[2], "at:", pca_plot_path))
    } else {
      print(paste("No significant genes found for comparison:", comp[1], "vs", comp[2]))
    }
    
    # Identify top 15 significant genes for PCA plotting
    topDEG <- resDEG[order(resDEG$padj, na.last = NA), ][1:15, ]
    topDEGs[[i]] <- topDEG
    exp_topDEGs <- t(log_dat[rownames(log_dat) %in% rownames(topDEG), ])
    
    # Save PCA plot for top 15 DEGs
    top_pca_plot_path <- paste0(RESULTS.FOLDER, "/PCA_top15degs_", comp[1], "_vs_", comp[2], ".jpeg")
    jpeg(top_pca_plot_path, units = "cm", width = 28, height = 28, res = 300)
    PlotPCA(exp_topDEGs, clin_data)
    dev.off()
    print(paste("PCA plot for top 15 DEGs saved for", comp[1], "vs", comp[2], "at:", top_pca_plot_path))
  }
  

  
  #

  ####################################### TASK 2 ########################################################
  #                                 Make Volcano plots for DEGs                                         # 
  
  
  ##### Create Volcano Plots #####
  PlotVolcano <- function(resDEG, label){
    # Vulcano plot
    deg <- resDEG
    deg$diffexpressed <- "No"
    deg$diffexpressed[deg$log2FoldChange > 0 & deg$padj < 0.5] <- "Slightly_up"
    deg$diffexpressed[deg$log2FoldChange > 0 & deg$padj < 0.1] <- "Up"
    deg$diffexpressed[deg$log2FoldChange > 0 & deg$padj < 0.05] <- "Significant_up"
    deg$diffexpressed[deg$log2FoldChange < 0 & deg$padj < 0.5] <- "Slightly_down"
    deg$diffexpressed[deg$log2FoldChange < 0 & deg$padj < 0.1] <- "Down"
    deg$diffexpressed[deg$log2FoldChange < 0 & deg$padj < 0.05] <- "Significant_down"
    deg$diffexpressed <- factor(deg$diffexpressed, levels = c("Significant_up", "Up", "Slightly_up", "No", "Slightly_down", "Down", "Significant_down"));
    deg$title <- label
    deg <- arrange(deg, padj)
    deg$X <- ""
    deg$X[1:15] <- rownames(deg)[1:15]
    
    ggplot(data=deg, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=X)) +
      geom_point() +
      theme_bw() +
      theme(axis.text = element_text(size=14),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # Enlarging and bolding the plot title
            strip.text = element_text(size = 16),
            plot.subtitle = element_text(size = 16, face = "bold"), # Enlarging and bolding the subtitle
            legend.position = "none") +
      ylab("-log10(adj. p-value)") +
      geom_text_repel(size=6) +
      scale_color_manual(values = c("#970009", "#CB181D", "#FB775A",  "#C0C0C0", "#6BAED6", "#084594", "#060578")) +
      geom_hline(yintercept=-log10(0.05), lty="dashed") +
      geom_hline(yintercept=-log10(0.10), lty="dotted") +
      geom_hline(yintercept=-log10(0.50), lty="dashed") +
      # facet_grid(. ~ title) +
      #ggtitle(label) +
      #labs () is my addition
      labs(
        title = paste("Volcano Plot:", comp[1], "vs", comp[2]),
        subtitle = "Differentiation of Significant Up/Downregulation",
        color = "P-value Thresholds"
      ) +
      annotate("text", label="adj. p-value < 0.05", x=-1.3, y=(-log10(0.05)+0.05), size=6) +
      annotate("text", label="adj. p-value < 0.10", x=-1.3, y=(-log10(0.1)+0.05), size=6) +
      annotate("text", label="adj. p-value < 0.50", x=-1.3, y=(-log10(0.5)+0.05), size=6)
  }
  
  for(i in 1:length(resDEGs)){
    comp <- comparisons[[i]]
    jpeg(paste0(RESULTS.FOLDER, "/Volcano_", comp[1], "_vs_", comp[2], ".jpeg"), units="cm", width=28, height=28, res=300)
    print(PlotVolcano(resDEGs[[i]], paste0(comp[1], " vs ", comp[2])))
    dev.off()
  }
  
  ##### Pathway analysis #####
  
  ### Kegg Pathways
  bhot_entrez <- read.csv(BHOT.ENTREZ.IDS,
                          sep = ";") # Import B-HOT Panel Entrez IDs
  data(kegg.gs) # Import KEGG pathway definitions
  # determine matching probes
  bhot_entrez <- bhot_entrez[!(duplicated(bhot_entrez$Entrez.Gene) | duplicated(bhot_entrez$Entrez.Gene, fromLast = TRUE)), ] # Remove probes with unique gene names
  present_gene_names <- intersect(rownames(log_dat), bhot_entrez$Original.BHOT.names)
  log_dat <- as.data.frame(log_dat[rownames(log_dat) %in% present_gene_names,])
  bhot_entrez <- bhot_entrez[bhot_entrez$Original.BHOT.names %in% present_gene_names,] # Retain only matching probes
  log_dat <- log_dat[bhot_entrez$Original.BHOT.names,] # Reorder to match
  rownames(log_dat) <- bhot_entrez$Entrez.Gene # Reassign ## problematic because of identical entrez ids which cannot be set as rownames
  rownames(clin_data) <- clin_data$File
  
  for(comp in comparisons){
    hn <- match(clin_data$File[clin_data$Rejection == comp[2]], colnames(log_dat)) # Column Index of reference class
    dcis <- match(clin_data$File[clin_data$Rejection == comp[1]], colnames(log_dat)) # Column Index of test class
    kegg_p_vals <- gage(log_dat, # Perform Pathway analysis
                        gsets = kegg.gs,
                        ref = hn,
                        samp = dcis,
                        compare = "unpaired")
    write.csv(kegg_p_vals$greater, file=paste0(RESULTS.FOLDER, "/Kegg_inc_pathways_", comp[1], "_vs_", comp[2], ".csv")) # Identify Pathways with increased expression
    write.csv(kegg_p_vals$less, file=paste0(RESULTS.FOLDER, "/Kegg_dec_pathways_", comp[1], "_vs_", comp[2], ".csv")) #  Identify Pathways with decreased expression
  }
  
  