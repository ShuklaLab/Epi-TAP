## Script for TCRa analysis
## Loading packages for environment----

library(dplyr)
library(ggplot2)
library(immunarch)
library(ggpubr)
library(RColorBrewer)
library(gtools)
library(ggchicklet)
library(viridis) 

## Setting working directory and import data ----
setwd("C:/TCR")

##importing TCR files
immdata = repLoad("C:/TCR/Input")

## Checking data & setting levels
immdata$meta
immdata$meta$Status <- factor(immdata$meta$Status, levels = c("Non-Hospitalised", "Hospitalised"))

## Sorting samples using meta data
desired_order <- immdata$meta$Sample
immdata$data <- immdata$data[desired_order]

## Mapping TCRs to known SARS-CoV-2 epitopes ----

##loading local vdjdb with dbLoad function - Latest dataset can be downloaded from https://vdjdb.cdr3.net/search - TRB can be used for beta chain analysis
db <- dbLoad("SARS_CDR3.tsv", "vdjdb", "HomoSapiens", "TRA",)

## Matching repertoire data to vdjdb using CDR3 sequences
res <- dbAnnotate(immdata$data, db, .data.col = c("CDR3.aa"), .db.col = c("CDR3"))

## renaming res columns
names(res)[names(res) == "CDR3.aa"] <- "CDR3"

## Combining results with db using CDR3 sequence - distinct only
merged_data <- left_join(res, db, by = "CDR3") %>%
  distinct(CDR3, .keep_all = TRUE)

## Selecting columns to add to results
res <- cbind(as.data.frame(merged_data)[c("Epitope", "Epitope gene")], res)


##export results as tables for quick overview - Not strictly necessary for pipeline
#write.csv(res, file="CDR3_epitope_mapping_all_samples.csv")
#write.csv(merged_data, file="TCRA_results_SARS_Downsampled_cohort.csv")

##Filtering immdata so individual samples only contain SARS-CoV-2 targeting TCRs ----

## Function to filter rows based on CDR3 

filter_rows <- function(sublist, sequence_list) {
  cdr3_aa <- sublist$CDR3.aa
  filtered_rows <- subset(sublist, cdr3_aa %in% sequence_list)
  return(filtered_rows) }

## Applying function for each sublist in immdata$data
for (sublist_name in desired_order) {
  immdata$data[[sublist_name]] <- filter_rows(immdata$data[[sublist_name]], res$CDR3) }

## Function to filter and append Epitope
filter_and_append_epitope <- function(sublist, res) {
  filtered_rows <- subset(sublist, CDR3.aa %in% res$CDR3)
  merged_data <- merge(filtered_rows, res[, c("CDR3", "Epitope")], 
  by.x = "CDR3.aa", by.y = "CDR3", all.x = TRUE)
  return(merged_data) }

## Applying the function for each sublist in immdata$data
for (sublist_name in desired_order) {
  immdata$data[[sublist_name]] <- filter_and_append_epitope(immdata$data[[sublist_name]], res)
}

## Creating directory for export
dir.create("1_SARS_CoV_2_TCRs", showWarnings = FALSE)

## Exporting individual results for each sample 
for (sublist_name in desired_order) {
  sublist <- immdata$data[[sublist_name]]
  output_file <- file.path("1_SARS_CoV_2_TCRs", paste0(sublist_name, ".csv"))
  write.csv(sublist, file = output_file, row.names = FALSE)
}

## TCR repertoire overview analysis----

## SARS-CoV-2 specific clone counts
clones <- repExplore(immdata$data, "clones")
status <- immdata$meta
clones <- left_join(clones, status, by = "Sample")
clones$Status <- factor(clones$Status, levels = c("Non-Hospitalised", "Hospitalised"))

clones_plot <- ggplot(clones, aes(x = Status, y = Clones, fill = Status)) +
  geom_bar(stat = "summary", fun = "mean", color = "black", linewidth = 0.5) +
  geom_point(aes(group = Sample), size = 2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
  scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
  labs(x = "Status", y = "Clone counts") +
  theme_pubr() +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_line(color = "gray", linewidth = 0.5)) +
  ylim(0, 6000)

print(clones_plot)

## immunarch function to amend plots

clones_plot %>% fixVis()

## Clonotype counts

Clonotype <- repExplore(immdata$data, "volume")
status <- immdata$meta
Clonotype <- left_join(Clonotype, status, by = "Sample")
clones$Status <- factor(clones$Status, levels = c("Non-Hospitalised", "Hospitalised"))

Clonotype_plot <- ggplot(Clonotype, aes(x = Status, y = Volume, fill = Status)) +
  geom_bar(stat = "summary", fun = "mean", color = "black", linewidth = 0.5) +
  geom_point(aes(group = Sample), size = 2) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
  scale_fill_manual(values = c("#377EB8", "#E41A1C")) +
  labs(x = "Status", y = "Unique Clonotypes") +
  theme_pubr() +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_line(color = "gray", linewidth = 0.5)) +
  ylim(0, 250)

print(Clonotype_plot)

## immunarch function to amend plots

Clonotype_plot %>% fixVis()

## Repertoire overlap

repOverlap(immdata$data,"jaccard") %>% vis(.by = "Sample", .meta = immdata$meta) %>% fixVis()


# CDR3 length distribution

repExplore(immdata$data, "len", .col = "aa", .coding = TRUE) %>% vis(.by = "Status", .meta = immdata$meta,
           .points = TRUE, .errorbars.off = TRUE, .test = FALSE, .errorbar.width = 0.0) -> CDR3_len

CDR3_len <- as.data.frame(CDR3_len$data)

# Splitting subgroups
CDR3_Hospitalised <- subset(CDR3_len, Group == "Hospitalised")
CDR3_Non_hospitalised <- subset(CDR3_len, Group == "Non-Hospitalised")

# Creating dataframe 
CDR3_H <- data.frame(CDR3_length = integer(),Value = integer(),Frequency = numeric())
CDR3_NH <- data.frame(CDR3_length = integer(),Value = integer(),Frequency = numeric())

# Loop through each row of CDR3_Non_Hospitalised - process to be repeated using CDR3_Hospitalised

for (i in 1:nrow(CDR3_Non_hospitalised)) {
  current_length <- CDR3_Non_hospitalised$Length[i]
  current_value <- CDR3_Non_hospitalised$Value[i]
  
  # Checking if current_length already exists in CDR3_NH - process to be repeated using CDR3_H
  if (current_length %in% CDR3_NH$CDR3_length) { 
    row_index <- which(CDR3_NH$CDR3_length == current_length)
    CDR3_NH$Value[row_index] <- CDR3_NH$Value[row_index] + current_value
    CDR3_NH$Frequency[row_index] <- CDR3_NH$Frequency[row_index] + 1 } 
    else { CDR3_NH <- rbind(CDR3_NH,data.frame(CDR3_length = current_length,
          Value = current_value,Frequency = 1))}}

# Calculating percentage frequency based on total of Value column - process to be repeated using CDR3_H

CDR3_NH$Frequency <- (CDR3_NH$Value / sum(CDR3_NH$Value)) *100


## CDR3 Length distribution plot - process to be repeated using CDR3_H with #E41A1C

CDR3_plot <- ggplot(CDR3_NH, aes(x = CDR3_length, y = Frequency)) +
  geom_bar(stat = "summary", fun = "mean", color = "black", fill = "#377EB8", size = 0.5, width = 1.0) +
  labs(x = "Length", y = "Value") +
  scale_fill_manual(values = "#377EB8") +
  ylim(0, 40) +
  scale_x_continuous(breaks = seq(7, 17, by = 1))

CDR3_plot %>% fixVis()


## Kmer analysis ----

## Splitting dataset by subgroup
non_hospitalised_data <- immdata$data[immdata$meta$Status == "Non-Hospitalised"]
hospitalised_data <- immdata$data[immdata$meta$Status == "Hospitalised"]

## Getting 6-mers - value can be changed to get kmers of different lengths
NH_kmers <- getKmers(non_hospitalised_data, 6)
H_kmers <- getKmers(hospitalised_data, 6)


## Custom palette for subgroups
H_palette <- c("#FFCCCC", "#FF9999", "#FF6666", "#FF3333", "#FF0000", "#CC0000", "#990000")
NH_palette <- c("#CCE5FF", "#99CCFF", "#66B2FF", "#3399FF", "#007FFF", "#0059B3", "#003366")

## Creating ggplot objects
NH_kmers <- vis(NH_kmers, .head = 20, .position = "fill") + scale_fill_manual(values = NH_palette)
H_kmers <- vis(H_kmers, .head = 20, .position = "fill") + scale_fill_manual(values = H_palette)

## Using FixVis function to amend plots
NH_kmers %>% fixVis()
H_kmers %>% fixVis()

## Exporting subgroup results for Immune viewer analysis----

## Create an empty data frame with required columns
NH_Immune_Viewer <- data.frame(
  `Read Count` = numeric(),
  Fraction = numeric(),
  `Clonal Sequence` = character(),
  `Clonal Sequence Quality` = numeric(),
  `CDR3 Min Quality` = numeric(),
  `CDR3 Sequence` = character(),
  `CDR3 Amino Acid Sequence` = character(),
  `Clonal Type` = character(),
  `Frame Shift` = character(),
  `Stop Codon` = character(),
  `Amino Acid Length` = integer(),
  `V segment` = character(),
  `all V hits` = character(),
  `D segment` = character(),
  `all D hits` = character(),
  `J segment` = character(),
  `all J hits` = character(),
  `C segment` = character(),
  `all C hits` = character(),
  stringsAsFactors = FALSE
)

# Assign the correct column names
colnames(NH_Immune_Viewer) <- c("Read Count", "Fraction", "Clonal Sequence", "Clonal Sequence Quality", 
                                "CDR3 Min Quality", "CDR3 Sequence", "CDR3 Amino Acid Sequence", "Clonal Type", "Frame Shift", 
                                "Stop Codon", "Amino Acid Length", "V segment", "all V hits", "D segment", "all D hits", 
                                "J segment", "all J hits", "C segment", "all C hits")
## Looping through sublists in non_hospitalised_data
for (sample_name in names(non_hospitalised_data)) {
  
  # Extracting sublist data
  sublist_df <- non_hospitalised_data[[sample_name]]
  
## Selecting and renaming the required columns
  renamed_df <- sublist_df %>%
    select(
      `Read Count` = Clones,
      `Clonal Sequence` = Sequence,
      `CDR3 Sequence` = CDR3.nt,
      `CDR3 Amino Acid Sequence` = CDR3.aa,
      `V segment` = V.name,
      `D segment` = D.name,
      `J segment` = J.name,
    Epitope) %>%
    mutate(Sample = sample_name) # Add sample name)
  
renamed_df <- rename_with(renamed_df, as.character) 

  
## Appending data to NH_Immune_Viewer
NH_Immune_Viewer <- bind_rows(NH_Immune_Viewer, renamed_df) 

## Trimming V and J segment data
NH_Immune_Viewer$`V segment` <- sub("\\*.*", "", NH_Immune_Viewer$`V segment`)
NH_Immune_Viewer$`J segment` <- sub("\\*.*", "", NH_Immune_Viewer$`J segment`) }

## Updating aa length and Fraction columns
NH_Immune_Viewer$`Amino Acid Length` <- nchar(NH_Immune_Viewer$`CDR3 Amino Acid Sequence`)
NH_Immune_Viewer$Fraction <- NH_Immune_Viewer$`Read Count` / sum(NH_Immune_Viewer$`Read Count`)

## Creating combined results dataframe for export
Combined_results <- NH_Immune_Viewer[, colSums(is.na(NH_Immune_Viewer)) == 0]

## Creating directory and exporting combined data 
dir.create("2_SARS_CoV_2_Combined_Results", showWarnings = FALSE)
write.csv(Combined_results,file.path("2_SARS_CoV_2_Combined_Results/NH_SARS_CoV_2_Combined_results.csv"), row.names = FALSE)

## Reformatting NH_Immune_Viewer in correct format
NH_Immune_Viewer <- NH_Immune_Viewer %>% select(-Sample, -Epitope)


## Creating directory and exporting Immune viewer data for analysis at:
## https://www.takarabio.com/products/next-generation-sequencing/bioinformatics-tools/cogent-ngs-immune-viewer
dir.create("3_Immune_Viewer_Input", showWarnings = FALSE)
write.csv(NH_Immune_Viewer,file.path("3_Immune_Viewer_Input/NH_SARS_CoV_2_CDR3_result.csv"), row.names = FALSE)



## SARS-CoV-2 epitopes bar chart----

## Creating table
NH_epitopes <- data.frame(Samples = integer(),Epitope = character(), Protein = character(),
                          Domains = character(), Clonotypes = integer(), Position = integer())


## import Data
epitopes <- read.csv("TCRa_H_Epitopes.csv")
epitopes$Position <- as.numeric(epitopes$Position)


## Colourblind friendly palette

colour_palette <- c("#362160FF", "#4143A6FF", "#4663D9FF", "#4681F7FF", "#3D9EFEFF", "#28BCEBFF", "#19D5CEFF", 
  "#1EE8B0FF", "#3DF58DFF", "#69FD66FF", "#96FE44FF", "#B6F735FF", "#D3E835FF",
  "#EBD339FF", "#FABA39FF", "#FE9C2DFF", "#F9791EFF", "#EE5711FF", "#DD3D08FF", 
  "#C52603FF", "#A81501FF", "#830702FF")

print(colour_palette)

Barchart <- ggplot(epitopes, aes(x = factor(Epitope, levels = Epitope[order(Position)]), y = Samples, fill = Protein)) +
  scale_fill_manual(breaks = c("NSP1","NSP2", "NSP3", "NSP4", "NSP5", "NSP6", "NSP8", "NSP9", "NSP10/16", "NSP12", 
  "NSP13","NSP14", "NSP15", "Spike", "ORF3a","Envelope", "Membrane","ORF6","ORF7a","ORF8", "ORF9b", "Nucleocapsid"), values = colour_palette) +
  geom_chicklet(radius = grid::unit(0.1, 'mm'), color = "black", size = 0.5) +
  labs(x = "Position along SARS-CoV-2 Proteome", y = "Individuals Expressing Clonotype", 
  title = "SARS-CoV-2 Epitopes Targeted by TCR Repertoires") +
  scale_y_continuous(breaks = seq(0, 7, by = 1),limits = c(0, 7)) +
  theme_pubr() +
  theme(axis.text.x = element_text(face = "bold")) +
  geom_text(aes(label = Clonotypes),  position = position_dodge(width = 0.01),
  vjust = +0.5, angle = 90, hjust = +1.0, color = "black", size = 4, fontface = "bold" )

print(Barchart)

## Adjust bar chart parameters 

Barchart %>% fixVis()

##Repertoire diversity barchart ----

##TRBV-J Distribution ----

Vgene <- read.csv("NH_V_gene.csv")

plot <- ggplot(Vgene, aes(x = Vgene, y = Frequency)) +
  geom_bar(stat = "summary", fun = "mean", color = "black", fill = "#377EB8", size = 0.5, width = 1.0) +
  labs(x = "TRBV Gene", y = "Frequency (%)") +
  scale_fill_manual(values = "#377EB8") +
  scale_y_continuous(breaks = seq(0, 25, by = 5),limits = c(0, 25))

print(plot)

plot %>% fixVis()

Jgene <- read.csv("NH_J_gene.csv")

plot <- ggplot(Jgene, aes(x = Jgene, y = Frequency)) +
  geom_bar(stat = "summary", fun = "mean", color = "black", fill = "#377EB8", size = 0.5, width = 1.0) +
  labs(x = "TRBJ Gene", y = "Frequency (%)") +
  scale_fill_manual(values = "#377EB8") +
  scale_y_continuous(breaks = seq(0, 25, by = 5),limits = c(0, 25))

print(plot)

plot %>% fixVis()