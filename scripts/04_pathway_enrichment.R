# pathway enrichment
library(dplyr)

# load annotated biomarkers
bio <- read.csv(
  "results/final_candidate_biomarkers_annotated.csv",
  stringsAsFactors = FALSE
)

# keep valid ENTREZ IDs only
bio_clean <- bio %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

# extract gene lists
entrez_genes <- bio_clean$ENTREZID
symbol_genes <- bio_clean$SYMBOL

length(entrez_genes)
symbol_genes


#######################
install.packages("enrichR")
library(enrichR)

dbs <- listEnrichrDbs()
head(dbs)

dbs_use <- c(
  "GO_Biological_Process_2023",
  "KEGG_2021_Human",
  "Reactome_2022"
)

enrich_res <- enrichr(symbol_genes, dbs_use)

go_bp <- enrich_res[["GO_Biological_Process_2023"]]

write.csv(
  go_bp,
  "results/GO_BP_enrichR_results.csv",
  row.names = FALSE
)

kegg <- enrich_res[["KEGG_2021_Human"]]

write.csv(
  kegg,
  "results/KEGG_enrichR_results.csv",
  row.names = FALSE
)

react <- enrich_res[["Reactome_2022"]]

write.csv(
  react,
  "results/Reactome_enrichR_results.csv",
  row.names = FALSE
)

# Visualization for GO, KEEG, Reactome
library(dplyr)

library(ggplot2)
go_bp <- enrich_res[["GO_Biological_Process_2023"]]
go_bp_top <- go_bp %>%
  dplyr::arrange(Adjusted.P.value) %>%
  head(15)

go_bp_top$Term <- factor(go_bp_top$Term, levels = rev(go_bp_top$Term))
#Bar plot
ggplot(go_bp_top, aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "GO Biological Process Enrichment",
    x = "",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal(base_size = 12)

#Dot plot
ggplot(go_bp_top,
       aes(x = -log10(Adjusted.P.value),
           y = Term,
           size = Overlap)) +
  geom_point(color = "darkred", alpha = 0.7) +
  labs(
    title = "GO BP Enrichment Dot Plot",
    x = "-log10(adjusted p-value)",
    y = "",
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 12)

#Prepare KEGG & Reactome results
kegg <- enrich_res[["KEGG_2021_Human"]]
reactome <- enrich_res[["Reactome_2022"]]

#Select top pathways
#KEGG
kegg_top <- kegg[order(kegg$Adjusted.P.value), ]
kegg_top <- head(kegg_top, 15)

kegg_top$Term <- factor(kegg_top$Term, levels = rev(kegg_top$Term))

#KEGG BAR PLOT
ggplot(kegg_top, aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +
  labs(
    title = "KEGG Pathway Enrichment",
    x = "",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal(base_size = 12)


#Reactome
reactome_filt <- reactome[
  reactome$Overlap >= 3 & reactome$Adjusted.P.value < 0.05,
]

react_top <- reactome_filt[order(reactome_filt$Adjusted.P.value), ]
react_top <- head(react_top, 10)

react_top$Term <- factor(react_top$Term, levels = rev(react_top$Term))

react_top$Term <- gsub(" R-HSA-[0-9]+", "", react_top$Term)

#REACTOME BAR PLOT

ggplot(react_top,
       aes(x = -log10(Adjusted.P.value),
           y = Term,
           size = Overlap)) +
  geom_point(color = "purple", alpha = 0.7) +
  labs(
    title = "Reactome Pathway Enrichment (Filtered)",
    x = "-log10(adjusted p-value)",
    y = "",
    size = "Gene count"
  ) +
  theme_minimal(base_size = 12)





#REACTOME BAR PLOT
ggplot(react_top, aes(x = Term, y = -log10(Adjusted.P.value))) +
  geom_bar(stat = "identity", fill = "purple") +
  coord_flip() +
  labs(
    title = "Reactome Pathway Enrichment",
    x = "",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal(base_size = 12)












