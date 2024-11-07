
required.packages <- c("Rtsne", "ggplot2", "RColorBrewer", "dplyr", "emmeans",  
                       "EmbedSOM", "tidyr", "coda", "scattermore",
                       "data.table", "ggrepel", "pheatmap", "parallelly",
                       "purrr", "lsa" )

for (req.package in required.packages){
  if(!requireNamespace(req.package, quietly=TRUE)){
    install.packages(req.package, repos='http://cran.us.r-project.org')
  }
}

bioconductor.packages <- c("ConsensusClusterPlus")

for (req.package in bioconductor.packages){
  if(!requireNamespace(req.package, quietly=TRUE)){
    BiocManager::install(req.package)
  }
}

invisible( lapply( c(required.packages, bioconductor.packages), library, character.only = TRUE ) )

flowcode.seed <- 20240915
flowcode.threads <- parallelly::availableCores() - 1
data.table::setDTthreads(flowcode.threads)

# import perCell data-----------
perCellData <- fread(file = list.files(pattern = "perCell_data.csv"))

non.procode.cells <- c( "no procode signal", "1 or 2 unexpected signals","3 or more unexpected signals")

perCellData$Procode_combination[perCellData$Procode_combination==""] <- "Untransduced"
perCellData$Id[perCellData$Id %in% non.procode.cells] <- "Other"

# add source and individual columns
perCellData <- perCellData %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5",
                                                         "X6", "mouse", "X8", "X9", "X10", "X11" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9, -X10, -X11)

# Correct the spelling in the perCellData data frame
perCellData$tissue <- gsub("Pancrease", "Pancreas", perCellData$tissue)

# List of mice for which the correction should be applied
mice_to_correct <- c("LSRB2.1a", "LSRB2.1b", "LSRB2.1c", "LSRB2.1d, LSRB2.1e") # Set1 mice

# Correct the tissue name in perCellData
perCellData <- perCellData %>%
  mutate(tissue = ifelse(tissue == "RT" & mouse %in% mice_to_correct, "Testes", tissue))

perCellData <- perCellData %>%
  mutate(tissue = ifelse(tissue == "RT", "FRT", tissue))

# optional: remove non-working guides here ###
library(stringr)

nrow(perCellData)

perCellData <- perCellData %>%
  filter(!str_detect(Id, "APC"))

nrow(perCellData)

any(str_detect(perCellData$Id, "APC"))



# read in percent matrix file-----------

pct_matrix <- read.csv(file = list.files(pattern = "pct_matrix.csv"))
pct_matrix <- pct_matrix %>% 
  separate(sample, remove = TRUE, sep = "[ _]", into = c("X1", "tissue", "X3", "X4", "X5",
                                                         "X6", "mouse", "X8", "X9", "X10", "X11" ) ) %>%
  select(-X1, -X3, -X4, -X5, -X6, -X8, -X9, -X10, -X11)

pct_matrix <- pct_matrix %>% select(-no.procode.signal, -X1.or.2.unexpected.signals, -X3.or.more.unexpected.signals ) 

# remove non-working guides here
pct_matrix <- pct_matrix %>% select(-APC)

# correct frequencies to be among transduced cells only
pct_matrix$total <- rowSums(pct_matrix[,-c(1,2)])

## Important: remove any rows without total = 0 ###
pct_matrix <- pct_matrix %>% filter(total != 0)

pct_matrix <- pct_matrix %>% mutate(across(where(is.numeric)) *100/ total)
anyNA(pct_matrix)

pct_matrix <- pct_matrix %>% 
  pivot_longer(cols = -c(tissue, mouse, total), names_to = "target_gene", values_to = "freq" )

# Correct the spelling in the pct_matrix data frame
pct_matrix$tissue <- gsub("Pancrease", "Pancreas", pct_matrix$tissue)

# Correct the tissue name in pct_matrix
pct_matrix <- pct_matrix %>%
  mutate(tissue = ifelse(tissue == "RT" & mouse %in% mice_to_correct, "Testes", tissue))

pct_matrix <- pct_matrix %>%
  mutate(tissue = ifelse(tissue == "RT", "FRT", tissue))

# calculate minimum detection frequency per sample from perCellData

tissue_counts <- perCellData %>%
  group_by( tissue, mouse, Id ) %>%
  filter( Id != "Other") %>%
  summarize( count = n()) %>%
  mutate(pct = round(count/sum(count)*100, 2)) 

tissue_mins <- tissue_counts %>%
  group_by( tissue, mouse ) %>%
  summarize( total = sum(count) ) %>%
  mutate( min_det = 100/total )

# use half-min frequency to replace zeros in percent matrix

tissue_mins$half_min <- tissue_mins$min_det/2
tissue_mins <- tissue_mins %>% unite(col = sample, tissue, mouse)

pct_matrix_corr <- pct_matrix
pct_matrix_corr <- pct_matrix_corr %>% unite(col = sample, tissue, mouse)

for (s in unique(pct_matrix_corr$sample)) {
  min_temp <- filter(tissue_mins, sample == s)
  pct_temp <- filter(pct_matrix_corr, sample == s)
  
  if (nrow(min_temp) == 1) {
    pct_temp$freq[pct_temp$freq == 0] <- min_temp$half_min
    pct_matrix_corr$freq[pct_matrix_corr$sample %in% unique(pct_temp$sample)] <- pct_temp$freq
  }
}

# Cell number cutoff: filter out sample with fewer than 20 cells

pct_matrix_corr <- left_join(pct_matrix_corr, tissue_mins, by = "sample")

pct_matrix_corr <- pct_matrix_corr %>% filter(total.y > 20)

unique(pct_matrix_corr[grep( "Brain", pct_matrix_corr$sample),]$sample)

pct_matrix_corr <- pct_matrix_corr %>% separate(sample, remove = TRUE, sep = "_", into = c("tissue", "mouse"))  %>%
  select(-total.x, -total.y, -min_det, -half_min)

anyNA(pct_matrix_corr)



## Simple tissue vs Spleen freq comparison----------------
graphs.dir <- paste0("./", flowcode.seed, "_graphs/")

if( !file.exists(graphs.dir) ){
  dir.create(graphs.dir)
}

comparison.dir <- "tissue_vs_Spleen/"
if( !file.exists( file.path(graphs.dir, comparison.dir) ) ){
  dir.create(file.path(graphs.dir, comparison.dir) )
}

reference_tissue <- "Spleen"

# set up collection of p.vals and log2FC
comparison.results <- data.frame()

for (t in unique(pct_matrix_corr$tissue)) {
  
  if(t != "Spleen"){
    tissue.dir <- paste0(graphs.dir, comparison.dir, t, "/")
    
    if( !file.exists(tissue.dir)){
      dir.create(tissue.dir)
    }
    
    toi = pct_matrix_corr %>% filter(tissue == t)
    
    spleen_toi = pct_matrix_corr %>% filter(tissue == reference_tissue)
    
    resdf = data.frame(target_gene = sort(unique(toi$target_gene)), log2FC = 0, pvalue = 1 )
    
    for (i in 1:nrow(resdf)) {
      g = resdf$target_gene[i]
      xvdf = (toi %>% filter(target_gene == g) %>% arrange(mouse))
      yvdf = spleen_toi %>% filter(target_gene == g)%>% filter(mouse %in% xvdf$mouse) %>% arrange(mouse)
      xv = (xvdf %>% filter(mouse %in% yvdf$mouse))$freq
      yv = yvdf$freq
      
      if (length(c(xv,yv)) >2){
        resdf$log2FC[i] = mean(log2((xv + 0.001) / (yv + 0.001)))
        
        res = t.test(x = xv,
                     y = yv,
                     paired = TRUE,
                     alternative = "two.sided")
        resdf$pvalue[i] = res$p.value
        
        ggplot(data.frame(toi = xv, spleen = yv, id = 1:length(xv)) %>% 
                 pivot_longer(cols = -id, names_to = "group", values_to = "freq")) +
          scale_x_discrete(name = "", labels = c(reference_tissue, t)) +
          scale_y_log10() +
          scale_fill_manual(values = c("darkgrey", "firebrick1")) +
          ggtitle(paste(g, ":", t, "vs", reference_tissue, "\n log2FC =", round(resdf$log2FC[i], 3), "; pvalue =", round(res$p.value, 3))) +
          geom_boxplot(aes(x = group, y = freq, fill = group), alpha = 0.2) +
          geom_point(aes(x = group, y = freq)) +
          geom_line(aes(x = group, y = freq, group = id)) +
          theme_bw() +
          theme(legend.position = "none")
        
        ggsave(paste0(tissue.dir, t, "_vs_", reference_tissue, "_", g, ".pdf"), height = 10, width = 10, units = "cm")
      }
    }
    
    write.csv(resdf, paste0(tissue.dir, t, "_vs_", reference_tissue, "_pvals.csv"))
    
    resdf$tissue <- rep(t)
    
    comparison.results <- rbind(comparison.results, resdf)
    
  }
  
}

# heatmap of gene impacts by tissue---------

save(comparison.results, file = "CD8_set1_comparison_results.Rds")

write.csv(comparison.results, file = "CD8_set1_comparison_results.csv")

# for any p>0.01, set log2FC to be 0
heatmap.input <- comparison.results
heatmap.input[heatmap.input$pvalue>0.01,]$log2FC <- 0
heatmap.input <- heatmap.input %>% select(-pvalue)

# reshape data into wide format: one row per gene, one column per tissue
heatmap.input <- heatmap.input %>% pivot_wider(names_from = tissue, values_from = log2FC)
heatmap.rownames <- heatmap.input$target_gene
heatmap.input <- heatmap.input %>% select(-target_gene)
heatmap.input <- as.matrix(heatmap.input)
rownames(heatmap.input) <- heatmap.rownames

palette.length <- 50
heatmap.colors <- colorRampPalette(c("blue", "white", "red"))(palette.length)
heatmap.breaks <- c(seq(min(heatmap.input), 0, length.out=ceiling(palette.length/2) + 1), 
                    seq(max(heatmap.input)/palette.length, max(heatmap.input), length.out=floor(palette.length/2)))


gene.impact.heatmap <- pheatmap(t(heatmap.input), color = heatmap.colors, 
                                breaks = heatmap.breaks,
                                cluster_rows = FALSE, cluster_cols = FALSE)

ggsave(gene.impact.heatmap, filename = paste0(graphs.dir, "Frequency_heatmap.png"), height = 10, 
       width = 20, units = "cm")

# volcano plots------------

maxmlog10pvalue <- max(-log10(comparison.results$pvalue))
maxlog2FC <- max(abs(comparison.results$log2FC))

comparison.results$dist0 <- sqrt(comparison.results$log2FC^2 + (-log10(comparison.results$pvalue) * maxlog2FC / maxmlog10pvalue)^2)

#comparison.results <- comparison.results %>% arrange(-dist0)
options(ggrepel.max.overlaps = Inf)

for (t in unique(comparison.results$tissue)) {
  
  if(t != "Spleen"){
    plot.figure <- paste0(t, "_vs_", reference_tissue, "_volcano.pdf" )
    
    resdf <- comparison.results %>% filter( tissue == t )
    resdf <- resdf %>% arrange(-dist0)
    
    ggplot(resdf) +
      ggtitle(paste(t, "vs", reference_tissue, "volcano")) +
      scale_x_continuous(limits = c(-maxlog2FC, maxlog2FC)) +
      scale_y_continuous(limits = c(0, maxmlog10pvalue))+
      geom_point(data = resdf, aes(x = log2FC, y = -log10(pvalue))) +
      geom_hline(yintercept = -log10(0.01), color = "firebrick1") +
      geom_text_repel(data = resdf %>% filter(log2FC < 0) %>% filter(pvalue < 0.01), 
                      aes(x = log2FC, y = -log10(pvalue), label = target_gene),
                      nudge_x = -0.5) +
      geom_text_repel(data = resdf %>% filter(log2FC > 0) %>% filter(pvalue < 0.01),
                      aes(x = log2FC, y = -log10(pvalue), label = target_gene),
                      nudge_x = 0.5) +
      theme_bw()
    
    ggsave(file.path(graphs.dir, comparison.dir, plot.figure), height = 10, width = 10, units = "cm")
  }
  
}




## Gene tissue profile analysis ####

tissue.profile.dir <- "./gene_tissue_profile/"

if( !file.exists( file.path(graphs.dir, tissue.profile.dir) ) ){
  dir.create(file.path(graphs.dir, tissue.profile.dir) )
}

for(g in unique(pct_matrix_corr$target_gene)){
  plot.figure <- paste0(g, "_gene_tissue_profile.pdf")
  
  ggplot(pct_matrix_corr %>% filter(target_gene == g))+
    ggtitle(g)+
    scale_y_log10()+
    geom_point(aes(x = tissue, y = freq))+
    geom_line(aes(x = tissue, y = freq, group = mouse, color = mouse))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(file.path(graphs.dir, tissue.profile.dir, plot.figure), height = 10, width = 20, units = "cm")
  
}
