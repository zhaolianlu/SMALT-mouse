

library("TarCA")
# Load necessary packages
library(ape)
library(tidyverse)
library(castor)
# Read the tree
tree_files <- list.files(pattern = "^68.*substituted_tumor_treev2\\.3\\.0\\.nwk$")
# Initialize an empty list to hold all Np dataframes
list_of_df <- list()
list_of_timing <- list()
print(tree_files)
# Loop over each tree file
for (nwkfile in tree_files) {
  # Read the tree
  tree <- read.tree(nwkfile)

  # Extract terminal names
  term_names <- tree$tip.label


  # Find parent node of "ref"
	ref_parent <- NULL
	for (edge in 1:(length(tree$edge) / 2)) {
  		if (tree$tip.label[tree$edge[edge, 2]] == "ref") {
    	ref_parent <- tree$edge[edge, 1]
    	break
  	}
	}

if (!is.null(ref_parent)) {
  print(ref_parent)
}
  # Annotate tips
  dfm <- data.frame(TipLabel = term_names)
  df <- dfm %>% mutate(TipAnn = str_split(TipLabel, "_", simplify = TRUE)[, 1])
  df <- df[which(df$TipAnn != 'ref'),]
  df$TipAnn <- gsub("Ade5-4|Ade5-3|Ade5-2|Ade5-1|Ade5-5", "Ade5", df$TipAnn)
  # Get subtree and calculate Np
  subtree = get_subtree_with_tips(tree, only_tips=df$TipLabel)$subtree
  subtree <- subtree %>% fun.MergePolytomy %>% .$tree
	temp_tree <- makeNodeLabel(subtree, method = "number", prefix = "Node__")
    
tmp.result <- Np_Estimator(
	Tree = temp_tree,
	Ann = df,
	Fileout = NULL,
	ReturnNp = TRUE
	)
df_Np <- tmp.result[["EffN"]]
df_Np$file <- nwkfile
phyfile <- gsub(".substituted_tumor_treev2.3.0.nwk", "", nwkfile)
# Read in the file
phy_file <- read.delim(phyfile,header = F, sep = " ",colClasses=c("character","character"),skip = 1,col.names = c("cellID","bi"))
print(head(term_names))
# Extract tips from the tree whose labels start with 'T'
tip_labels_T <- term_names[grep("^(?!N)", term_names, perl=TRUE)]

# Extract sequence names from the phy_file using strsplit and spaces
print(head(phy_file))
# Display the number of sequences before subsetting
cat("Number of sequences before subsetting:", length(phy_file$cellID), "\n")
# Subset the phy file sequences based on the extracted tip labels
phy_file <- phy_file[which(phy_file$cellID %in% tip_labels_T),]  # +1 to adjust for the header row in phy_file

# Display the number of sequences after subsetting (subtract 1 to exclude the header)
cat("Number of sequences after subsetting:", nrow(phy_file) - 1, "\n")

# Extract cell types and sequences, excluding the header
celltypes <- gsub("\\s+.+$", "", phy_file$cellID)
sequences <- phy_file$bi
# Function to sum the numbers in a sequence
sum_sequence <- function(seq) {
  return(sum(as.numeric(unlist(strsplit(seq, "")))))
}

# Get the sum for each sequence
sums <- sapply(sequences, sum_sequence)
print("OK")
print(head(sums))
# Combine cell types and sums
df_distances_tips <- data.frame(node = celltypes, mutation_numbers = sums)
PureNode <- tmp.result[['PureNode']]
PureNode_type <- PureNode %>% select(c('Parent','TipAnn'))
PureNode_type <- unique(PureNode_type)
PureNode_distance <- merge(PureNode, df_distances_tips, by.x=c('TipLabel'),by.y=c('node'))

distance_mean <- PureNode_distance %>%
  group_by(Parent) %>%
  summarise(Mean_MB = mean(mutation_numbers, na.rm = TRUE))
    
distance_median <- PureNode_distance %>%
  group_by(Parent) %>%
  summarise(Median_MB = median(mutation_numbers, na.rm = TRUE))

distance_mean_type <- merge(distance_mean, PureNode_type, by=c('Parent'))
distance_median_type <- merge(distance_median, PureNode_type, by=c('Parent'))
distance_type <- merge(distance_mean_type, distance_median_type, by=c('Parent','TipAnn'))
mean_M0 <- 3.79
median_M0 <- 3
distance_type$Estimated_node_timing_Mean <- (34*7+20)*(3.4-distance_type$Mean_MB/mean_M0)/2.4
distance_type$Estimated_node_timing_Median <- (34*7+20)*(3.4-distance_type$Mean_MB/median_M0)/2.4

print(head(distance_type))

  # Append to list
  list_of_df[[nwkfile]] <- df_Np
      list_of_timing[[nwkfile]] <- distance_type

}

# Combine all dataframes into a single dataframe
df_all <- bind_rows(list_of_df)
timing_all <- bind_rows(list_of_timing)

# Write combined dataframe to file
write.table(df_all, "68_combined_trimed_Np.tsv", sep="\t", row.names=FALSE, na = "NA",quote=F)

write.table(timing_all, "68_combined_trimed_timing.tsv", sep="\t", row.names=FALSE, na = "NA",quote=F)



