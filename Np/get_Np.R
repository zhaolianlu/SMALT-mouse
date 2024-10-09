###This script is used to perform analysis in 'Estimating the number of progenitor cells seeding each mouse neoplasm' of methods
args = commandArgs(trailingOnly = TRUE)
library("TarCA")
# Load necessary packages
library(ape)
library(tidyverse)
library(castor)
#### Get Clade detail
fun.GetCladeSizeDetail <- function(PPPureNode2Organ){
    PPPureNode2Organ %>% 
    group_by(TipAnn,Count) %>% summarise(CCC=n()) %>% 
    arrange(TipAnn,Count) %>% 
    group_by %>% mutate(BinFreq=paste0(Count," (",CCC,")")) %>% 
    group_by(TipAnn) %>% summarise(CladeSize=paste(BinFreq,collapse = ", "))
}

fun.GetEffN <- function(PPPureNode2Organ){
    PPPureNode2Organ %>% 
    mutate(C_2inN=Count*(Count-1)/2) %>% 
    group_by(TipAnn) %>% summarise(I_num=sum(C_2inN),Total=sum(Count)) %>% 
    mutate(I_deno=Total*(Total-1)/2) %>% 
    mutate(I_index=I_num/I_deno) %>% 
    arrange(I_index) %>% 
    mutate(EffN=1/I_index) %>% 
    select(TipAnn,Total,EffN) %>% 
    arrange(-EffN)
}
# Read all tree files
tree_files <- list.files(pattern = ".*substituted_tumor_treev2\\.3\\.0\\.nwk$")

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
#sample <- sub("_filtered_re.phy.treefile$", "", nwkfile)
print(nwkfile)
# Extract unique identifiers along with T pattern from file names
sample <- sub("(_T[0-9]*).*", "\\1", nwkfile)
sample <- sub("_sampling.*", "", sample)
print(sample)
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
#df_Np <- tmp.EffN
df_Np$file <- nwkfile
df_Np <- df_Np[which(df_Np$TipAnn!='N'),]
cell_num <- nrow(df[which(df$TipAnn!="N"),])
df_Np[which(df_Np$TipAnn!='N'),]$TipAnn <- sample
df_Np$Cell_num <- cell_num
df_Np$Type <- 'Pruned_3'
  # Append to list
  list_of_df[[nwkfile]] <- df_Np

}

# Combine all dataframes into a single dataframe
df_all <- bind_rows(list_of_df)
timing_all <- bind_rows(list_of_timing)
filename <- paste('CRC_combined_trimed_3.0_Np.tsv',sep='')
# Write combined dataframe to file
write.table(df_all, filename, sep="\t", row.names=FALSE, na = "NA",quote=F)




