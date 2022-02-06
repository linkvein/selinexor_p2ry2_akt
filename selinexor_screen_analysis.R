# Wood Lab 2022

##############################################
# selinexor sensitizer screen analysis
# corresponds to Fig. 1b
##############################################

sensitizer_screen_data = read.csv("sensitizer_screen.csv")

# subsetting normalized seli datasets
seli1_norm = sensitizer_screen_data[,7] / sum(sensitizer_screen_data[,7])
seli2_norm = sensitizer_screen_data[,8] / sum(sensitizer_screen_data[,8])

# calculating normalized vehicle dataset (combining two vehicle pools for a more representative reference point)
veh_norm = (sensitizer_screen_data[,5] + sensitizer_screen_data[,6]) / (sum(sensitizer_screen_data[,5]) + sum(sensitizer_screen_data[,6]))

# calculating construct-wise selinexor depletion scores
seliveh_1 = seli1_norm/veh_norm
seliveh_2 = seli2_norm/veh_norm

gene_list = unique(gsub("_.*","",sensitizer_screen_data[,1]))

# collapsing construct-level scores down to gene-level scores
seliveh_genes_1 = matrix(seliveh_1, nrow=5)
colnames(seliveh_genes_1) = gene_list
seliveh_genes_2 = matrix(seliveh_2, nrow=5)
colnames(seliveh_genes_2) = gene_list

# log2 transform; this constitutes Table S1
log_seliveh_1 = log(colMeans(seliveh_genes_1, na.rm = TRUE),2)
log_seliveh_2 = log(colMeans(seliveh_genes_2, na.rm = TRUE),2)
seliveh_avg = rowMeans(cbind(log_seliveh_1, log_seliveh_2))

# plot showing replicate consistency
plot(log_seliveh_1,log_seliveh_2)

# 'snake' plot of genewise ordinals
plot(c(1:length(seliveh_avg)),seliveh_avg[order(seliveh_avg)])

##############################################
# selinexor FACS screen analysis
# corresponds to Fig. 3d
##############################################

facs_screen_data = read.csv("facs_screen.csv")
construct_ID = facs_screen_data[,2]

# subsetting normalized top and bottom datasets
seli_bot1_norm = facs_screen_data[,3] / sum(facs_screen_data[,3])
seli_top1_norm = facs_screen_data[,4] / sum(facs_screen_data[,4])
seli_bot2_norm = facs_screen_data[,5] / sum(facs_screen_data[,5])
seli_top2_norm = facs_screen_data[,6] / sum(facs_screen_data[,6])

# calculating FACS screen gene score (FSGS): ratio of bottom/top in seli-treate samples
FSGS_1 = log(seli_bot1_norm/seli_top1_norm,2)
FSGS_2 = log(seli_bot2_norm/seli_top2_norm,2)
FSGS_table = cbind(FSGS_1, FSGS_2)
rownames(FSGS_table) = construct_ID

# excluding constructs with fewer than 100 counts in untreated samples
min_table = cbind(facs_screen_data[,c(2,7:10)],as.numeric(apply(facs_screen_data[,c(7:10)],1,min)))
FSGS_table2 = FSGS_table[-which(min_table[,6] < 100),]

# identifying genes with fewer than two eligible constructs
# 737 genes with only 1 construct were excluded. 1030 w/2 constructs, 3431 w/3 constructs, 12634 w/4 constructs.
ineligible_genes = names(which(table(gsub("_.*","",rownames(FSGS_table2))) < 2))

# identifying genes that are essential
tkov3_log2_essential_genes = read.csv("tkov3_log2_essential_genes.csv")

# finalizing genes to be excluded
excluded_genes = unique(c(tkov3_log2_essential_genes[,1], ineligible_genes))

# collapsing constructs to genes by averaging
FSGS_frame = as.data.frame(FSGS_table2)
FSGS_frame = cbind(FSGS_frame, gsub("_.*","",rownames(FSGS_frame)))
colnames(FSGS_frame)[3] <- "gene"

gene_1 = aggregate(FSGS_frame$FSGS_1, list(FSGS_frame$gene), FUN=mean)
gene_2 = aggregate(FSGS_frame$FSGS_2, list(FSGS_frame$gene), FUN=mean)
gene_table = cbind(gene_1,gene_2[,2])
colnames(gene_table) = c("Gene","FSGS_1","FSGS_2")

# this constitutes Table S3
final_gene_table = gene_table[-which(gene_table[,1] %in% excluded_genes),]

# plot showing replicate consistency
plot(final_gene_table[,2],final_gene_table[,3])







