#load following packages
library(ape)
library(dplyr)
library(ggplot2)
library(tidytree)
library(treeio)
library(ggtree)
library(tidyverse)


#import nwk as ggtree object
all_tree <- read.tree("snps.strain_list.fa.treefile.nwk")
#reduce large edge lengths for better viewing
#all_tree$edge.length[all_tree$edge.length  > 0.18  ]  <- 0.1

#import associated data as ggtree object
general_info <- read.csv("gdw_annotator.csv",header = TRUE)

#import virulence/persistence data
heatmap_info <- read.csv("virpersis.tsv",sep = "\t",header = TRUE)

#tip label names to annotation
lb <- general_info$SAMPLE
anno <- general_info$SNP
d <- tibble(label=lb,label2=anno)

#link info to the tree with %<+%
#create ggtree and colour tips by sample type
p1 <- ggtree(all_tree) %<+% 
  general_info

#rotate so outgroup is on top - need to figure out which node to rotate on using
node_tree <- p1 + geom_text(aes(label=node))
p1 <- p1 %>% ggtree::rotate(18)

#change label
p2_change_label <- p1 %<+%
  d +
  geom_tiplab(aes(label=label2,
                  colour=SAMP_TYPE),
              size=2) +
  scale_color_manual(values = c(Food = 'blue', Human = 'red', Environmental = 'dark green', Outgroup = 'black', other = 'black')) +
  #ggtitle("CC217 Listeria \nOfficial Sensitive\n04.10.2024") +
  geom_treescale(y = 2) +
  theme(plot.title = element_text(colour = 'red'))

#create dataframe containing tip label and value
#heat_data <- data_frame(label=lb,gdw_check=all_info$ERROR_SMP_NON_GDW)
#for gheatmap, the row name must be the tip label name as in (all_tree$tip.label) AND NOT A COLUMN NAME!
#delete molis column in parse_genefinder output
heatmap_info <- heatmap_info %>% remove_rownames %>% column_to_rownames(var = 'sample')

p3_gheatmap <- p1 %<+%
  d +
  geom_tiplab(aes(label=label2,
                  colour=SAMP_TYPE),
              size=3.5,
              #offset = 0.01,
              align=TRUE,
              family='mono',
              linetype="dotted",
              linesize = .3) +
  scale_color_manual(values = c(Food = 'blue', Human = 'red', Environmental = 'dark green', Outgroup = 'black'), name="Sample Type") +
  ggtitle("CC217") +
  geom_treescale(y = 2) +
  theme(plot.title = element_text(colour = 'red')) +
  coord_cartesian(clip = 'off') 
#run tree first then run gheatmap second

p4 <- gheatmap(p3_gheatmap, heatmap_info,color = 'black', offset=0.12, width=0.5, font.size=3.5,colnames_angle=+90, hjust=0, colnames_position = "top") +
  scale_fill_manual(breaks = c("D","ND","U"),
                    values = c("firebrick","black","steelblue"), name="Detection State") +
  theme(legend.position=c(.1, 0.5),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.key = element_rect())


pdf("annotated_tree_pdf", width = 20, height = 20) 
p4
dev.off()
#save plot with ggsave("vpProfile.pdf", device = "pdf", path = "~/Downloads/phylogeny_requests/")