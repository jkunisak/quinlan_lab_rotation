######################################################
##### Generate plots to compare CCR calculations #####
######################################################
## Load necessary packages
library(ggplot2)
library(data.table)
library(reshape2)

## Read in chr20 data 
jason_CCR = fread("~/git/quinlan_lab_rotation/new_CCR_no_sift_polyphen_filter_no_vaf_filter.txt")
jim_CCR = fread("~/git/quinlan_lab_rotation/chr20_jim_ccr.txt")
jason_CCR[which(start > 62073700 & end < 62073900)]

## Generate an ID that cats chromosome, start, and end positions
jason_CCR$CCR_ID = paste(jason_CCR$`#chrom`, paste(jason_CCR$start, jason_CCR$end, sep = "-"), sep=":")
jim_CCR$CCR_ID = paste("chr", paste(jim_CCR$`#chrom`, paste(jim_CCR$start, jim_CCR$end, sep = "-"), sep=":"), sep="")
jim_CCR$window_size = jim_CCR$end - jim_CCR$start

## Subset the CCR window datasets based on relevant columns
jason_CCR_subset = jason_CCR[,c("CCR_ID", "window_size", "gene")]
jason_CCR_subset$maker = "jason"
jim_CCR_subset = jim_CCR[,c("CCR_ID", "window_size", "gene")]
jim_CCR_subset$maker = "jim"

## First get the intersection
combined_CCR = merge(jason_CCR_subset, jim_CCR_subset, by="CCR_ID")

## Get the CCRs in Jason but NOT in Jim
jason_CCR_set = setdiff(jason_CCR_subset$CCR_ID, combined_CCR$CCR_ID)
unique_jason_CCR = jason_CCR_subset[jason_CCR_subset$CCR_ID %in% jason_CCR_set,]
colnames(unique_jason_CCR) = c("CCR_ID", "window_size.x", "maker.x", "gene.x")
unique_jason_CCR$window_size.y = 0
unique_jason_CCR$maker.y = "jim"

## Get the CCRs in Jim but NOT in Jason
jim_CCR_set = setdiff(jim_CCR_subset$CCR_ID, combined_CCR$CCR_ID)
unique_jim_CCR = jim_CCR_subset[jim_CCR_subset$CCR_ID %in% jim_CCR_set]
colnames(unique_jim_CCR) = c("CCR_ID", "window_size.y", "gene.y", "maker.y")
unique_jim_CCR$window_size.x = 0
unique_jim_CCR$maker.x = "jason"
unique_jim_CCR = unique_jim_CCR[,c("CCR_ID", "window_size.x", "maker.x", "window_size.y", "maker.y", "gene.y")]
unique_jason_CCR[window_size.y > 20,]


## Plot a scatter plot
ggplot(data = combined_CCR, aes(window_size.x, window_size.y)) + geom_point() +
  labs(x = "Jason CCR Length", y = "Jim CCR Length") +
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_smooth(method="lm", formula=y~x, se = FALSE, fullrange=FALSE) +
  theme(panel.border = element_blank()) +
  expand_limits(x=0, y=0) #+
#scale_x_discrete(limits=seq(0, max(plot_data$jason_CCR), by = 10), breaks = seq(0, max(plot_data$jason_CCR), by = 10), expand = c(0,0)) +
#scale_y_discrete(limits=seq(0, max(plot_data$jim_CCR), by = 10), breaks = seq(0, max(plot_data$jim_CCR), by = 10), expand = c(0,0))

