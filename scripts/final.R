library(tidyverse)

## Programmer: 2.4-2.5
#scripts:
  #run_tophat.qsub (align paired sample reads to mouse genome mm9)
  #run_cufflinks.qsub (count aligned reads / normalize)
  #run_cuffdiff.qsub (identify differential expression in genes)

fpkm_vals <- read.table("~/BF528/final-project/data/P0_1_cufflinks/genes.fpkm_tracking") %>% #counts in FPKM
  tibble ()

colnames(fpkm_vals) <- fpkm_vals[1,] #colnames are in row 1
fpkm_vals <- fpkm_vals %>% filter(tss_id != "tss_id") %>% #remove header row
  mutate(across(FPKM : FPKM_conf_hi, ~as.numeric(.x))) %>%
  filter(FPKM < 500, FPKM > 0)

#number of genes: 20,487 (FPKM >0)
#20,321 on histogram (had to remove high outliers for visualization)

#plot histogram
ggplot(fpkm_vals) + geom_histogram(aes(FPKM), fill = "darkslateblue", bins = 40) +
  labs(x = "Alignments (FPKM)", y = "Count",
       title = "Distribution of Alignments (FPKM)")


## Analyst: 2.6

diff_exp <- read.table("~/BF528/final-project/data/cuffdiff_out/gene_exp.diff") %>%
  tibble()
colnames(diff_exp) <- diff_exp[1,] #colnames are in row 1
diff_exp <- diff_exp %>% filter(test_id != "test_id") %>% #remove header row
  arrange(q_value) %>% #order by ascending q value
  mutate(across(value_1:q_value, ~as.numeric(.x))) #make columns numeric

#2.6.1 table:
table1 <- diff_exp %>%
  filter(q_value == min(q_value), !is.nan(test_stat)) %>%
  arrange(desc(abs(`log2(fold_change)`))) %>%
  filter(row_number() < 11) #top 10 differentially expressed

write_csv(table1, file = "~/BF528/final-project/plots/2-6-1Table.csv")

#2.6.2 hist:
ggplot(diff_exp) +
  geom_histogram(aes(`log2(fold_change)`), fill = "darkslateblue", bins = 100) +
  labs(x = "Log2 Fold Change", y = "Count",
       title = "Log 2 Fold Change of Differentially Expressed Genes")

#2.6.3 significant df
diff_exp_sig <- diff_exp %>%
  filter(significant == "yes")
#q value falls between 0.003303 and 0.0497 ~ q_value < 0.05
#p_value falls between 0.00005 and 0.0016 ~ p_value < 0.005

#2.6.4 hist:
ggplot(diff_exp_sig) +
  geom_histogram(aes(`log2(fold_change)`), fill = "darkslateblue", bins = 100) +
  labs(x = "Log2 Fold Change", y = "Count",
       title = "Log 2 Fold Change of Significant Differentially Expressed Genes")

#2.6.5 up and down dfs
up_df <- diff_exp_sig %>%
  filter(`log2(fold_change)` > 0)
#2757 upregulated genes

down_df <- diff_exp_sig %>%
  filter(`log2(fold_change)` < 0)
#2431 downregulated genes

#write up and downregulated genes to csv for DAVID
lapply(up_df %>%
         separate_rows(gene, sep = ",") %>%
         pull(gene), write, "~/BF528/final-project/data/up_regulated.txt", append=TRUE)
lapply(down_df %>%
         separate_rows(gene, sep = ",") %>%
         pull(gene), write, "~/BF528/final-project/data/down_regulated.txt", append=TRUE)
