#LOAD REQUIRED PACKAGES##############
library(readxl)
library(tidyverse)
library(glue)
library(dplyr)
library(ggtext)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(broom)
library(purrr)

#SET ENVIRONMENT#####################################################
#Set working directory
setwd("~/Desktop/mothur_mixed_colonization_analysis/otu_analysis/")
#Set seed for all graphing
set.seed(20220212)

#READ IN METADATA#####################################################
metadata<-read_excel(path="raw_data/case_control_metadata.xlsx")
mixed_status<-read_excel(path="raw_data/mixed_status.xlsx") %>%
  filter(pr_num != 24765)
mixed_metadata<-left_join(mixed_status, metadata, by = "pr_num") %>%
  select(group, pr_num, case_control, mixed)
#ASSIGNING TAG PREFIXES#####################################################
tag_prefix<-"final_mixed_col.opti_mcc"
graph_prefix<-"mixed_colonization_status_otu"

#DEFINING AESTETICS#####################################################
group_name<-"Group"
group_breaks<-c("Yes", "No")
group_colors<-c("red", "black")
group_shapes<-c(16,15)
group_labels<-c("Mixed", "Mono")

#DEFINING SAMPLE SIZE#####################################################
#Write sample size function
run_samplesize<-function(tag){
  samplesize<-glue('./mothur "#count.groups(shared={tag}.shared, inputdir=raw_data)"')
  system(samplesize)
}

#Run samples size function
run_samplesize(tag_prefix)

#CREATING SUBSAMPLE SHARE FILE#####################################################
#Assign sample size REMEMBER TO CUSTOMIZE FOR EACH PROJECT
subsample_size<-8363

#Create subsample function
run_subsample<-function(tag, size){
  subsample<-glue('./mothur "#sub.sample(shared={tag}.shared, size={size}, inputdir=raw_data, outputdir=raw_data)"')
  system(subsample)
}

#Running subsample
run_subsample(tag_prefix, subsample_size)

#ALPHA DIVERSITY#####################################################
#Write alpha diversity function
run_alpha<-function(tag){
  alpha<-glue('./mothur "#summary.single(shared={tag}.0.03.subsample.shared, subsample=F, inputdir=raw_data, outputdir=raw_data)"')
  system(alpha)
}

#Running alpha diversity function
run_alpha(tag_prefix)

#Create alpha diversity metadata file
alpha_diversity <- read_tsv(glue("raw_data/{tag_prefix}.0.03.subsample.groups.summary")) %>%
  mutate(invsimpson = 1/simpson)
metadata_alpha <- inner_join(mixed_metadata, alpha_diversity, by="group") %>%
  na.omit()

#Test for differences in alpha diversity
alpha_pairwise_invsimpson <- pairwise.t.test(metadata_alpha$invsimpson, g=metadata_alpha$mixed, p.adjust.method = "none")
alpha_pairwise_shannon <- pairwise.t.test(metadata_alpha$shannon, g=metadata_alpha$mixed, p.adjust.method = "none")
alpha_pairwise_chao <- pairwise.t.test(metadata_alpha$chao, g=metadata_alpha$mixed, p.adjust.method = "none")

#Graph inverse simpson
metadata_alpha %>%
  ggplot(aes(x=mixed, y=invsimpson, shape=mixed, color=mixed, fill=mixed)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, size = 0.75, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(47, 47)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=48.5), fill = NA, label.color = NA, label="*p* = 0.93" ,aes(x=x, y=y), inherit.aes=FALSE, size=4) +
  labs(x=NULL, 
       y="Inverse Simpson") +
  scale_x_discrete(breaks=group_breaks,
                   labels=group_labels) +
  scale_color_manual(breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_shape_manual(breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_fill_manual(breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("graphs/{graph_prefix}_inverse_simpson.pdf"), height = 10, width = 6, unit = "cm")

#Graph shannon
metadata_alpha %>%
  ggplot(aes(x=mixed, y=shannon, shape=mixed, color=mixed, fill=mixed)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, size = 0.75, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(5.5, 5.5)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=5.75), fill = NA, label.color = NA, label="*p* = 0.50" ,aes(x=x, y=y), inherit.aes=FALSE, size=4) +
  labs(x=NULL, 
       y="Shannon") +
  scale_x_discrete(breaks=group_breaks,
                   labels=group_labels) +
  scale_color_manual(breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_shape_manual(breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_fill_manual(breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("graphs/{graph_prefix}_shannon.pdf"), height = 10, width = 6, unit = "cm")

#Graph chao
metadata_alpha %>%
  ggplot(aes(x=mixed, y=chao, shape=mixed, color=mixed, fill=mixed)) +
  geom_boxplot(show.legend = FALSE, width = 0.5, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0) + 
  geom_jitter(show.legend = FALSE, width=0.25, size = 0.75, alpha = 0.5) +
  geom_line(data=tibble(x=c(1,2), y=c(850, 850)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1.5, y=900), fill = NA, label.color = NA, label="*p* = 0.48" ,aes(x=x, y=y), inherit.aes=FALSE, size=4) +
  labs(x=NULL, 
       y="Chao") +
  scale_x_discrete(breaks=group_breaks,
                   labels=group_labels) +
  scale_color_manual(breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_shape_manual(breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_fill_manual(breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(axis.line = element_line(size = 0.5),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_blank())
ggsave(glue("graphs/{graph_prefix}_chao.pdf"), height = 10, width = 6, unit = "cm")

#BETA DIVERSITY#####################################################
otu_counts<-read_tsv(file=glue("raw_data/{tag_prefix}.0.03.subsample.shared")) %>%
  select(-label, -numOtus) %>%
  rename(group = Group) %>%
  left_join(., mixed_metadata, by = "group") %>%
  na.omit() %>%
  select(-pr_num, -case_control, -mixed) %>%
  pivot_longer(-group, names_to="otu", values_to = "count")
# set up the input data for the distance matrix calculation
dist_data_all<-otu_counts%>%
  pivot_wider(names_from = otu, values_from = count)%>%
  column_to_rownames("group")
# generate distance matrices
otu_dist_all=avgdist(dist_data_all, sample = subsample_size, dmethod = "robust.aitchison")
# amova comparisons
adonis2(otu_dist_all ~ mixed_metadata$mixed, permutations = 1000)
# generate pcoa inputs
otu_all_pcoa=cmdscale(otu_dist_all, eig=TRUE, add=TRUE)
otu_all_pcoa_positions<-otu_all_pcoa$points
colnames(otu_all_pcoa_positions) = c("axis_1", "axis_2")
#axis variability calculation
otu_all_pcoa$eig/sum(otu_all_pcoa$eig)*100
#Graph case control PCoA
otu_all_pcoa_positions %>%
  as_tibble(rownames = "group") %>%
  left_join(., mixed_metadata, by = "group") %>%
  ggplot(aes(x=axis_1, y=axis_2, color=mixed, fill=mixed, shape=mixed)) +
  labs(x=glue("PCoA Axis 1 (8.7%)"), 
       y=glue("PCoA Axis 2 (6.1%)"),
       subtitle = glue("Adjusted *p*-value = 0.34"))+
 stat_ellipse(level=0.95,
               geom="polygon",
               alpha=0.1,
               show.legend=FALSE) +
  geom_point(size = 1.5, alpha = 0.5) +
  coord_fixed(xlim=c(-15, 15), ylim = c(-15, 15)) +
  scale_color_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_shape_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_fill_manual(name=group_name,
                    breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) +
  theme_classic() +
  theme(legend.key.height = unit(0.2, "cm"),
        legend.position = "bottom",
        legend.background = element_rect(color="black", linewidth=0.4),
        legend.margin = margin(t=3, b=2, r=2, l=2),
        legend.title.align = 0.5,
        legend.title = element_text(face="bold", size=12),
        legend.text = element_markdown(size=10),
        axis.line = element_line(linewidth = 0.5),
        axis.title = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 10),
        plot.subtitle = element_markdown(size = 8))
ggsave(glue("graphs/{graph_prefix}_pcoa.pdf"), width = 12, height = 12, unit = "cm")
#READ AND JOIN OTU AND TAXONOMY FILES#####################################################
#Read in taxonomy and re-format taxon names
taxonomy<-read_tsv(file=glue("raw_data/{tag_prefix}.0.03.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";") %>%
  mutate(pretty_otu = str_replace(string=otu,
                                  pattern="otu",
                                  replacement = "otu"),
         genus = str_replace(string=genus,
                             pattern="(.*)",
                             replacement="*\\1*"),
         genus = str_replace(string=genus,
                             pattern="\\*(.*)_unclassified\\*",
                             replacement="Unclassified<br>*\\1*"),
         taxon = glue("{pretty_otu} {genus}")) %>%
  select(otu, taxon)

#Read in complete taxonomy
complete_taxonomy<-read_tsv(file=glue("raw_data/{tag_prefix}.0.03.cons.taxonomy")) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("kingdom", "phylum", "class", "order", "family", "genus"),
           sep=";")

#Create relative abundance df
otu_rel_abun<-inner_join(mixed_metadata, otu_counts, by="group") %>%
  inner_join(., taxonomy, by = "otu") %>%
  group_by(group) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) 

#Create relative abundance df sorted by taxonomic levels
otu_rel_abun_taxon_level<-inner_join(mixed_metadata, otu_counts, by="group") %>%
  inner_join(., complete_taxonomy, by = "otu") %>%
  group_by(group) %>%
  mutate(rel_abund = 100*(count / sum(count))) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("kingdom", "phylum", "class", "order", "family", "genus", "otu"),
               names_to="level",
               values_to="taxon") 
#===================================================================
# SIG OTUs
#===================================================================
#Build relative abundance data frame
#Determine significantly different asvs
# identify the high abundance asvs to avoid killing the stats
otu_pool <- otu_rel_abun %>%
  group_by(taxon, mixed) %>%
  summarize(mean=mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarize(pool = max(mean) > 1, ### <-feel free to change this cutoff
            mean = max(mean),
            .groups="drop")
# pool the low abundance asvs into "other"
otu_high_rel_abun <- inner_join(otu_rel_abun, otu_pool, by="taxon") %>%
  mutate(taxon = case_when(pool == FALSE ~ "Other", 
                           pool == TRUE ~ taxon)) %>%
  group_by(group, taxon, mixed) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = max(mean),
            .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc=FALSE))
# identify any significant asvs
sig_diff_otus <- otu_high_rel_abun%>%
  nest(data = -taxon) %>%
  mutate(experiment_tests = map(.x=data, 
                                ~pairwise.t.test(.x$rel_abund, .x$mixed) %>%
                                  tidy())) %>%
  unnest(experiment_tests) %>%
  mutate(raw_p_value = p.adjust(p.value, method ="none")) %>%
  mutate(adj_p_value = p.adjust(p.value, method ="fdr")) %>% ### <-feel free to change your correction method
  select(taxon, data, raw_p_value, adj_p_value) 

otu_rel_abun %>%
  filter(taxon == "Otu00001 *Klebsiella*") %>%
  ggplot(aes(y=rel_abund, x=taxon, shape=mixed, fill=mixed, color = mixed)) +
  geom_boxplot(width = 0.8, size=0.5, alpha=0.25, outlier.shape = NA, coef= 0, fatten = 0.75) + 
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), size = 0.75, alpha=0.5) + 
  geom_line(data=tibble(x=c(0.75,1.25), y=c(105, 105)), aes(x=x, y=y), inherit.aes=FALSE) +
  geom_richtext(data=tibble(x=1, y=110), fill = NA, label.color = NA, label="*p* = 0.14" ,aes(x=x, y=y), inherit.aes=FALSE, size=3) +
  labs(x=NULL, 
       y="Relative Abundance (%)") +
  scale_y_continuous(limits=c(0, 110),
                     breaks=c(0, 25, 50, 75, 100),
                     labels=c(0, 25, 50, 75, 100)) +
  scale_shape_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_shapes,
                     labels=group_labels) +
  scale_color_manual(name=group_name,
                     breaks=group_breaks,
                     values=group_colors,
                     labels=group_labels) +
  scale_fill_manual(name=group_name,
                    breaks=group_breaks,
                    values=group_colors,
                    labels=group_labels) + 
  theme_classic() +
  theme(legend.key.height = unit(0.75, "cm"),
        legend.position = "bottom",
        legend.background = element_rect(color="black", size=0.6),
        legend.margin = margin(t=1, b=2, r=2, l=2),
        legend.title.align = 0.5,
        legend.title = element_markdown(face="bold", size=12),
        legend.text = element_markdown(size=10),
        axis.title = element_markdown(size = 12, face = "bold"),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_markdown(color = "black", size = 10),
        axis.text.y = element_markdown(color = "black",size = 10))
ggsave(glue("graphs/{graph_prefix}_Kp.pdf"), width = 7, height = 12, unit = "cm")