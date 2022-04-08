
library(tidyverse)
library(ggalluvial)
library(sf)
library(gganimate)
library(cowplot)

#### General Data ####
west_hex_maps <- parlitools::west_hex_map 
eight_colors <- c('#2CBDFE','#47DBCD',
                  '#F3A0F2','#9D2EC5',
                  '#ECD4F5', '#F9D094',
                  '#661D98','#F5B14C')

model_data <- readxl::read_excel("Clustering_Results_2.xlsx")
X_gb_with_snp <- read_csv("X_gb_with_snp.csv")

clusters <- list.files("All Model Results") %>% 
  map_dfr(~read_csv(paste("All Model Results",.x,sep="/")) %>% 
            mutate(Model = str_remove(.x,".csv"),
                   `Model ID` = str_remove(Model, "model_") %>% as.numeric(),
                   Cluster = as.character(Cluster))) %>% 
  left_join(model_data, by = "Model ID")

hac_umap_snp_false <- list.files("HAC_UMAP_SNP_False") %>% 
  map_dfr(~read_csv(paste("HAC_UMAP_SNP_False",.x,sep="/")) %>% 
            mutate(Model = str_remove(.x,".csv"),
                   `Model ID` = str_remove(Model, "HAC_Clustering_") %>% as.numeric(),
                   Cluster = as.character(Cluster))) 

hac_pca_snp_true <- list.files("HAC_PCA_SNP_True") %>% 
  map_dfr(~read_csv(paste("HAC_PCA_SNP_True",.x,sep="/")) %>% 
            mutate(Model = str_remove(.x,".csv"),
                   `Model ID` = str_remove(Model, "HAC_Clustering_") %>% as.numeric(),
                   Cluster = as.character(Cluster))) 


#### All ####

all_pca <- clusters %>% 
  filter(`Dimensionality Reduction` == "PCA") %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  mutate(Label = paste(`Algorithm`,`Feature Selection`,`Metric`, sep="\n")) %>% 
  st_as_sf() %>% 
  {ggplot(., aes(fill = Cluster)) + 
      facet_wrap(~Label, ncol = 6) + 
      geom_sf(color = "gray85", size = rel(0.1)) + 
      theme_classic() + 
      scale_fill_viridis_d(option = "A", direction = -1) +
      labs(title = "PCA") + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = rel(0.7)),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")}

all_umap <- clusters %>% 
  filter(`Dimensionality Reduction` == "UMAP") %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  mutate(Label = paste(`Algorithm`,`Feature Selection`,`Metric`, sep="\n")) %>% 
  st_as_sf() %>% 
  {ggplot(., aes(fill = Cluster)) + 
      facet_wrap(~Label, ncol = 6) + 
      geom_sf(color = "gray85", size = rel(0.1)) + 
      theme_classic() + 
      scale_fill_viridis_d(option = "A", direction = -1) +
      labs(title = "UMAP") + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = rel(0.7)),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none")}

plot_grid(all_pca, all_umap, nrow = 1)


#### K-Means ####

clusters %>% 
  filter(Algorithm == "K-Means", Model %in% c("model_03","model_35","model_17","model_49")) %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  mutate(Label = paste(`Dimensionality Reduction`,`Feature Selection`,`Metric`,sep="\n")) %>% 
  st_as_sf() %>% 
  {ggplot(., aes(fill = as.factor(Cluster))) + 
      facet_wrap(~Label, nrow = 1) + 
      geom_sf(color = "gray85", size = rel(0.1), show.legend = FALSE) + 
      theme_classic() + 
      scale_fill_viridis_d(option = "A", direction = -1) + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank())} 

test_kmeans <- clusters %>% 
  select(ID:Model,-Colors) %>% 
  spread(Model, Cluster) %>% 
  count(model_03, model_35, model_17, model_49) %>% 
  arrange(desc(n)) %>% 
  group_by(model_17) %>% 
  mutate(Aux_Order = max(n)) %>% 
  ungroup %>% 
  mutate(model_17 = reorder(model_17, Aux_Order))

test_kmeans %>% 
  ggplot(aes(axis1 = model_35, axis2 = model_03, axis3 = model_17, axis4 = model_49, 
             y = n)) + 
  geom_alluvium(aes(fill = model_49)) + 
  geom_stratum() + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

#### HAC ####

clusters %>% 
  filter(`Algorithm` == "HAC") %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  mutate(Label = paste(`Dimensionality Reduction`, `Feature Selection`,`Metric`,`Model ID`,sep="\n")) %>% 
  st_as_sf() %>% 
  {ggplot(., aes(fill = Cluster)) + 
      facet_wrap(~Label) + 
      geom_sf(color = "gray85", size = rel(0.1)) + 
      theme_classic() + 
      scale_fill_viridis_d(option = "A", direction = -1) +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = "none")}


# UMAP SNP False
set.seed(412365)
aux_hac <- hac_umap_snp_false %>% 
  select(-Colors,-`Model ID`) %>% 
  spread(Model,Cluster) %>% 
  count(across(starts_with("HAC_Clustering"))) %>% 
  arrange(desc(n)) %>% 
  mutate(Col_12 = sample(viridis::inferno(12,direction = -1)),
         Cluster_12 = LETTERS[1+as.integer(HAC_Clustering_12)]) %>% 
  gather(Model,Cluster,starts_with("HAC_Clustering")) %>% 
  mutate(K = str_sub(Model, start = -2) %>% as.integer()) %>% 
  select(Model,K,n,Cluster,Cluster_12,Col_12)

hac_colors_df <- aux_hac %>% 
  select(-K) %>% 
  spread(Model,Cluster) %>% 
  # Color 11
  group_by(HAC_Clustering_11) %>% 
  mutate(Col_11 = sample(Col_12,1, prob = n)) %>% 
  # Color 10
  group_by(HAC_Clustering_10) %>% 
  mutate(Col_10 = sample(Col_11,1, prob = n)) %>% 
  # Color 9
  group_by(HAC_Clustering_09) %>% 
  mutate(Col_09 = sample(Col_10,1, prob = n)) %>% 
  # Color 8
  group_by(HAC_Clustering_08) %>% 
  mutate(Col_08 = sample(Col_09,1, prob = n)) %>% 
  # Color 7
  group_by(HAC_Clustering_07) %>% 
  mutate(Col_07 = sample(Col_08,1, prob = n)) %>% 
  # Color 6
  group_by(HAC_Clustering_06) %>% 
  mutate(Col_06 = sample(Col_07,1, prob = n)) %>% 
  # Color 5
  group_by(HAC_Clustering_05) %>%
  mutate(Col_05 = sample(Col_06,1, prob = n)) %>%
  # Color 4
  group_by(HAC_Clustering_04) %>%
  mutate(Col_04 = sample(Col_05,1, prob = n)) %>%
  # Color 3
  group_by(HAC_Clustering_03) %>%
  mutate(Col_03 = sample(Col_04,1, prob = n)) %>%
  # Color 2
  group_by(HAC_Clustering_02) %>%
  mutate(Col_02 = sample(Col_03,1, prob = n)) %>%
  ungroup() %>% 
  select(Cluster_12,starts_with("Col"),everything()) %>% 
  gather(Model,Cluster,starts_with("HAC_Clustering")) %>% 
  gather(Aux_Color, Color, starts_with("Col")) %>% 
  filter(str_sub(Model,start = -2) == str_sub(Aux_Color,start=-2)) %>% 
  mutate(K = str_sub(Aux_Color,start=-2))

test_hac <- hac_umap_snp_false %>% 
  left_join(hac_colors_df) %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  st_as_sf() 

for(k in 2:12){
  test_hac %>% 
    filter(Model == paste("HAC_Clustering",str_pad(k,2,"left",0),sep="_")) %>% 
    {ggplot(.,aes(fill = Color, group = K)) + 
        geom_sf(color = "gray85", size = rel(0.1)) +
        scale_fill_identity() +
        theme_void()} %>% 
    plot()
}


plot <- test_hac %>% 
  ggplot(aes(fill = Color, group = K)) + 
  geom_sf(color = "gray85", size = rel(0.1)) +
  transition_states(K, transition_length = 1, state_length = 4, wrap = FALSE) + 
  scale_fill_identity() + 
  theme_void() + 
  ggtitle("HAC with K = {closest_state} for the UMAP without %SNP dataset")

animate(plot, nframes = 200, rewind = TRUE)

# PCA SNP True
set.seed(412365)
aux_hac <- hac_pca_snp_true %>% 
  select(-Colors,-`Model ID`) %>% 
  spread(Model,Cluster) %>% 
  count(across(starts_with("HAC_Clustering"))) %>% 
  arrange(desc(n)) %>% 
  mutate(Col_12 = sample(viridis::inferno(12,direction = -1)),
         Cluster_12 = LETTERS[1+as.integer(HAC_Clustering_12)]) %>% 
  gather(Model,Cluster,starts_with("HAC_Clustering")) %>% 
  mutate(K = str_sub(Model, start = -2) %>% as.integer()) %>% 
  select(Model,K,n,Cluster,Cluster_12,Col_12)

hac_colors_df <- aux_hac %>% 
  select(-K) %>% 
  spread(Model,Cluster) %>% 
  # Color 11
  group_by(HAC_Clustering_11) %>% 
  mutate(Col_11 = sample(Col_12,1, prob = n)) %>% 
  # Color 10
  group_by(HAC_Clustering_10) %>% 
  mutate(Col_10 = sample(Col_11,1, prob = n)) %>% 
  # Color 9
  group_by(HAC_Clustering_09) %>% 
  mutate(Col_09 = sample(Col_10,1, prob = n)) %>% 
  # Color 8
  group_by(HAC_Clustering_08) %>% 
  mutate(Col_08 = sample(Col_09,1, prob = n)) %>% 
  # Color 7
  group_by(HAC_Clustering_07) %>% 
  mutate(Col_07 = sample(Col_08,1, prob = n)) %>% 
  # Color 6
  group_by(HAC_Clustering_06) %>% 
  mutate(Col_06 = sample(Col_07,1, prob = n)) %>% 
  # Color 5
  group_by(HAC_Clustering_05) %>%
  mutate(Col_05 = sample(Col_06,1, prob = n)) %>%
  # Color 4
  group_by(HAC_Clustering_04) %>%
  mutate(Col_04 = sample(Col_05,1, prob = n)) %>%
  # Color 3
  group_by(HAC_Clustering_03) %>%
  mutate(Col_03 = sample(Col_04,1, prob = n)) %>%
  # Color 2
  group_by(HAC_Clustering_02) %>%
  mutate(Col_02 = sample(Col_03,1, prob = n)) %>%
  ungroup() %>% 
  select(Cluster_12,starts_with("Col"),everything()) %>% 
  gather(Model,Cluster,starts_with("HAC_Clustering")) %>% 
  gather(Aux_Color, Color, starts_with("Col")) %>% 
  filter(str_sub(Model,start = -2) == str_sub(Aux_Color,start=-2)) %>% 
  mutate(K = str_sub(Aux_Color,start=-2))

test_hac <- hac_pca_snp_true %>% 
  left_join(hac_colors_df) %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  st_as_sf() 

for(k in 2:12){
  test_hac %>% 
    filter(Model == paste("HAC_Clustering",str_pad(k,2,"left",0),sep="_")) %>% 
    {ggplot(.,aes(fill = Color, group = K)) + 
        geom_sf(color = "gray85", size = rel(0.1)) +
        scale_fill_identity() +
        theme_void()} %>% 
    plot()
}


plot <- test_hac %>% 
  ggplot(aes(fill = Color, group = K)) + 
  geom_sf(color = "gray85", size = rel(0.1)) +
  transition_states(K, transition_length = 1, state_length = 4, wrap = FALSE) + 
  scale_fill_identity() + 
  theme_void() + 
  ggtitle("HAC with K = {closest_state} for the PCA with %SNP dataset")

animate(plot, nframes = 200, rewind = TRUE)

#### DBScan ####

aux_eight_colors <- c("steelblue4",rev(eight_colors))
aux_eight_colors[4] <- "firebrick"
clusters %>% 
  filter(Model %in% c("model_11","model_42","model_58")) %>% 
  mutate(aux_index = as.integer(Cluster) + 2,
         aux_color = map_chr(aux_index,~aux_eight_colors[.x]),
         Cluster_Color = case_when(Model == "model_42" ~ 
                                     if_else(Cluster == "-1", "steelblue4", eight_colors[8]),
                                   Model == "model_11" ~ 
                                     if_else(Cluster == "1", eight_colors[2], eight_colors[8]),
                                   T ~ aux_color)) %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  mutate(Label = paste(`Dimensionality Reduction`,`Feature Selection`,`Metric`,sep="\n")) %>% 
  st_as_sf() %>% 
  {ggplot(., aes(fill = Cluster_Color)) + 
      facet_wrap(~Label, nrow = 1) + 
      geom_sf(color = "gray85", size = rel(0.1), show.legend = FALSE) + 
      theme_classic() + 
      scale_fill_identity() + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank())} 

#### GMM #####

clusters %>% filter(Algorithm == "GMM") %>% 
  mutate(Weight_Tier = cut(Weights, c(0, seq(0.6,1, by = 0.05))), 
         Label = paste(`Dimensionality Reduction`,`Feature Selection`,`Metric`,`Model ID`, sep="\n")) %>% 
  filter(Weights <= 0.9) %>% 
  count(Label, Weight_Tier) %>% 
  ggplot(aes(x = Weight_Tier, y = n, color = n, label = n)) + 
  facet_wrap(~Label) + 
  geom_segment(aes(xend = Weight_Tier, y = 0, yend = n)) + 
  geom_point() + 
  geom_label(aes(label = n)) + 
  theme_classic() + 
  coord_flip()

clusters %>% 
  filter(Model %in% c("model_45","model_13")) %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  mutate(Weight_Tier = cut(Weights, c(0, seq(0.5, 0.9, by = 0.1), 0.95, 1)), 
         Label = paste(`Dimensionality Reduction`,`Feature Selection`,`Metric`, sep=" + ")) %>% 
  st_as_sf() %>% 
  {ggplot(., aes(fill = Weight_Tier)) + 
      facet_wrap(~Label) + 
      geom_sf(color = "gray85", size = rel(0.1)) + 
      theme_classic() + 
      scale_fill_viridis_d(option = "D") +
      guides(fill = guide_legend(title = "Cluster Probability",nrow = 1)) + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.text = element_text(size = rel(1)),
            strip.background = element_rect(color = "transparent"),
            panel.background = element_rect(fill = "gray25"),
            legend.position = "bottom")}

#### Different Models Comparison ####

test_1 <- clusters %>% 
  filter(Model %in% c("model_48","model_16","model_07","model_04")) %>% 
  select(ID:Model,-Colors) %>% 
  spread(Model, Cluster) %>% 
  count(model_48,model_16,model_07,model_04) %>% 
  arrange(desc(n)) %>% 
  group_by(model_16) %>% 
  mutate(Aux_Order = max(n)) %>% 
  ungroup %>% 
  mutate(model_16 = reorder(model_16, Aux_Order))

test_1 %>% 
  ggplot(aes(axis1 = model_48, axis2 = model_07, axis3 = model_04, axis4 = model_16, 
             y = n)) + 
  geom_alluvium(aes(fill = model_16)) + 
  geom_stratum() + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) + 
  scale_x_continuous(breaks = 1:4,
                     labels = paste("model",c(48,07,04,16),sep = "_")) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

cluster_names <- clusters %>% 
  filter(Model %in% c("model_48","model_16","model_07","model_04")) %>% 
  select(-Colors) %>% 
  mutate(Cluster = str_pad(Cluster,2,"left",0),
         Cluster_Name = paste(Model,Cluster,sep="_"),
         Cluster_Name = case_when(Cluster_Name %in% c("model_16_02",
                                                      "model_04_05",
                                                      "model_07_04",
                                                      "model_48_01") ~ "A",
                                  Cluster_Name %in% c("model_16_08",
                                                      "model_04_01",
                                                      "model_07_05") ~ "B",
                                  Cluster_Name %in% c("model_16_03",
                                                      "model_04_02",
                                                      "model_07_03") ~ "C",
                                  Cluster_Name %in% c("model_16_04",
                                                      "model_04_03",
                                                      "model_07_02",
                                                      "model_48_00") ~ "D",
                                  Cluster_Name %in% c("model_16_06",
                                                      "model_04_00",
                                                      "model_07_01") ~ "E",
                                  Cluster_Name %in% c("model_16_07",
                                                      "model_04_04",
                                                      "model_07_00") ~ "F",
                                  Cluster_Name %in% c("model_16_00",
                                                      "model_04_06") ~ "G",
                                  Cluster_Name %in% c("model_16_05",
                                                      "model_04_07") ~ "H",
                                  Cluster_Name == "model_16_01" ~ "I",
                                  Cluster_Name == "model_16_10" ~ "J",
                                  Cluster_Name == "model_16_09" ~ "K"),
         Color = map_chr(Cluster_Name, ~ sankey_palette_random[which(LETTERS == .x)])) 

test_2 <- cluster_names %>% 
  select(ID:Model,Cluster_Name,-Cluster) %>% 
  spread(Model, Cluster_Name) %>% 
  count(model_48,model_16,model_07,model_04) %>% 
  arrange(desc(n)) %>% 
  group_by(model_16) %>% 
  mutate(Aux_Order = max(n)) %>% 
  ungroup %>% 
  mutate(model_16 = reorder(model_16, Aux_Order))

sankey_palette <- c("#340052","#4E007A","#7400b8", "#6930c3", "#5e60ce", 
                    "#5390d9", "#4ea8de", "#48bfe3", 
                    "#56cfe1", "#64dfdf", "#ADFFE8")
set.seed(272602)
sankey_palette_random <- sample(viridis::inferno(12,direction = -1))
sankey_palette_guide <- cluster_names %>% 
  distinct(Cluster_Name,Color) %>% 
  {set_names(pull(.,Color),pull(.,Cluster_Name))}
sankey_palette_guide[3] <- "#48BFE3"

sankey_clusters <- test_2 %>% 
  ggplot(aes(axis1 = model_16, axis2 = model_04, axis3 = model_07, axis4 = model_48, 
             y = n)) + 
  geom_alluvium(aes(fill = model_16)) + 
  geom_stratum() + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) + 
  scale_x_continuous(breaks = 1:4,
                     labels = c("UMAP + GMM",
                                "UMAP + K-Means",
                                "UMAP + HAC",
                                "PCA + GMM")) +
  scale_fill_manual(values = sankey_palette_guide) +
  labs(title = "Sankey Diagram comparing flows under different models including %SNP") + 
  theme_classic() + 
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.25), hjust = 0.5),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


comparison_map <- cluster_names %>% 
  left_join(west_hex_maps, by = c("ID"="gss_code")) %>% 
  mutate(Aux = row_number()) %>% 
  mutate(Model = ordered(Model, 
                         levels = c("model_16","model_04","model_07","model_48"),
                         labels = c("UMAP + GMM",
                                    "UMAP + K-Means",
                                    "UMAP + HAC",
                                    "PCA + GMM"))) %>% 
  st_as_sf() %>% 
  ggplot(aes(fill=Cluster_Name)) + 
  facet_wrap(~Model) + 
  geom_sf(color = "gray85", size = rel(0.1)) + 
  scale_fill_manual(values = sankey_palette_guide) + 
  theme_classic() + 
  theme(legend.position = "none")


plot_grid(sankey_clusters, comparison_map)

#### Final Model Heat Map ####

cluster_names %>% 
  filter(Model == "model_07") %>% 
  left_join(X_gb_with_snp) %>% 
  select(ID,Cluster_Name,
         Pop_Density,`%Unemployment`,`2019_Wage`,HousePrice,
         `%Heavy Industry & Manufacturing`,`%White`, `%NoQuals`,`%Level4+`,
         `con%`,`lab%`,`snp%`,`%LeaveVote`,`60+%`,`MbpsSpeed`,`%Sciences`) %>% 
  group_by(Cluster_Name) %>% 
  summarise_if(is.numeric,mean) %>% 
  mutate_if(is.numeric,.funs = ~scales::rescale(.,to=c(0,1))) %>% 
  gather(Variable,Score,-Cluster_Name) %>% 
  ggplot(aes(y=Variable,x=Cluster_Name, fill = Score)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "chocolate2", high = "steelblue4", mid = "gray95", 
                       midpoint = 0.5) + 
  theme_classic()
