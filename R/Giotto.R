# Rscript Giotto.R [matrix rds]

args<-commandArgs(T)

#setwd("/data")
library(Giotto)
library(Seurat)


#Consider to install these (optional) packages to run all possible Giotto commands for spatial analyses:  
# scran MAST smfishHmrf trendsceek SPARK multinet RTriangle FactoMiner

###### create object
# from csv files
#expression_matrix = readExprMatrix("mouseE9.5_hvg3000_genexcell.csv")
#cell_locations = data.table::fread("mouseE9.5_bin50_sample_info.csv")
# change rownames of cell location
#cell_locations$V1<-paste0("X", cell_locations$V1)

SeuObj<-readRDS(args[1])
expression_matrix<-SeuObj@assays$RNA@counts
cell_locations<-SeuObj@meta.data[,c("coord_x","coord_y")]

cell_locations$cell_id<-rownames(cell_locations)
rownames(cell_locations)<-NULL

#  coord_x coord_y cell_id
#1   19880    9378  170675
#2   19876    9349  170785

#expression_matrix[1:3, 1:3]
#head(cell_locations)

GioObj = createGiottoObject(raw_exprs = expression_matrix, spatial_locs = cell_locations)
GioObj

# filter
GioObj <- filterGiotto(gobject = GioObj, 
                                 expression_threshold = 1, 
                                 gene_det_in_min_cells = 10, 
                                 min_det_genes_per_cell = 5)

# normalise
GioObj <- normalizeGiotto(gobject = GioObj, scalefactor = 6000, verbose = T)

# add stats
GioObj <- addStatistics(gobject = GioObj)
# account for technical confounders
GioObj <- adjustGiottoMatrix(gobject = GioObj, 
                                   expression_values = c('normalized'),
                                   covariate_columns = c('nr_genes', 'total_expr'))

head(GioObj@cell_metadata)

###### dim reduction (HVG)
GioObj <- calculateHVG(gobject = GioObj, save_plot=TRUE)
# run PCA
GioObj <- runPCA(gobject = GioObj)
# umap
GioObj <- runUMAP(GioObj, dimensions_to_use = 1:10) #n_components=3, n_threads=4
# tsne
#GioObj <- runtSNE(GioObj, dimensions_to_use = 1:5)
#plotTSNE(gobject = GioObj)

# identify most informative principal components
screePlot(GioObj, ncp = 20, save_plot=T)
#jackstrawPlot(GioObj, ncp = 20)
plotPCA(GioObj, save_plot=T)
plotUMAP(gobject = GioObj, save_plot=T)


###### clustering
GioObj <- createNearestNetwork(gobject = GioObj, dimensions_to_use = 1:10, k = 8)
# leiden
GioObj = doLeidenCluster(GioObj, resolution=0.4, name = 'leiden_clus') # resolution=0.25, n_iterations=200,name = 'leiden_0.25.200'
# louvain
GioObj = doLouvainCluster(GioObj, name = 'louvain_clus') # works
# kmeans
#GioObj = doKmeans(GioObj, centers = 4, name = 'kmeans_clus')
# hierarchical
#GioObj = doHclust(GioObj, k = 4, name = 'hier_clus')


# calculate cluster similarities to see how individual clusters are correlated
cluster_similarities = getClusterSimilarity(GioObj, cluster_column = 'leiden_clus')
# merge similar clusters based on correlation and size parameters
GioObj = mergeClusters(GioObj, cluster_column = 'leiden_clus', new_cluster_name = 'leiden_clus_m',
                                        min_cor_score = 0.7, force_min_group_size = 4) #max_group_size=100, max_sim_clusters=10

# visualize
cell_metadata=pDataDT(GioObj)
head(cell_metadata)

#      cell_ID nr_genes perc_genes total_expr louvain_clus
#   1:  170675       94 0.75794227   531.7494            0
#   2:  170785       10 0.08063216    91.2457            2

write.table(cell_metadata, file="giotto_cellmeta.txt", quote=F)

# plots

# clustering 
#plotUMAP_2D(GioObj, cell_color = 'leiden_clus', point_size = 3)
plotUMAP_2D(GioObj, cell_color = 'leiden_clus', show_NN_network = T, point_size = 3, save_param = list(save_name = '4_a_UMAP_leiden'), save_plot=T)
plotUMAP_2D(GioObj, cell_color = 'louvain_clus', point_size = 3, save_plot=T)
#plotUMAP_2D(GioObj, cell_color = 'kmeans_clus', point_size = 3)
#plotUMAP_2D(GioObj, cell_color = 'hier_clus', point_size = 3)

# merged and spatial 
plotUMAP_2D(GioObj, cell_color='leiden_clus_m', point_size=2.5, show_NN_network=T, edge_alpha=0.05, save_param = list(save_name = '4_a_UMAP_leiden_m'), save_plot=T)
spatDimPlot(gobject = GioObj, cell_color = 'leiden_clus_m', spat_point_shape = 'voronoi', save_param = list(save_name = '4_a_DimPlot_leiden_m_voronoi'), save_plot=T)
spatDimPlot(gobject = GioObj, cell_color = 'leiden_clus_m', spat_point_size=1.5, save_param = list(save_name = '4_a_DimPlot_leiden_m'), save_plot=T)


# dendrogram splits
# install.packages('dendextend')
#splits = getDendrogramSplits(GioObj, cluster_column = 'merged_cluster')

###### marker genes
# three methods: gini, scran, and mast
gini_markers = findMarkers_one_vs_all(gobject = GioObj, method = 'gini', expression_values = 'normalized', cluster_column = 'leiden_clus')

write.table(gini_markers, file="gini_markers_leiden.txt", quote=F)
# get top 6 genes per cluster and visualize with heatmap
topgenes_gini2 = gini_markers[, head(.SD, 6), by = 'cluster']
plotMetaDataHeatmap(GioObj, selected_genes = topgenes_gini2$genes,
                    metadata_cols = c('leiden_clus'), save_plot=T)

####### cell types in tissue
# need manual annotation


####### spatial grid
GioObj <- createSpatialGrid(gobject = GioObj, sdimx_stepsize = 500, sdimy_stepsize = 500, minimum_padding = 50)

spatPlot(gobject = GioObj, show_grid = T, cell_color='leiden_clus_m', grid_color='red', point_size = 1.5, save_param = c(save_name = '8_a_grid'), save_plot=T)

###### cell networks and genes
# install.packages('geometry')
plotStatDelaunayNetwork(gobject= GioObj, method = 'delaunayn_geometry',
                        maximum_distance = 100, show_plot = F, save_plot = T)

# Create Spatial Network using Delaunay geometry  or kNN
GioObj = createSpatialNetwork(gobject = GioObj, delaunay_method = 'delaunayn_geometry', 
                               minimum_k = 2, maximum_distance_delaunay = 100)

# visualize spatial networks
highexp_ids = cell_metadata[total_expr>=100]$cell_ID
subGioObj = subsetGiotto(GioObj, cell_ids = highexp_ids)

spatPlot(gobject = subGioObj, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network',
         point_size = 1.5, cell_color = 'leiden_clus_m',
         save_param = list(save_name = '9_Delaunay_network_spatPlot_highExp'), save_plot=T)

# spatially coherent expressed genes using binSpect
spat_genes = binSpect(GioObj, expression_values = "scaled", bin_method="kmeans", spatial_network_name = "Delaunay_network") #kNN_network

spatGenePlot(GioObj, expression_values = "scaled", genes = spat_genes[1:4]$genes, point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5, cow_n_col = 2,
             save_param = list(save_name = '10_a_spatialgenes_km'), save_plot=T)


# spatial gene co-expr modules
ext_spatial_genes = spat_genes[1:1000]$genes
spat_cor_netw_DT = detectSpatialCorGenes(GioObj, method = 'network',
			spatial_network_name = 'Delaunay_network', subset_genes = ext_spatial_genes)

spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = 8)

heatmSpatialCorGenes(GioObj, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', save_param=c(save_name = '10_b_spatialCoExpression_heatmap', base_height = 6, base_width = 8, units = 'cm'), save_plot=T)

# rank spatial correlated clusters
netw_ranks = rankSpatialCorGroups(GioObj, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                  save_param = c(save_name = '10_c_spatialcoexpression_rank',
                                                 base_height = 3, base_width = 5))

## metagenes
cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT,
	use_clus_name = 'spat_netw_clus', show_top_genes = 1)

# create metagenes from cluster modules and visualize
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$gene_ID
GioObj = createMetagenes(GioObj, gene_clusters = cluster_genes, name = 'cluster_metagene')

spatCellPlot(GioObj, spat_enr_names = 'cluster_metagene',
	cell_annotation_values = netw_ranks$clusters, point_size = 1.5, cow_n_col = 3, save_param = list(save_name = '10_c_spatialMetagenes'), save_plot=T)

######## spatial domains HMRF
#hmrf_folder = paste0('./','11_HMRF/')
#if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

# do HMRF with different betas
#HMRF_spat_genes = doHMRF(GioObj, expression_values="scaled", spatial_genes=names(cluster_genes),
#	spatial_network_name = "Delaunay_network", zscore = "none", k = 8,
#	betas = c(28,2,3),
#	output_folder = paste0(hmrf_folder, "/", 'Spatial_genes/'))

#Error: package or namespace load failed for ‘dplyr’ in dyn.load(file, DLLpath = DLLpath, ...):
# unable to load shared object '/usr/local/lib/R/site-library/rlang/libs/rlang.so':
#  /usr/local/lib/R/site-library/rlang/libs/rlang.so: undefined symbol: R_ActiveBindingFunction 

## add HMRF of interest to giotto object
#GioObj = addHMRF(gobject = GioObj, HMRFoutput = HMRF_spatial_genes,
#                  k = 9, betas_to_add = c(28), hmrf_name = 'HMRF_2')

## visualize
#spatPlot(gobject = GioObj, cell_color = 'HMRF_2_k8_b.28', point_size = 3, coord_fix_ratio = 1, 
#         save_param = c(save_name = '11_HMRF_2_k8_b.28', base_height = 3, base_width = 9, save_format = 'pdf'))


######## cell neighborhood: 
cell_proximities = cellProximityEnrichment(gobject = GioObj, cluster_column = 'leiden_clus_m',
                                           spatial_network_name = 'Delaunay_network', adjust_method = 'fdr',
                                           number_of_simulations = 2000)

## barplot
cellProximityBarplot(gobject = GioObj,
                     CPscore = cell_proximities, min_orig_ints = 5, min_sim_ints = 5, 
                     save_param = c(save_name = '12_a_barplot_cell_cell_enrichment'), save_plot=T)
## heatmap
cellProximityHeatmap(gobject = GioObj, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),
                     save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in'), save_plot=T)
## network
cellProximityNetwork(gobject = GioObj, CPscore = cell_proximities, remove_self_edges = T,
                     only_show_enrichment_edges = T,
                     save_param = c(save_name = '12_c_network_cell_cell_enrichment'), save_plot=T)

## network with self-edges
cellProximityNetwork(gobject = GioObj, CPscore = cell_proximities, remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F, rescale_edge_weights = T, node_size = 8,
                     edge_weight_range_depletion = c(1, 2), edge_weight_range_enrichment = c(2,5),
                     save_param = c(save_name = '12_d_network_cell_cell_enrichment_self',
                                    base_height = 5, base_width = 5, save_format = 'pdf'), save_plot=T)
