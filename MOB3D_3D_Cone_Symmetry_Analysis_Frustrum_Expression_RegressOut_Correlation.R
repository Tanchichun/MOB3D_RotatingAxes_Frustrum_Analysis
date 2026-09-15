args = commandArgs(trailingOnly=TRUE)
library(SingleCellExperiment)
library(ggplot2)
library(tidyverse)
library(Seurat)

##### Collecting 3D Coordinate ####
# IRIS_Domain = readRDS(paste0("~/data/Collaboration/AlexF/MOB/IRIS_Leiden_classifyDV_Analysis/DomainLabel/IRIS_30Section_",numCluster,"Domain.rds"))
IRIS_Domain = read.csv(paste0("~/data/Collaboration/AlexF/MOB/IRIS_Leiden_classifyDV_Analysis/DomainLabel/Fleischmann_leiden_clusters.csv.gz"))
coord3d = read.csv("~/data/Collaboration/AlexF/MOB_NatureRev/symmetrized_spots_align15.csv")

##### Combine IRIS Results and Visium 3D Coordinate and D-V layers
coord3d$slice = unlist(lapply(coord3d$folder,function(x){strsplit(x,split = "_")[[1]][1]}))
coord3d$slice =  unlist(lapply(coord3d$slice,function(x){
  temp = strsplit(x,split = "-")[[1]]
  return(paste0("Slide",temp[2],"_",temp[3]))}))

### for Leiden Label File
IRIS_Domain.slide =  unlist(lapply(IRIS_Domain$zone,function(x){strsplit(x,split = "_ob")[[1]][1]}))
IRIS_Domain.section =  unlist(lapply(IRIS_Domain$zone,function(x){strsplit(x,split = "_ob")[[1]][2]}))
IRIS_Domain$Slice = paste0(IRIS_Domain.slide,"_",IRIS_Domain.section)


# IRIS_Domain$slice_barcode = paste0(IRIS_Domain$Slice,"_",IRIS_Domain$spotName)
IRIS_Domain.spotName = paste0(unlist(lapply(IRIS_Domain$barcode,function(x){strsplit(x,split = "-")[[1]][1]})),"-",
                              (unlist(lapply(IRIS_Domain$barcode,function(x){strsplit(x,split = "-")[[1]][2]}))) )
IRIS_Domain$spotName = IRIS_Domain.spotName
IRIS_Domain$slice_barcode = paste0(IRIS_Domain$Slice,"_",IRIS_Domain$spotName)

coord3d$slice_barcode = paste0(coord3d$slice,"_",coord3d$barcode,"-1")
IRIS_Domain$iris_domain = paste0("Leiden Domain",IRIS_Domain$leiden_1_2)

coord3d_IRIS = merge(coord3d,
                     IRIS_Domain,
                     by.x = c("slice_barcode"),
                     by.y = c("slice_barcode"))



###### Select Sister Glomeruli #####
select_glom = readxl::read_xlsx("~/data/Collaboration/AlexF/MOB_NatureRev/Glomerulus_list_two_column.xlsx")


###### Laoding Count Data #####
##### Filter Low Quality Genes
count_all = readRDS("~/data/Collaboration/AlexF/MOB_NatureRev/allsections_count.rds")
count_all = count_all[rowSums(count_all > 0) > 5,] ### drop genes with no greater than 5 spots expressed.


##### Normalized Gene '
# x = CreateSeuratObject(counts = count_all, assay = "RNA")
# ##### Normalization
# x = NormalizeData(x)
# normcount_all = x@assays$RNA$data
# saveRDS(normcount_all,"~/data/Collaboration/AlexF/MOB_NatureRev/allsections_logtransform.rds")
normcount_all = readRDS("~/data/Collaboration/AlexF/MOB_NatureRev/allsections_logtransform.rds")



##### Loading Frustum Structure ######
Rb = as.numeric(args)/2 ### radius of bottom circle (mitral layer spot)
# ###### Approach 0: All Capturing Spot
# covered_spots_all = readRDS("~/data/Collaboration/AlexF/MOB_NatureRev/Leiden_Cluster/DomainMap3D/Rt50_Rb200Extend_100_SelectedGlomerlus/frustum_covered_spots.rds")
###### Approach 1: Spots Filtered by Number of Markers Expressing
covered_spots_greater10 = readRDS(paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Leiden_Cluster/DomainMap3D/Rt50_Rb",Rb,"Extend_100_SpotReAlign_MitralMarker/frustum_covered_spots.rds"))
###### Approach 2: Spots Filtered by Number of Markers Expressing
covered_spots_geometry = readRDS(paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Leiden_Cluster/DomainMap3D/Rt50_Rb",Rb,"Extend_100_SpotReAlign_Geometry//frustum_covered_spots.rds"))
###### Approach 3: Cylinder Around Glomeruli
# covered_spots_sphere =  readRDS(paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Leiden_Cluster/DomainMap3D/Sphere_R",Rb,"_SelectedGlomerulus/sphere_covered_spots.rds"))


# ##### Investigate the overlap between two filtering spots method
# mitral_layer = coord3d_IRIS %>%
#   filter(iris_domain == "Leiden Domain6")
# 
# library(ggVennDiagram)
# glomerulus_list = names(covered_spots_all)
# 
# 
# venn_list = list()
# frustum_spots_all = NULL
# for(ig in 1:length(glomerulus_list)){
#   print(glomerulus_list[ig])
#   all_spot = covered_spots_all[[glomerulus_list[ig]]]$slice_barcode
#   marker_spot = covered_spots_greater10[[glomerulus_list[ig]]]$slice_barcode
#   geomtry_spot = covered_spots_geometry[[glomerulus_list[ig]]]$slice_barcode
#   
#   #### Record Spots membership for each frustum 
#   spot_df = data.frame(barcode = c(all_spot,marker_spot,geomtry_spot),
#                        type = c(rep("All Spots",length(all_spot)),
#                                 rep("Marker-filtered Spots",length(marker_spot)),
#                                 rep("Geometry-filtered Spots",length(geomtry_spot))
#                                 )
#                        )
#   spot_df$glomeruli = glomerulus_list[ig]
#   frustum_spots_all = rbind(frustum_spots_all,spot_df)
#   
#   
#   ##### Measure the overlap between filtered spot
#   all_spot_name = paste0("All Spots\n(",length(all_spot),")")
#   marker_spot_name = paste0("Marker-filtered Spots\n(",length(marker_spot),")")
#   geomtry_spot_name = paste0("Geometry-filtered Spots\n(",length(geomtry_spot),")")
#   
#   venn_sets <- list(
#     all_spot_name = all_spot,
#     marker_spot_name = marker_spot,
#     geomtry_spot_name = geomtry_spot
#   )
#   
#   names(venn_sets) = c(all_spot_name,marker_spot_name,geomtry_spot_name)
#  
#   p.venn <- ggVennDiagram(
#     venn_sets,
#     label = "count",            
#     label_alpha = 0,
#     label_size = 10,
#     set_size = 10,
#     # shape_id = "201",
#     category.names = names(venn_sets)
#   ) +
#     scale_x_discrete(
#       expand = expansion(mult = c(0.3, 0.3))   # 10% headroom
#     ) +
#     scale_y_discrete(
#       expand = expansion(mult = c(0.2, 0.2))   # 10% headroom
#     ) +
#     scale_fill_gradientn(colours = c("#fec89a","#ffd7ba","#fffcf2")) + # keep fills simple; customize below
#     theme(
#       legend.position = "none"
#     )
#   
#   
#   mitral_layer$marker_spot = mitral_layer$slice_barcode %in% marker_spot
#   mitral_layer$geomtry_spot = mitral_layer$slice_barcode %in% geomtry_spot
#   
#   mitral_layer = mitral_layer %>% arrange(marker_spot)
#   p.marker_spot = ggplot(mitral_layer,aes(x = x.x,y = y.x,z = z.x,color = marker_spot)) +
#     geom_point(aes(alpha = marker_spot),size =2) +
#     scale_color_manual(values = c("#adb5bd","#e63946")) +
#     labs(color = "Marker-filtered Spot") +
#     theme(panel.background = element_blank(),
#           panel.grid = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           axis.title = element_blank(),
#           legend.position = "bottom",
#           legend.title = element_text(size = 30),
#           legend.text = element_text(size = 30)) +
#     guides(alpha = "none")
#   
#   mitral_layer = mitral_layer %>% arrange(geomtry_spot)
#   p.geometry_spot = ggplot(mitral_layer ,aes(x = x.x,y = y.x,z = z.x,color = geomtry_spot)) +
#     geom_point(aes(alpha = geomtry_spot),size =2) +
#     scale_color_manual(values = c("#adb5bd","#e63946")) +
#     labs(color = "Geometry-filtered Spot") +
#     theme(panel.background = element_blank(),
#           panel.grid = element_blank(),
#           axis.ticks = element_blank(),
#           axis.text = element_blank(),
#           axis.title = element_blank(),
#           legend.position = "bottom",
#           legend.title = element_text(size = 30),
#           legend.text = element_text(size = 30)) +
#     guides(alpha = "none")
#   
#   
#   output_file = paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Check_Frustum_Overlap/",glomerulus_list[ig],".pdf")
#   
#   out_dir = dirname(output_file)
#   if (!dir.exists(out_dir)) {
#     dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#   }
#   
#   pdf(output_file,width = 30,height = 8)
#   print(cowplot::plot_grid( p.venn, p.marker_spot,p.geometry_spot,nrow = 1))
#   dev.off()
#   
#   
#   p.venn = p.venn + labs(title = glomerulus_list[ig]) + theme(plot.title = element_text(size = 30))
#   venn_list[[glomerulus_list[ig]]] = ggplotify::as.grob(p.venn)
# }
# 
# 
# write.table(frustum_spots_all,
#             "~/data/Collaboration/AlexF/MOB_NatureRev/frustum_spots_all.txt")
# 
# 
# 
# output_file = paste0("~/data/Collaboration/AlexF/MOB_NatureRev/frustum_spots_overlap_venn.png")
# 
# out_dir = dirname(output_file)
# if (!dir.exists(out_dir)) {
#   dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# }
# 
# png(output_file,width = 3000,height = 4000)
# print(cowplot::plot_grid( plotlist = venn_list,nrow = 8))
# dev.off()



#### mitral layer marker
# marker_list = read_csv("~/data/Collaboration/AlexF/Mitral_cell_layer_genes.csv",col_names = F)
all_gene = rownames(count_all)
marker_select = c("allgene")


# spots_filter = c("MitralMarker","Geometry","Sphere")
spots_filter = c("MitralMarker","Geometry")
glomeruli = names(covered_spots_greater10)




##### identify sister glomerulis for all glomerulis pairs
#### Glomeruli Location
glomerulus_location = read.csv("~/data/Collaboration/AlexF/MOB_NatureRev/geometrically_confirmed_correct.csv")

#### Glomeruli Domain
glomerulus_domain = read.csv("~/data/Collaboration/AlexF/MOB_NatureRev/OCAM_OMAC_glomeruli_ID.csv")

#### Identify Sister Glomerulis
select_glom = readxl::read_xls("~/data/Collaboration/AlexF/MOB_NatureRev/glom_for_correlation_480.xls")
select_glom = as.data.frame(select_glom)
colnames(select_glom) = c("Glomerulus_1","Glomerulus_2")

select_glom = select_glom[-104,] ###### remove Olfr237.ps1_259 and Olfr237.ps1_260
select_glom = select_glom[-171,] ###### remove Olfr715b_668 and Olfr715b_665

glomeruli_pair = expand.grid(c(select_glom$Glomerulus_1,select_glom$Glomerulus_2),
                             c(select_glom$Glomerulus_1,select_glom$Glomerulus_2))
colnames(glomeruli_pair) = c("Glomerulus_1","Glomerulus_2")
glomeruli_pair = glomeruli_pair[glomeruli_pair$Glomerulus_1 != glomeruli_pair$Glomerulus_2,]
glomeruli_pair$Glomerulus_1 = as.character(glomeruli_pair$Glomerulus_1)
glomeruli_pair$Glomerulus_2 = as.character(glomeruli_pair$Glomerulus_2)

## 1. Build canonical keys (sorted within each pair) for sister pairs
sister_keys <- select_glom %>%
  mutate(
    key1 = pmin(Glomerulus_1, Glomerulus_2),
    key2 = pmax(Glomerulus_1, Glomerulus_2)
  ) %>%
  transmute(key1, key2, is_sister = TRUE) %>%
  distinct()

## 2. Build canonical keys for all pairs and join
glomeruli_pair <- glomeruli_pair %>%
  mutate(
    key1 = pmin(Glomerulus_1, Glomerulus_2),
    key2 = pmax(Glomerulus_1, Glomerulus_2)
  ) %>%
  left_join(sister_keys, by = c("key1", "key2")) %>%
  mutate(
    is_sister = replace_na(is_sister, FALSE)
  ) %>%
  select(Glomerulus_1, Glomerulus_2, is_sister)


###### glomeruli gemotry info
glomerulus_location = glomerulus_location %>%
  mutate(slice_barcode = paste0(gene,"_",X))
rownames(glomerulus_location) = glomerulus_location$slice_barcode

glomerulus_location$domain.ML = sub(".*-", "", glomerulus_location$domain)
glomeruli_pair$Glomerulus_2_ML = glomerulus_location[glomeruli_pair$Glomerulus_2,]$domain.ML
glomeruli_pair$Glomerulus_1_ML = glomerulus_location[glomeruli_pair$Glomerulus_1,]$domain.ML

glomeruli_pair$Glomerulus_1.x = glomerulus_location[glomeruli_pair$Glomerulus_1,]$x
glomeruli_pair$Glomerulus_1.y = glomerulus_location[glomeruli_pair$Glomerulus_1,]$y
glomeruli_pair$Glomerulus_1.z = glomerulus_location[glomeruli_pair$Glomerulus_1,]$z

glomeruli_pair$Glomerulus_2.x = glomerulus_location[glomeruli_pair$Glomerulus_2,]$x
glomeruli_pair$Glomerulus_2.y = glomerulus_location[glomeruli_pair$Glomerulus_2,]$y
glomeruli_pair$Glomerulus_2.z = glomerulus_location[glomeruli_pair$Glomerulus_2,]$z


###### glomeruli domain info
glomerulus_domain = glomerulus_domain %>%
  mutate(slice_barcode = paste0(gene,"_",X))
glomerulus_domain = glomerulus_domain %>%
  filter(!duplicated(slice_barcode))
rownames(glomerulus_domain) = glomerulus_domain$slice_barcode



glomeruli_pair$Glomerulus_1.OMACs_domain = glomerulus_domain[glomeruli_pair$Glomerulus_1,]$OMACs_domain
glomeruli_pair$Glomerulus_1.OCAM_domain = glomerulus_domain[glomeruli_pair$Glomerulus_1,]$OCAM_domain
glomeruli_pair$Glomerulus_1.Domain = glomerulus_domain[glomeruli_pair$Glomerulus_1,]$Domain

glomeruli_pair$Glomerulus_2.OMACs_domain = glomerulus_domain[glomeruli_pair$Glomerulus_2,]$OMACs_domain
glomeruli_pair$Glomerulus_2.OCAM_domain =glomerulus_domain[glomeruli_pair$Glomerulus_2,]$OCAM_domain
glomeruli_pair$Glomerulus_2.Domain = glomerulus_domain[glomeruli_pair$Glomerulus_2,]$Domain


glomeruli_pair$OMAC_pair = paste0(glomeruli_pair$Glomerulus_1.OMACs_domain,"-",glomeruli_pair$Glomerulus_2.OMACs_domain)
glomeruli_pair$OCAM_pair = paste0(glomeruli_pair$Glomerulus_1.OCAM_domain,"-",glomeruli_pair$Glomerulus_2.OCAM_domain)
glomeruli_pair$Domain_pair = paste0(glomeruli_pair$Glomerulus_1.Domain,"-",glomeruli_pair$Glomerulus_2.Domain)


#### Identify Most Anterior and Most Dorsal Glomeruli
glomerulus_selected_locations = glomerulus_location %>%
  filter(slice_barcode %in% c(select_glom$Glomerulus_1,
                             select_glom$Glomerulus_2)) 


# smallest z per medial/lateral
min_z <- glomerulus_selected_locations %>%
  group_by(domain.ML) %>%
  slice_min(order_by = z, n = 1, with_ties = FALSE) %>%
  select(domain.ML, slice_barcode, z)

# smallest y per medial/lateral
min_y <- glomerulus_selected_locations %>%
  group_by(domain.ML) %>%
  slice_min(order_by = y, n = 1, with_ties = FALSE) %>%
  select(domain.ML, slice_barcode, y)


lateral_anterior_glomeruli = "Olfr726_372"
medial_anterior_glomeruli = "Olfr1281_524"
lateral_dorsal_glomeruli = "Olfr1309_521"
medial_dorsal_glomeruli = "Olfr1228_121"

##### Measuring Correlation between Frustrum's Expression
for(imarker in 1:length(marker_select)){ ##### For Each Approach of Adjusting Genes for "Regress-Out" analysis
  
  ####### Select Adjusting Marker list #####
  if(marker_select[imarker] == "allgene"){
    marker_gene = all_gene
  }
  

  

  for(ispotfilter in 1:length(spots_filter)){ ##### For Each Approach of Caputering Spots
    start_time <- Sys.time()
    cat(spots_filter[ispotfilter],"...\n")

    ####### Select Spots Filtering methods #####
    if(spots_filter[ispotfilter] == "Sphere"){
      covered_spots = covered_spots_sphere
    }
    if(spots_filter[ispotfilter] == "MitralMarker"){
      covered_spots = covered_spots_greater10 
    }
    if(spots_filter[ispotfilter] == "Geometry"){
      covered_spots = covered_spots_geometry 
    }
    
    ####### identified all cover spots ######
    all_cover_spots = do.call(rbind,covered_spots)
    all_cover_barcode = unique(all_cover_spots$slice_barcode)
    
    
    ####### Subset to Covered Normalized Expression ######
    normcount_covered = normcount_all[,all_cover_barcode]
    
    
    
    #### No Processed on the Gene Expression #####
    normcount_raw = normcount_covered
    
    
    #### mitral layer marker Mean out ##### 
    normcount_marker_mean = colMeans(normcount_covered[marker_gene,])
    normcount_meanout = normcount_covered
    normcount_meanout = sweep(normcount_meanout, 2, normcount_marker_mean, FUN = "-") 
    
    #### mitral layer marker Mean out - By Gene ##### 
    normcount_marker_mean = rowMeans(normcount_covered[marker_gene,])
    normcount_rowout = normcount_covered
    normcount_rowout = sweep(normcount_rowout, 1 , normcount_marker_mean, FUN = "-") 
    
    
    
    #### mitral layer marker Regress out #######
    normcount_marker_mean = colMeans(normcount_covered[marker_gene,])
    normcount_regressout = normcount_covered
    normcount_regressout = as.matrix(normcount_regressout)
    for(i in 1:nrow(normcount_regressout)){
      y = normcount_regressout[i,]
      x = normcount_marker_mean
      lm.fit = lm(y ~ x)
      residuals <- resid(lm.fit)
      normcount_regressout[i,] = residuals
    }
    
    #### mitral layer marker Regress out - By Gene #######
    normcount_marker_mean = rowMeans(normcount_covered[marker_gene,])
    normcount_regressrow = normcount_covered
    normcount_regressrow = as.matrix(normcount_regressout)
    for(j in 1:ncol(normcount_regressrow)){
      y = normcount_regressrow[,j]
      x = normcount_marker_mean
      lm.fit = lm(y ~ x)
      residuals <- resid(lm.fit)
      normcount_regressrow[,j] = residuals
    }
    
    
    
    
    #### mitral layer marker Regress PC out ######
    normcount_marker = normcount_covered[marker_gene,]
    marker_pc = irlba::irlba( normcount_marker)
    marker_pc = marker_pc$v
    
    normcount_regressoutPC = normcount_covered
    normcount_regressoutPC = as.matrix(normcount_regressoutPC)
    for(i in 1:nrow(normcount_regressout)){
      y = normcount_regressoutPC[i,]
      x = marker_pc
      lm.fit = lm(y ~ x)
      residuals <- resid(lm.fit)
      normcount_regressoutPC[i,] = residuals
    }
    
    end_time <- Sys.time()
    
    cat("Finished All Expression Adjustment...\n",
        as.numeric(difftime(end_time, start_time, units = "secs")),
        "seconds\n")
    
    ######### Calculate Correlation
    p_list = list()
    adjust_type = c("Raw","MeanOut","RowOut","RegressMeanOut","RegressRowOut","RegressPCOut")
    for(itype in 1:length(adjust_type)){ ##### For each expression adjustment method
      start_time <- Sys.time()
      # cat(adjust_type[itype],"...\n")
      ##### Identifying Adjustment for Expression #####
      if(adjust_type[itype] == "Raw"){
        expr =  normcount_raw
      }
      if(adjust_type[itype] == "MeanOut"){
        expr = normcount_meanout
      }
      if(adjust_type[itype] == "RowOut"){
        expr = normcount_rowout
      }
      if(adjust_type[itype] == "RegressMeanOut"){
        expr =  normcount_regressout
      }
      if(adjust_type[itype] == "RegressRowOut"){
        expr =  normcount_regressrow
      }
      if(adjust_type[itype] == "RegressPCOut"){
        expr = normcount_regressoutPC
      }
      
      
      
      
      ##### Average Expression Apporach #####
      frustrum_all = matrix(NA,nrow = length(covered_spots),
                            ncol = nrow(expr))
      
      rownames(frustrum_all) = names(covered_spots)
      colnames(frustrum_all) = rownames(expr)
      
      
      for(iglom in 1:nrow(frustrum_all)){
        # print(rownames(frustrum_all)[iglom])
        glom = rownames(frustrum_all)[iglom]
        spots = covered_spots[[ glom]]
        
        if(nrow(spots) < 10){
          next
        }
        frustrum_all[glom,] = rowMeans(expr[,spots$slice_barcode])
      }
      
    
      
      
      ##### measure the correlation between glomerulis' mean expression profiles by anterior and dorsalin medial and lateral ####
      pearson_corr= lapply(1:nrow(glomeruli_pair),function(i){
        cor(frustrum_all[glomeruli_pair$Glomerulus_1[i],],frustrum_all[glomeruli_pair$Glomerulus_2[i],])
      })
      pearson_corr = unlist(pearson_corr)
      glomeruli_pair$corr = pearson_corr
      
      saveRDS(glomeruli_pair,paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Frustum_Analysis_Result_Summary_SpotReAlign//Mean_Expression_Correlation/Rb_specific/Rb",Rb,"_",
                                    spots_filter[ispotfilter],"_",
                                    adjust_type[itype],"_Correlation_Results.rds"))
      
      end_time <- Sys.time()
      
      cat(adjust_type[itype],"...\n",
          as.numeric(difftime(end_time, start_time, units = "secs")),
          "seconds\n")

      
      # p.lateral.ap = glomeruli_pair %>%
      #   mutate(
      #          `A-P Position` = Glomerulus_2.z,
      #          `D-V Position` = Glomerulus_2.y) %>%
      #   filter(Glomerulus_2_ML == "lateral" &
      #            Glomerulus_1 == lateral_anterior_glomeruli) %>%
      #   ggplot(aes(x = `A-P Position`,y = corr,color = Glomerulus_2_ML)) +
      #   geom_point(size = 2) +
      #   labs(y = "Frustum Correlation",
      #        color = "M-L Domain",
      #        title =  paste0("Most anterior lateral glomerulus:\n",
      #                        lateral_anterior_glomeruli)) +
      #   scale_color_manual(values = c("green")) +
      #   theme(panel.background =element_blank(),
      #         panel.grid = element_blank(),
      #         plot.title = element_text(size = 20),
      #         axis.title.y = element_text(size = 15),
      #         axis.title.x = element_text(size = 15),
      #         axis.text.y = element_text(size = 15),
      #         axis.text.x = element_text(size = 15),
      #         legend.title = element_text(size = 15),
      #         legend.text = element_text(size = 15))
      # 
      # p.lateral.dv = glomeruli_pair %>%
      #   mutate(
      #          `A-P Position` = Glomerulus_2.z,
      #          `D-V Position` = Glomerulus_2.y) %>%
      #   filter(Glomerulus_2_ML == "lateral" &
      #            Glomerulus_1 == lateral_dorsal_glomeruli) %>%
      #   ggplot(aes(x = `D-V Position`,y = corr,color = Glomerulus_2_ML)) +
      #   geom_point(size = 2) +
      #   labs(y = "Frustum Correlation",
      #        color = "M-L Domain",
      #        title =  paste0("Most dorsal lateral glomerulus:\n",
      #                        lateral_dorsal_glomeruli)) +
      #   scale_color_manual(values = c("green")) +
      #   theme(panel.background =element_blank(),
      #         panel.grid = element_blank(),
      #         plot.title = element_text(size = 20),
      #         axis.title.y = element_text(size = 15),
      #         axis.title.x = element_text(size = 15),
      #         axis.text.y = element_text(size = 15),
      #         axis.text.x = element_text(size = 15),
      #         legend.title = element_text(size = 15),
      #         legend.text = element_text(size = 15))
      # 
      # p.medial.ap = glomeruli_pair %>%
      #   mutate(
      #          `A-P Position` = Glomerulus_2.z,
      #          `D-V Position` = Glomerulus_2.y) %>%
      #   filter(Glomerulus_2_ML == "medial"&
      #            Glomerulus_1 == medial_anterior_glomeruli) %>%
      #   ggplot(aes(x = `A-P Position`,y = corr,color = Glomerulus_2_ML)) +
      #   geom_point(size = 2) +
      #   labs(y = "Frustum Correlation",
      #        color = "M-L Domain",
      #        title =   paste0("Most anterior medial glomerulus:\n",
      #                         medial_anterior_glomeruli)) +
      #   scale_color_manual(values = c("blue")) +
      #   theme(panel.background =element_blank(),
      #         panel.grid = element_blank(),
      #         plot.title = element_text(size = 20),
      #         axis.title.y = element_text(size = 15),
      #         axis.title.x = element_text(size = 15),
      #         axis.text.y = element_text(size = 15),
      #         axis.text.x = element_text(size = 15),
      #         legend.title = element_text(size = 15),
      #         legend.text = element_text(size = 15))
      # 
      # p.medial.dv = glomeruli_pair %>%
      #   mutate(
      #          `A-P Position` = Glomerulus_2.z,
      #          `D-V Position` = Glomerulus_2.y) %>%
      #   filter(Glomerulus_2_ML == "medial"&
      #            Glomerulus_1 == medial_dorsal_glomeruli) %>%
      #   ggplot(aes(x = `D-V Position`,y = corr,color = Glomerulus_2_ML)) +
      #   geom_point(size = 2) +
      #   labs(y = "Frustum Correlation",
      #        color = "M-L Domain",
      #        title =   paste0("Most dorsal medial glomerulus:\n",
      #                         medial_dorsal_glomeruli)) +
      #   scale_color_manual(values = c("blue")) +
      #   theme(panel.background =element_blank(),
      #         panel.grid = element_blank(),
      #         plot.title = element_text(size = 20),
      #         axis.title.y = element_text(size = 15),
      #         axis.title.x = element_text(size = 15),
      #         axis.text.y = element_text(size = 15),
      #         axis.text.x = element_text(size = 15),
      #         legend.title = element_text(size = 15),
      #         legend.text = element_text(size = 15))
      # 
      # 
      #   plot_file = paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Frustum_Analysis_Result_Summary_SpotReAlign//Mean_Expression_Correlation/",
      #                      spots_filter[ispotfilter],"_",
      #                      adjust_type[itype],"_MostAnterior_MostDorsal_Comparisons.pdf")
      # 
      #   out_dir <- dirname(plot_file)
      #   if (!dir.exists(out_dir)) {
      #     dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      #   }
      # 
      #   pdf(plot_file,width = 10,height = 6)
      #   print(
      #     cowplot::plot_grid(p.lateral.ap,p.lateral.dv,
      #                        p.medial.ap,p.medial.dv,nrow = 2)
      #   )
      #   dev.off()
      #   
      #   
      #   
      #   
      #   ##### measure the correlation between glomerulis' by OMAC, OCAM, Domain ####
      #   p.OCAM = glomeruli_pair %>%
      #     # filter(OCAM_pair %in% c("No-No",
      #     #                         "No-Yes",
      #     #                         "Yes-Yes")) %>%
      #     ggplot(aes(x = OCAM_pair, y = corr,fill = OCAM_pair)) +
      #     geom_boxplot( size = 1
      #     ) +
      #     labs(y = "Pearson Correlation Coefficient of\n Mean Expression Profiles",
      #          x = "OCAM",
      #          title = paste0(adjust_type[itype]," Expression\n Covered Spot:",spots_filter[ispotfilter])) +
      #     theme(panel.background =element_blank(),
      #           panel.grid = element_blank(),
      #           plot.title = element_text(size = 20,face = "bold"),
      #           axis.title.y = element_text(size = 15),
      #           axis.title.x = element_blank(),
      #           axis.text.y = element_text(size = 15),
      #           axis.text.x = element_text(angle = 45,size = 12,hjust =  1),
      #           legend.title = element_text(size = 15),
      #           legend.text = element_text(size = 15)) +
      #     guides( color = guide_legend(
      #       override.aes   = list(size = 6)
      #     ))
      #   
      #   p.OMAC = glomeruli_pair %>%
      #     # filter(OMAC_pair %in% c("No-No",
      #     #                         "No-Yes",
      #     #                         "Yes-Yes")) %>%
      #     ggplot(aes(x = OMAC_pair, y = corr,fill = OMAC_pair)) +
      #     geom_boxplot( size = 1
      #     ) +
      #     labs(y = "Pearson Correlation Coefficient of\n Mean Expression Profiles",
      #          x = "OMAC",
      #          title = paste0(adjust_type[itype]," Expression\n Covered Spot:",spots_filter[ispotfilter])) +
      #     theme(panel.background =element_blank(),
      #           panel.grid = element_blank(),
      #           plot.title = element_text(size = 20,face = "bold"),
      #           axis.title.y = element_text(size = 15),
      #           axis.title.x = element_blank(),
      #           axis.text.y = element_text(size = 15),
      #           axis.text.x = element_text(angle = 45,size = 12,hjust =  1),
      #           legend.title = element_text(size = 15),
      #           legend.text = element_text(size = 15)) +
      #     guides( color = guide_legend(
      #       override.aes   = list(size = 6)
      #     ))
      #   
      #   
      #   p.Domain = glomeruli_pair %>%
      #     # filter(Domain_pair %in% c("classI-classI",
      #     #                         "classI-classII",
      #     #                         "classII-classII")) %>%
      #     ggplot(aes(x = Domain_pair, y = corr,fill = Domain_pair)) +
      #     geom_boxplot( size = 1
      #     ) +
      #     labs(y = "Pearson Correlation Coefficient of\n Mean Expression Profiles",
      #          x = "Domain",
      #          title = paste0(adjust_type[itype]," Expression\n Covered Spot:",spots_filter[ispotfilter])) +
      #     theme(panel.background =element_blank(),
      #           panel.grid = element_blank(),
      #           plot.title = element_text(size = 20,face = "bold"),
      #           axis.title.y = element_text(size = 15),
      #           axis.title.x = element_blank(),
      #           axis.text.y = element_text(size = 15),
      #           axis.text.x = element_text(angle = 45,size = 12,hjust =  1),
      #           legend.title = element_text(size = 15),
      #           legend.text = element_text(size = 15)) +
      #     guides( color = guide_legend(
      #       override.aes   = list(size = 6)
      #     ))
      #   
      #   
      #   plot_file = paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Frustum_Analysis_Result_Summary_SpotReAlign//Mean_Expression_Correlation/",
      #                      spots_filter[ispotfilter],"_",
      #                      adjust_type[itype],"_OCAM_OMAC_Domain_Comparisons.pdf")
      #   
      #   out_dir <- dirname(plot_file)
      #   if (!dir.exists(out_dir)) {
      #     dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      #   }
      #   
      #   pdf(plot_file,width = 16,height = 4)
      #   print(
      #     cowplot::plot_grid(p.OCAM,p.OMAC,
      #                        p.Domain,nrow = 1)
      #   )
      #   dev.off()
      # 
        # ##### Figure: For Each Glomeruli, Visualize correaltion by A-P,D-V,M-L #####
        # cat("Per Glomeruli Correlation Visualization...\n")
        # for(ig in 1:length(glomeruli)){
        #   focus_glomeruli = glomeruli[ig]
        #   sub_pair = glomeruli_pair %>%
        #     filter(Glomerulus_1 == focus_glomeruli)
        #   
        #   y_min =  min(sub_pair$corr,na.rm = T)
        #   y_max =  max(sub_pair$corr,na.rm = T)
        #   
        #   p.lateral.ap = sub_pair %>%
        #     mutate(point_color = case_when(is_sister == T ~ "Sister",
        #                                    is_sister == F ~ Glomerulus_2_ML),
        #            `A-P Position` = Glomerulus_2.z,
        #            `D-V Position` = Glomerulus_2.y) %>%
        #     filter(Glomerulus_2_ML == "lateral") %>%
        #     ggplot(aes(x = `A-P Position`,y = corr,color = point_color)) +
        #     geom_point(size = 2) +
        #     labs(y = "Frustum Correlation",
        #          color = "M-L Domain",
        #          title =  focus_glomeruli) +
        #     scale_y_continuous(limits = c( y_min, y_max)) +
        #     scale_color_manual(values = c("green","red")) +
        #     theme(panel.background =element_blank(),
        #           panel.grid = element_blank(),
        #           plot.title = element_text(size = 20),
        #           axis.title.y = element_text(size = 15),
        #           axis.title.x = element_text(size = 15),
        #           axis.text.y = element_text(size = 15),
        #           axis.text.x = element_text(size = 15),
        #           legend.title = element_text(size = 15),
        #           legend.text = element_text(size = 15))
        #   
        #   p.lateral.dv = sub_pair %>%
        #     mutate(point_color = case_when(is_sister == T ~ "Sister",
        #                                    is_sister == F ~ Glomerulus_2_ML),
        #            `A-P Position` = Glomerulus_2.z,
        #            `D-V Position` = Glomerulus_2.y) %>%
        #     filter(Glomerulus_2_ML == "lateral") %>%
        #     ggplot(aes(x = `D-V Position`,y = corr,color = point_color)) +
        #     geom_point(size = 2) +
        #     labs(y = "Frustum Correlation",
        #          color = "M-L Domain",
        #          title =  focus_glomeruli) +
        #     scale_y_continuous(limits = c( y_min, y_max)) +
        #     scale_color_manual(values = c("green","red")) +
        #     theme(panel.background =element_blank(),
        #           panel.grid = element_blank(),
        #           plot.title = element_text(size = 20),
        #           axis.title.y = element_text(size = 15),
        #           axis.title.x = element_text(size = 15),
        #           axis.text.y = element_text(size = 15),
        #           axis.text.x = element_text(size = 15),
        #           legend.title = element_text(size = 15),
        #           legend.text = element_text(size = 15))
        #   
        #   p.medial.ap = sub_pair %>%
        #     mutate(point_color = case_when(is_sister == T ~ "Sister",
        #                                    is_sister == F ~ Glomerulus_2_ML),
        #            `A-P Position` = Glomerulus_2.z,
        #            `D-V Position` = Glomerulus_2.y) %>%
        #     filter(Glomerulus_2_ML == "medial") %>%
        #     ggplot(aes(x = `A-P Position`,y = corr,color = point_color)) +
        #     geom_point(size = 2) +
        #     labs(y = "Frustum Correlation",
        #          color = "M-L Domain",
        #          title =  focus_glomeruli) +
        #     scale_y_continuous(limits = c( y_min, y_max)) +
        #     scale_color_manual(values = c("blue","red")) +
        #     theme(panel.background =element_blank(),
        #           panel.grid = element_blank(),
        #           plot.title = element_text(size = 20),
        #           axis.title.y = element_text(size = 15),
        #           axis.title.x = element_text(size = 15),
        #           axis.text.y = element_text(size = 15),
        #           axis.text.x = element_text(size = 15),
        #           legend.title = element_text(size = 15),
        #           legend.text = element_text(size = 15))
        #   
        #   p.medial.dv = sub_pair %>%
        #     mutate(point_color = case_when(is_sister == T ~ "Sister",
        #                                    is_sister == F ~ Glomerulus_2_ML),
        #            `A-P Position` = Glomerulus_2.z,
        #            `D-V Position` = Glomerulus_2.y) %>%
        #     filter(Glomerulus_2_ML == "medial") %>%
        #     ggplot(aes(x = `D-V Position`,y = corr,color = point_color)) +
        #     geom_point(size = 2) +
        #     labs(y = "Frustum Correlation",
        #          color = "M-L Domain",
        #          title =  focus_glomeruli) +
        #     scale_y_continuous(limits = c( y_min, y_max)) +
        #     scale_color_manual(values = c("blue","red")) +
        #     theme(panel.background =element_blank(),
        #           panel.grid = element_blank(),
        #           plot.title = element_text(size = 20),
        #           axis.title.y = element_text(size = 15),
        #           axis.title.x = element_text(size = 15),
        #           axis.text.y = element_text(size = 15),
        #           axis.text.x = element_text(size = 15),
        #           legend.title = element_text(size = 15),
        #           legend.text = element_text(size = 15))
        #   
        #   
        #   plot_file = paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Frustum_Analysis_Result_Summary/Mean_Expression_Correlation/",
        #                      spots_filter[ispotfilter],"/",
        #                      adjust_type[itype],"/",
        #                      focus_glomeruli,".pdf")
        #   
        #   out_dir <- dirname(plot_file)
        #   if (!dir.exists(out_dir)) {
        #     dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        #   }
        #   
        #   pdf(plot_file,width = 10,height = 6)
        #   print(
        #     cowplot::plot_grid(p.lateral.ap,p.lateral.dv,
        #                        p.medial.ap,p.medial.dv,nrow = 2)
        #   )
        #   dev.off()
        #   
        # 
        # }
        
        
        
        
        
        # ##### Figure: Correlation Displayed by Medial Lateral #####
        # p.noML = glomeruli_pair %>%
        #   arrange(is_sister) %>%
        #   ggplot(aes(x = Glomerulus_1, y = corr,color = is_sister)) +
        #   geom_point( size = 2
        #   ) +
        #   labs(y = "Pearson Correlation Coefficient of\n Mean Expression Profiles",
        #        x = "Selected Glomeruli",
        #        color = "Sister Glomeruli",
        #        title = paste0(adjust_type[itype]," Expression\n Covered Spot:",spots_filter[ispotfilter])) +
        #   scale_color_manual(values = c("#f2e8cf","#c1121f")) +
        #   theme(panel.background =element_blank(),
        #         panel.grid = element_blank(),
        #         plot.title = element_text(size = 20,face = "bold"),
        #         axis.title.y = element_text(size = 15),
        #         axis.title.x = element_blank(),
        #         axis.text.y = element_text(size = 15),
        #         axis.text.x = element_text(angle = 45,size = 12,hjust =  1),
        #         legend.title = element_text(size = 15),
        #         legend.text = element_text(size = 15)) +
        #   guides( color = guide_legend(
        #     override.aes   = list(size = 6)
        #   ))
        # 
        # p.ML = glomeruli_pair %>%
        #   mutate(point_color = case_when(is_sister == T ~ "Sister",
        #                                  is_sister == F ~ Glomerulus_2_ML)) %>%
        #   arrange(point_color) %>%
        #   ggplot(aes(x = Glomerulus_1, y = corr,color = point_color,group = Glomerulus_2_ML,alpha = is_sister)) +
        #   geom_point(size = 2,
        #     position = position_dodge(width = 1)
        #   ) +
        #   labs(y = "Pearson Correlation Coefficient of\n Mean Expression Profiles",
        #        x = "Selected Glomeruli",
        #        color = "M-L Domain",
        #        title = paste0(adjust_type[itype]," Expression\n Covered Spot:",spots_filter[ispotfilter])) +
        #   scale_color_manual(values = c("green","blue","red")) +
        #   theme(panel.background =element_blank(),
        #         panel.grid = element_blank(),
        #         plot.title = element_text(size = 20,face = "bold"),
        #         axis.title.y = element_text(size = 15),
        #         axis.title.x = element_blank(),
        #         axis.text.y = element_text(size = 15),
        #         axis.text.x = element_text(angle = 45,size = 12,hjust =  1),
        #         legend.title = element_text(size = 15),
        #         legend.text = element_text(size = 15)) +
        #   guides( color = guide_legend(
        #     override.aes   = list(size = 6)
        #   ),
        #   alpha = "none")
        # 
        # p.sister = ggplot(glomeruli_pair,aes(x = is_sister, y = corr,fill = is_sister)) +
        #   geom_boxplot(
        #   ) +
        #   geom_jitter(alpha = 0.2
        #   ) +
        #   labs(y = "Pearson Correlation Coefficient of\n Mean Expression Profiles",
        #        x = "Sister Glomeruli",
        #        fill = "Sister Glomeruli",
        #        title = paste0(adjust_type[itype]," Expression\n Covered Spot:",spots_filter[ispotfilter])) +
        #   scale_fill_manual(values = c("#f2e8cf","#c1121f")) +
        #   theme(panel.background =element_blank(),
        #         panel.grid = element_blank(),
        #         plot.title = element_text(size = 20,face = "bold"),
        #         axis.title.y = element_text(size = 15),
        #         axis.title.x = element_blank(),
        #         axis.title = element_text(size = 20),
        #         axis.text.y = element_text(size = 20),
        #         axis.text.x = element_blank(),
        #         legend.position = "right",
        #         legend.title = element_text(size = 15),
        #         legend.text = element_text(size = 15)) +
        #   guides( color = guide_legend(
        #     override.aes   = list(size = 6)
        #   ))
        # 
        # 
        # p_combined = cowplot::plot_grid(p.noML,p.ML,p.sister,nrow = 1)
        # p_combined = ggplotify::as.grob(p_combined)
        # p_list[[ adjust_type[itype] ]] = ggplotify::as.grob(    p_combined )
      
    }  ##### For each expression adjustment method
    
    # plot_file = paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Frustum_Analysis_Result_Summary_SpotReAlign//Mean_Expression_Correlation/Covering_Spots_",
    #                    spots_filter[ispotfilter],".pdf")
    # 
    # pdf(plot_file,width = 35,height = 25)
    # print(
    #   cowplot::plot_grid(plotlist = p_list,nrow = length(adjust_type))
    # )
    # dev.off()


  }   ##### For Each Approach of Caputering Spots

  
}  ##### For Each Approach of Adjusting Genes for "Regress-Out" analysis
