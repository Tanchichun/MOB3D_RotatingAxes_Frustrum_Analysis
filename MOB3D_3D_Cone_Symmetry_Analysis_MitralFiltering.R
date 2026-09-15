args = commandArgs(trailingOnly=TRUE)
library(SingleCellExperiment)
library(ggplot2)
library(tidyverse)
library(Seurat)

#### Function: change the first letter to capital
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
#### Load All MOB Slices in Matrix 
countMat = readRDS("~/data/Collaboration/AlexF/MOB/AllSection_NormCount.rds")


######## IRIS domain and 3D Coordiante
##### IRIS Results and 3D Coordinate
spatial_countMat_list = readRDS("~/data/Collaboration/AlexF/MOB/spatial_countMat_list.RDS")
spatial_location_list = readRDS("~/data/Collaboration/AlexF/MOB/spatial_location_list.RDS")


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

coord3d  = coord3d %>%
  filter(slice %in% names(spatial_countMat_list))

section_slide = coord3d %>%
  dplyr::select(slice,section) %>%
  unique() %>%
  dplyr::filter(slice %in% names(spatial_countMat_list)) %>%
  dplyr::arrange(section)

section_s = section_slide$section
visium_s = section_slide$slice

coord3d = coord3d %>%
  filter(section %in% section_s)
IRIS_Domain = IRIS_Domain %>%
  filter(Slice %in% visium_s)

# IRIS_Domain$slice_barcode = paste0(IRIS_Domain$Slice,"_",IRIS_Domain$spotName)
IRIS_Domain.spotName = paste0(unlist(lapply(IRIS_Domain$barcode,function(x){strsplit(x,split = "-")[[1]][1]})),"-",
                              (unlist(lapply(IRIS_Domain$barcode,function(x){strsplit(x,split = "-")[[1]][2]}))) )
IRIS_Domain$spotName = IRIS_Domain.spotName
IRIS_Domain$slice_barcode = paste0(IRIS_Domain$Slice,"_",IRIS_Domain$spotName)

coord3d$slice_barcode = paste0(coord3d$slice,"_",coord3d$barcode,"-1")
# IRIS_Domain$iris_domain = paste0("IRIS Domain",IRIS_Domain$IRIS_domain)
IRIS_Domain$iris_domain = paste0("Leiden Domain",IRIS_Domain$leiden_1_2)

coord3d_IRIS = merge(coord3d,
                     IRIS_Domain,
                     by.x = c("slice_barcode"),
                     by.y = c("slice_barcode"))


# #### Slide29_2 Main Figure Focus
# count = spatial_countMat_list[["Slide29_2"]]
# location = spatial_location_list[["Slide29_2"]]
# saveRDS(count,"~/data/Collaboration/AlexF/MOB_Slide29_2_Count.rds")
# saveRDS(location,"~/data/Collaboration/AlexF/MOB_Slide29_2_Location.rds")

Rt = 100/2 ### radius of top circle (glomeruli)
Rb = as.numeric(args)/2 ### radius of bottom circle (mitral layer spot)
extend_height <- 100

radius_dir = paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Leiden_Cluster/DomainMap3D/Rt",Rt,"_Rb",Rb,"Extend_",extend_height,"_SpotReAlign/")
if (!dir.exists(radius_dir)) dir.create(radius_dir, recursive = TRUE)


#### Glomerulus Location
glomerulus_location = read.csv("~/data/Collaboration/AlexF/MOB_NatureRev/geometrically_confirmed_correct.csv")
glomerulus_forMerge =  glomerulus_location %>%
  mutate(slice_barcode = paste0(gene,"_",X),
         iris_domain = "Glomeruli",
         Slice = NA) %>%
  select(slice_barcode,x,y,z,iris_domain,Slice)

#### 3D Visualization
library(rgl)
library(gtools)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(plotly)

coord3d_IRIS_mitral = coord3d_IRIS %>%
  filter(iris_domain == "Leiden Domain6") %>%
  select(slice_barcode,x.x,y.x,z.x,iris_domain,Slice) %>%
  rename(
    x = x.x,
    y = y.x,
    z = z.x
  )

coord3d_IRIS_mitral_and_glomer = rbind(coord3d_IRIS_mitral,
                                       glomerulus_forMerge)




##### Frustum Definition 
Rt = 100/2 ### radius of top circle (glomeruli)
Rb = as.numeric(args)/2 ### radius of bottom circle (mitral layer spot)
extend_height <- 100

#### Glomerulus Location
glomerulus_location = read.csv("~/data/Collaboration/AlexF/MOB_NatureRev/geometrically_confirmed_correct.csv")
glomerulus_forMerge =  glomerulus_location %>%
  mutate(slice_barcode = paste0(gene,"_",X),
         iris_domain = "Glomeruli",
         Slice = NA) %>%
  select(slice_barcode,x,y,z,iris_domain,Slice)


coord3d_IRIS_mitral = coord3d_IRIS %>%
  filter(iris_domain == "Leiden Domain6") %>%
  select(slice_barcode,x.x,y.x,z.x,iris_domain,Slice) %>%
  rename(
    x = x.x,
    y = y.x,
    z = z.x
  )

coord3d_IRIS_mitral_and_glomer = rbind(coord3d_IRIS_mitral,
                                       glomerulus_forMerge)



# ###### Select Sister Glomeruli
# select_glom = readxl::read_xlsx("~/data/Collaboration/AlexF/MOB_NatureRev/Glomerulus_list_two_column.xlsx")



##### Symmetric Anlaysis
covered_df = readRDS(paste0(radius_dir,"/frustum_covered_spots.rds")) 
count_all = readRDS("~/data/Collaboration/AlexF/MOB_NatureRev/allsections_count.rds")

##### Filter Low Quality Genes
count_all = count_all[rowSums(count_all > 0) > 5,] ### drop genes with no greater than 5 spots expressed.


##### All Genes
x = CreateSeuratObject(counts = count_all, assay = "RNA")
##### Normalization
x = NormalizeData(x)
normcount_all = x@assays$RNA$data


#### mitral layer marker
marker_list = read_csv("~/data/Collaboration/AlexF/Mitral_cell_layer_genes.csv",col_names = F)
marker_gene = marker_list$X1

normcount_markrer = normcount_all[marker_gene,]





###### Spots Filtering 1: filtering spots with a certain number of marker genes expressing ######
radius_dir = paste0("~/data/Collaboration/AlexF/MOB_NatureRev/Leiden_Cluster/DomainMap3D/Rt",Rt,"_Rb",Rb,"Extend_",extend_height,"_SpotReAlign_MitralMarker/")
if (!dir.exists(radius_dir)) dir.create(radius_dir, recursive = TRUE)


##### Check the number of spots expressing marker genes
num_df = NULL
for(iglom in 1:length(covered_df)){
  print(names(covered_df)[iglom])
  glom = names(covered_df)[iglom]
  covered_spots = covered_df[[ glom]]
  
  normcount_covered_spots = normcount_markrer[,covered_spots$slice_barcode]
  marker_exp_ornot =   1 * (normcount_covered_spots == 0)
  spot_exp_sum = colSums(marker_exp_ornot)
  
  temp = data.frame(spot = names(spot_exp_sum),
                    exp_num =   spot_exp_sum)
  temp$glomeruli  =   glom 
  
  num_df = rbind(num_df,temp)
}



p_spot_num_markers = num_df %>%
  ggplot(aes(x = exp_num)) +
  geom_histogram() + 
  labs(x = "The Number of Marker Genes Expressing in Spots",
         y = "The Number of Spots") +
  facet_wrap(~ glomeruli,nrow = 5) +
  
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 20))

pdf(paste0(radius_dir,"/selected_glomerulis_Num_of_Markers_Express_in_Spots.pdf"),width =16,height = 8)
p_spot_num_markers
dev.off()



num_df = num_df %>%
  mutate(num_markers_factor = case_when(exp_num > 40 & exp_num <= 50 ~ "41~50 Markers",
                                        exp_num > 30 & exp_num <= 40 ~ "31~40 Markers",
                                        exp_num > 20 & exp_num <= 30 ~ "21~30 Markers",
                                        exp_num > 10 & exp_num <= 20 ~ "11~20 Markers",
                                        exp_num <= 10 ~ "No More Than 10 Markers"))

prop_df <- num_df %>%
  group_by(glomeruli, num_markers_factor) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(glomeruli) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

num_df$num_markers_factor = factor(num_df$num_markers_factor,
                                   levels = c("No More Than 10 Markers","11~20 Markers",
                                              "21~30 Markers","31~40 Markers","41~50 Markers")
                                   )


p_glom_prop_spot = ggplot(num_df, aes(x = glomeruli, fill = num_markers_factor)) +
  geom_bar(position = "fill") +  # stacks to 1 (i.e., proportions)
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Glomerulus",
    y = "Proportion of spots",
    fill = "Number of Markers Expressing"
  ) +
  ggsci::scale_fill_jco() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 20,angle = 45,hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title =  element_text(size = 20),
        legend.text =  element_text(size = 20),
        legend.position = "bottom"
        )

pdf(paste0(radius_dir,"/selected_glomerulis_Marker_Expressed_Level_Per_Spot.pdf"),width =20,height = 10)
p_glom_prop_spot 
dev.off()


##### Covered Spot After Filtering 
covered_spots_select = num_df %>%
  filter(!num_markers_factor %in% c("No More Than 10 Markers"
  ))


covered_df_filtered = covered_df
for(iglom in 1:length(covered_df)){
  print(names(covered_df)[iglom])
  glom = names(covered_df)[iglom]
  covered_spots = covered_df[[ glom]]
  filter_spots = intersect(covered_spots$slice_barcode,covered_spots_select$spot)
  print(length(filter_spots))
  
  covered_df_filtered[[iglom]] = covered_df_filtered[[iglom]] %>%
    filter(slice_barcode %in% filter_spots)
}

saveRDS(covered_df_filtered,paste0(radius_dir,"/frustum_covered_spots.rds")) 












result_df = read.table(paste0(radius_dir,"/frustum_top_bottom_center.txt"))
mitral_df = coord3d_IRIS_mitral
in_idx_list = readRDS(file.path(radius_dir, "frustum_covered_spots.rds"))
all_spots <- coord3d_IRIS %>%
  rename(
    x = x.x,
    y = y.x,
    z = z.x
  )


##### Visualization for covered spots to check algoritm work
unit <- function(v) v / sqrt(sum(v*v))
cross3 <- function(u, v) c(
  u[2]*v[3] - u[3]*v[2],
  u[3]*v[1] - u[1]*v[3],
  u[1]*v[2] - u[2]*v[1]
)

# Build an orthonormal frame {u1, u2} perpendicular to axis a
perp_frame <- function(a) {
  a <- unit(a)
  ref <- if (abs(a[1]) < 0.9) c(1,0,0) else c(0,1,0)
  u1 <- unit(cross3(a, ref))
  u2 <- cross3(a, u1)  # already unit because a ⟂ u1 and both unit
  list(u1 = u1, u2 = u2, a = a)
}


make_frustum_mesh <- function(upper_center, lower_center, Rt, Rb,
                              n = 24, close_caps = TRUE) {
  a <- lower_center - upper_center
  h <- sqrt(sum(a*a))
  if (h == 0) stop("Upper and lower centers coincide (zero height).")
  frm <- perp_frame(a); u1 <- frm$u1; u2 <- frm$u2
  
  # param angles
  theta <- seq(0, 2*pi, length.out = n + 1)[-(n+1)]
  
  # top & bottom rings
  top_pts <- t(vapply(theta, function(t)
    upper_center + Rt * (cos(t)*u1 + sin(t)*u2), numeric(3)))
  bot_pts <- t(vapply(theta, function(t)
    lower_center + Rb * (cos(t)*u1 + sin(t)*u2), numeric(3)))
  
  # vertices: top ring (1..n), bottom ring (n+1..2n), optional centers
  verts <- rbind(top_pts, bot_pts)
  idx_top <- seq_len(n)
  idx_bot <- n + seq_len(n)
  
  if (close_caps) {
    verts <- rbind(verts, upper_center, lower_center)
    idx_top_center <- nrow(verts) - 1
    idx_bot_center <- nrow(verts)
  }
  
  # side triangles (two per segment)
  side_tris <- do.call(rbind, lapply(1:n, function(k) {
    k2 <- if (k < n) k + 1 else 1
    tk  <- idx_top[k];  tk2 <- idx_top[k2]
    bk  <- idx_bot[k];  bk2 <- idx_bot[k2]
    rbind(
      c(tk,  bk,  bk2),   # tri 1
      c(tk,  bk2, tk2)    # tri 2
    )
  }))
  
  # caps (optional): triangle fans
  top_tris <- bot_tris <- NULL
  if (close_caps) {
    top_tris <- do.call(rbind, lapply(1:n, function(k) {
      k2 <- if (k < n) k + 1 else 1
      c(idx_top_center, idx_top[k2], idx_top[k])
    }))
    bot_tris <- do.call(rbind, lapply(1:n, function(k) {
      k2 <- if (k < n) k + 1 else 1
      c(idx_bot_center, idx_bot[k], idx_bot[k2])
    }))
  }
  
  # convert to 0-based for plotly::add_trace(type="mesh3d")
  to0 <- function(M) list(i = M[,1]-1, j = M[,2]-1, k = M[,3]-1)
  
  out <- list(
    x = verts[,1], y = verts[,2], z = verts[,3],
    side = to0(side_tris),
    top  = if (!is.null(top_tris)) to0(top_tris) else NULL,
    bot  = if (!is.null(bot_tris)) to0(bot_tris) else NULL
  )
  out
}

# Output dir for HTML files
out_dir <- paste0(radius_dir,"frustum_3d_covered_spot")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Colors
col_spots_all   <- "#7f7f7f"  # mitral spots (all gray)
col_spots_cover <- "#ff7f0e"                 # covered spots (orange)
col_glom        <- "#1f77b4"                 # top center (blue)
col_mitral      <- "#111111"                 # nearest mitral spot (black)
col_mesh        <- "#FFFF99"    # frustum surface (translucent)
col_link        <- "#d62728"                 # link line (red)
col_base        <- "#e41a1c"                 # NEW: extended base center (red-ish)                # link line (red)

# Optional: downsample mitral_df for speed if huge
mitral_plot_df <- mitral_df
# set.seed(1); if (nrow(mitral_plot_df) > 8000) mitral_plot_df <- mitral_plot_df[sample(nrow(mitral_plot_df), 8000), ]

file_safe <- function(s) gsub("[^A-Za-z0-9._-]", "_", s)

for (i in 1:nrow(result_df)) {
  gid <- if ("glomerulus_id" %in% names(result_df)) result_df$glomerulus_id[i] else paste0("g", i)
  sid <- if ("nearest_spot_id" %in% names(result_df)) result_df$nearest_spot_id[i] else paste0("s", i)
  
  # Centers
  up  <- c(result_df$g_x[i], result_df$g_y[i], result_df$g_z[i])  # glomerulus (top)
  low <- c(result_df$b_x[i], result_df$b_y[i], result_df$b_z[i])  # extended base center B (CHANGED)
  
  # Frustum mesh (single) – now G → B
  fr <- make_frustum_mesh(up, low, Rt = Rt, Rb = Rb, n = 36, close_caps = TRUE)
  
  # ALL mitral-layer spots (as-is)
  mitral_all <- mitral_df
  
  # Covered all_spots for this frustum
  idx_in <- in_idx_list[[i]]
  bar_in <- as.character(idx_in$slice_barcode)
  covered_df <- subset(all_spots, slice_barcode %in% bar_in)
  
  # Nearest mitral spot S (same as before)
  if ("nearest_spot_id" %in% names(result_df) && "slice_barcode" %in% names(mitral_df)) {
    bottom_df <- subset(mitral_df, slice_barcode == result_df$nearest_spot_id[i])
  } else {
    tol <- 1e-6
    bottom_df <- subset(
      mitral_df,
      abs(x - result_df$s_x[i]) <= tol &
        abs(y - result_df$s_y[i]) <= tol &
        abs(z - result_df$s_z[i]) <= tol
    )
  }
  
  # 🔴 NEW: extended base center B as its own tiny df
  base_df <- data.frame(
    x = result_df$b_x[i],
    y = result_df$b_y[i],
    z = result_df$b_z[i]
  )
  
  # apex–base link (G → B, not G → S)
  seg_x <- c(up[1], low[1], NA)
  seg_y <- c(up[2], low[2], NA)
  seg_z <- c(up[3], low[3], NA)
  
  p <- plot_ly() |>
    # (1) ALL mitral-layer spots
    add_markers(
      data = mitral_all,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      inherit = FALSE,
      name = paste0("Mitral spots (all, n=", nrow(mitral_all), ")"),
      marker = list(size = 2.5, opacity = 0.5, color = col_spots_all),
      hoverinfo = "text",
      text = ~slice_barcode,
      showlegend = TRUE
    ) |>
    # (2) Covered spots
    add_markers(
      data = covered_df,
      x = ~x, y = ~y, z = ~z,
      name = paste0("Covered spots (", nrow(covered_df), ")"),
      type = "scatter3d", mode = "markers",
      marker = list(size = 4, opacity = 0.95, color = col_spots_cover),
      hoverinfo = "text",
      text = ~paste0(slice_barcode, "<br>Slice: ", Slice)
    ) |>
    # (3) Frustum mesh
    add_trace(
      type = "mesh3d",
      x = fr$x, y = fr$y, z = fr$z,
      i = c(fr$side$i, if (!is.null(fr$top)) fr$top$i, if (!is.null(fr$bot)) fr$bot$i),
      j = c(fr$side$j, if (!is.null(fr$top)) fr$top$j, if (!is.null(fr$bot)) fr$bot$j),
      k = c(fr$side$k, if (!is.null(fr$top)) fr$top$k, if (!is.null(fr$bot)) fr$bot$k),
      color = col_mesh, opacity = 0.40, showscale = FALSE,
      name = "Frustum"
    ) |>
    # (4) Glomerulus (top)
    add_markers(
      x = up[1], y = up[2], z = up[3],
      type = "scatter3d", mode = "markers",
      marker = list(size = 8, color = col_glom),
      name = "Glomerulus (top)"
    ) |>
    # (5) Nearest mitral spot S (black)
    add_markers(
      data = bottom_df,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      marker = list(size = 8, color = col_mitral),
      name = "Nearest mitral spot"
    ) |>
    # (6) 🔴 Extended base center B
    add_markers(
      data = base_df,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers",
      marker = list(size = 8, color = col_base),
      name = "Extended base center"
    ) |>
    # (7) Link G → B
    add_trace(
      x = seg_x, y = seg_y, z = seg_z,
      type = "scatter3d", mode = "lines",
      line = list(width = 4, color = col_link),
      hoverinfo = "none", showlegend = FALSE
    ) |>
    layout(
      title = paste0("Frustum ", gid, " → ", sid, "  (Rt=", Rt, ", Rb=", Rb, ", +", extend_height, "µm)"),
      legend = list(orientation = "h", y = -0.05),
      scene = list(
        xaxis = list(title = "x"),
        yaxis = list(title = "y"),
        zaxis = list(title = "z"),
        aspectmode = "data"
      ),
      margin = list(l=0, r=0, b=0, t=30)
    )
  
  fname <- file.path(out_dir, sprintf("frustum_%s_to_%s.html",
                                      file_safe(gid), file_safe(sid)))
  htmlwidgets::saveWidget(as_widget(p), file = fname, selfcontained = TRUE)
}



##### Covered Spots Analysis ######
covered_df = readRDS(paste0(radius_dir,"/frustum_covered_spots.rds"))

#### Number of covering spots per cone
covered_count = data.frame(count = unlist(lapply(covered_df,nrow)))

p = ggplot(covered_count,aes(x = count)) +
  geom_histogram() +
  labs(title = "Number of Capturing Spots Per Frustum", x= "Number of Capturing Spots") +
  theme(plot.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

pdf(paste0(radius_dir,"/number_of_covering_spots_histogram.pdf"),width = 8,height = 8)
p
dev.off()


#### Number of covering slices per cone
covered_count = data.frame(count = unlist(lapply(covered_df,
                                                 function(x){
                                                   length(unique(x$Slice))
                                                 })))

p = ggplot(covered_count,aes(x = count)) +
  geom_histogram() +
  labs(title = "Number of Capturing Slices Per Frustum", x= "Number of Capturing Slices") +
  theme(plot.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15))

pdf(paste0(radius_dir,"/number_of_covering_slices_histogram.pdf"),width = 8,height = 8)
p
dev.off()



#### Overlap between Frustum
# ---- 1) Normalize to character vectors of barcodes per frustum ----
get_barcodes <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.data.frame(x)) {
    if (!"slice_barcode" %in% names(x)) return(character(0))
    return(as.character(x$slice_barcode))
  }
  if (is.character(x)) return(x)
  # if it’s integer indices into all_spots, convert if you prefer (needs all_spots)
  stop("Unsupported element in covered_df ; expected data.frame with 'slice_barcode' or character vector.")
}

sets <- lapply(covered_df , get_barcodes)

# ---- 2) Build a sparse incidence matrix: rows = frusta, cols = unique barcodes ----
all_barcodes <- unique(unlist(sets))
if (length(all_barcodes) == 0L) stop("No barcodes found in covered_df .")

barcode_index <- setNames(seq_along(all_barcodes), all_barcodes)

# i = row indices (frustum id), j = column indices (barcode id)
i_idx <- integer(0)
j_idx <- integer(0)
for (i in seq_along(sets)) {
  if (length(sets[[i]]) == 0L) next
  j <- barcode_index[ sets[[i]] ]
  j <- j[!is.na(j)]
  if (length(j)) {
    i_idx <- c(i_idx, rep.int(i, length(j)))
    j_idx <- c(j_idx, j)
  }
}

# sparse matrix (frusta x barcodes) with 1 where barcode is covered by frustum
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
library(Matrix)
M <- sparseMatrix(i = i_idx, j = j_idx, x = 1L,
                  dims = c(length(sets), length(all_barcodes)),
                  dimnames = list(NULL, all_barcodes))

# ---- 3) Pairwise intersections & unions, then Jaccard ----
# intersections = M %*% t(M) (tcrossprod) -> counts of shared barcodes
intersections <- tcrossprod(M)          # (n_frusta x n_frusta), sparse

# sizes per frustum
sizes <- Matrix::rowSums(M)             # vector length n_frusta

# unions = |A| + |B| - |A∩B|
# convert intersections to "dense" only at the end for plotting
unions <- outer(as.numeric(sizes), as.numeric(sizes), "+") - as.matrix(intersections)

# Jaccard = intersection / union; define 0 where union == 0
J <- as.matrix(intersections)
J[unions > 0] <- J[unions > 0] / unions[unions > 0]
J[unions == 0] <- 0

# ---- 4) Optional: add labels from result_df (e.g., glomerulus_id) ----
if (exists("result_df") && "glomerulus_id" %in% names(result_df)) {
  rownames(J) <- colnames(J) <- as.character(result_df$glomerulus_id)
} else {
  rownames(J) <- colnames(J) <- paste0("F", seq_len(nrow(J)))
}

# ---- 5) Draw a heatmap ----
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(pheatmap)

# For large matrices, clustering can be slow; set cluster_* = FALSE to keep order
pheatmap(J,
         color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(200),
         breaks = seq(0, 1, length.out = 201),
         legend = TRUE,
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Jaccard similarity of covered barcodes across frusta")

# ---- 6) Optional: save to file ----
# ggsave won't capture pheatmap directly; use png()/dev.off()
pdf(paste0(radius_dir,"/jaccard_frusta_heatmap.pdf"), width = 12, height = 10)
pheatmap(J,
         color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(200),
         breaks = seq(0, 1, length.out = 201),
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         main = "Jaccard similarity of covered barcodes across frusta")
dev.off()






# ---- deps ----
# install.packages("ggplot2")  # if needed
library(ggplot2)

# ---- helpers ----
unit   <- function(v) v / sqrt(sum(v*v))
cross3 <- function(u, v) c(
  u[2]*v[3] - u[3]*v[2],
  u[3]*v[1] - u[1]*v[3],
  u[1]*v[2] - u[2]*v[1]
)

# orthonormal frame in plane ⟂ to axis a
perp_frame <- function(a) {
  a <- unit(a)
  axis <- which.min(abs(a))
  ref  <- c(0,0,0); ref[axis] <- 1
  u1 <- unit(cross3(a, ref))
  u2 <- unit(cross3(a, u1))
  list(u1=u1, u2=u2, a=a)
}

# build top/bottom ring points in 3D then drop z → XY
frustum_rings_xy <- function(upper_center, lower_center, Rt=Rt, Rb=Rb, n=64) {
  a <- lower_center - upper_center
  h <- sqrt(sum(a*a))
  if (!is.finite(h) || h < 1e-12) stop("Degenerate frustum (height ~ 0).")
  f <- perp_frame(a); u1 <- f$u1; u2 <- f$u2
  theta <- seq(0, 2*pi, length.out = n + 1)[-(n+1)]
  top3 <- t(vapply(theta, function(t) upper_center + Rt*(cos(t)*u1 + sin(t)*u2), numeric(3)))
  bot3 <- t(vapply(theta, function(t) lower_center + Rb*(cos(t)*u1 + sin(t)*u2), numeric(3)))
  list(top_xy = top3[,1:2, drop=FALSE], bot_xy = bot3[,1:2, drop=FALSE])
}


# safe filename helper
file_safe <- function(s) {
  s <- as.character(s); s <- iconv(s, to = "ASCII//TRANSLIT")
  s <- gsub("[^A-Za-z0-9._-]", "_", s); s <- gsub("_+", "_", s)
  substr(s, 1, 150)
}

# ---- parameters ----
n_ring <- 16

# ---- output dir ----
out_dir <- paste0(radius_dir,"frustum_xy_with_overlap_mitral_layer")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# input data
all_spots = coord3d_IRIS %>%
  rename(
    x = x.x,
    y = y.x,
    z = z.x
  )
in_idx_list = covered_df
# ---- loop over frusta ----
for (i in seq_len(nrow(result_df))) {
  # centers
  up  <- c(result_df$g_x[i], result_df$g_y[i], result_df$g_z[i])  # glomerulus (top center)
  low <- c(result_df$b_x[i], result_df$b_y[i], result_df$b_z[i])  # extended base center
  
  # 2D rings
  rings <- try(frustum_rings_xy(up, low, Rt=Rt, Rb=Rb, n=n_ring), silent = TRUE)
  if (inherits(rings, "try-error")) { message("Skip ", i, " (degenerate)"); next }
  top_df <- data.frame(x=rings$top_xy[,1], y=rings$top_xy[,2])
  bot_df <- data.frame(x=rings$bot_xy[,1], y=rings$bot_xy[,2])
  
  # ---- captured spots (use in_idx_list[[i]] directly) ----
  # Expect a 2-col df with 'slice_barcode' and 'Slice'
  cap_keys <- in_idx_list[[i]]
  if (!is.null(cap_keys) && nrow(cap_keys)) {
    # join to get XY for captured barcodes
    m <- match(cap_keys$slice_barcode, all_spots$slice_barcode)
    m <- m[!is.na(m)]
    covered_xy <- all_spots[m, c("slice_barcode","Slice","x","y")]
  } else {
    covered_xy <- all_spots[0, c("slice_barcode","Slice","x","y")]
  }
  
  # ---- overlap-by-slice: all spots whose Slice is in captured slices ----
  overlap_xy <- if (nrow(covered_xy)) {
    subset(all_spots,
           Slice %in% unique(covered_xy$Slice),
           select = c(slice_barcode, Slice, x, y, iris_domain))
  } else {
    all_spots[0, c("slice_barcode","Slice","x","y","iris_domain")]
  }
  
  # centers in 2D
  up2  <- up[1:2]; low2 <- low[1:2]
  
  # local bounds
  x_all <- c(top_df$x, bot_df$x, overlap_xy$x)
  y_all <- c(top_df$y, bot_df$y, overlap_xy$y)
  if (!length(x_all)) { x_all <- c(up2[1], low2[1]); y_all <- c(up2[2], low2[2]) }
  pad <- 40
  xlim_local <- range(x_all, na.rm=TRUE) + c(-pad, pad)
  ylim_local <- range(y_all, na.rm=TRUE) + c(-pad, pad)
  
  gid <- if ("glomerulus_id" %in% names(result_df)) result_df$glomerulus_id[i] else paste0("g", i)
  sid <- if ("nearest_spot_id" %in% names(result_df)) result_df$nearest_spot_id[i] else paste0("s", i)
  
  # Extended base center B (already using b_x, b_y)
  bottom_df <- data.frame(
    slice_barcode = result_df$nearest_spot_id[i],
    Slice         = result_df$nearest_spot_slice[i],
    x             = result_df$b_x[i],
    y             = result_df$b_y[i]
  )
  
  # 🔸 NEW: nearest mitral spot S (original, not extended)
  nearest_df <- data.frame(
    slice_barcode = result_df$nearest_spot_id[i],
    Slice         = result_df$nearest_spot_slice[i],
    x             = result_df$s_x[i],
    y             = result_df$s_y[i]
  )
  
  seg_df <- data.frame(
    x1 = top_df$x, y1 = top_df$y,
    x2 = bot_df$x, y2 = bot_df$y,
    k  = seq_len(nrow(top_df))
  )
  
  # include slice of nearest spot in the levels (in case it wasn't in covered_xy)
  slice_levels <- sort(unique(c(covered_xy$Slice, result_df$nearest_spot_slice[i])))
  
  # (Optional) thin the segments to avoid clutter:
  k_step <- max(1, floor(nrow(seg_df) / 24))   # ~24 spokes
  seg_df <- seg_df[seq(1, nrow(seg_df), by = k_step), , drop = FALSE]
  
  overlap_xy <- overlap_xy %>% filter(iris_domain == "Leiden Domain6")
  
  p <- ggplot() +
    # (S) spokes
    geom_segment(data = seg_df,
                 aes(x = x1, y = y1, xend = x2, yend = y2, group = k),
                 linewidth = 0.35, alpha = 0.35, color = "#d62728") +
    geom_segment(aes(x = up2[1], y = up2[2], xend = low2[1], yend = low2[2]),
                 linewidth = 0.8, color = "#d62728", linetype = "22") +
    # (BG) overlapping-slice spots
    geom_point(
      data = overlap_xy,
      aes(x = x, y = y, color = Slice),
      size = 1, alpha = 1, show.legend = TRUE
    ) +
    # (R) frustum rings
    geom_path(data = bot_df, aes(x = x, y = y, group = 1),
              linewidth = 0.6, color = "#7f7f7f") +
    geom_path(data = top_df, aes(x = x, y = y, group = 1),
              linewidth = 0.6, color = "#1f77b4") +
    # (Centers) top center as blue dot
    geom_point(aes(x = up2[1], y = up2[2]),
               color = "#1f77b4", size = 1.8, show.legend = FALSE) +
    # ⭐ Extended base center B (star)
    geom_point(
      data = bottom_df,
      aes(x = x, y = y, color = Slice),
      fill = "#FF073A",color = "#FF073A",
      shape = 8, size = 2.4, stroke = 1.0, show.legend = FALSE
    ) +
    # 🔸 NEW: nearest mitral spot S (orange filled circle)
    geom_point(
      data = nearest_df,
      aes(x = x, y = y, color = Slice),
      shape = 21, color = "#39FF14",fill = "#39FF14", size = 2.6, stroke = 0.8, show.legend = FALSE
    ) +
    coord_equal(xlim = xlim_local, ylim = ylim_local, expand = FALSE) +
    scale_color_discrete(limits = slice_levels, name = "Overlapping slices") +
    labs(
      title = sprintf("Frustum XY: %s → %s (Rt=%d, Rb=%d)", gid, sid, Rt, Rb),
      subtitle = "Blue = glomerulus; ⭐ = extended base; orange dot = nearest mitral spot; \nRed = frustum axis; colored = overlapping slices",
      x = "x", y = "y"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right") +
    guides(color = guide_legend(ncol = 2, override.aes = list(alpha = 0.8, size = 1.5)))
  
  # save
  fname <- file.path(out_dir, sprintf("frustum_%s_to_%s.png",
                                      file_safe(gid), file_safe(sid)))
  ggsave(fname, p, width = 12, height = 8, dpi = 300)
}

message("Done. PNGs in: ", normalizePath(out_dir))
