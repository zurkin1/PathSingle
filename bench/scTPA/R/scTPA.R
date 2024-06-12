start_call = Sys.time()

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("bigstatsr"))
suppressPackageStartupMessages(library("optparse"))
#suppressPackageStartupMessages(library("scImpute"))
#suppressPackageStartupMessages(library("Rmagic"))
#suppressPackageStartupMessages(library("scran"))
#suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("pheatmap"))


options(digits=5)

homo_paths = c('kegg','reactome','biocarta','smpdb','humancyc',
               'nci','panther','pharmgkb','acsn2','rb',
               'c2.cgp','c2.cp',
               'c4.cgn','c4.cm',
               'c5.bp','c5.mf','c5.cc',
               'c6.all','c7.all','h.all')

mus_paths = c('kegg','reactome','pathbank','smpdb',
              'c5.bp','c5.mf','c5.cc','other')

option_list = list(
  make_option(c("-f", "--file"),
              type = "character",
              default = "NULL",
              dest = 'file',
              help = "gene expression profile, genes X cells",
              metavar = "file"),

              make_option(c("--cellType"),
              type = "character",
              default = "NULL",
              dest = 'cellType',
              help = "cell type file. First column is cell name (same as the colnames of gene expression profile), second column is cell type. No header names.[default= %default]",
              metavar = "cellType"),

              make_option(c("--work_dir"),
              type = "character",
              default = "./",
              dest = 'work_dir',
              help = "Workshop direction. [default= %default]",
              metavar = "work_dir"),

              make_option(c("--normalize"),
              type = "character",
              default = 'none',
              dest = 'normalize_method',
              help = "methods used for normalization. one of log, CLR, RC, scran[default= %default]",
              metavar = "normalize_method"),

              make_option(c("--min_cells"),
              type = "integer",
              default = 3,
              dest = 'min_cells',
              help = "genes must be in a minimum number of cells. Used for filtering genes[default= %default]",
              metavar = "min_cells"),

              make_option(c("--min_features"),
              type = "integer",
              default = 200,
              dest = 'min_features',
              help = "cells must have at least the minimum number of genes. Used for filtering cells[default= %default]",
              metavar = "min_features"),

              make_option(c("--species"),
              type = "character",
              default = "homo",
              dest = 'species',
              help = "species[default= %default]",
              metavar = "species"),

              make_option(c("--imputation"),
              type = "character",
              default = "none",
              dest = 'imputation',
              help = "Imputation method[default= %default]",
              metavar = "imputation"),

              make_option(c("--data_type"),
              type = "character",
              default = 'TPM',
              dest = 'data_type',
              help = "data type of gene expression profile TPM or count [default= %default]",
              metavar = "file"),

              make_option(c("--pathway_database"),
              type = "character",
              default = "kegg",
              dest = 'pathway_database',
              help = "pathway database, detials see https://github.com/sulab-wmu/scTPA[default= %default]",
              metavar = "pathway_database"),

              make_option(c("--user_pathway"),
              type = "character",
              default = "NULL",
              dest = "user_pathway",
              help = "user defined pathway file，only for gmt format[default = %default]",
              metavar = "user_pathway"),

              make_option(c("--pas_method"),
              type = "character",
              default = "gsva",
              dest = 'pas_method',
              help = "method for calculating PAS. gsva, ssgsea, zscore or plage [default= %default]",
              metavar = "pas_method"),

              make_option(c("--para_size"),
              type = "integer",
              default = 4,
              dest = 'para_size',
              help = "number of kernels used for parallel[default= %default]",
              metavar = "para_size"),

              make_option(c("--cluster_method"),
              type = "character",
              default = "seurat",
              dest = 'cluster_method',
              help = "clustering method. seurat, hclust, simlr, kmedoids, kmeans or dbscan [default= %default]",
              metavar = "cluster_method"),

              make_option(c("--seurat_dims"),
              type = "integer",
              default = 8,
              dest = 'seurat_dims',
              help = "dimensions used in Seurat clustering[default= %default]",
              metavar = "seurat_dims"),

              make_option(c("--seurat_resolution"),
              type = "double",
              default = 0.5,
              dest = 'seurat_resolution',
              help = "resolution used for Seurat clustering[default= %default]",
              metavar = "seurat_resolution"),

              make_option(c("--k_cluster"),
              type = "integer",
              default = 5,
              dest = 'k_cluster',
              help = "number of clusters, useless if clustering method is Seurat or dbscan[default= %default]",
              metavar = "k_cluster"),

              make_option(c("--min_pts"),
              type = "integer",
              default = 3,
              dest = 'min_pts',
              help = "parameter in DBSCAN[default= %default]",
              metavar = "min_pts"),

              make_option(c("--dims"),
              type = "integer",
              default = 20,
              dest = 'dims',
              help = "number of PCA dimensionas used for TSNE or UMAP[default= %default]"),

              make_option(c("--marker_method"),
              type = "character",
              default = "wilcox",
              dest = "marker_method",
              help = "method of finding siginificant markers[default= %default]",
              metavar="find_maker_method"),

              make_option(c("--logFC_thre"),
              type = "double",
              default = 0.25,
              dest = "logFC_thre",
              help = "threshold of logFC (Detail see Seurat)[default= %default]",
              metavar="threshold_logFC"),

              make_option(c("--min_pct"),
              type = "double",
              default = 0.1,
              dest = "min_pct",
              help = "only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations.[default= %default]",
              metavar="min_pct"),

              make_option(c("--pic_type"),
              type = "character",
              default = 'png',
              dest = "pic_type",
              help = "type of picture, png or pdf [default= %default]",
              metavar="pic_type"),

              make_option(c("-o","--out_dir"),
              type = "character",
              default = "NULL",
              dest = "out_dir",
              help = "output folder[default= %default]",
              metavar="out_dir")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$file <- ".\\test\\expression.csv"
opt$cellType <- ".\\test\\cell_type.csv"
opt$work_dir <- "."
opt$species <- "homo"
opt$pathway_database <- "kegg"
opt$para_size <- 1
opt$out_dir <- ".\\results\\"

pipeline_cal = function(expr_path,
                        cell_type_path,
                        data_type = c('TPM','count'),
                        normalize_method = c('log', 'CLR', 'RC', 'none','scran'),
                        species = c('homo', 'mus'),
                        min_cells = 3,
                        min_features = 200,
                        imputation = c("none", "scImpute", "magic"),
                        pathway_database,
                        pas_method = c('ssgsea','gsva','plage','zscore'),
                        user_pathway = NULL,
                        para_size = 4,
                        cluster_method,
                        shown_method = c('tsne','umap'),
                        seurat_dims = 8,
                        seurat_resolution = 0.5,
                        k_cluster = 5,
                        min_pts = 3,
                        dims = 20,
                        marker_method = c('wilcox','bimod','MAST','t','LR'),
                        logFC_thre = 0.25,
                        min.pct = 0.1,
                        pic_type = c('png','pdf'),
                        out_dir,
                        work_dir
)
{
  src_dir = file.path(work_dir,'src')
  pathway_dir = file.path(work_dir,'data/pathways')
  r_dir = file.path(work_dir,'R')

  source(file.path(r_dir,'visualization.R'))
  source(file.path(r_dir,'PAS.R'))
  source(file.path(r_dir,'function.R'))
  source(file.path(r_dir,'markers.R'))
  source(file.path(r_dir,'clustering.R'))

  out_dir = paste0(out_dir, '/')
  out_dir = gsub('/+','/',out_dir)

  min_cells = max(2, min_cells)
  idType='symbol'

  start_time = Sys.time()
  cat("start... \n")

  ##################################################
  #####   step 1. load expression profile   ########
  ##################################################

  if(is.null(expr_path)){
    stop("Please specific input gene expression file path through --file or -f")
  }

  if(is.null(out_dir)){
    stop("Please specific output folder through --out_dir")
  }

  if(!file_test("-d", out_dir)){
    dir.create(out_dir,recursive=TRUE)
  }

   ret = load_expr(expr_path,
                   cell_type_path,
                   imputation,
                   data_type,
                   species,
                   out_dir,
                   para_size,
                   normalize_method = normalize_method,
                   min_cells = min_cells,
                   min_features = min_features,
                   pas_method = pas_method,
                   work_dir = work_dir)

  expr = ret$expr
  cell_type = ret$cell_type
  genes_name = ret$genes_name
  cells_name = ret$cells_name
  #print("####################6##############")
  #print(genes_name[1:10])

  rm(ret, min_cells, min_features, normalize_method, data_type)
  gc()
  genes_name = toupper(genes_name)


  ##################################################
  ########    step 2. select pathways     ##########
  ##################################################

  if(species == 'homo'){
    if(pathway_database %in% homo_paths){
      pathway_dir = file.path(pathway_dir, species, paste0(pathway_database,'.rds'))
      path_list = readRDS(pathway_dir)
      #weight_list = 'none'
    }else{
      stop("Do not have this pathway database for homo species,
           please ensure your specific pathway database is legal")
    }

  }else if(species == 'mus'){
    if(pathway_database %in% mus_paths){
      pathway_dir = file.path(pathway_dir, species, paste0(pathway_database,'.rds'))
      path_list = readRDS(pathway_dir)
      #weight_list = 'none'
    }else{
      stop("Do not have this pathway database for mus musculus,
           please ensure your specific pathway database is legal")
    }
  }else{
    stop("The --species must be homo or mus")
  }


  ##################################################
  ####  step-optional. user defined pathways   #####
  ##################################################

  if(! user_pathway == 'NULL'){
    cat(paste("Loading user defined pathway from:",user_pathway))
    cat('\n')
    user_pathway = load_gmt_data(user_pathway)
    aa = gsub('-|,|/| |\\(|\\)','_',names(user_pathway))
    aa = gsub('_+','_',aa)
    aa = gsub('_\\b','',aa)
    aa = gsub(',','',aa)
    names(user_pathway) = paste0('user_', names(user_pathway))
    path_list = c(path_list, user_pathway)
    #weight_list = 'none'
  }

  cat("pathway preparation successful\t")

  gc()


  ##################################################
  ########    step 3. calculate PAS       ##########
  ##################################################
  pas_method <- c('plage')
  cluster_mat = PAS(expr,
                    path_list,
                    genes_name = genes_name,
                    cells_name = cells_name,
                    src_folder = src_dir,
                    method = pas_method,
                    parallel.sz = para_size,
                    verbose=TRUE)


  pas_dir = file.path(out_dir,'pas.csv')
  write.csv(cluster_mat, file = pas_dir, quote=F)

  cat("Calculating pathway activity score successful.\n")
  rm(pas_dir)
  gc()

  #print(cluster_mat[1:3,1:3])


  ##################################################
  ########    step 4. cluster idents      ##########
  ##################################################

  if(cell_type == 'none'){
    cat("Start clustering...\n")
    cluster_idents = cal_cluster(cluster_mat,
                                 method=cluster_method,
                                 seurat_dims = seurat_dims,
                                 seurat_resolution = seurat_resolution,
                                 k_cluster = k_cluster,
                                 min_pts = min_pts)
    cluster_idents = as.numeric(cluster_idents)

    ### ensure the minimum cluster number is 1
    if(min(cluster_idents) == 0){
        cluster_idents = cluster_idents + 1
    }

    cluster_idents = factor(cluster_idents,levels = sort(unique(cluster_idents)))
    #cluster_idents = paste0('C',cluster_idents)
    names(cluster_idents) = colnames(cluster_mat)
    #print(cluster_idents[1:10])
    cluster_idents = sort(cluster_idents)
    #print(cluster_idents[1:10])
    cluster_mat = cluster_mat[,names(cluster_idents)]
    cluster_idents = factor(paste0('C',cluster_idents),levels = paste0('C',levels(cluster_idents)))
    names(cluster_idents) = colnames(cluster_mat)
    #print(cluster_idents[1:10])
    #pas = pas[,names(cluster_idents)]
    cat("Clustering cells successful\n")
  }else{
    cat("known cell type","\n")
    cluster_idents = factor(cell_type,levels = sort(unique(cell_type)))
    names(cluster_idents) = colnames(cluster_mat)
    #print(cluster_idents[1:10])
    cluster_idents = sort(cluster_idents)
    #print(cluster_idents[1:10])
    cluster_mat = cluster_mat[,names(cluster_idents)]
  }





  ##################################################
  ########    step 5. visulization        ##########
  ##################################################

  cells = c()
  if(ncol(cluster_mat) > 10000){
    ## extract some cells if the mat is too big
    print("extract subset cells")
    pst = 10000/ncol(cluster_mat)
    #print(dim(mat))
    #mat_new = matrix(nrow=nrow(mat))
    for(clust in unique(cluster_idents)){
      #print(clust)
      curr_cells = which(cluster_idents == clust)
      #print(length(curr_cells))
      curr_cells = sample(curr_cells, round(length(curr_cells)*pst))
      cells = c(cells, curr_cells)
      #curr_mat = mat[,curr_cells]
      #print(dim(curr_mat))
      #mat_new = cbind(mat_new, curr_mat)
      #print(dim(mat_new))
    }
    cat("extract subset of original expression profile, dimensions of new matrix is ",length(cells),"\n")
    #mat = mat_new[,2:ncol(mat_new)]
    #cat("extract subset of original expression profile, dimensions of new matrix is ",paste(dim(mat),collapse=','),"\n")
    #rm(mat_new)
    #gc()
  }else{
    cells = 1:ncol(cluster_mat)
  }
  #print("#################   cells    #######################")
  #print(cells[1:10])
  #print(length(cells))


  options(bitmapType='cairo')
  plotPAS_heatmap(mat = cluster_mat[,cells],
                  cluster_idents = cluster_idents[cells],
                  out_path = out_dir,
                  pic_type = pic_type)



  obSeurat_pas = prepa_seuratOb(cluster_mat,
                                cluster_idents,
                                seurat_dims,
                                out_dir)
  rm(cluster_mat)
  gc()

  plot_cluster(obSeurat = obSeurat_pas, out_dir=out_dir,pic_type = pic_type)

  cat("Ploting cluster picture successful\n")



  ##################################################
  ########    step 6. finding markers      #########
  ##################################################



  errMarkers = FALSE
  if(user_pathway == 'NULL'){
    parseAllMarkers(obSeurat_pas,
                    NULL,
                    method = marker_method,
                    logfc_thre = logFC_thre,
                    min.pct=min.pct,
                    out_dir = out_dir,
                    para_size = para_size,
                    expr = expr,
					genes_name = genes_name,
                    cells_name = cells_name,
                    cell_type = cluster_idents,
                    path_list = path_list,
                    pic_type = pic_type,
                    cells = cells)
  }else{
    parseAllMarkers(obSeurat_pas,
                    names(user_pathway),
                    method = marker_method,
                    logfc_thre = logFC_thre,
                    min.pct = min.pct,
                    out_dir = out_dir,
                    para_size = para_size,
                    expr = expr,
					genes_name = genes_name,
                    cells_name = cells_name,
                    cell_type = cluster_idents,
                    path_list = path_list,
                    pic_type = pic_type,
                    cells = cells)
  }

  cat("Finding marker pathways successful\n")
}








#print(opt$pas_method)

start_time = Sys.time()
pipeline_cal(expr_path = opt$file,
             cell_type = opt$cellType,
             normalize_method = opt$normalize_method,
             min_cells = opt$min_cells,
             min_features = opt$min_features,
             data_type = opt$data_type,
             imputation = opt$imputation,
             species = opt$species,
             pathway_database = opt$pathway_database,
             pas_method = opt$pas_method,
             user_pathway = opt$user_pathway,
             para_size = opt$para_size,
             cluster_method = opt$cluster_method,
             seurat_dims = opt$seurat_dims,
             seurat_resolution = opt$seurat_resolution,
             k_cluster = opt$k_cluster,
             min_pts = opt$min_pts,
             dims = opt$dims,
             marker_method = opt$marker_method,
             pic_type = opt$pic_type,
             logFC_thre = opt$logFC_thre,
             min.pct = opt$min_pct,
             out_dir = opt$out_dir,
             work_dir = opt$work_dir
)
print(Sys.time()-start_call)