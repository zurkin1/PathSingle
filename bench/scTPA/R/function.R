
#' basic functions for scTPA
#'
#' import scImpute
#' import Rmagic


## get all genes in pathways
#' @param paths: pathway names in vector format
#' @param path_list: pathways in List format
#' @param expr_genes: gene names in expression profile
#' @return genes in pathways, every item represent genes in certain pathway

#library(pryr)
library(data.table)

get_geneLis = function(paths,
                       path_lis,
                       expr_genes){
  res = sapply(paths,FUN = function(x){
    genes = path_lis[[x]]
    genes = intersect(genes,expr_genes)
    paste(genes,collapse = ' ')
  })
  as.character(res)
}

#' convert Ensembl ID to gene name
#' convert files are download from Ensembl Gene Browser
#' DATABASE: Ensembl Genes 99
#' homo/Human Genes (GRCh38.p13)
#' mus/Mouse Genes (GRCm38.p6)
#' @param expr: expression profile
#' @param species: species
#' @param idConver_path: path of file unsed for convert
#' @param thre_pst: minimum persent threshold for convert
#' @return converted matrix
conveEnsembl = function(expr,      # data.table format
                        genes_name,
                        species,
                        idConver_path,
                        thre_pst = 0.1){
  cat(file.path(idConver_path, paste0(species, '.csv')),"\n")
  ids = read.csv(file.path(idConver_path, paste0(species, '.csv')), row.names=1)
  #print(head(ids),quote=F)
  match_index = match(genes_name, rownames(ids))
  matched_genes_index = na.omit(match_index)
  if(length(matched_genes_index)> thre_pst*nrow(expr)){
    #print(dim(ids))
    ids = ids[matched_genes_index,]
    #print(length(ids),quote=F)
    expr = as_FBM(expr[!is.na(match_index)])
    #print(head(ids),quote=F)
    genes_name = as.character(ids)
  }else{
    expr = as_FBM(expr)
  }
  list(expr = expr, genes_name=genes_name)
}

#' get file type
#' @param filename: file path
#' @return file type
getFiletype = function(filename){
  filename_lis = tolower(strsplit(filename, split='[.]')[[1]])
  return(filename_lis[length(filename_lis)])
}



#' loading gene expression profile ane cell type file
#' @param expr_file: expression profile path
#' @param cell_type_path: cell type file path, default is `none`.
#' @param imputation: imputation method, `scImpute`, `magic` or `none`.
#' scImpute: scImpute(0.0.9), magic: Rmagic(2.0.3)
#' @param data_type: data type of expression profile. `count` or `tpm`.
#' @param out_dir: output direction.
#' Just useful if imputation method is `scImpute`.
#' @param para_size: number of cores in parallel.
#' Just useful if imputation method is `scImpute`.
#' @param normalize_method: normalization method.
#' `log`,`CLR`,`RC`,`scran`,`sctransform`,or `none`
#' @param min_cells: genes must be in a minimum number of cells. Default is 3
#' @param min_features: cells must have at least the minimum number of genes. Default is 200.
#' @param pas_method: method for calculating pathway activity score.
#' `gsva`,`ssgsea`,`zscore` or `plage`. Default is `gsva`
#' @param work_dir: workshop direction.
#' `workshop direction`/data/genome store gene length file and ID convert file
#' @return a List for expression matrix and cell type
load_expr = function(expr_file,
                     cell_type_path,
                     imputation,
                     data_type,
                     species,
                     out_dir,
                     para_size,
                     normalize_method,
                     min_cells,
                     min_features,
                     pas_method,
                     work_dir)
{
  # load gene expression matrix
  # file format:
  #            rownames: gene symbol;
  #            colnames: sample (cell) ID;
  #            sep: comma for '.csv' or tab or '.txt'


  data_dir = file.path(work_dir,'data/genome')

  ###################################################
  ######       read expression profile        #######
  ###################################################

  expr_path_type = getFiletype(expr_file)
  if(expr_path_type == 'csv'){
    #expr = read.csv(expr_file, stringsAsFactors = F,row.names=NULL,header=T,check.names=FALSE)
    expr = data.table::fread(expr_file,header=TRUE,check.names=FALSE)
    genes_name = as.character(as.data.frame(expr[,1])[,1])
    cells_name = colnames(expr)[2:ncol(expr)]
    if(sum(duplicated(genes_name) > 0)){
      #expr = bigstatsr::as_FBM(expr[!duplicated(genes_name), 2:ncol(expr)])
      expr = expr[!duplicated(genes_name), 2:ncol(expr)]
      genes_name = genes_name[!duplicated(genes_name)]
    }else{
      expr = expr[, 2:ncol(expr)]
      #expr = as_FBM(expr[, 2:ncol(expr)])
    }

  }else if(expr_path_type == 'txt'){
    #expr = read.table(expr_file,sep='\t',stringsAsFactors = F,row.names=NULL,header=T,check.names=FALSE)
    expr = data.table::fread(expr_file,header=TRUE,check.names=FALSE,sep='\t')
    genes_name = as.character(as.data.frame(expr[,1])[,1])
    cells_name = colnames(expr)[2:ncol(expr)]
    if(sum(duplicated(genes_name) > 0)){
      #expr = bigstatsr::as_FBM(expr[!duplicated(genes_name), 2:ncol(expr)])
      expr = expr[!duplicated(genes_name), 2:ncol(expr)]
      genes_name = genes_name[!duplicated(genes_name)]
    }else{
      expr = expr[, 2:ncol(expr)]
      #expr = as_FBM(expr[, 2:ncol(expr)])
    }

  }else if(expr_path_type == 'rds'){
    expr = readRDS(expr_file)
    expr = na.omit(expr)
    genes_name = rownames(expr)
    cells_name = colnames(expr)
  }else{
    stop("Format of file must be 'csv', 'txt', or 'rds'! Please check it!",quote=F)
  }

  genes_name = gsub(" +","",genes_name)
  cells_name = gsub(" +","",cells_name)
  cat(paste("Reading expr file success, dimensions of matrix are ",
            paste(dim(expr),collapse=',')),"\n")
  gc()
  #print(mem_used())
  #print("###################2######################")
  #print(genes_name[1:10])



  ###################################################
  ########       read cell type file        #########
  ###################################################

  if(! cell_type_path == 'NULL'){
    cat("Reading cell type file...\n")
    #print(cell_type_path)
    cellLabel_path_type = getFiletype(cell_type_path)
    if(cellLabel_path_type == 'csv'){
      cell_type = read.csv(cell_type_path, stringsAsFactors = F,row.names=1,header=F,check.names=FALSE)
    }else if(cellLabel_path_type == 'txt'){
      cell_type = read.table(cell_type_path,sep='\t',stringsAsFactors = F,row.names=1,header=F,check.names=FALSE)
    }
    rownames(cell_type) = gsub(" +","",rownames(cell_type))
    ### check cells' name
    if(length(intersect(cells_name, rownames(cell_type))) > 0){
      if(length(setdiff(cells_name, rownames(cell_type))) > 0){
        warning("Some cells in expression matrix are not included in first line of cell type file. Throwing out them ...")
        cells_name = intersect(cells_name,rownames(cell_type))
        expr = expr[, cells_name]
        cell_type = cell_type[cells_name,]
        cat("After throwing out cells not included in cell type file, dimensions of expression matrix are ",
            paste(dim(expr),collapse=','),"\n")
      }
    }else{
      stop("None of the gene names in the expression profile appear in the cell type file,Please check it!",quote=F)
    }

  }else{
    cat("Do not specific cell type file\n")
    cell_type = 'none'
  }

  ### convert expression matrix
  #expr = as.matrix(expr)
  #rownames(expr) = gsub(" +","",rownames(expr))

  gc()


  ###################################################
  ##############     filtering         ##############
  ###################################################

  ### filtering cells
  if (min_features > 0) {
    print("start filtering cells")
    nfeatures = colSums(expr)
    index = as.integer(which(x = nfeatures >= min_features))
    if(length(index) < ncol(expr)){
      expr = expr[, ..index]
      cells_name = cells_name[index]
    }
    rm(nfeatures, index, min_features)
    #nfeatures = bigstatsr::big_apply(x, a.FUN = function(X, ind) {colSums(X[, ind])}, a.combine = 'c')
    #nfeatures = bigstatsr::big_colstats(expr, ncores = 5)$sum
    #print("?????????????????////")
    #expr = as_FBM(expr[, which(x = nfeatures >= min_features)])
    cat("After filtering cells, dimensions of expression matrix are ",
        paste(dim(expr),collapse=','),'\n')
  }else{
    cat("Do not filter cells","\n")
  }

  if (min_cells > 0) {
    print("start filtering genes")
    num_cells = rowSums(x = expr)
    index = as.integer(which(x = num_cells >= min_cells))
    if(length(index) < nrow(expr)){
      expr = expr[index]
      genes_name = genes_name[index]
    }
    rm(num_cells, index, min_cells)
    #big_rowSums = function(X, ind){rowSums(X[ind,])}
    #num_cells = bigstatsr::big_parallelize(expr, p.FUN = big_rowSums, p.combine = 'c',ind = rows_along(expr), ncores = 5)
    #num_cells = Matrix::rowSums(x = expr > 0)
    #expr = as_FBM(expr[which(x = num_cells >= min_cells), ])
    #genes_name = genes_name[which(x = num_cells >= min_cells)]
    cat("After filtering genes, dimensions of expression matrix are ",
        paste(dim(expr),collapse=','),"\n")
  }else{
    cat("Do not filter genes","\n")
  }
  #print("###################3###################")
  #print(genes_name[1:10])
  gc()
  ### convert character matrix to numeric matrix
  #if(pas_method == 'gsva' && imputation != 'none' && class(expr[1,1]) == 'character'){
  #  warning("Character type for expression matrix is unuseful when pas_method is gsva. Casting data type to numeric...")
  #  genes = rownames(expr)
  #  expr = apply(expr,2,as.numeric)
  #  rownames(expr) = genes
  #}


  ## convert Ensembl ID to gene name
  ## convert files are download from Ensembl Gene Browser
  ## DATABASE: Ensembl Genes 99
  ## homo/Human Genes (GRCh38.p13)
  ## mus/Mouse Genes (GRCm38.p6)
  idConver_path = file.path(data_dir, 'ID_convert')
  ret = conveEnsembl(expr = expr,
                     genes_name = genes_name,
                     species = species,
                     idConver_path = idConver_path)
  expr = ret$expr
  genes_name = ret$genes_name
  rm(ret)
  gc()
  #print("######################4######################")
  #print(genes_name[1:10])

  #print(mem_used())

  ###################################################
  ##############     normalization     ##############
  ###################################################

  if(normalize_method == 'log'){
    expr = as_FBM(NormalizeData(expr[],normalization.method = 'LogNormalize',verbose=0))
  }else if(normalize_method == 'CLR'){
    expr = as_FBM(NormalizeData(expr[],normalization.method = 'CLR',verbose=0))
  }else if(normalize_method == 'RC'){
    expr = as_FBM(NormalizeData(expr[],normalization.method = 'RC',verbose=0))
  }else if(normalize_method == 'scran'){
    #genes = rownames(expr)
    #cells = colnames(expr)
    #expr = apply(expr,2,function(x) {storage.mode(x) <- 'integer'; x})
    sc = SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = expr[])
    )
    rm(expr)
    gc()
    clusters = scran::quickCluster(sc)
    gc()
    sc = scran::computeSumFactors(sc, clusters=clusters)
    rm(cluster)
    gc()
    sc = scater::normalize(sc)
    expr = as_FBM(sc@assays$data$logcounts)
    rm(sc)
    gc()
    #rownames(expr) = genes
    #colnames(expr) = cells
  }else if(normalize_method == 'sctransform'){
    expr_ob = Seurat::CreateSeuratObject(counts=expr[])
    rm(expr)
    gc()
    expr_ob = Seurat::SCTransform(expr_ob,verbose=FALSE)
    expr = as_FBM(expr_ob@assays$SCT@data)
  }else{
    cat("Do not procedure normalization","\n")
  }
  #print("after normalization")
  #print(mem_used())



  ###################################################
  ###############      imputation     ##############
  ###################################################

  if(imputation == 'scImpute'){
    library(scImpute)
    geneLength_path = file.path(data_dir,'all_symbol_length.txt')
    gene_len = read.table(geneLength_path,stringsAsFactors = F,sep="\t",header=F,row.names=1)
    match_tmp = match(rownames(gene_len),genes_name)
    if (length(na.omit(tmp)) < nrow(expr)){
      warning("check the length file, some genes in expression profile do not have length. Throwing them ...")
      cat(length(setdiff(genes_name,rownames(gene_len))),"\n")
      expr = expr[!is.na(match_tmp),]
      genes_name = genes_name[!is.na(match_tmp)]
      cat("After throwing out cells not included in cell type file, dimensions of expression matrix are ",
          paste(dim(expr),collapse=','),"\n")
      gene_len = gene_len[na.omit(match_tmp),]
      expr = expr[]
      rownames(expr) = genes_name
      colnames(expr) = cells_name
      write.csv(expr,file = file.path(out_dir, 'expression.csv'),quote=F)
      expr_file = file.path(out_dir, 'expression.csv')
    }else{
      gene_len = gene_len[genes_name, ]
      write.csv(expr,file = file.path(out_dir, 'expression.csv'),quote=F)
      expr_file = file.path(out_dir, 'expression.csv')
    }

    if(cell_type == 'none'){
      cat("Start imputing missing value of expression profile without specific cell type using scImpute 0.0.9 ...","\n")
      scimpute(count_path = expr_file,
               type = data_type,
               infile = "csv",
               outfile = "csv",
               genelen = gene_len,
               out_dir = out_dir,
               drop_thre = 0.5,
               labeled = F,
               Kcluster = 8,
               ncores = para_size)
      expr = as_FBM(read.csv(file.path(out_dir,'scimpute_count.csv'), stringsAsFactors = F,row.names=1))
    }else{
      cat("Start imputing missing value of expression profile with cell type using scImpute 0.0.9 ...","\n")
      scimpute(count_path = expr_file,
               type = data_type,
               infile = "csv",
               outfile = "csv",
               genelen = gene_len,
               out_dir = out_dir,
               drop_thre = 0.5,
               labeled = F,
               Kcluster = 8,
               ncores = para_size)
      expr = as_FBM(read.csv(file.path(out_dir,'scimpute_count.csv'), stringsAsFactors = F,row.names=1))
    }
  }else if(imputation == 'magic'){
    library(Rmagic)
    cat("Start imputing missing value of expression profile using Rmagic 2.0.3 ...","\n")
    expr = as_FBM(magic(expr[])$results)
  }else{
    cat("Do not impute missing value of expression profile.","\n")
  }


  if(length(cell_type != 'none') > 0){
    cell_type = cell_type[cells_name,]      ## convert a dataframe to vector
  }
  #print("#####################5###############")
  #print(genes_name[1:10])
  #print(dim(expr))
  #print(class(expr))
  #print(length(genes_name))
  #rint(class(genes_name))
  #print(length(cells_name))
  #print(cell_type[1:5])
  return (list(expr = expr,cell_type=cell_type,genes_name=genes_name,cells_name=cells_name))



}


load_gmt_data = function(gmt_file_path){
  # load pathway gene sets' gmt file
  # file format:
  #            first index: pathway's name/ID
  #            second index: pathway's url or others, it dosen't matter
  #            third to all: gene symbols in pathway
  #            sep: \t
  # return a list
  tmp = readLines(gmt_file_path)
  gsets = list()
  for(i in 1:length(tmp)){
    t = strsplit(tmp[i],'\t')[[1]]
    genes = t[3:length(t)]
    genes = genes[which(genes != "")]
    gsets[[t[1]]] = genes
  }
  return (gsets)
}




### deserted by zyr on 2020.04.13
calcu_PAS = function(expr,
                     dataSets,
                     method='ssgsea',
                     topo=TRUE){
  # calculate pathway activity score for every gene set
  # expr: gene expression matrix throught normalization
  # dataSets: list format, genes in pathways
  # method: method used for calculate pathway activity score for every sample.
  #         By default is sssgsea and other option are zscore, gsva, plage.
  # topo: use topological information as weight. Default: TRUE.
  if(topo){
    weight = NA
  }else{
    weight = NA
  }
  score = PAS(as.matrix(expr), dataSets, method = method, weight = weight)
  score
}

#' clustering cells
#' @param mat: matrix for clustering procedure. genes X cells
#' @param method: clustering method.
#' `seurat`, `hclust`, `dbscan`, `simlr`, `kmeans` or `kmedoids`
#' @param resolution: for Seurat FindClusters
#' @param dimensions: for Seurat FindNeighbors
#' @param k_cluster: number of clusters
#' Just useful if clustering method is `kmeans`, `kmedoids`, `simlr` or `hclust`
#' @param min_pts: number of minimum points in the eps region
#' Just useful if clustering method is `dbscan`
#' @return a vector of clustering type
cal_cluster = function(mat,
                       method='seurat',
                       seurat_dims,
                       seurat_resolution,
                       k_cluster,
                       min_pts){
  ##return a character
  res_cluster = clustering(mat,
                           method=method,
                           seurat_dims = seurat_dims,
                           seurat_resolution = seurat_resolution,
                           k_cluster,
                           min_pts)
  return(res_cluster)
}



#' preparation seurat object
#' writing tsne-2D.csv, tsne-3D.csv, umap-2D.csv and umap-3D.csv
#' @param mat: gene expression profile
#' @param cluster_idents: clustering idents in factor format.
#' @param dims: dimensions for Seurat RunTSNE or RunUMAP
#' @param out_dir: output direction
#' @return a seurat object containing cell type idents
prepa_seuratOb = function(mat,
                          cluster_idents,
                          dims,
                          out_dir){
  ### define number of PCs and number of PCA dimensions used for TSNE/UMAP
  number_pc = min(50, ncol(mat)-5)
  cat("number of pcs: ",number_pc,"\n")
  dims = min(dims, number_pc-2)
  cell_name_frame = data.frame(colnames(mat))
  colnames(cell_name_frame) = "cellName"
  ###### creat seurat object
  obSeurat = CreateSeuratObject(counts=mat,min.cells = 0,min.features = 0)
  obSeurat = FindVariableFeatures(obSeurat, verbose = FALSE)
  obSeurat = ScaleData(obSeurat, verbose = FALSE)
  obSeurat = RunPCA(obSeurat,verbose=F,seed.use = 42,npcs = number_pc)
  ###### locate idents
  Idents(obSeurat) = cluster_idents

  if((ncol(mat)-2)%/%3 < 1){
    cat("The number of cells is too small\n")
  }
  perplexity = min(30, (ncol(mat)-2)%/%3)

  obSeurat = RunTSNE(obSeurat,dims=1:dims,seed.use = 42, dim.embed = 3, perplexity = perplexity)
  tsne_seurat = as.data.frame(Embeddings(obSeurat[['tsne']]))
  tsne_seurat[,'cell_type'] = Idents(obSeurat)
  #print(dim(cell_name_frame))
  #print(dim(tsne_seurat))
  tsne_seurat = cbind(cell_name_frame, tsne_seurat)
  #rownames(tsne_seurat) = tsne_seurat[,1]
  write.csv(format(tsne_seurat, digits = 3), file.path(out_dir,'tsne_3D.csv'),quote=F, row.names=F)

  ###### write 2D-tsne
  obSeurat = RunTSNE(obSeurat,dims=1:dims,seed.use = 42, perplexity = perplexity)
  tsne_seurat = as.data.frame(Embeddings(obSeurat[['tsne']]))
  tsne_seurat[,'cell_type'] = Idents(obSeurat)
  tsne_seurat = cbind(cell_name_frame, tsne_seurat)
  #rownames(tsne_seurat) = tsne_seurat[,1]
  write.csv(format(tsne_seurat, digits = 3), file.path(out_dir,'tsne_2D.csv'),quote=F, row.names=F)
  rm(tsne_seurat)
  gc()

  ###### write 3D-umap
  obSeurat = RunUMAP(obSeurat,dims=1:dims,seed.use = 42, umap.method = "umap-learn",metric = "correlation", n.components=3,verbose=FALSE)
  umap_seurat = as.data.frame(Embeddings(obSeurat[['umap']]))
  umap_seurat[,'cell_type'] = Idents(obSeurat)
  umap_seurat = cbind(cell_name_frame, umap_seurat)
  #rownames(umap_seurat) = umap_seurat[,1]
  write.csv(format(umap_seurat, digits = 3), file.path(out_dir,'umap_3D.csv'), quote=F, row.names=F)

  ###### write 2D-umap
  obSeurat = RunUMAP(obSeurat,dims=1:dims,seed.use = 42, umap.method = "umap-learn",metric = "correlation",verbose=FALSE)
  umap_seurat = as.data.frame(Embeddings(obSeurat[['umap']]))
  umap_seurat[,'cell_type'] = Idents(obSeurat)
  umap_seurat = cbind(cell_name_frame, umap_seurat)
  #rownames(umap_seurat) = umap_seurat[,1]
  write.csv(format(umap_seurat, digits = 3), file.path(out_dir,'umap_2D.csv'), quote=F, row.names=F)
  ##### write cell_type
  write.table(data.frame(Idents(obSeurat)),
              file.path(out_dir,'cell_type.csv'),sep=',',
              quote = F, col.names = F)
  return(obSeurat)
}


#' deserted by zyr
plot_cluster = function(obSeurat,
                        method = c('tsne','umap'),
                        out_dir,
                        pic_type
){
  #print(out_dir)
  #print(method)
  options(bitmapType='cairo')
  method = match.arg(method)
  #if(method == 'tsne'){
    #print(out_dir)
    #print(dim(obSeurat))
    p = DimPlot(object = obSeurat, group.by = 'ident', reduction='tsne',label=T)
    if(pic_type == 'png'){
      png(file = file.path(out_dir,'tsne.png'))
      print(p)
      dev.off()
    }else if(pic_type == 'pdf'){
      pdf(file = file.path(out_dir, 'tsne.pdf'))
      print(p)
      dev.off()
    }else{
      stop("--pic_type must be png or pdf")
    }

  #}else if(method =='umap'){

    p = DimPlot(object = obSeurat, group.by = 'ident', reduction='umap', label=T)
    if(pic_type == 'png'){
      png(file = file.path(out_dir,'umap.png'))
      print(p)
      dev.off()
    }else if(pic_type == 'pdf'){
      pdf(file = file.path(out_dir, 'umap.pdf'))
      print(p)
      dev.off()
    }else{
      stop("--pic_type must be png or pdf")
    }
  #}
}


# deserted by zyr
findMarker = function(obSeurat,
                      method
)
{
  markers = FindAllMarkers(obSeurat, test.use=method,
                           min.cells.feature = 3,
                           min.cells.group = 3)
  # return all markers, include adjusted p-Value > 0.05
  return(markers)
}



#' Finding all markers
#' Ploting all markers,including scatter, violin and heatmap plot.
#'
#'
#'
#' @return a heatmap plot and a series of folder
#' Scene of significantly activity pathway: Heatmap.pdf/Heatmap.png
#' Creating folder for top 50 markers.
#' Every marker folder contain three picture.
#' marker_name.fet.pdf(png): Highlight scatter plot
#' marker_name.vil.pdf(png): violin plot
#' marker_name.het.pdf(png): heatmap plot
parseAllMarkers = function(obSeurat,
                           user_pathways,
                           method,
                           logfc_thre,
                           min.pct,
                           out_dir,
                           para_size,
                           expr,
						   genes_name,
						   cells_name,
                           cell_type,
                           path_list,
                           pic_type,
                           cells

){
  all_markers = FindAllMarkers_me(obSeurat,
                                  test.use = method,
                                  logfc.threshold = logfc_thre,
                                  min.pct = min.pct,
                                  only.pos = T,
                                  min.cells.feature = 3,
                                  min.cells.group = 3,
                                  para_size = para_size)
  #print("from function::::::")
  #print(dim(all_markers))
  if(nrow(all_markers) == 0){
    stop("No differential pathways are identified.")
    #return(data.frame())
  }

  ## subset pathways
  #all_markers = subset(all_markers, select=-c(pct.1,pct.2))
  #print("after subset:")
  #print(dim(all_markers))

  ###### write excel
  path_new_names = gsub('-','_',all_markers$gene)
  pathway_frame = as.data.frame(path_new_names)
  colnames(pathway_frame) = "pathways"
  all_markers_write = cbind(pathway_frame, all_markers[,1:(ncol(all_markers)-1)])
  genes = get_geneLis(paths = path_new_names,
                      path_lis = path_list,
                      expr = genes_name)
  print(genes[1])
  #print(length(genes))
  all_markers_write['geneList'] = genes
  #print(dim(all_markers_write))

  write.csv(format(all_markers_write, digits = 3), file.path(out_dir,'pas_markers.csv'),quote=F,row.names=F)
  cat("write pas_markers success\n")
  #print(head(all_markers))
  n_top = 10
  tops = all_markers %>% group_by(cluster) %>% top_n(n = n_top, wt = avg_logFC)
  #print("groupby success")

  #while( nrow(tops) > 50 ){
  #  n_top = n_top - 1
  #  tops = all_markers %>% group_by(cluster) %>% top_n(n = n_top, wt = avg_logFC)
  #}

  #print(dim(tops))
  p_heatmap = DoHeatmap_me(obSeurat,
                           features = unique(tops$gene),
                           size=6, angle = 90, group.bar.height = 0.01,
                           label=FALSE)
  #png(file = file.path(out_dir, 'Heatmap.png'),width = 1200, height = 800)
  #print(p_heatmap)
  #dev.off()

  if(pic_type == 'png'){
    png(file.path(out_dir, 'Heatmap.png'))
    print(p_heatmap)
    dev.off()
  }else if(pic_type == 'pdf'){
    pdf(file.path(out_dir, 'Heatmap.pdf'))
    print(p_heatmap)
    dev.off()
  }else{
    stop("--pic_type must be png or pdf")
  }

  rm(p_heatmap)
  gc()

  #p_bubble = bubble_plot(obSeurat,as.data.frame(tops))
  #png(file = file.path(out_dir, 'Bubble.png'),width = 1200, height = 800)
  #print(p_bubble)
  #dev.off()

  #pdf(file.path(out_dir, 'Bubble.pdf'),width = 16,height = 10)
  #print(p_bubble)
  #dev.off()

  cat("plot success\n")


  #if(!is.null(user_pathways)){
  #  feature_pathways = unique(tops$gene, user_pathways)
  #}else{
  #  feature_pathways = unique(tops$gene)
  #}


  #expr = expr[,colnames(obSeurat)]
  #print(head(cell_type))
  expr = scale(expr[,cells])
  rownames(expr) = genes_name
  colnames(expr) = cells_name[cells]
  cell_type = as.data.frame(Idents(obSeurat)[cells])
  #print("after seurat ...")
  #print(head(cell_type))
  #print(dim(cell_type))
  #print(dim(expr))
  gc()
  print("zscore of expression success, dimensions: ")
  print(dim(expr))
  print(expr[1:5,1:5])

  x_angle=0



  mclapp = get('mclapply', envir = getNamespace('parallel'))
  options(mc.cores = para_size)



  feature_pathways = unique(tops$gene)

  mclapp(1:length(feature_pathways),function(i){
    feature_name = gsub('-','_',feature_pathways[i])
    path_dir = file.path(out_dir, feature_name)
    dir.create(path_dir)
    #print(feature_name)
    #pic_name_vil = paste0(feature_name,'.vil.png')
    #pic_name_fet = paste0(feature_name,'.fet.png')
    #pic_name_het = paste0(feature_name,'.het.png')
    #mat_name_het = paste0(feature_name,'.csv')
    feature_vio(obSeurat,
                feature_pathways[i],
                path_dir,
                feature_name,
                x_angle,
                pic_type)
    feature_highlight(obSeurat,
                      feature_pathways[i],
                      path_dir,
                      feature_name,
                      pic_type)
    #print(expr[1:5,1:5])
    #print(path_list[[feature_name]])
	#print(cell_type)
	#print(feature_name)
    genes_heatmap(expr,
                  path_list[[feature_name]],
                  cell_type,
                  path_dir,
                  feature_name,
                  pic_type)
  })
  #for(i in 1:length(feature_pathways)){
  #print(feature_pathways[i])
  #print(names(path_list)[1:30])
  #print(feature_pathways[i])
  #print(length(path_list[[feature_pathways[i]]]))
  # feature_name = gsub('-','_',feature_pathways[i])
  #feature_name = gsub(' ','_',feature_pathways[i])
  #path_dir = file.path(out_dir, feature_name)
  #dir.create(path_dir)
  #pic_name_vil = paste0(feature_name,'.vil.png')
  #pic_name_fet = paste0(feature_name,'.fet.png')
  #pic_name_het = paste0(feature_name,'.het.png')
  #mat_name_het = paste0(feature_name,'.csv')


  #print(expr[1:3,1:3])
  #print(path_list[[feature_name]][1:5])
  # feature_vio(obSeurat, feature_pathways[i], file.path(path_dir, pic_name_vil))
  # feature_highlight(obSeurat, feature_pathways[i], file.path(path_dir, pic_name_fet))
  # genes_heatmap(expr, path_list[[feature_name]], cell_type, file.path(path_dir, pic_name_het), file.path(path_dir, mat_name_het))
  #}

  #feature_frame = as.data.frame(feature_pathways)
  #write.table(feature_frame, file = file.path(out_dir, 'pathway_names.csv'), col.names=F, row.names=F, quote=F)
  #return(all_markers)
}

# deserted
plot_feature = function(obSeurat, feature, reduction){
  FeaturePlot(obSeurat, features = feature, reduction = reduction)
}


#load_user_pathway()


#' overview of highly variable pathway activity score
plotPAS_heatmap = function(mat,
                           cluster_idents,
                           out_path,
                           pic_type,
                           max_path = 500){
  #print(cluster_idents[1:5])
  colors = hue_pal()(length(unique(cluster_idents)))
  #print(length(cluster_idents))
  #print(colors)
  names(colors) = levels(cluster_idents)
  #print(colors)
  #print(file.path(out_path,'pas_color.csv'))
  write.csv(as.data.frame(colors),file = file.path(out_path,'pas_color.csv'),quote=F)
  cat("Writing color success\n")
  #names(cluster_idents) = colnames(mat)
  cluster_idents = as.data.frame(cluster_idents)
  #print(head(cluster_idents))
  colnames(cluster_idents) = "cell types"
  #mat = mat[,rownames(cluster_idents)]

  ####filter pathways
  if(dim(mat)[1] > max_path){
    #pathway_sd = apply(mat,1, sd)
    #mat = mat[names(sort(pathway_sd,decreasing = T)[0:max_path]),]
    # edited by zyr on 2020.04.07
    mat_obj = CreateAssayObject(mat)
    mat_obj = FindVariableFeatures(mat_obj, nfeatures=mat_path)
    mat = mat[mat_obj@var.features,]
  }


  p = pheatmap(mat,
               color = colorRampPalette(c("forestgreen", "white", "orange"))(100),
               #color = "PiYG",
               show_rownames = F,
               show_colnames = F,
               fontsize = 20,
               annotation_names_col = F,
               annotation_colors = list("cell types" = colors),
               cluster_cols = F,
               annotation_col = cluster_idents)

  if(pic_type == 'png'){
    png(file.path(out_path, 'pas_heatmap.png'))
    print(p)
    dev.off()
  }else if( pic_type == 'pdf'){
    pdf(file.path(out_path, 'pas_heatmap.pdf'))
    print(p)
    dev.off()
  }else{
    stop("--pic_type must be png or pdf")
  }

}




