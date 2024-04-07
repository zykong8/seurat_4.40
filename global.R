

UpdateObj <- function(pcas, res, obj){

  if (pcas == 15 & res == 0.6){
    return(obj)
  }else{
    
    if (pcas == 15 & res != 0.6){
      obj <- FindClusters(obj, resolution = res)
    }else{
      obj <- FindNeighbors(obj, dims = 1:pcas)
      obj <- FindClusters(obj, resolution = res)
      obj <- RunUMAP(obj, dims = 1:pcas)
      obj <- RunTSNE(obj, dims = 1:pcas)
    }
    
    return(obj)
    
  }
}

