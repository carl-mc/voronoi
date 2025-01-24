# SAMPLE VORONOI CELLS WITHIN AREA #####

# (function) Sample regular voronoi cells within entries of SpatialPolygonsDataFrame
sample_vorcells <- function(spdf, res, size, size.km2 = F, seed = 1, iter.max = 100, sample.type = "nonaligned",
                            ncore = 1){
  require(rgeos)
  require(raster)
  if(ncore == 1){
    # Sampe cells per entry in spdf
    result <- do.call(rbind,
                      lapply(c(1:length(spdf)), function(i){
                        print(i)
                        t <- kmean_spat(shape = spdf[i,], res = res, size = size, 
                                        size.km2 = size.km2, seed = seed, iter.max = iter.max, sample.type = sample.type)
                        t$unit.id <- i
                        sp::proj4string(t) <- sp::proj4string(spdf)
                        row.names(t) <- paste0(i, ".", row.names(t))
                        return(t)
                      }))
  } else if(ncore > 1) {
    
    # Libraries
    require(foreach)
    require(doParallel)
    
    # Setup cluster
    try(registerDoSEQ())
    cl <- makeCluster(getOption("cl.cores", ncore), outfile = "error.txt")
    clusterExport(cl, list("voronoi_poly", "kmean_spat"), envir = .GlobalEnv)
    registerDoParallel(cl)
    
    # Compute tesselation
    result <- foreach(i = c(1:length(spdf)), .packages = c("deldir","raster","rgeos","geosphere","sp"),
                      .options.multicore = list(preschedule = FALSE)) %dopar% {
                        t <- kmean_spat(shape = spdf[i,], res = res, size = size, 
                                        size.km2 = size.km2, seed = seed, iter.max = iter.max, sample.type = sample.type)
                        t$unit.id <- i
                        proj4string(t) <- proj4string(spdf)
                        row.names(t) <- paste0(i, ".", row.names(t))
                        t
                      }
    result <- do.call(rbind, result)
    
    # Close cluster
    stopCluster(cl)
    rm(cl)
    
  } else {
    stop("Please enter a valid number of cores > 0")
  }
  
  # Clean result
  result <- rgeos::gBuffer(result, byid = T, width = 0)
  
  # Return
  return(result)
}

# (function) Spatial Voronoi Polygons
voronoi_poly <- function(vpoints, extent = NULL){
  require(deldir)
  require(sp)
  
  # ... tesselation
  z = deldir(data.frame(vpoints@coords), rw = as.vector(extent))
  w = tile.list(z)
  
  # ... make polygons
  polys = vector(mode='list', length=length(w))
  for (i in seq(along=polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(w[[i]]$ptNum))
  }
  return(SpatialPolygonsDataFrame(SpatialPolygons(polys),
                                  data = data.frame(cell.id = z$ind.orig), 
                                  match.ID = F))
}

# (function) Use kmeans-clustering for optimized voronoi tesselation of a delimited plane
kmean_spat <- function(shape, res, size, size.km2 = F,
                       seed = 1, iter.max = 100, sample.type = "nonaligned"){
  # Packages
  require(raster)
  require(sp)
  require(deldir)
  require(geosphere)
  require(rgeos)
  set.seed(seed)
  
  # ... make coords
  # ... --- raster
  base.r <- raster::raster(extent(shape), res = res)
  shape$is.shape <- 1
  base.r <- raster::rasterize(shape, base.r, field = "is.shape", background = NA)
  
  # ... --- points
  points <- na.omit(as.data.frame(base.r, xy = T))
  coords <- points[,1:2]
  
  # ... how many clusters?
  if(size.km2){
    require(geosphere)
    N <- max(c(1,round((areaPolygon(shape)/1000^2) / (size))))
  } else {
    N <- max(c(1,round(gArea(shape) / (size))))
  }
  
  # ... return if N == 1
  if(N == 1 | nrow(coords) == 0){
    tess <- shape
    tess@data <- data.frame(cell.id = 1)
    return(tess)
  }
  
  # ... if N > points
  if(N >= nrow(coords)){
    # ... tesselate
    tess <- try(voronoi_poly(SpatialPoints(coords),
                             extent = extent(shape)))
    if(class(tess) == "try-error"){print(paste("ERROR IN VORONOI TESSELATION"))}
    
  } else {
    # ... regular sampling of centers
    set.seed(seed)
    m <- 0
    while(T & m < iter.max){
      m = m + 1
      centers <- sp::spsample(shape, n = N, type = sample.type)@coords
      if(nrow(centers) == N){
        # ... kmeans of coords
        # set.seed(seed)
        km <- stats::kmeans(x = coords, centers = centers, algorithm = "Lloyd", iter.max = iter.max)
        
        # ... check result
        if(length(unique(km$cluster)) == N){
          break
        }
      }
    }
    
    
    # ... tesselate
    tess <- try(voronoi_poly(SpatialPoints(na.omit(km$centers)),
                             extent = extent(shape)))
    if(class(tess) == "try-error"){
      print(paste("ERROR IN VORONOI TESSELATION"))
      return(NULL)
    }
  }
  
  # ... crop
  if(!is.null(shape)){
    if(length(shape) > 1){
      shape <- gUnaryUnion(shape, id = rep(1, length(shape)))
    }
    proj4string(shape) <- proj4string(tess)
    tess <- SpatialPolygonsDataFrame(gIntersection(tess, shape, byid = TRUE), 
                                     tess@data, match.ID = FALSE) # raster::crop(tess, shape)
  }
  
  # Return
  return(tess)
}


# (function) Node all admin shapes, to get constitutive parts that are never split
get_poly_const_parts <- function(poly, size.cutoff = .001, base.shape = NULL){
  require(rgeos)
  require(sp)
  
  # Base sspacefillr# Base shape
  if(is.null(base.shape)){
    base.shape <- rgeos::gUnaryUnion(poly, id = rep(1, length(poly)))
  }
  
  # # Delete small holes
  # base.shape <- remove.holes(base.shape, max.size = 1)
  
  # # Intersect with all shapes
  # poly.sup <- poly[gIntersects(poly,base.shape, byid = T)[1,],]
  # poly.sup <- raster::crop(poly.sup, base.shape)
  
  # Delete below size cutoff
  poly.sup <- poly[rgeos::gArea(poly, byid = T) > size.cutoff,]
  
  # As lines (with buffer)
  poly.lines <- as(poly.sup, "SpatialLines")
  poly.lines <- rgeos::gBuffer(poly.lines, width = size.cutoff^2)
  poly.lines <- rgeos::gUnaryUnion(poly.lines, 
                           id = rep(1, length(poly.lines)), 
                           checkValidity = 2L)
  
  # Split base polygon
  base.split <- rgeos::gDifference(base.shape, poly.lines)     
  
  # Get holes
  holes <- lapply(seq_along(base.split@polygons[[1]]@Polygons), function(p){
    if(base.split@polygons[[1]]@Polygons[[p]]@hole){
      sp::Polygons(list(base.split@polygons[[1]]@Polygons[[p]]), ID = as.character(p))
    } else {
      NULL
    }
  })
  holes[sapply(holes, is.null)] <- NULL
  holes <- sp::SpatialPolygons(holes)
  
  # Delete below size cutoff
  holes <- holes[gArea(holes, byid = T) > size.cutoff,]
  
  # Valid polygons
  base.split <- lapply(1:length(base.split@polygons[[1]]@Polygons), function(p){
    if(base.split@polygons[[1]]@Polygons[[p]]@hole){
      NULL
    } else {
      sp::Polygons(list(base.split@polygons[[1]]@Polygons[[p]]), ID = as.character(p))
    }
  })
  base.split[sapply(base.split, is.null)] <- NULL
  base.split <- sp::SpatialPolygons(base.split)
  
  # Delete below size cutoff
  base.split <- base.split[rgeos::gArea(base.split, byid = T) > size.cutoff,]
  
  # Correct holes
  if(length(holes) > 0){
    save.rn <- row.names(base.split)
    for(h in seq_along(holes)){
      this.h <- holes[h, ]
      this.p.id <- try(which(rgeos::gContains(base.split, this.h, byid = T)), silent = T)
      if(class(this.p.id) == "try-error" | length(this.p.id) == 0 | length(this.p.id) > 1){
        next
      } else {
        base.split@polygons[[this.p.id]] <- rgeos::gDifference(base.split[this.p.id,], this.h)@polygons[[1]] 
      }
    }
    row.names(base.split) <- save.rn
  }
  
  # Return
  final <- sp::SpatialPolygonsDataFrame(base.split, 
                                    data.frame(const.id = 1:length(base.split)),
                                    match.ID = F)
  return(final)
}


#Remove holes in polygon function
remove.holes <-  function(Poly, max.size = NULL){
  require(sp)
  
  # get holes
  holes <- lapply(Poly@polygons, function(x) sapply(slot(x, "Polygons"), slot,
                                                    "hole"))
  
  # subset by size
  size <- lapply(Poly@polygons, function(x) sapply(slot(x, "Polygons"), slot,
                                                   "area"))
  
  # Delete
  if(is.null(max.size)){
    res <- lapply(1:length(Poly@polygons), function(i) slot(Poly@polygons[[i]],
                                                            "Polygons")[!holes[[i]]])
  } else {
    res <- lapply(1:length(Poly@polygons), function(i) slot(Poly@polygons[[i]],
                                                            "Polygons")[!holes[[i]] | size[[i]] > max.size ])
  }
  
  # Finish
  IDs <- row.names(Poly)
  Polyfill <- SpatialPolygons(lapply(1:length(res), function(i)
    Polygons(res[[i]], ID=IDs[i])), proj4string=CRS(proj4string(Poly)))
  return(Polyfill)
}

