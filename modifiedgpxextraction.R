## modified version of https://github.com/cran/plotKML/blob/master/R/readGPX.R
track.readGPX.element <- function(file) {
  # element = "metadata", "wpt", "rte", "trk"
  
  ret <- xmlTreeParse(file, useInternalNodes = TRUE)
  # top structure: 
  top <- xmlRoot(ret)
  
  # check if there is any content:
    # tracks:
      ret <- NULL
      nu <- which(names(top) %in% element)
      for(c in seq_along(nu)){
        lst <- which(names(top[[nu[c]]]) %in% "trkseg")
        nm <- names(top[[nu[c]]][[lst[1]]][[1]])
        ret[[c]] <- list(NULL)
        for(i in seq_along(lst)) {
          trkpt <- top[[nu[c]]][[lst[i]]]
          ret[[c]][[i]] <- data.frame(NULL)
          ## get columns (http://www.topografix.com/GPX/1/1/#type_wptType)
          lon <- as.numeric(unlist(xmlSApply(trkpt, xmlGetAttr, "lon")))
          lat <- as.numeric(unlist(xmlSApply(trkpt, xmlGetAttr, "lat")))
          ret[[c]][[i]][1:length(lon),"lon"] <- lon
          ret[[c]][[i]][1:length(lat),"lat"] <- lat

          if(!length(nm)==0){
            for(j in 1:length(nm)){
              xm <- as.character(sapply(sapply(xmlChildren(trkpt), function(x) x[[nm[[j]]]]), xmlValue))
              ret[[c]][[i]][1:length(xm), nm[[j]]] <- xm 
            }
          } 
        }
        names(ret[[c]]) <- xmlValue(top[[nu[c]]][["name"]])
      }   
  
  return(ret)
}

