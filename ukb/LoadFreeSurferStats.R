
LoadFreeSurferStats <- function(fsdir, 
                                sublistdir,
                                missing_gwr){
  
  # This function is used to directly read in FreeSurfer .stats files to R and format them for ComBat harmonization. 
  # 
  # input:
  #   - fsdir:
  #     character string of the directory where FreeSurfer output files are stored. This directory should contain
  #     the full sample of sub-* directories that contain the stats/ directories.
  #
  #   - sublistdir:
  #     character string of the directory where a participant ID list can be found. Participant ID list should
  #     be saved as 'sublist.csv'.
  #
  #   - missing_gwr:
  #     TRUE/FALSE as to whether gray-white signal intensity ratio phenotypes are unavilable. Default is FALSE.
  #
  #
  # output:
  #   - The output of this function will be a data.frame with the output of the entire sample's .stats files. Each row
  #     will have a unique participant and each column will have a different morphometric estimate.
  #
  #
  # example:
  # 
  # LoadFreeSurferStats(fsdir = '/Users/ew198/Documents/brainpace/data/dns/freesurfer_stats/',
  #                     sublistdir = '/Users/ew198/Documents/brainpace/data/dns/')
  #
  #
  # Written by Ethan Whitman (ethan.whitman@duke.edu)
  
  # setting default settings
  if (missing(missing_gwr)){
    missing_gwr <- FALSE
  }

  library(ggplot2)
  library(gridExtra)
  library(ggseg)
  
  # find first participant label
  
# profvis({
    
    sublist <- read.csv(paste0(sublistdir, 'sublist.csv'))$ID
    example_sub <- sublist[1]
    
    # load and format data from freesurfer to set up empty object for everyone
    
    aseg <- read_freesurfer_stats(paste0(fsdir, example_sub, '/stats/aseg.stats'))
    lh_aparc <- read_freesurfer_stats(paste0(fsdir, example_sub, '/stats/lh.aparc.stats'))
    rh_aparc <- read_freesurfer_stats(paste0(fsdir, example_sub, '/stats/rh.aparc.stats'))
    if (missing_gwr == FALSE){
      lh_wg <- read_freesurfer_stats(paste0(fsdir, example_sub, '/stats/lh.w-g.pct.stats'))
      rh_wg <- read_freesurfer_stats(paste0(fsdir, example_sub, '/stats/rh.w-g.pct.stats'))
    }

    # label hemispheres
    
    lh_aparc$label <- paste0(lh_aparc$label, '_left')
    rh_aparc$label <- paste0(rh_aparc$label, '_right')
    if (missing_gwr == FALSE){
      lh_wg$label <- paste0(lh_wg$label, '_left')
      rh_wg$label <- paste0(rh_wg$label, '_right')
    }

    # label phenotypes 
    aseg_temp <- t(aseg[,c('label', 'Volume_mm3')])
    aseg_temp[1,] <- gsub("^(Right|Left)-(.*)$", "\\2_\\1", aseg_temp[1,])
    aseg_temp[1,] <- gsub("3rd-Ventricle", "X3rd.Ventricle", aseg_temp[1,])
    aseg_temp[1,] <- gsub("4th-Ventricle", "X4th.Ventricle", aseg_temp[1,])
    aseg_temp[1,] <- gsub("5th-Ventricle", "X5th.Ventricle", aseg_temp[1,])
    
    surfarea_temp <- cbind(t(lh_aparc[,c('label', 'SurfArea')]),
                           t(rh_aparc[,c('label', 'SurfArea')]))
    surfarea_temp[1,] <- paste0('SA_',surfarea_temp[1,])
    
    thickavg_temp <-  cbind(t(lh_aparc[,c('label', 'ThickAvg')]),
                            t(rh_aparc[,c('label', 'ThickAvg')]))
    thickavg_temp[1,] <- paste0('CT_',thickavg_temp[1,])
    
    grayvol_temp <-  cbind(t(lh_aparc[,c('label', 'GrayVol')]),
                           t(rh_aparc[,c('label', 'GrayVol')]))
    grayvol_temp[1,] <- paste0('GMV_',grayvol_temp[1,])
    
    if (missing_gwr == FALSE){
      wg_temp <-  cbind(t(lh_wg[,c('label', 'Mean')]),
                      t(rh_wg[,c('label', 'Mean')]))
      wg_temp[1,] <- paste0('GWR_',wg_temp[1,])
    }

    # combine into one object
    if (missing_gwr == FALSE){
      labels <- cbind('ID', aseg_temp, surfarea_temp, thickavg_temp, grayvol_temp, wg_temp)[1,]
    } else if (missing_gwr == TRUE){
      labels <- cbind('ID', aseg_temp, surfarea_temp, thickavg_temp, grayvol_temp)[1,]
    }
    
    # data <- data.frame(matrix(ncol=length(labels), nrow = 0))
    data <- data.frame(matrix(ncol=length(labels), nrow = length(sublist)))
    colnames(data) <- labels
    
    # load everyone at once
    x <- 1
    asegs <- vector(mode = 'list', length = length(sublist))
    lh_aparcs <- vector(mode = 'list', length = length(sublist))
    rh_aparcs <- vector(mode = 'list', length = length(sublist))
    if (missing_gwr == FALSE){
      lh_wgs <- vector(mode = 'list', length = length(sublist))
      rh_wgs <- vector(mode = 'list', length = length(sublist))
    }
    for (sub in sublist){
      asegs[[x]] <- read.table(paste0(fsdir, sub, '/stats/aseg.stats'))
      lh_aparcs[[x]] <- read.table(paste0(fsdir, sub, '/stats/lh.aparc.stats'))
      rh_aparcs[[x]] <- read.table(paste0(fsdir, sub, '/stats/rh.aparc.stats'))
      if (missing_gwr == FALSE){
        lh_wgs[[x]] <- read.table(paste0(fsdir, sub, '/stats/lh.w-g.pct.stats'))
        rh_wgs[[x]] <- read.table(paste0(fsdir, sub, '/stats/rh.w-g.pct.stats'))
      }
      x <- x+1
    }
    
    # load labels
    aseg_txt <- readLines(file(paste0(fsdir, sub, '/stats/aseg.stats')))
    aseg_colheaders <- c(strsplit(aseg_txt[grep('# ColHeaders', aseg_txt)], "\\s+")[[1]])
    lh_aparc_txt <- readLines(file(paste0(fsdir, sub, '/stats/lh.aparc.stats')))
    lh_aparc_colheaders <- c(strsplit(lh_aparc_txt[grep('# ColHeaders', lh_aparc_txt)], "\\s+")[[1]])
    rh_aparc_txt <- readLines(file(paste0(fsdir, sub, '/stats/rh.aparc.stats')))
    rh_aparc_colheaders <- c(strsplit(rh_aparc_txt[grep('# ColHeaders', rh_aparc_txt)], "\\s+")[[1]])
    if (missing_gwr == FALSE){
      lh_wg_txt <- readLines(file(paste0(fsdir, sub, '/stats/lh.w-g.pct.stats')))
      lh_wg_colheaders <- c(strsplit(lh_wg_txt[grep('# ColHeaders', lh_wg_txt)], "\\s+")[[1]])
      rh_wg_txt <- readLines(file(paste0(fsdir, sub, '/stats/rh.w-g.pct.stats')))
      rh_wg_colheaders <- c(strsplit(rh_wg_txt[grep('# ColHeaders', rh_wg_txt)], "\\s+")[[1]])
    }
    
    # rename 
    y <- 1
    for (sub in sublist){
      aseg <- asegs[[y]]
      colnames(aseg) <- aseg_colheaders[3:length(aseg_colheaders)]
      
      lh_aparc <- lh_aparcs[[y]]
      colnames(lh_aparc) <- lh_aparc_colheaders[3:length(lh_aparc_colheaders)]
      
      rh_aparc <- rh_aparcs[[y]]
      colnames(rh_aparc) <- rh_aparc_colheaders[3:length(rh_aparc_colheaders)]    
      
      if (missing_gwr == FALSE){
        lh_wg <- lh_wgs[[y]]
        colnames(lh_wg) <- lh_wg_colheaders[3:length(lh_wg_colheaders)]
      
        rh_wg <- rh_wgs[[y]]
        colnames(rh_wg) <- rh_wg_colheaders[3:length(rh_wg_colheaders)] 
      }
      
      # label hemispheres
      
      lh_aparc$StructName <- paste0(lh_aparc$StructName, '_left')
      rh_aparc$StructName <- paste0(rh_aparc$StructName, '_right')
      if (missing_gwr == FALSE){
        lh_wg$StructName <- paste0(lh_wg$StructName, '_left')
        rh_wg$StructName <- paste0(rh_wg$StructName, '_right')
      }
      
      # label phenotypes 
      aseg_temp <- t(aseg[,c('StructName', 'Volume_mm3')])
      aseg_temp[1,] <- gsub("^(Right|Left)-(.*)$", "\\2_\\1", aseg_temp[1,])
      aseg_temp[1,] <- gsub("3rd-Ventricle", "X3rd.Ventricle", aseg_temp[1,])
      aseg_temp[1,] <- gsub("4th-Ventricle", "X4th.Ventricle", aseg_temp[1,])
      aseg_temp[1,] <- gsub("5th-Ventricle", "X5th.Ventricle", aseg_temp[1,])
      
      
      surfarea_temp <- cbind(t(lh_aparc[,c('StructName', 'SurfArea')]),
                             t(rh_aparc[,c('StructName', 'SurfArea')]))
      surfarea_temp[1,] <- paste0('SA_',surfarea_temp[1,])
      
      thickavg_temp <-  cbind(t(lh_aparc[,c('StructName', 'ThickAvg')]),
                              t(rh_aparc[,c('StructName', 'ThickAvg')]))
      thickavg_temp[1,] <- paste0('CT_',thickavg_temp[1,])
      
      grayvol_temp <-  cbind(t(lh_aparc[,c('StructName', 'GrayVol')]),
                             t(rh_aparc[,c('StructName', 'GrayVol')]))
      grayvol_temp[1,] <- paste0('GMV_',grayvol_temp[1,])
      
      if (missing_gwr == FALSE){
        wg_temp <-  cbind(t(lh_wg[,c('StructName', 'Mean')]),
                         t(rh_wg[,c('StructName', 'Mean')]))
        wg_temp[1,] <- paste0('GWR_',wg_temp[1,])
      }
      
      if (missing_gwr == FALSE){
        data_temp <- cbind(c('ID',paste0(sub)), aseg_temp, surfarea_temp, thickavg_temp, grayvol_temp, wg_temp)
      } else if (missing_gwr == TRUE){
        data_temp <- cbind(c('ID',paste0(sub)), aseg_temp, surfarea_temp, thickavg_temp, grayvol_temp)
      }
      
      colnames(data_temp) <- data_temp[1,]
      rownames(data_temp) <- c('label', 'value')
      
      data[y,] <- data_temp[2,]
      y <- y+1
    }
    
    # merge in IDs
   # data <- merge(sublist, data, by = 'ID')
  
  # check for missing data
  
  if (sum(colSums(is.na(data))[colSums(is.na(data)) != 0]) == 0){
    print('looks like you have no missing data. yay!')
  } else if (sum(colSums(is.na(data))[colSums(is.na(data)) != 0]) > 0){
    
    cat("\033[31mlooks like you have some missing data. here is the N missing values from your data:\033[0m\n")
    print(colSums(is.na(data))[colSums(is.na(data)) != 0])
    
    cat("\033[31mComBat cannot handle missing data. we recommend excluding any
    participants with missing ROIs or covariates, or you may want to run some kind of imputation procedure. 
    
    If you have certain ROIs with high incidences of missingness, you could
    consider just excluding those ROIs. If you do this (or have missing ROIs for some other reason), 
    we will automatically impute with data from the Dunedin Study in our DunedinPACNI estimates.
    ***** BEAR IN MIND: imputing ROI data from the Dunedin Study will worsen the accuracy of 
    DunedinPACNI estimates.\033[0m\n")
  }
  
# })
  
  return(data)
 
}

