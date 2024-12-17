library(pacman)
p_load("ggplot2", "tidyverse", "doParallel", "foreach")
options(pillar.sigfig = 5)


# # # #
# Function to extract chromatogram from a raw file
# # # #

extractChromato = function(rawFile, curr_mol_data, hdata, rt_window, mztol_ppm){
  # Extract the scan matching the precursor mass
  current_mz = curr_mol_data$Frag
  
  # if MS2 available, filter scan based on precursor m/z
  if(unique(curr_mol_data$MS == 2)){
    hdata = hdata %>% filter(StartTime <= rt_window$up & StartTime >= rt_window$down &
                               hdata$polarity == unique(curr_mol_data$Polarity) & 
                               hdata$msLevel == unique(curr_mol_data$MS))
    
    hprecmz = unique(hdata$precursorMass)
    precmz_diff = abs(hprecmz - curr_mol_data$Precursor)
    select_scan = hdata$scanType[which(hdata$precursorMass == hprecmz[which.min(precmz_diff)])] %>% unique()
  }else{
    select_scan = hdata$scanType[which(hdata$polarity == unique(curr_mol_data$Polarity) & 
                                         hdata$msLevel == unique(curr_mol_data$MS))] %>% unique()
  }
  
  chromato = tibble()
  # Extract chromatogram for given mz windows
  if(length(select_scan) == 1){
    # Try to extract chromatogram and catch error if there is no data
    tryCatch({
      chromato = readChromatogram(rawFile, mass=current_mz, tol=mztol_ppm, filter=select_scan, type="xic")[[1]]
      chromato = tibble(intensity = chromato$intensities, rtime = chromato$times) %>%
        filter(rtime <= rt_window$up & rtime >= rt_window$down)
    }, error=function(er){
    }
    )
  } 
  return(chromato)
}

# # # #
# Function to extract chromatograms and save a ggplot for the TargetPeak as input to Deep Learning model
# # # #

getChromato = function(rawFile, hdata, mol_data, mztol_ppm=10, plot_name="", rt_window_min=0.5, 
                       rt_mol=NULL, rt_ISTD = NULL, correct_RT = F){
  
  # RT window
  if(is.null(rt_mol)) rt_mol = unique(mol_data$RT)
  rt_window = list(down = rt_mol - rt_window_min, up = rt_mol + rt_window_min)
  
  # Select the closest ISTD to calculate devRRT
  closest_ISTD = which.min(abs(rt_mol - rt_ISTD$th))
  
  # If option TRUE, correct the RT based on the ISTD if they were detected
  if(correct_RT & length(rt_ISTD) > 0){
    # Correct the RT only if the shift is above 0.2s
    if(abs(rt_ISTD$obs[closest_ISTD] - rt_ISTD$th[closest_ISTD]) > 0.2){
      # Calibrate chromato based on closest ISTD ratio
      ratio_ISTD = rt_ISTD$obs[closest_ISTD]/rt_ISTD$th[closest_ISTD]
      rt_window = lapply(rt_window, function(x) x*ratio_ISTD)
    }
  }
  
  chromato_tab = tibble(intensity = numeric(), rt = numeric(), mz = numeric(),
                        typePeak = character())
  
  # For each MZ, extract the chromatogram
  # End the function at the TargetPeak for ISTD
  for(id_mz in 1:nrow(mol_data)){
    
    chromato = extractChromato(rawFile, mol_data[id_mz,], hdata, rt_window, mztol_ppm)
    
    # If there is a chromatogram
    if(nrow(chromato)){
      chromato_tab = add_row(chromato_tab, tibble(intensity=chromato$intensity, 
                                                  rt=chromato$rtime, mz=mol_data[id_mz,]$Frag,
                                                  typePeak = mol_data[id_mz,]$Workflow))
    }
    
    # If there is some peak (intensity > 0)
    if(sum(chromato$intensity > 0)){
      # Pick peaks using cardidates package
      peaksInfo = peakwindow(chromato$rtime, chromato$intensity)
      peaksRT = peaksInfo$peaks$x
      
      # It's an ISTD, get the RT and end the function
      if(is.null(rt_ISTD) & mol_data$Workflow[id_mz] == "TargetPeak"){
        ISTD_peak = round(peaksRT[which.min(abs(peaksRT-rt_mol))], 2)
        ISTD_RT = list(obs = ISTD_peak, th=rt_mol)
        return(ISTD_RT)
      }
      # If there is no peak but it's an ISTD
    }else if(is.null(rt_ISTD)){
      return(NULL)
    }
  } #end for
  
  
  # Select the confirming peak to use based on the max of intensity 
  conf_to_use = chromato_tab %>% filter(typePeak=="Confirming") %>% 
    filter(intensity==max(intensity))
  
  # If there is multiple line, keep the lower m/z
  if(nrow(conf_to_use) > 1) conf_to_use = conf_to_use[which.min(conf_to_use$mz),]
  
  chromato_tab$typePeak[chromato_tab$mz == conf_to_use$mz] = "Confirming used"
  
  # Write JPEG input for the IA-assisted tool if it's the Target Peak
  # Or if it's the confirming peak to use
  chromato_to_save = chromato_tab %>% filter(typePeak=="TargetPeak" | typePeak=="Confirming used") %>%
    group_by(mz) %>% group_map(function(chromato, ...){
      
      plotnameEdited = ifelse(unique(chromato$typePeak) == "TargetPeak", plot_name, 
                              paste0(plot_name, "_confirming"))
      peakinfo_name = ifelse(unique(chromato$typePeak) == "TargetPeak", "peakinfo.csv", 
                           "peakinfo_confirming.csv")
      
      if(sum(chromato$intensity > 0)){
        
        # Pick peaks using cardidates package
        peaksInfo = peakwindow(chromato$rt, chromato$intensity)
        peaksRT = peaksInfo$peaks$x
        
        # If rt_ISTD is given, so the compound is not an ISTD, calculate the RRT deviation
        if(!is.null(rt_ISTD)){
          
          # If there is no peak to extract (ex: peak is not entirely in the RT window), set devRRT to 1
          if(!length(peaksRT)){
            devRRT=1
            peaksArea=0
            nbPoints=0
            
            # Else, there is at least one peak
          }else{
            
            # Calculate RRT deviation and keep the closest
            devRRT = abs(peaksRT/rt_ISTD$obs[closest_ISTD] - rt_mol/rt_ISTD$th[closest_ISTD])
            selected_peak = which.min(devRRT)
            devRRT = round(devRRT[selected_peak], 5)
            
            # Get intensity at the apex of the peak of interest
            peaksArea = round(sum(chromato$intensity[which(peaksInfo$peakid == selected_peak)]), 2)
            # Get the number of points around the peak of interest
            nbPoints = sum(chromato$intensity[which(peaksInfo$peakid == selected_peak)] > 0)
            
          }
          
          readr::write_csv(tibble(plotname=paste0(plotnameEdited, ".jpeg"), area=peaksArea, 
                                  points=nbPoints, devRRT=devRRT), 
                           file.path(dirname(plotnameEdited), peakinfo_name), append = T)
          
          # Smooth intensity with Savitzky-Golay algorithm using quadratic filter (2)
          chromato$intensity_smoothed = pracma::savgol(chromato$intensity,5,2,0)
          chromato$intensity_smoothed[chromato$intensity_smoothed < 0.1] = 0
          #chromato$intensity_smoothed = rev(chromato$intensity_smoothed)
          
          ggp = ggplot(chromato) + geom_line(aes(rt, intensity_smoothed)) +
            #geom_line(aes(rt, intensity)) +
            theme_classic() +
            theme(axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank())
          
          ggsave(filename = basename(paste0(plotnameEdited, ".jpeg")), device = "jpeg",
                 plot = ggp, path = dirname(plotnameEdited), units = "in", height = 2, width = 2)
        }
      }else if(sum(chromato$intensity == 0) & !is.null(rt_ISTD)){
        # Save an empty plot if the sum of intensity is equal to zero
        
        readr::write_csv(tibble(plotname=paste0(plotnameEdited, ".jpeg"), area=0, 
                                points=0, devRRT=1), 
                         file.path(dirname(plotnameEdited), peakinfo_name), append = T)
        
        ggp = ggplot() + geom_hline(yintercept = 0) +
          coord_cartesian(ylim = c(0,1)) +
          theme_classic() +
          theme(axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank())
        
        ggsave(filename = basename(paste0(plotnameEdited, ".jpeg")), device = "jpeg",
               plot = ggp, path = dirname(plotnameEdited), units = "in", height = 2, width = 2)
      }
      
    })#end for
  
  if(!is.null(rt_ISTD)){
    return(chromato_tab)
  }
  
}


# # # #
# Function to randomly shift RT
# # # #

rt_shift = function(rt){
  shift = sample(seq(-0.05,0.05, by=0.02), size=1)
  return(rt+shift)
}


# # # #
# Function to get a smoothed figure in the PDF report for a molecule in a sample
# Need : 
#   - extracted chromatogram
#   - information about the compound selected
#   - the name of the current sample
# Return :
#   - An arranged plot with stacked chromato and split transitions for the PDF report
# # # #
plotReport = function(chromato_data, mol_data, current_file){
  # Build ggplot of chromatogram and save it
  chromato_data_plot = merge(mol_data, chromato_data, by.x="Frag", by.y="mz") %>% 
    mutate(facet_title = paste(Frag, typePeak)) %>% 
    select(rt, intensity, Frag, facet_title) %>% 
    # Smooth intensity with Savitzky-Golay algorithm using quadratic filter (2)
    group_by(Frag) %>% group_modify(function(frag_tab, ...){
      frag_tab$intensity_smoothed = pracma::savgol(frag_tab$intensity,5,2,0)
      frag_tab$intensity_smoothed[frag_tab$intensity_smoothed < 0.1] = 0
      frag_tab
    })
  
  chromato_data_plot$intensity[is.na(chromato_data_plot$intensity)] = 0
  
  facet_ggp = ggplot(data=chromato_data_plot, aes(x=rt, y=intensity_smoothed)) + 
    facet_wrap(~facet_title, scales = "free", ncol=1)+
    geom_line(alpha=0.8)+
    scale_color_brewer(palette="Set1")+
    theme_half_open()+
    scale_y_continuous(expand=c(0,0)) +
    theme(strip.background = element_blank(), axis.text = element_text(size=10), 
          strip.text = element_text(size=10), axis.title = element_text(size=10))
  
  # Plot all ions
  full_ggp = ggplot(data=chromato_data_plot, aes(x=rt,y=intensity_smoothed, color=as.factor(Frag)))+
    geom_line(alpha=0.8)+
    scale_color_brewer(palette="Set1")+
    theme_half_open()+
    scale_y_continuous(expand=c(0,0)) +
    theme(legend.position = 'bottom',legend.text = element_text(size=10), legend.title = element_blank(),
          axis.text = element_text(size=10), axis.title = element_text(size=10))
  
  return(plot_grid(plotlist = list(facet_ggp, full_ggp), ncol=2) %>%
           annotate_figure(top = text_grob(
             paste(tools::file_path_sans_ext(basename(current_file)), unique(mol_data$Compound), sep="\n")
             , size = 12)))
}


# # # #
# Function to write CSV input for DL model
# Need : 
#   - path to raw files
#   - path to CSV containing targets information
#   - output directory
#   - mz error in ppm
#   - number of cores to use
#   - allow correction of RT based on ISTD
# Return :
#   - CSV file as input for Deep Learning model
#   - Sample PDF report with each molecules
# # # #
formatCSV_DLinput = function(pathRAW, moleculesCSV, out_dir, mz_ppm = 10, 
                             cores=20, correct_RT = F){
  
  # Create directory if necessary
  list_outdir = paste0(out_dir, c("", "/Rapport"))
  dir_to_create = which(!dir.exists(list_outdir))
  lapply(list_outdir[dir_to_create], dir.create)
  
  
  # Get each RAW files in /datas/DL_LCH
  path_to_data = list.files(pathRAW, full.names = T, pattern="raw")
  
  molecules = readr::read_csv2(moleculesCSV, show_col_types = FALSE)
  molecules_ISTD = molecules %>% filter(ISTD == T) %>% group_by(Compound)
  molecules_bycomp = molecules %>% group_by(Compound)
  
  # Init devRRT csv
  readr::write_csv(tibble(plotname=character(), area=integer(), points=integer(), devRRT=integer()), 
                   file.path(out_dir, "peakinfo.csv"))
  readr::write_csv(tibble(plotname=character(), area=integer(), points=integer(), devRRT=integer()), 
                   file.path(out_dir, "peakinfo_confirming.csv"))
  
  cl = makeCluster(cores)
  registerDoParallel(cl)
  
  start = Sys.time()
  
  # For each mzML file, search for molecules in screening
  chromatogram = foreach(current_file=path_to_data, .packages = c("MSnbase", "rawrr", "tibble", "dplyr", "Rdisop", 
                                                                  "stringr", "cowplot", "ggpubr", "ggplot2", "cardidates"), 
                         .export = c("getChromato", "rt_shift", "plotReport", "extractChromato"))%dopar%{
                           
                           hdata = readIndex(current_file)
                           
                           hdata$polarity = ifelse(grepl("\\+", hdata$scanType), 1, 0)
                           hdata$msLevel = ifelse(hdata$MSOrder=="Ms2", 2, 1)
                           
                           # Extract ISTD retention time
                           rt_ISTD = list(obs=c(), th=c())
                           for(ISTD in unique(molecules_ISTD$Compound)){
                             mol_ISTD = molecules_ISTD %>% filter(Compound==ISTD)
                             curr_rt_ISTD = getChromato(rawFile = current_file,
                                                      hdata = hdata,
                                                      mol_data = mol_ISTD,
                                                      mztol_ppm = mz_ppm,
                                                      rt_window_min = 0.5)
                             
                             rt_ISTD$obs = c(rt_ISTD$obs, curr_rt_ISTD$obs)
                             rt_ISTD$th = c(rt_ISTD$th, curr_rt_ISTD$th)
                           }
                           
                           
                           # For each molecule to search, extract chromatogram and save target peak plot in the out_dir
                           plots.list = molecules_bycomp %>% group_map(function(mol_data, comp){
                             # shifted_rt = rt_shift(unique(mol_data["RT"])$RT)
                             chromato_data = getChromato(rawFile = current_file,
                                                         hdata = hdata,
                                                         mol_data=mol_data,
                                                         mztol_ppm = mz_ppm,
                                                         plot_name = paste0(out_dir, "/", tools::file_path_sans_ext(basename(current_file)), "_comp=", comp),
                                                         # rt_mol = shifted_rt*60,
                                                         rt_window_min = 0.5, 
                                                         rt_ISTD = rt_ISTD, correct_RT = correct_RT)
                             
                             if(nrow(chromato_data)){
                               mol_data$Compound = comp
                               return(plotReport(chromato_data, mol_data, current_file))
                             }
                             
                           })
                           
                           plots = gridExtra::marrangeGrob(plots.list[lengths(plots.list) > 0], nrow=2, ncol=1)
                           # To save to file, here on A4 paper
                           ggsave(paste0(out_dir, "/Rapport/", tools::file_path_sans_ext(basename(current_file)), ".pdf"), 
                                  plots, width = 21, height = 29.7, units = "cm")
                         }
  
  stopCluster(cl)
  
  stop=Sys.time()
  
  print(stop-start)
  
}

###
# LAUNCH THE MAIN FUNCTION TO GENERATE INPUTS
# 
###

main(pathRAW = "/rawData/", 
     moleculesCSV = "my_method.csv",
     out_dir = "/outputs/",
     mz_ppm=10, correct_RT = T)

