source("FJX_Subroutines.R")
library(tidyverse)
library(fuzzyjoin)
library(janitor)
library(purrr)

# ---addX_FJX_sss.f  -- with user supplied subroutine that supplies X-section x Q-yield
# ---generates fast-JX 18-bin X-sections  revised and updated v73 (mprather,4/2015)
# Rewritten for R by J.F. Brewer, 4/22/22

# FTBL = 'FJX_Dir/XACET_LinkX.dat'
# NB = 100
# NS = 40000
# NJ = 18

# SRB = matrix(nrow = 15, ncol = NJ)
# WBIN = array(NB+1)
# FBIN = array(NB)
# ABIN = FBIN

# FFBIN = array(NJ)
# AABIN = FFBIN
# IJX = array(NB)

# W = array(NS)
# F = W

# IBINJ = array(0, dim = NS)

# MM = array(3)
# 
# TTT = array(3)
# PPP = TTT
addx <- function(FTBL, FJX_Dir){
  ###### These are variable names that actually make sense ######
  PratmoBins <- read_table(paste0(FJX_Dir, "/wavel-bins.dat"), col_names = c("PratmoBin", "Lambda_Start", "Lambda_End"), skip = 1, n_max = 78)
  PratmoBins$Lambda_End[78] = 2100 #Fix this assignment so everything else goes in bin 78
  # bin_assignment_percent <- read_csv("Rewrite_in_R/Bin_Assignment.csv", col_names = T, skip = 1)
  JXBin_Num <- read_table(paste0(FJX_Dir,"/wavel-bins.dat"), skip = 90, col_names = c("PratmoBin", "FJX_Bin"))
  SUSIM_Spectra <- read_table(paste0(FJX_Dir,"/solar-p05nm-UCI.dat"), skip = 2, col_names = c("SUSIM_Wavelength", "Flux"))
  
  # Here we read in the variables used - details are included in the Sample file
  # I didn't figure out how to read in the 'notes' column - not convinced it matters.
  # But units are in hPa for pressure and 1e19 for M
  PT_Combo_Used <- read_table(FTBL, skip = 3, col_names = c("Variable", "13", "5", "0"),
                              n_max = 3) %>% # Read in the crucial 3 lines for the altitude values
    transpose_df(.) %>% # flop rows and columns for better R formatting
    row_to_names(row_number = 1) %>% # First row is the Header/column name
    mutate_if(is.character,as.numeric) %>% # Convert remaining values to doubles
    dplyr::rename(`Altitude (km)` = Variable) %>% # Rename the Altitude vector
    mutate(`[M]` = `[M]`*1e19) # Convert M to correct units
  
  ###### Here they are in Prather-ese #####
  # WBIN <- PratmoBins$Lambda_Start
  SRB <- read_csv(paste0(FJX_Dir,"/Bin_Assignment.csv"), col_names = T, skip = 1) %>%
    dplyr::select(-`SR:`) # Opacity distribution function bins
  # IJX <- JXBin_Num
  # W = SUSIM_Spectra$Wavelength
  # Flux = SUSIM_Spectra$Flux
  
  # Here is where the temperature-Pressure interpolation happens. 
  # Temperature-Pressure is iterated using K, which runs through the three temps, pressures, and M values
  # Which Prather terms, respectively, XT, XP, and XM
  
  Integrated_XSQY <- X_MeVK(SUSIM_Spectra, FTBL)
  
  ##### Now assign each bin from I = 1:77 to each p05 nm microbin J from the SUSIM (1:40000) not already handled #####
  # Fuzzy Join allows us to two columns that are not equal according to some match function.
  # In this case, we join the PratmoBins array (Prather's binning from the PrAtmo model)
  # to the SUSIM wavelengths. SUSIM is in 0.05 nm bins, while PrAtmo uses 1 nm bins.
  # Therefore, I join each 1 nm bin as PratmoBin$labmda_start < SUSIM_Wavelength <= PratmoBin$Lambda_End
  
  Pratmo_Summary <- Do_PrAtmo_Binning(Integrated_XSQY = Integrated_XSQY, 
                                      SUSIM_Spectra = SUSIM_Spectra,
                                      PratmoBins = PratmoBins)
  
  ####### Then he converts from the 77-bin pratmo to the 18-bin fast-JX #####

  Final_FJX <- Do_FJX_Binning(Pratmo_Summary = Pratmo_Summary,
                              JXBin_Num = JXBin_Num,
                              SRB = SRB) %>%
    left_join(., PT_Combo_Used) # Add the scenario metadata back into the output
  
  return(Final_FJX)
}


