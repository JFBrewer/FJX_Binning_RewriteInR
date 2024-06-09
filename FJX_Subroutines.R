X_MeVK <- function(SUSIM_Spectra, FTBL){
  # -----------------------------------------------------------------------
  # WW = wavelength (nm) to calc Xsection for
  # XT = temerature (K) for interpolation
  # XP = pressure (hPa) can be used in Stern-Volmer formula if need be
  # XM = air density (#/cm3), ditto
  # XNEW = cross section (cm2) as a function of WW and XT (and XP, XM)
  # INIT = initialization:
  # if INIT.eq.0 .then reads in any tables and sets Xsect name to TITLNEW
  # 
  # inputs: WW, XT, XP, XM
  # outputs: XNEW, MM, TTT, PPP, ISX, ISP, TITLENEW, TITLTBL
  # inout: INIT
  # -----------------------------------------------------------------------
  
  # ---include the Temperatures (p, M) that you want to interpolate to and give tables.
  # ---general format for generating all Xsections
  
  # When you run this in "INIT" form, it just reads in a bunch of stuff and itializes variables
  # Run it in non-INIT form, and it returns an XNEW that corresponds to the species info (XS, QY)
  # as well as the local conditions (T, P, M)
  
  # In its new vectorized form, it just reads in and then calculates the appropriate XNEW in one go
  
  #### Read in the input.dat file #####
  # The example in Prather's code is 'XMeVK_JPL11X.dat'
  # This code handles similarly formatted data for now
  Species_name <- read_csv(FTBL, col_names = "Species Name",  n_max = 1) 
  Flags <- read_delim(FTBL, skip = 1, n_max = 2, col_names = c("Flag", "Target Value", "Info"),
                      del = " = ") %>%
    mutate(Use = case_when(Flag == `Target Value` ~ T, 
                           Flag != `Target Value` ~ F))
  Strat_Flag = Flags$Use[1]
  Pres_Flag = Flags$Use[2]
  
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
  
  # For whatever reason, this is where Prather wants you to read in the 'ACTUAL TEMPERATURES FOR INTERPOLATION'
  # This line, the last before the format statement, MUST correspond directly to the calculated values below
  # In contrast, the upper temperature statement (on line 6 of FTBL) can be whatever you please and it will get interpolated
  Temperature_Limits <- as.numeric(read_table(FTBL, skip = 12, col_names = c("13", "5", "0"), n_max = 1))
  
  Max_To_Read <- read_table(FTBL, skip = 13, n_max = 1, col_names = c("MaxN", "FormatC"))
  
  Combined_XSQY <- read_table(FTBL, skip = 14, 
                              col_names = c("Referenced_Wavelength", "Cross-Section", "13", "5", "0"),
                              n_max = Max_To_Read$MaxN) %>%
    pivot_longer(cols = c(`13`, `5`, `0`), names_to = "Altitude (km)", values_to = "QY") %>%
    mutate_if(is.character,as.numeric) %>% # Convert remaining values to doubles
    left_join(., PT_Combo_Used) %>% # If needed for nesting, do this
    # This is how you nest by temperature
    nest(Pressure_Case = c(Referenced_Wavelength, `Cross-Section`, QY)) %>%
    # The interpolators need to know about the Lambda, XS, and QY of the next wavelength step for interpolation
    mutate(Pressure_Case = map(Pressure_Case, ~.x %>%
                                 mutate(Next_WL = lead(Referenced_Wavelength),
                                        Next_XS = lead(`Cross-Section`),
                                        Next_QY = lead(QY)))) %>%
    unnest(Pressure_Case) #%>%
    # filter(!is.na(Next_WL)) # Filter out the very tail end case so it ends
  # # This is how you nest by wavelength
  # nest(Species_Data = c(`Altitude (km)`, QY, `p/hPa`, `[M]`, `T/K`, `Cross-Section`)) %>%
  # mutate(Next_WL = lead(Referenced_Wavelength),
  #        Next_Species_Data = lead(Species_Data))
  
  # The distinction between Temperatures_Used and Temperature_Limits is made in the Prather code so we maintain it here
  # It exists to allow you to set these things differently
  Temperatures_Used = unique(Combined_XSQY$`T/K`)
  
  SUSIM_Flux <- SUSIM_Spectra %>%
    filter(SUSIM_Wavelength >= min(Combined_XSQY$Referenced_Wavelength),
           SUSIM_Wavelength < max(Combined_XSQY$Referenced_Wavelength)) %>%
    fuzzy_left_join(., Combined_XSQY, 
                    by = c('SUSIM_Wavelength' = 'Referenced_Wavelength', 
                           'SUSIM_Wavelength' = 'Next_WL'),
                    match_fun = list(`>=`, `<`))
  
  # --Note that Xsects need to be interpolated for wavelength, but keyed to
  #   Sturm-Volmer and atmos lapse rate.- evaluted for 3 specific values in table
  #   can allow for interp if need be - Locate T-p interp numbers for all wavels
  
  XNEW_Output <- Calculate_XNew(Photovalue_df = SUSIM_Flux, 
                                Temperature_Limits = Temperature_Limits,
                                Temperatures_Used = Temperatures_Used) %>%
    left_join(., PT_Combo_Used)
  
  return(XNEW_Output)
  
}

Calculate_XNew <- function(Photovalue_df, Temperature_Limits, Temperatures_Used){
  # Photovalue_df is the SUSIM_FLUX given by the X_MeVK subroutine
  # Temperature_Limits correspond to Prather's QT3
  
  # QY, P, M, TK, and XS are the values given in the stored file.
  # Quantum Yield, Pressure, [M], Temperature in Kelvin, and the cross-sections
  # If for some reason, you were going to ask the subroutine to interpolate to a point outside this range,
  # Prather includes code here to force your temperatures back to within the range given.
  
  # I don't really see why you'd want this behavior. For now, I'm adding in a warning that will never trigger
  
  ##### The three referenced Lambda values are added to make this vectorizable
  # Referenced_Lambda refers to the wavelengths read-in in XMeVK_JPL11X.dat
  # Given_Lambda is the SUSIM Lambda
  # Referenced_Lambda_Difference is the difference between the wavelengths read-in in XMeVK_JPL11X.dat
  # This is not very robust right now, and currently ONLY WORKS IF THE SCHEME INCREMENTS REGULARLY
  # However, in his code, Prather has a check to make sure that diff can never be greater than 1 anyway
  # So I'm not suuuuper concerned about that
  
  # Of this section, Prather says:
  # --Note that Xsects need to be interpolated for wavelength, but keyed to
  #   Sturm-Volmer and atmos lapse rate.- evaluted for 3 specific values in table
  #   can allow for interp if need be - Locate T-p interp numbers for all wavels
  
  # Note that if you're using the referenced temperatures (here 220, 272, and 295),
  # the Temp interpolator QFACT will be (0, 1, and 1), respectively
  Wavelength_Interpolated_df <- Photovalue_df %>%
    rowwise() %>%
    # The following variable names are preserved mostly from Prather's code
    ##### These steps set up for what comes next ######
  mutate(QFACT = Calc_TempInterpolator(`T/K`, Temperature_Limits), # Temperature interpolator
         FW = Calc_LambdaDiff(Flux_Wavelength = SUSIM_Wavelength,
                              Bin_Lower_Lambda = Referenced_Wavelength,
                              Bin_Upper_Lamda = Next_WL), # Fractional wavelength
         ##### These steps perform the wavelength interpolation to put the XS and QY on the same WL as the Flux ######
         # Linearly interpolate the XS term for the fractional wavelength
         XS_Interpolated = QYXS_Interpolator(ThisValue = `Cross-Section`, NextValue = Next_XS, FW = FW),
         # Linearly interpolate the QY term for the fractional wavelength
         QY_Interpolated = QYXS_Interpolator(ThisValue = QY, NextValue = Next_QY, FW = FW)
  )
  
  ########## Now, we need to change paradigms ######## 
  # instead of leading and interpolating the columns by wavelength, we have to swap dimensionality to temperature
  # Step 1 is to get rid of the superfluous columns for neatness
  # We are already in SUSIM-Wavelength space, so we don't need the references anymore
  
  # This is causing me problems to code and also it seems pointless
  # The only reason to have this is if you for some reason want to feed in different QY values than the temperatures you wanna use

  if(sum(Temperature_Limits != Temperatures_Used) != 0){
    print("Hey, your temperature values in the Reference table don't match! Time for some INTERPOLATIN'")
    
    Temperature_Interpolated_df <- Wavelength_Interpolated_df %>%
      dplyr::select(-Next_WL, -Next_XS, -Next_QY, -`Cross-Section`, -QY, -Referenced_Wavelength, -FW, 
                    -`p/hPa`, -`[M]`, -`Altitude (km)`) %>%
      nest(wavelength_case = c(`T/K`, QFACT, XS_Interpolated, QY_Interpolated))
    
    # Now interpolate!
    # Prather says "interpolate X-section vs. Wavelength"
    # I initially assumed this interpolated across P-T regimes
    # But in fact that happened above.
    Fully_Interpolated <- Temperature_Interpolated_df %>%
      rowwise() %>%
      mutate(SV_Correction = SternVollmer_Interpolator(Wavelength_Case = wavelength_case,
                                                       Temperature_Limits = Temperature_Limits, 
                                                       Temperatures_Used = Temperatures_Used))
    
    # Clean up
    Final_Output <- Fully_Interpolated %>%
      dplyr::select(-wavelength_case) %>%
      unnest(SV_Correction) %>%
      dplyr::select(-This_QY, -Next_QY, -QFACT, -XS_Interpolated, -QY_Interpolated, -Temp_Corrected_QY)
  } else {
    # If your Temperature_Used and Temperature_Limits values are the same, the time-consuming interpolation step is entirely unnecessary.
    Final_Output <- Wavelength_Interpolated_df %>%
      mutate(Integrated_QYXS = XS_Interpolated*QY_Interpolated,
             XNEW = Integrated_QYXS*1e-20) %>%
      dplyr::select(SUSIM_Wavelength, Flux, `T/K`, Integrated_QYXS, XNEW)
  }
  
  return(Final_Output)
}

transpose_df <- function(df) {
  # This is a custom function to transpose a tibble while preserving names
  # It's necessary so I can put things in a useful tibble format
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}

Calc_TempInterpolator <- function(TK, Temperature_Limits){
  # If for some reason, you were going to ask this fitting function to interpolate to a Temp outside this range,
  # Prather includes code here to force your temperatures back to within the range given.
  
  if(is.na(TK)){
    return(NA)
  } else {
    
    # I don't really see why you'd want this behavior. For now, I'm adding in a warning that will never trigger
    if(TK < min(Temperature_Limits)) {
      print(paste0("You've tried to interpolate to a temperature below the acceptable range. Resetting to ",
                   min(Temperature_Limits)))
      Use_TK = min(Temperature_Limits)
    } else if(TK > max(Temperature_Limits)) {
      print(paste0("You've tried to interpolate to a temperature above the acceptable range. Resetting to ",
                   max(Temperature_Limits)))
      Use_TK = max(Temperature_Limits)
    } else{
      Use_TK = TK
    }
    
    # Here we allow for interpolation within temperatures above the midpoint of Temperature_Limits
    # (272 in the case of MeVK)
    if(Use_TK > Temperature_Limits[2]){
      IT = 2
    } else {
      IT = 1
    }
    
    # Now calculate QFACT for use in the code
    QFACT = (Use_TK - Temperature_Limits[IT])/(Temperature_Limits[IT+1] - Temperature_Limits[IT])
    return(QFACT)
  }
}

Calc_LambdaDiff <- function(Flux_Wavelength, Bin_Lower_Lambda, Bin_Upper_Lamda){
  # This interpolator allows you to correct for flux wavelengths of smaller resolution than your bin sizing
  # Lambda_Diff is the distance between the SUSIM wavelength and the lower limit of the reference bin,
  # Normalized by the width of the reference bin. It corresponds to FW in Prather's code
  if(is.na(Flux_Wavelength) | is.na(Bin_Upper_Lamda) | is.na(Bin_Lower_Lambda)){
    return(NaN)
  } else {
    # However, despite allowing for bin sizing larger than 1, Prather then immediately forces max bin size to 1
    Bin_Size = Bin_Upper_Lamda - Bin_Lower_Lambda
    LambdaDiff = (Flux_Wavelength - Bin_Lower_Lambda)/Bin_Size
    if(Bin_Size > 1){
      print("Max bin size allowed is 1 nm. This result may not perform as expected.")
    }
    LambdaDiff = min(1,max(0, LambdaDiff))
    return(LambdaDiff)
  }
}

QYXS_Interpolator <- function(ThisValue, NextValue, FW){
  # Simple linear interpolation of the quantum yield and cross-section by the LambdaDiff term calculated in the previous step
  if(is.na(ThisValue) | is.na(FW)){
    return(NaN)
  } else if(is.na(NextValue)) {
    #Final interpolation ends on 0 to preserve behavior
    Interpolated_Value = 0 
  } else {
    Interpolated_Value = ThisValue + FW*(NextValue - ThisValue)
  }
  return(Interpolated_Value)
}

SternVollmer_Interpolator <- function(Wavelength_Case, Temperature_Limits, Temperatures_Used){
  # This is the final step of this interpolation, and it's messy because I'm trying to vectorize it.
  # Basically, we've shifted into wavelength space for the interpolation.
  # The actual interpolation happens in the very last step, everything else is prep.
  
  # First, while we went to all that trouble to calculate an interpolated QY, we actually have to change our thinking
  # Prather has created a case-handling system that allows him to set Temperature limits that don't correspond to this actually used
  # It's unclear to me why you'd do this, since presumably the QY temperature cases would correspond directly to the temperatures used
  # However, in order to be consistent with his code, I've maintained the behavior.
  
  # In fact, if Temperatures_Used = Temperature_Limits, all of this is pointless.
  
  # First, unpack the list
  # Wavelength_Case <- Wavelength_Case[[1]]
  
  # print(paste0("This is a test: ", Wavelength_Case$`T/K`[1]))
  
  ##### Now set the case values for This_QY and Next_QY #####
  # Case 1 is the value for This_QY and Next_QY to be used if the Temperature used is greater than Temperature_Limits[2]
  ThisQY_Case_1 <- filter(Wavelength_Case, `T/K` == Temperatures_Used[2]) %>%
    dplyr::select(QY_Interpolated) %>%
    as.numeric(.)
  NextQY_Case_1 <- filter(Wavelength_Case, `T/K` == Temperatures_Used[3]) %>%
    dplyr::select(QY_Interpolated) %>%
    as.numeric(.)
  
  # Case 2 is the value for This_QY and Next_QY to be used if the Temperature used is less than or equal to Temperature_Limits[2]
  ThisQY_Case_2 <- filter(Wavelength_Case, `T/K` == Temperatures_Used[1]) %>%
    dplyr::select(QY_Interpolated) %>%
    as.numeric(.)
  NextQY_Case_2 <- filter(Wavelength_Case, `T/K` == Temperatures_Used[2]) %>%
    dplyr::select(QY_Interpolated) %>%
    as.numeric(.)
  
  
  ##### Now apply these cases #####
  Case_Logic_For_Temperature <- Wavelength_Case %>%
    mutate(This_QY = case_when(`T/K` > Temperature_Limits[2] ~ ThisQY_Case_1,
                               `T/K` <= Temperature_Limits[2] ~ ThisQY_Case_2),
           Next_QY = case_when(`T/K` > Temperature_Limits[2] ~ NextQY_Case_1,
                               `T/K` <= Temperature_Limits[2] ~ NextQY_Case_2))
  
  ####### Perform the Temperature interpolation HERE ######
  # This step corrects for Temperature effects and the QY differences!
  Interpolated_Integrated_QYXS <- Case_Logic_For_Temperature %>%
    # rowwise() %>%
    mutate(Temp_Corrected_QY = This_QY + QFACT*(Next_QY-This_QY),
           Integrated_QYXS = XS_Interpolated*Temp_Corrected_QY,
           XNEW = Integrated_QYXS*1e-20)
  
  return(list(Interpolated_Integrated_QYXS))
}

Do_FJX_Binning <- function(Pratmo_Summary, JXBin_Num, SRB = SRB){
  # This function used to have more roles but now it just does the application of the
  # Opacity distribution function, which contains all of the relevant stuff anyway!
  
  OpacityDistOutput <- Pratmo_Summary %>% 
    group_by(`Altitude (km)`) %>%
    nest(Pratmo_Binning = c(PratmoBin, PratmoFlux, PratmoJ)) %>%
    mutate(FJX_Binning = map(Pratmo_Binning, ~apply_OpacityDistFn(PratmoNest = .,
                                                      JXBin_Num_in = JXBin_Num,
                                                      SRB_in = SRB)))
  
  # and then for good measure, we get this, for all J 1-18
  # if FFBIN[J] > 0 {AABIN[J] = AABIN[J]/FFBIN[J]}
  
  return(OpacityDistOutput)
}

apply_OpacityDistFn <- function(JXBin_Num_in, PratmoNest, SRB_in){
  # This starts with an iterator I which runs from 16-on, which corresponds to wavelengths 202+
  # This binning is then correlated to a JX_bin_No, given in our JXBin_Num (his IJX) using the J iterator
  # FBIN is collected by FFBIN, ABIN is collected by AABIN, but MODIFIED BY FBIN
  
  # FFBIN[J] = FFBIN[J] + FBIN[I]
  # AABIN[J] = AABIN[J] + FBIN[I]*ABIN[I]
  
  # Add this in order to account for the extra 4 bins that don't get mention in Prather's dictionary
  head_of_FJX_Binning <- tibble(FJX_Bin = seq(1,4), FJX_Flux = rep(0, 4), FJX_J = rep(0, 4))
  
  FJX_Binning <- left_join(PratmoNest, JXBin_Num_in) %>% # Use Prather's dictionary to match FJX, PrAtmo bins
    filter(!is.na(FJX_Bin)) %>% # This seems unnecessary but what the hey
    ungroup() %>%
    dplyr::select(FJX_Bin, PratmoFlux, PratmoJ) %>% # ditch the pratmo bin, we don't care anymore
    group_by(FJX_Bin) %>% # condense and sum on the FJX bins!
    summarise(FJX_Flux = sum(PratmoFlux, na.rm = T),
              FJX_J = sum(PratmoJ*PratmoFlux, na.rm = T)) %>%
    bind_rows(head_of_FJX_Binning, .) # Assemble full FJX Binning
  
  # For PratmoBins 1-15, we do something else, with two iterators, I and J
  # I runs from 1-15, accounting for each PratmoBin
  # J runs from 1-18
  
  # FFBIN[J] = FFBIN[J] + FBIN[I]*SRB[I,J]
  # AABIN[J] = AABIN[J] + FBIN[I]*ABIN[I]*SRB[I,J]
  # This is irrelevant for tropospheric species but important for stratospheric, I guess?
  
  # Declare the relevant FFBIN (Flux) and AABIN (J-Value) arrays
  FFBIN = numeric(length = max(JXBin_Num_in$FJX_Bin))
  AABIN = numeric(length = max(JXBin_Num_in$FJX_Bin))
  SRB = as.data.frame(SRB_in)
  
  # Here we apply the Opacity Distribution Function from Prather's code
  for(ii in 1:15){
    for(jj in 1:max(JXBin_Num_in$FJX_Bin)){
      ###### Master equation for flux:  FFBIN[J] = FFBIN[J] + FBIN[I]*SRB[I,J] ####
      # First calculate the flux
      FluxThisStep = PratmoNest$PratmoFlux[PratmoNest$PratmoBin == ii]*SRB[jj, ii]
      # Now sum across the relevant binning
      FFBIN[jj] = sum(FFBIN[jj], FluxThisStep,  na.rm = T)
      
      ###### Master Equation for J-Value: AABIN[J] = AABIN[J] + FBIN[I]*ABIN[I]*SRB[I,J] ####
      # First calculate the J-Value
      JThisStep = PratmoNest$PratmoJ[PratmoNest$PratmoBin == ii]*
        PratmoNest$PratmoFlux[PratmoNest$PratmoBin == ii]*SRB[jj, ii]
      # Now sum across the relevant binning
      AABIN[jj] = sum(AABIN[jj],  JThisStep, na.rm = T)
    }
  }
  
  # Add the condensed Opacity Distibutions for flux and J-value to the relevant
  # columns from the prior FJX-binning
  FJX_Binning$FJX_Flux = FJX_Binning$FJX_Flux + FFBIN
  FJX_Binning$FJX_J = FJX_Binning$FJX_J + AABIN
  
  FJX_out <- FJX_Binning %>%
    rowwise() %>%
    mutate(FJX_J = case_when(FJX_Flux > 0 ~ FJX_J/FJX_Flux))
  # The last step is once again the wavelength weighting
  
  return(FJX_out)
}

Do_PrAtmo_Binning <- function(Integrated_XSQY, SUSIM_Spectra, PratmoBins){
  ##### Now assign each bin from I = 1:77 to each p05 nm microbin J from the SUSIM (1:40000) not already handled #####
  # Fuzzy Join allows us to two columns that are not equal according to some match function.
  # In this case, we join the PratmoBins array (Prather's binning from the PrAtmo model)
  # to the SUSIM wavelengths. SUSIM is in 0.05 nm bins, while PrAtmo uses 1 nm bins.
  # Therefore, I join each 1 nm bin as PratmoBin$labmda_start <= SUSIM_Wavelength < PratmoBin$Lambda_End
  Altitudes <- tibble(`Altitude (km)` = c(0, 5, 13))
  
  Shorter_WLs <- SUSIM_Spectra %>%
    filter(SUSIM_Wavelength < min(Integrated_XSQY$SUSIM_Wavelength)) %>%
    fuzzy_left_join(., PratmoBins,
                    by = c('SUSIM_Wavelength' = 'Lambda_Start', 'SUSIM_Wavelength' = 'Lambda_End'),
                    match_fun = list(`>`, `<=`)) %>%
    filter(!is.na(PratmoBin)) %>% # remove wavelengths for which no PratmoBin is assigned
    crossing(., Altitudes) # Add in the altitude cases to make future steps easier
  
  Shorter_WL_Sub16 <- Shorter_WLs %>%
    filter(PratmoBin <= 15) # Separate bins 1-15 for spectral oscillator function treatment
  
  Shorter_WL_16Plus <- Shorter_WLs %>%
    filter(PratmoBin > 15) # Separate bins 16+ for standard XS-QY interpolation
  
  # This incorporates all wavelengths longer than those for which we have XS/QY data
  # It's only really relevant to match prather's output and contributes no new information.
  Longer_WLs <- SUSIM_Spectra %>%
    filter(SUSIM_Wavelength > max(Integrated_XSQY$SUSIM_Wavelength)) %>%
    fuzzy_left_join(., PratmoBins,
                    by = c('SUSIM_Wavelength' = 'Lambda_Start', 'SUSIM_Wavelength' = 'Lambda_End'),
                    match_fun = list(`>`, `<=`)) %>%
    filter(!is.na(PratmoBin)) %>%
    crossing(., Altitudes)
  
  # At this point, Prather's code loops over the total number of SUSIM Bins (40,000), NS_
  # His iterator for this is J
  # He then sets a second iterator I, which is the Prather Bin number associated with the SUSIM bin
  # If I is >0 (i.e., if PratmoBin isn't NA), then we call the subroutine X_MeVK
  
  # X_MeVK subroutine returns an integrated XS-QY term XNEW on the 0.05 nm SUSIM wavelength binning.
  # Prather then converts to Pratmo binning by summing up the flux and J-Value (FBIN and ABIN)
  # FBIN = FBIN[I] + F[J]  - summed flux per PrAtmo bin
  # ABIN = ABIN[I] + F[J]*XNEW - summed J-value per PrAtmo bin
  
  # This is then further modified: if FBIN[I] > 0, then NB77 = I and ABIN[I] = ABIN[I]/FBIN[I]
  # This is the wavelength normalization
  
  # Prather then writes out the 77-binned UCI data for pratmo - that is, ABIN
  
  Pratmo_Binning <- Integrated_XSQY %>%
    fuzzy_left_join(., PratmoBins, 
                    by = c('SUSIM_Wavelength' = 'Lambda_Start', 'SUSIM_Wavelength' = 'Lambda_End'),
                    match_fun = list(`>`, `<=`))
  
  # We treat the first 15 bins differently than the rest! Only relevant in the very short-UV (so, Stratosphere, p much)
  Sub16s <- Pratmo_Binning %>% 
    filter(PratmoBin <= 15) %>%
    bind_rows(., Shorter_WL_Sub16) %>%
    mutate(J_Value = Flux*XNEW) %>%
    group_by(PratmoBin, `Altitude (km)`) %>%
    summarise(PratmoFlux = sum(Flux),
              PratmoJ = sum(J_Value)) %>%
    mutate(PratmoJ = case_when(PratmoFlux > 0 ~ PratmoJ/PratmoFlux,
                               PratmoFlux <= 0 ~ PratmoJ))
  
  # normalize by wavelength weighting
  
  Pratmo_Summary <- Pratmo_Binning %>%
    filter(PratmoBin >= 15) %>%
    bind_rows(Shorter_WL_16Plus, .) %>% # Account for the shorter WLs, I guess?
    bind_rows(., Longer_WLs) %>%
    mutate(J_Value = Flux*XNEW) %>%
    group_by(PratmoBin, `Altitude (km)`) %>%
    summarise(PratmoFlux = sum(Flux),
              PratmoJ = sum(J_Value)) %>%
    mutate(PratmoJ = case_when(PratmoFlux > 0 ~ PratmoJ/PratmoFlux,
                               PratmoFlux <= 0 ~ PratmoJ)) %>%
    bind_rows(Sub16s, .)
  # normalize by wavelength weighting
  
  return(Pratmo_Summary)
}
