# HEADER --------------------------------------------------------------
# SCRIPT NAME: mfafunctions.R
#
# DATE CREATED: 19 June 2020
#
# AUTHOR: Joshua Walters
#
# DESCRIPTION: This script contains functions for extracting data
# from XML files for the microplate feeder assay.
#
#
# NOTES
# - libraries required:
#                       lubridate
#                       tidyverse
#                       xml2
#                       abind

#' parse_time: converts date-time to lubridate format for easier
#' time calculations
#'
#' @param date_time_string A character string containing the
#' date and time as exported from Molecular Devices iD5 plate reader.
#' i.e. HH:MM am/pm M/D/YYYY
#'
#' @keywords time, date
#'
#' @export
#'
#' @examples
#' ReadTime = "2:54 PM 7/2/2020"
#' time = parse_time(ReadTime)
#' print(time)
#'
#' Output: 2020-07-02 14:54:00
parse_time = function(date_time_string) {
  time = readr::parse_datetime(date_time_string,"%I%.%M %p %m/%d/%Y")
  time = lubridate::force_tz(time, tzone = "US/Eastern")
  return(time)
}

#' parse_IDs: converts a vector of sample IDS into a matrix of
#' factors for each respective ID.
#'
#' @param id A vector containing character strings of sample IDs
#' with each factor separated by a character. sep The character
#' separating factors within the sample ID. Default sep = "_".
#'
#' @keywords IDs, samples, factors
#'
#' @export
#'
#' @examples
#' sampleID = c("DGRP_0120_M_MUSH",
#'              "DGRP_5430_F_BLIP")
#' factors = parse_IDs(sampleID)
#' print(factors)
#'
#' Output:
#'      [,1]                [,2]    [,3]    [,4] [,5]
#' [1,] "DGRP_0120_M_MUSH", "DGRP", "0120", "M", "MUSH"
#  [2,] "DGRP_5430_F_BLIP", "DGRP", "5430", "F", "BLIP"
parse_IDs = function(id, sep = "_") {
  id.list = data.frame(NULL)

  # Iterates over IDs in id.list to split and combine the IDs and
  # individual factors
  for (i in seq(1,length(id))) {
    temp = unlist(strsplit(id[i],split = sep))
    temp = append(temp, id[i], after=0)
    id.list = rbind(id.list, temp)
  }

  colnames(id.list) = c("ID", "Project", "Line", "Sex", "Assay")
  return(id.list)
}


#' import_plates: generates a full readout of all plate reads in an
#' xml document. Calls parse_IDs and parse_time functions to reformat
#' the xml attributes.
#'
#'
#' @param imp_xml An internal xml document imported via xml2::read_xml().
#'
#' @keywords import, plates, read information
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#' all_reads = import_plates(imp_xml)
#' print(all_reads)
#'
#' Output:
#' ID                  Project Line   Sex  Assay    ReadTime
#' HARD_0420_M_BLZR    HARD    0420   M    BLZR     2020-07-02 11:54:00
#' TIPS_8008_F_SPEC    TIPS    8008   F    SPEC     2020-07-03 08:49:00
#' TORK_TAPE_F_SPEC    TORK    TAPE   F    SPEC     2020-07-02 14:10:00
#' HARD_BAPE_M_SPEC    HARD    BAPE   M    SPEC     2020-07-03 10:53:00
import_plates = function(imp_xml) {
  # Identify all read plates
  plates = xml_children(imp_xml)

  # Retrieve read information from individual plates
  read.info = as.data.frame(xml_attrs(xml_children(plates)))

  # Format read information for easier handling
  read.info = apply(read.info,2, unlist)
  colnames(read.info) = NULL
  read.info = as.data.frame(t(read.info))

  # Transpose SampleIDs and append to read.info dataframe
  read.info = cbind(parse_IDs(read.info$Name), select(read.info, -c("Name", "InstrumentInfo")))

  # Parse read time for subsequent calculations
  read.info$ReadTime = parse_time(read.info$ReadTime)

  return(as.data.frame(read.info))
}


#' sample_list: generates a simplified sample list containing unique IDs
#' and optional duration between reads.
#'
#'
#' @param imp_xml An internal xml document imported via xml2::read_xml().
#'
#' calc_duration A logical variable indicating whether to return
#' estimated consumption duration. Default calc_duration = FALSE.
#'
#' @keywords consumption time, elapsed time
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' # Get all reads for double checking
#' all_reads = import_plates(imp_xml)
#'
#' print(all_reads)
#'
#' Output:
#' ID                  Project  Line   Sex  Assay   ReadTime
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    2020-07-02 11:54:00
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    2020-07-03 08:49:00
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    2020-07-02 14:10:00
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    2020-07-03 10:53:00
#'
#' # Calculate consumption times from matching IDs
#' consumption_times = sample_list(imp_xml, calc_duration = TRUE)
#'
#' print(consumption_times)
#'
#' ID                  Project  Line   Sex  Assay   Duration(hr)
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    20.9166666666667
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    20.7166666666667
#'
sample_list = function(imp_xml, calc_duration = FALSE) {
  # Import data from internal xml document
  all_reads = import_plates(imp_xml)

  # Identify rows with matching sample IDs
  # Create dataframe to store time calculations for each ID
  summary.data = data.frame(ID = unique(all_reads$ID))

  if (calc_duration == TRUE) {
    comp_df = con_time(all_reads, summary.data)
  } else {
    comp_df = parse_IDs(summary.data$ID)
  }

  return(comp_df)
}


#' con_time: estimates the consumption time.
#'
#'
#' @param
#'
#' @keywords consumption time, elapsed time
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' # Get all reads for double checking
#' all_reads = import_plates(imp_xml)
#'
#' print(all_reads)
#'
#' Output:
#' ID                  Project  Line   Sex  Assay   ReadTime
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    2020-07-02 11:54:00
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    2020-07-03 08:49:00
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    2020-07-02 14:10:00
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    2020-07-03 10:53:00
#'
#' # Calculate consumption times from matching IDs
#' consumption_times = con_time(imp_xml)
#'
#' print(consumption_times)
#'
#' ID                  Project  Line   Sex  Assay   Duration(hr)
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    20.9166666666667
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    20.7166666666667
#'
con_time = function(all_reads, summary.data) {
  # Expand the logical vectors to make them pairwise comparisons
  exp_matr = get_logical_comps(all_reads, summary.data)

  # Calculate the time interval for matching IDs in logical matrix
  for (j in seq(8,ncol(exp_matr))){
    pairs = exp_matr[,j]

    temp = exp_matr[pairs,]

    temp = arrange(temp, desc(ReadTime))

    temp_time = lubridate::interval(start = temp[2,"ReadTime"], end = temp[1,"ReadTime"])

    temp_duration = lubridate::int_length(temp_time)/3600

    temp_vec = data.frame(exp_matr[min(which(pairs == TRUE)),1:6],
                          temp_duration)

    temp_vec$ReadNo = paste(exp_matr[pairs,"ReadNo"],collapse = " vs ")

    colnames(temp_vec) = c("ID", "Project", "Line", "Sex", "Assay", "Comparison", "Duration(hr)")

    if (j == 8) {
      comp_df = temp_vec
    } else {
      comp_df = rbind(comp_df, temp_vec)
    }
  }


  return(comp_df)
}


#' get_logical_comps:
#'
#'
#' @param
#'
#' @keywords
#'
#' @export
#'
#' @examples
get_logical_comps = function(all_reads, summary.data) {
  # Generate logical matrix showing which rows share the same ID
  indices = sapply(all_reads$ID,function(x) {x == summary.data[,"ID"]})
  colnames(indices) = NULL
  indices = t(indices)

  read_ind = replicate(nrow(indices),0)

  read_ind = as.vector(apply(indices, 2, function(x) { read_ind[x] = seq(1, sum(x))}))

  # Split logical vector into independent vectors each with only
  # a single value
  # Iterate over each full logical vector
  for (i in seq(1,ncol(indices))){
    vec = indices[,i]
    init = TRUE

    while (sum(vec) > 0){
      new_vec = logical(length(vec))

      loc = min(which(vec == TRUE))
      vec[loc] = FALSE
      new_vec[loc] = TRUE

      if (init){
        temp_matr = as.matrix(new_vec)
        init = FALSE
      } else {
        temp_matr = cbind(temp_matr, as.matrix(new_vec))
      }
    }

    if (ncol(temp_matr) < 2){
      temp_pairs = temp_matr
    } else{
      pairs = combn(ncol(temp_matr),2)
      temp_pairs = apply(pairs, 2, function(x) {rowSums(temp_matr[,x])})
    }


    if (i == 1) {
      exp_matr = as.matrix(temp_pairs)
    } else {
      exp_matr = cbind(exp_matr, temp_pairs)
    }

  }

  exp_matr = sapply(as.data.frame(exp_matr), as.logical)

  all_reads = add_column(all_reads,ReadNo = read_ind, .after = "Assay")

  exp_matr = cbind(all_reads, exp_matr)

  return(exp_matr)
}



#' make_1536: re-formats the plates from 384 with 4 sub-areas to 1536
#'
#'
#' @param imp_xml An internal xml document imported via xml2::read_xml().
#' plate_node The location of the node containing the plate reads that
#' are being re-formatted. Note: For clarification concerning navigating
#' xml documents, search for XPATH tutorials.
#'
#' @keywords format, 1536, 384
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#' plate_node = "/Experiment/PlateSections[2]/PlateSection"
#'
#' # Get a 32x48 matrix of absorbance readings
#' abs_1536 = make_1536(imp_xml, plate_node)
#'
#' print(abs_1536)
#'
#' Output:
#'      [,1]   [,2]   [,3]   [,4]   ...[,48]
#' [1,] 0.2774 0.1596 0.1390 0.4679
#' [2,] 0.1667 0.1481 0.3128 0.1542
#' [3,] 0.1526 0.2950 2.2933 0.1468
#' [4,] 0.1509 0.1556 0.1655 0.2997
#' ...
#' [32,]
make_1536 = function(imp_xml,plate_node) {
  # Create empty matrix to hold absorbance readings
  abs = matrix(data = NA, nrow = 0, ncol = 48)

  # Iterate over the plate by rows
  for (i in 1:16) {
    # Iterate over each row by columns
    for (j in 1:24) {
      # Reinitialize temp_row variable at the beginning of each loop
      if(j==1) {temp_row = matrix(data = NA, nrow = 2, ncol = 0)}

      # Read values from the "RawData" node
      temp = xml_text(xml_find_first(imp_xml,str_glue(plate_node,"/Wavelengths/Wavelength/Wells/Well[{j+(24*(i-1))}]/RawData")))

      # Process data from node (split string, convert to numeric, unlist)
      temp = strsplit(temp, split = " ")
      temp = unlist(lapply(temp, as.numeric))

      # Shape readings into a 2x2 matrix
      temp_matrix = matrix(data = temp, nrow = 2, ncol = 2, byrow = TRUE)

      # Append readings via cbind to the temp_row variable
      temp_row = cbind(temp_row,temp_matrix)
    }

    # Append temp_row data to abs variable to build the matrix conting absorbance readings
    abs = rbind(abs, temp_row)
  }
  return(abs)
}


#' read_plates: creates a 3 dimensional array containing the plate reads
#' for each sample in a [plate row, plate column, plate ID] format
#'
#'
#' @param imp_xml An internal xml document imported via xml2::read_xml().
#' cutoff The absorbance threshold for identifying wells of interest.
#' Full wells tend to have absorbances of ~2.2. Default cutoff = 1.0.
#' blank A logical variable indicating whether to substract the
#' background absorbance from each reading. Default = TRUE.
#'
#' @keywords read, absorbance, cutoff, threshold, blank
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' abs_reads = read_plates(imp_xml)
#'
#' print(abs_reads)
#'
#' Output:
#'      [,1]   [,2]   [,3]   [,4]   ...[,12]
#' [1,] 0.2774 0.1596 0.1390 0.4679
#' [2,] 0.1667 0.1481 0.3128 0.1542
#' [3,] 0.1526 0.2950 2.2933 0.1468
#' [4,] 0.1509 0.1556 0.1655 0.2997
#' ...
#' [8,]
read_plates = function(imp_xml, cutoff = 1, blank = TRUE) {
  # Generate consumption dataframe using sample_list function and pass
  # the sample IDs to the names variable
  names = sample_list(imp_xml)$ID

  ## Find the nodes in the xml document for each sample ID
  # Use sapply to return list of nodes matching IDs in names variable
  nodes = lapply(
    names,
    # Function substitutes each ID into the xpath string
    # at {x} using str_glue, and searches (xml_find_all()) the data
    # variable for any plates (PlateSection) nodes with
    # the attribute "Name" matching the ID. The node paths
    # are returned using xml_path().
    function(x) {
      xml_path(
        xml_find_all(
          imp_xml,
          str_glue("//PlateSection[@Name=\"{x}\"]")
        )
      )
    }
  )



  # Iterate over the number of unique IDs
  for (i in seq(1, length(names))){

    # Iterate over the nodesets for each ID indicated by i
    for (j in seq(1,length(nodes[[i]]))){

      # Retrieve and reorganize plate at the node using make_1536
      abs = make_1536(imp_xml,nodes[[i]][j])

      # Find wells and blank the plate before binding
      abs = find_wells(abs, cutoff = cutoff, blank = blank)

      # Initialize or add to the array holding the plates and the indices
      if (i == 1 & j == 1){
        sample_plates = array(abs, dim = c(nrow(abs), ncol(abs), 1))
        mat_index = paste(names[i], j)
      } else {
        sample_plates = abind(sample_plates,abs, along=3)
        mat_index = rbind(mat_index,paste(names[i], j))
      }
    }
  }

  # Set the indices for the matrices in the sample plates array
  dimnames(sample_plates)[[3]] = mat_index
  return(sample_plates)
}


#' find_wells: blank the plates and identify wells of interest based
#' on the specified absorbance cutoff
#'
#' @param readings A matrix containing the absorbance readings.
#' cutoff The absorbance threshold for identifying wells of interest.
#' Full wells tend to have absorbances of ~2.2. Default cutoff = 1.0.
#' bin A logical variable indicating whether to zero the wells below
#' the threshold value. Default bin = TRUE.
#' blank A logical variable indicating whether to substract the
#' background absorbance from each reading. Default = TRUE.
#'
#' @keywords read, absorbance, cutoff, threshold, blank
#'
#' @export
find_wells = function(readings,cutoff = 1.0, bin = TRUE, blank=TRUE) {
  # Creates logical matrix: TRUE = above cutoff, FALSE = below cutoff
  well_loc = readings > cutoff

  # Calculate the background as the median of readings below the cutoff
  if (blank == TRUE) {readings = readings - median(readings[!well_loc])}

  # Sets wells below cutoff to zero based on bin argument
  if (bin == TRUE) {readings[!well_loc] = 0}

  # Sums across rows/columns, then converts vector to logical
  # Gives logical vector of rows/columns containing readings (1) or not containing readings (0)
  keep.row = as.logical(apply(well_loc, 1, sum))
  keep.col = as.logical(apply(well_loc, 2, sum))

  # Use logical vectors to keep only those rows/columns containing values
  readings = readings[keep.row,]
  readings = readings[,keep.col]

  return(readings)
}


#' get_stats: apply a standard set of statistical tests to the data
#'
#' @param con_list A list of consumption values belonging to a single
#' group for analysis.
#' remove_outliers A logical variables specifying whether to remove
#' outliers from the data via IQR.
#'
#' @keywords statistics, normality, descriptive
#'
#' @export
get_stats = function(x, remove_outliers = FALSE) {
  if (remove_outliers == FALSE) {
    norm_test = shapiro.test(x)
    ttest = t.test(x, mu = 0, alternative = "two.sided")
    q = quantile(x)
    upper = q[4]+1.5*IQR(x)
    lower = q[2]-1.5*IQR(x)
    outlier_range = paste0("(",round(lower,3), ")-(", round(upper,3), ")")
    outlier_loc = x > upper | x < lower
    outlier_count = length(x[outlier_loc])
    inlier_count = length(x[!outlier_loc])
    stats = rbind(median(x), mean(x), var(x), sd(x), sd(x)/sqrt(length(x)), sd(x)/mean(x), norm_test$p.value, ttest$p.value, outlier_range, paste0(outlier_count,"/",inlier_count))
<<<<<<< HEAD
    rownames(stats) = c("Median", "Mean", "Var", "SD", "SEM", "Coef. of Variation", "Normality (p>0.05)", "t test vs 0", "Outlier Range (1.5IQR)", "Outliers (Out/In)")
=======
    rownames(stats) = c("Median", "Mean", "Var", "SD", "SEM", "Coef. of Variation", "Normality (p>0.05)", "t test vs 0", "Outlier Range (± 1.5IQR)", "Outliers (Out/In)")
>>>>>>> 0392cae8185bf4674cd8bf5d234f028709473a20
    colnames(stats) = deparse(substitute(x))
  } else if (remove_outliers == TRUE) {
    q = quantile(x)
    upper = q[4]+1.5*IQR(x)
    lower = q[2]-1.5*IQR(x)
    outlier_range = paste0("(",round(lower,3), ")-(", round(upper,3), ")")
    outlier_loc = x > upper | x < lower
    outlier_count = length(x[outlier_loc])
    inlier_count = length(x[!outlier_loc])

    y = x[!outlier_loc]

    norm_test = shapiro.test(y)
    ttest = t.test(y, mu = 0, alternative = "two.sided")
    stats = rbind(median(y), mean(y), var(y), sd(y), sd(y)/sqrt(length(y)), sd(y)/mean(y), norm_test$p.value, ttest$p.value, outlier_range, paste0(outlier_count,"/",inlier_count))
<<<<<<< HEAD
    rownames(stats) = c("Median", "Mean", "Var", "SD", "SEM", "Coef. of Variation", "Normality (p>0.05)", "t test vs 0", "Outlier Range (? 1.5IQR)", "Outliers (Out/In)")
=======
    rownames(stats) = c("Median", "Mean", "Var", "SD", "SEM", "Coef. of Variation", "Normality (p>0.05)", "t test vs 0", "Outlier Range (± 1.5IQR)", "Outliers (Out/In)")
>>>>>>> 0392cae8185bf4674cd8bf5d234f028709473a20
    colnames(stats) = deparse(substitute(x))
  }
  return(stats)
}


#' consumption_from_array: calculate consumption from pairs of plates
#'
#'
#' @param plate_array The 3 dimensional matrix of absorbances generated
#' by the read_plates function.
#'
#' vol The initial volume dispensed into the wells of the microplate.
#' Default vol = 10.
#'
#' matrix_type A character variable indicating the type of built-in
#' layout used to organize the 96-well culture plate. Default = "cols12".
#' See collate_samples() function for details about layouts.
#'
#' For custom matrix layouts, matrix_type = "custom".
#' layout_matrix A matrix containing a user-provided layour matrix for
#' sample in the 96-well culture plate. Default layout_matrix = NULL.
#' matrix_type must be matrix_type = "custom".
#'
#' @keywords read, absorbance, layout, consumption
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' abs_reads = read_plates(imp_xml)
#'
#' con_data = consumption_from_array(abs_reads)
#'
#' Output:
#' Line   Sex  Rep    Cons(ÂµL)
#' 0000   M    1      0.11549404
#' 0000   M    2      0.00803428
#' 0000   M    3      0.08090379
#' 0000   M    4      0.19400821
#' 0000   M    5      0.06189788
#' 0000   M    6      -0.07962333
consumption_from_array = function(plate_array, vol = 10, matrix_type = "cols12", layout_matrix = NULL) {
  # Search dimnames for matching IDs
  samIDs = unlist(dimnames(plate_array))
  samIDs = as.data.frame(cbind(samIDs, t(as.data.frame(str_split(samIDs, pattern = " ")))))
  rownames(samIDs) = NULL
  colnames(samIDs) = c("SampleID", "RootID", "Read No")

  # Logical dataframe of plates to compare
  samIDs = cbind(samIDs,sapply(unique(samIDs$RootID), function(x) {x == samIDs$RootID}))

  # Get sampleID from corresponding logical columns
  for (i in seq(4, ncol(samIDs),by=1)){
    # Use sampleID to call the appropriate plates for calculating from plate_array
    con_array = plate_array[,,samIDs[samIDs[,i],1]]

    consumption = consumption_formula(vol, con_array[,,1], con_array[,,2])

    ID = unique(samIDs[samIDs[,i],2])

    temp_data = collate_samples(consumption,matrix_type = matrix_type, layout_matrix = layout_matrix, ID = ID)

    if (i == 4) {
      data_list = temp_data
    } else {
      data_list = rbind(data_list, temp_data)
    }
  }

  colnames(data_list) = c("Line", "Sex", "Rep", "Cons(ÂµL)")

  return(data_list)
}


#' consumption_from_xml: calculate consumption from pairs of plates
#'
#'
#' @param imp_xml An internal xml document imported via xml2::read_xml().
#'
#' vol The initial volume dispensed into the wells of the microplate.
#' Default vol = 10.
#'
#' matrix_type A character variable indicating the type of built-in
#' layout used to organize the 96-well culture plate. Default = "cols12".
#' See collate_samples() function for details about layouts.
#'
#' For custom matrix layouts, matrix_type = "custom".
#' layout_matrix A matrix containing a user-provided layour matrix for
#' sample in the 96-well culture plate. Default layout_matrix = NULL.
#' matrix_type must be matrix_type = "custom".
#'
#' @keywords read, absorbance, layout, consumption
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' abs_reads = read_plates(imp_xml)
#'
#' con_data = consumption_from_array(abs_reads)
#'
#' Output:
#' Line   Sex  Rep    Cons(?L)
#' 0000   M    1      0.11549404
#' 0000   M    2      0.00803428
#' 0000   M    3      0.08090379
#' 0000   M    4      0.19400821
#' 0000   M    5      0.06189788
#' 0000   M    6      -0.07962333
consumption_from_xml = function(imp_xml, Blank = TRUE, vol = 10, matrix_type = "cols12", layout_matrix = NULL) {
  plate_array = read_plates(imp_xml, Blank)

  # Search dimnames for matching IDs
  samIDs = unlist(dimnames(plate_array))
  samIDs = as.data.frame(cbind(samIDs, t(as.data.frame(str_split(samIDs, pattern = " ")))))
  rownames(samIDs) = NULL
  colnames(samIDs) = c("SampleID", "RootID", "Read No")

  # Logical dataframe of plates to compare
  samIDs = cbind(samIDs,sapply(unique(samIDs$RootID), function(x) {x == samIDs$RootID}))

  # Get sampleID from corresponding logical columns
  for (i in seq(4, ncol(samIDs),by=1)){
    # Use sampleID to call the appropriate plates for calculating from plate_array
    con_array = plate_array[,,samIDs[samIDs[,i],1]]

    consumption = consumption_formula(vol, con_array[,,1], con_array[,,2])

    ID = unique(samIDs[samIDs[,i],2])

    temp_data = collate_samples(consumption,matrix_type = matrix_type, layout_matrix = layout_matrix, ID = ID)

    if (i == 4) {
      data_list = temp_data
    } else {
      data_list = rbind(data_list, temp_data)
    }
  }

  colnames(data_list) = c("Line", "Sex", "Rep", "Cons(ÂµL)")

  return(data_list)
}


#' consumption_formula: formula used to calculate consumption
#'
#'
#' @param vol The initial volume dispensed into the wells of the microplate.
#' Default vol = 10.
#'
#' plate1 The pre-exposure plate in the calculation.
#'
#' plate2 The post-exposure plate in the calculation.
#'
#' @keywords formula, calculate
#'
#' @export
consumption_formula = function(vol, plate1, plate2) {
  # plate1 - pre-exposure plate reading
  # plate2 - post-exposure plate reading
  # Positive = decreased absorbance
  # Negative = increased absorbance
  delta_vol = vol*((plate1-plate2)/plate1)
  return(delta_vol)
}


#' collate_samples: takes a matrix of consumption values and re-formats
#' it to give a dataframe in Tidy format containing Line, Group, Rep, and
#' consumption values. Sample replicates are numbered down columns;
#' e.g. A1, B1, C1, D1 -> 1, 2, 3, 4
#'
#'
#' @param con_matrix A matrix containing consumption values
#' in a wellplate format.
#'
#' vol The initial volume dispensed into the wells of the microplate.
#' Default vol = 10.
#'
#' matrix_type A character variable indicating the type of built-in
#' layout used to organize the 96-well culture plate. Default = "cols12".
#' See collate_samples() function for details about layouts.
#'
#' For custom matrix layouts, matrix_type = "custom".
#' layout_matrix A matrix containing a user-provided layour matrix for
#' sample in the 96-well culture plate. Default layout_matrix = NULL.
#' matrix_type must be matrix_type = "custom".
#'
#' ID A character string containing the sample ID for the plate.
#'
#' @keywords
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' abs_reads = read_plates(imp_xml)
#'
#' con_data = consumption_from_array(abs_reads)
#'
#' Output:
#' Line   Group   Rep    Cons
#' 0000   M       1      0.11549404
#' 0000   M       2      0.00803428
#' 0000   M       3      0.08090379
#' 0000   F       1      0.19400821
#' 0000   G       1      0.06189788
#' 0000   G       2     -0.07962330
collate_samples = function(con_matrix,
                           matrix_type = "cols12",
                           layout_matrix = NULL,
                           ID = "PROJ_LINE_S_ASSY") {

  if (matrix_type == "cols12") {
    des_matrix = matrix(c("C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F"),
                        nrow = 8,
                        ncol = 12,
                        byrow = TRUE)
  } else if (matrix_type == "row1") {
    des_matrix = matrix(c("C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
                          "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E",
                          "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F",
                          "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F",
                          "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F",
                          "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F",
                          "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F",
                          "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F"),
                        nrow = 8,
                        ncol = 12,
                        byrow = TRUE)
  } else if (matrix_type == "diag") {
    des_matrix = matrix(c("C", "M", "M", "M", "M", "C", "C", "F", "F", "F", "F", "C",
                          "M", "C", "M", "M", "M", "M", "F", "F", "F", "F", "C", "F",
                          "M", "M", "C", "M", "M", "M", "F", "F", "F", "C", "F", "F",
                          "M", "M", "M", "C", "M", "M", "F", "F", "C", "F", "F", "F",
                          "M", "M", "M", "M", "C", "M", "F", "C", "F", "F", "F", "F",
                          "M", "M", "M", "M", "M", "C", "C", "F", "F", "F", "F", "F",
                          "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F",
                          "M", "C", "M", "M", "M", "M", "F", "F", "F", "F", "C", "F"),
                        nrow = 8,
                        ncol = 12,
                        byrow = TRUE)
  } else if (matrix_type == "custom") {
    des_matrix = layout_matrix
  }

  # Identify unique group identifiers in the design matrix
  groups = unique(as.vector(des_matrix))

  # Isolate the line number from the ID argument
  Line = parse_IDs(ID)$Line

  # Iterate over the groups variable to collate the data
  ## Iterate over the number of groups
  for (i in seq(1,length(groups))) {
    # Generate a logical matric indicating the location of the group
    loc = des_matrix == groups[i]

    # Create or grow the dataframe collecting the collated data
    if (i == 1) {
      collated_df = data.frame(Line = Line,
                               Group = groups[i],
                               Rep = seq(1,length(des_matrix[loc])),
                               Cons = con_matrix[loc])
    } else {
      temp = data.frame(Line = Line,
                        Group = groups[i],
                        Rep = seq(1,length(des_matrix[loc])),
                        Cons = con_matrix[loc])

      collated_df = rbind(collated_df, temp)
    }
  }

  return(collated_df)
}
