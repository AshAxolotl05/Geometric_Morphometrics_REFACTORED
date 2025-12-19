library('r2r')
library('Morpho')
library('geomorph')
library('roxygen2')

# A collection of functions based on the Morpho and geomorph packages that handles most of the logic necessary for performing Generalized Procrustes Analysis.
# These are designed to work for a macaque only, human only, or macaque-human cross species analysis, but could be modified to work with other NHP in the future.

#' Makes a hashmap to store the sulci for each species(macaques, humans, and both combined)
#' + surface(gray matter/pial as 'pi' and white matter as 'wm') combination.
#' sulci are referred to by the same abbreviations as in the tracing spreadsheets.
#'
makeSulci = function() {
  map = hashmap()

  pial = hashmap()
  mPi = list('ce', 'arc', 'pr', 'la', 'st', 'lu', 'po', 'of', 'ip', 'di', 'cb', 'io')
  hPi = list('ar', 'cb', 'ce', 'hr', 'if', 'ip', 'it', 'la', 'of', 'pc', 'po', 'pt', 'sf', 'st')
  cPi = list('ce', 'la', 'st', 'po', 'of', 'ip', 'cb')

  pial[['macaque']] = mPi
  pial[['human']] = hPi
  pial[['combined']] = cPi

  white = hashmap()
  mWh = list('ce', 'arc', 'pr', 'la', 'st', 'lu', 'co', 'ci', 'po', 'of', 'ca', 'ip', 'di', 'io')
  hWh = list('ar', 'ca', 'ce', 'ci', 'co', 'hr', 'if', 'ip', 'it', 'la', 'of', 'pc', 'po', 'pt', 'sf', 'st')
  cWh = list('ca', 'ce', 'ci', 'co','ip', 'la', 'of', 'po', 'st')

  white[['macaque']] = mWh
  white[['human']] = hWh
  white[['combined']] = cWh

  map[['pi']] = pial
  map[['wh']] = white

  return(map)
}

#' Returns the sulci that have been traced for a given surface (wm or pi) and species (macaque, human, or both 'combined') combination formatted as a list
#'
#' @param surface A string representing the surface, either 'pi' for gray matter or 'wm' for white matter
#' @param species A string representing the species, either 'macaque', 'human', or 'combined' for sulci traced in both
#'
getSulci = function(surface, species) {
  return(SULCI[[surface]][[species]])
}

#' Resamples a single sulcus to a given number of equidistant points along the curve
#'
#' @param sulcus A string abbreviation for the sulcus to be resampled. Use the same abbreviation as in the tracing spreadsheets.
#' @param subject The subject number (ex 1, 12, 25) to reference. For humans, this information is saved in an excel spreadsheet alongside other ways to reference teh same subject
#' @param species The species (human or macaque) to reference.
#' @param hemi The hemisphere('lh' for left, 'rh' for right, 'central' for the median longitudinal sulcus) to reference.
#' @param time A string representation of the time point to reference. Should be in the format '1month' '2weeks' '12months' etc.
#' @param df The dataframe to use for looking up tracing points. This should have been read in from combined_tracings.csv or something with the same format.
#' @param numpts The number of points to resample the sulcus to.
resample = function(sulcus, subject, species, hemi, time, df, numpts) {
  temp = df[(df$sulcus == sulcus & df$subject == subject & df$species == species & df$hemisphere == hemi & df$time == time), c('x', 'y','z')]

  # exclude empty sulci
  if (nrow(temp) == 0) {
   stop(paste(subject, time, sulcus, 'DNE', sep=' '))
  }

  # resample curve
  temp = data.frame(resampleCurve(temp, n=300))
  curve = cbind(temp$X1, temp$X2, temp$X3)

  # get starting point
  start = curve[1,]

  temp = data.frame(digit.curves(start, curve, (numpts - 2), closed = FALSE)) #semi lm + start + end

  # reintroduce metadata
  temp$species = species
  temp$subject = subject
  temp$hemisphere = hemi
  temp$time = time
  temp$surface = surface
  temp$sulcus = sulcus

  return(temp)
}

#' Handles all resampling for a given surface (White or Gray matter)
#'
#' @param df The dataframe to use for looking up tracing points. This should have been read in from combined_tracings.csv or something with the same format.
#' @param surface The surface ('pi' for gray matter or 'wm' for white matter) to investigate.
#' @param species The species (human, macaque, or combined) to reference.
#' @param sulcalLengths A hashmap or hasmaps formatted so that the first set of keys correspond to surface and the second set of keys correspond to sulcus
#' that stores the number of points to resample each sulcus to. Null by default, in which case every sulcus will be resampled to 100 points.
resampleCurves = function(df, surface, species, sulcalLengths=NULL) {
  #create data frame to store resampled points
  resampled = data.frame(
  x = numeric(),
  y = numeric(),
  z = numeric(),
  species = character(),
  subject = character(),
  surface = character(),
  time = character(),
  sulcus = character(),
  hemisphere = character()
)

  # select the relevant sulci and species
  sulci = getSulci(surface, species)
  if(species == 'combined') {
    speciesList = list('macaque', 'human')
  } else { speciesList = list(species)}

  # if a reference hashmap was provided, take note. Otherwise all sulci are resampled to 100 points
  REFERENCE = !is.null(sulcalLengths)
  numpoints = 100

  for (species in speciesList) {
    for (subject in unique(df[df$species == species, 'subject'])) {
      for (time in unique(df[df$species == species & df$subject == subject, 'time'])) {
        for (hemi in list('lh', 'rh')) {
          for (sulcus in sulci) {
            #handle intitial vs second resampling
            if(REFERENCE) {
              numpoints = sulcalLengths[[hemi]][[sulcus]]
            }

            #resample
            resampled = rbind(resample(sulcus, subject, species, hemi, time, df, numpoints), resampled)

          }
        }

        if (surface == 'pi') {
          for (sulcus in list('pl', 'ol', 'fl')) {

            #handle intitial vs second resampling
            if(REFERENCE) {
              numpoints = sulcalLengths[['central']][[sulcus]]
            }

            #resample
            resampled = rbind(resample(sulcus, subject, species, 'central', time, df, numpoints), resampled)
          }
        }
      }
    }
  }

  # renamed resampled points
  colnames(resampled)[colnames(resampled) == "X1"] ="x"
  colnames(resampled)[colnames(resampled) == "X2"] ="y"
  colnames(resampled)[colnames(resampled) == "X3"] ="z"

  resampled$specimen = paste(resampled$subject, resampled$time, resampled$species, sep='_')
  return(resampled)
}

#' Excludes a given list of outliers from the analysis
#'
#' @param df The dataframe to exclude specimens from. This should have been read in from combined_tracings.csv or something with the same format. It must have a specimen column.
#' @param outliers A list of specimen names in the format [subject]_[time]_[species] (1_12months_macaque) to exclude.
excludeOutliers= function(df, outliers){
  filteredDf = df
  filteredDf$specimen = paste(filteredDf$subject, filteredDf$time, filteredDf$species, sep='_')

  #remove outliers as defined by plotOutliers
  for (outlier in outliers) {
    filteredDf = subset(filteredDf, filteredDf$specimen != outlier)
  }

  return(filteredDf)
}

#' Basic euclidean distance function for calculating distance between two points in 3D
#'
#' @param x1 First x coordinate
#' @param y1 First y coordinate
#' @param z1 First z coordinate
#' @param x2 Second x coordinate
#' @param y2 Second y coordinate
#' @param z2 Second z coordinate
euclidDist = function(x1, y1, z1, x2, y2, z2) {
  return(sqrt((x2-x1)**2 + (y2- y1)**2 + (z2-z1)**2))
}

#' Calculates sulcal length by summing the euclidean distance between every point along the curve.
#'
#' @param i1 Row in the matrix of coordinates that coordesponds to the first landmark of this sulcus
#' @param i2 Row in the matrix of coordinates that coordesponds to the last landmark of this sulcus
#' @param referenceMatrix A 2D matrix of all sulcal coordinates (x, y, z) to use as a reference.
curveDist = function(i1, i2, referenceMatrix) {
  sum = 0
  for(i in seq(i1, i2 - 1)) {
    sum = sum + euclidDist(referenceMatrix[i,][1], referenceMatrix[i,][2], referenceMatrix[i,][3],
                           referenceMatrix[i+1,][1], referenceMatrix[i+1,][2], referenceMatrix[i+1,][3])
  }
  return(sum)
}

#' Create a reference hashmap of sulcal lengths
#'
#' @param surface The string surface ('pi for gray matter or 'wm' for white matter) to investigate.
#' @param species The string species ('human' 'macaque', or 'combined') to investigate.
#' @param referenceMatrix A 2D matrix of all sulcal coordinates (x, y, z) to use as a reference.
findSulcalLengths = function (surface, species, referenceMatrix) {
  sulci = getSulci(surface, species)

  #Make hasmap containing all sulci and their lengths in the consensus configuration
  sulcalLengths = hashmap()

  # Handle pial vs white matter
  if (surface == 'pi') {
    i = 1

    #handle median longitudinal sulcus
    sulcMap = hashmap()
    for (sulcus in list('fl', 'ol', 'pl')) {
      length =  curveDist(i, i+99, referenceMatrix)
      sulcMap[[sulcus]] = length
      i = i + 100
    }
    sulcalLengths[['central']] = sulcMap
    i = 301 # keep track of landmark number
  } else {
    i = 1 # keep track of landmark number
  }

  #handle left and right hemispheres
  for (hemi in list('lh', 'rh')) {
    sulcMap = hashmap()

    for (sulcus in sort(unlist(sulci))) {
      length =  curveDist(i, i+99, referenceMatrix)
      i = i + 100
      sulcMap[[sulcus]] = length
    }
    sulcalLengths[[hemi]] = sulcMap
  }

  #translate length into number of points based on central sulcus length
  ratio = 100 / sulcalLengths[['lh']][['ce']]

  # calculate number of points to resample each sulcus to
  for (hemi in keys(sulcalLengths)) {
    for (sulcus in keys(sulcalLengths[[hemi]])) {
      sulcalLengths[[hemi]][[sulcus]] = round(sulcalLengths[[hemi]][[sulcus]] * ratio)
    }
  }

  return(sulcalLengths)
}

#' Computes the total number of landmarks from a reference hashmap of sulcal lengths.
#'
#' @param lengths A hashmap of hashmaps formatted so that the first set of keys correspond to surface and the second set of keys correspond to sulcus
#' that stores the number of points each  sulcus was resampled to.
sumlm = function(lengths) {
  sum = 0

  #iterate over map
  for(key in keys(lengths)) {
    for (length in keys(lengths[[key]])) {
      sum = sum + lengths[[key]][[length]]
    }
  }

  return(sum)
}

#' Transforms the dataframe into a p x k x n array compatible with geomorph functions,
#' where p is the number of landmarks, k is the number of dimensions (2 or 3), and n is the number of specimens
#'
#' @param df A dataframe containing all landmark data properly labeled with sulcus, surface, time, subject, hemisphere and subject.
#' It must have a specimen column in the format [subject]_[time]_[species] (1_12months_macaque).
#' @param num Total number of landmarks.
transformToArray = function(df, num) {
  # order df; this is crucial to measuring sulcal lengths later
  df = df[order(df$subject, df$time, df$species, df$hemisphere, df$sulcus),]

  #format array
  data = arrayspecs(cbind(df$x, df$y, df$z), num, 3)

  # add specimen ids
  ids = c()
  species = c()
  subject = c()

  for(i in seq(1, nrow(df), by = num)) {
    ids = c(ids, df[i, 'specimen'])
    species = c(species, df[i, 'species'])
    subject = c(subject, paste(df[i, 'subject'], df[i, 'species'], sep='_'))
  }

  dimnames(data) = list(
    Landmark = 1:num ,  # names for landmarks
    Dimension = c("X", "Y", "Z"),       # names for dimensions
    Specimen = ids      # names for specimens
  )


  return(data)
}

#' Performs Generalized Procrustes Superimposition, minimizing bending energy and utilizing sliding semi landmarks.
#' Returns an object of class gpagen from the geomorph package.
#'
#' @param num The number of landmarks.
#' @param data A p x k x n array compatible with geomorph functions, where p is the number of landmarks, k is the number of dimensions (2 or 3), and n is the number of specimens
#' @param sulcalLengths A hashmap or hasmaps formatted so that the first set of keys correspond to surface and the second set of keys correspond to sulcus
#' that stores the number of points to resample each sulcus to. Null by default, in which case it assumes every sulcus was resampled to 100 points.
procrustes = function(num, data, sulcalLengths=NULL) {
  #handle semilandmarks (num = number of landmarks)
  semiLm = c() # landmark to be slid
  slideStart = c() # previous lm
  slideEnd = c()  # next lm

  # If given no reference for sulcus lengths, there is a "real" lm at the start and end of every 100 pt curve
  if (is.null(sulcalLengths)) {
    #there should be a 'real' lm every 100 lms
    for (i in 1:num) {
      if (i %% 100 != 0 & i %% 100 != 1) {
        semiLm = c(semiLm, i)
        slideStart = c(slideStart, i - 1)
        slideEnd = c(slideEnd, i + 1)
      }
    }
  } else { # OTHERWISE use sulcalLengths as reference
    # get a list of 'real' landmarks based on lengths
    reallm = c()
    i = 0

    for(hemi in sort(unlist(keys(sulcalLengths)))) {
      for (sulcus in sort(unlist(keys(sulcalLengths[[hemi]])))) {
        reallm = c(reallm, i + 1, i + sulcalLengths[[hemi]][[sulcus]])
        i = i + sulcalLengths[[hemi]][[sulcus]]
      }
    }

    # slide all landmarks that are not 'real' landmarks between the previous and next landmark
    for (i in 1:num) {
          if (!(i %in% reallm)) {
            semiLm = c(semiLm, i)
            slideStart = c(slideStart, i - 1)
            slideEnd = c(slideEnd, i + 1)
        }
    }
  }

  sliders <- cbind(slideStart, semiLm, slideEnd)
  gpa = gpagen(data, curves = sliders, surfaces = NULL, PrinAxes = TRUE, max.iter = NULL, ProcD = FALSE, Proj = TRUE, print.progress = FALSE)
  return(gpa)
}

SULCI = makeSulci()