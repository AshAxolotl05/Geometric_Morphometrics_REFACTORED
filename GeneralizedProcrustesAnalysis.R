library('r2r')

#Resamples a single sulcus to a given number of points
resample = function(sulcus, subject, species, hemi, time, df, numpts) {
  temp = df[(df$sulcus == sulcus & df$subject == subject & df$species == species & df$hemisphere == hemi & df$time == time), c('x', 'y','z')]

  # exclude empty sulci
  if (nrow(temp) == 0) {
   stop(paste(subject, time, sulcus, 'DNE', sep=' '))
  }

  # resample curve
  temp = data.frame(resampleCurve(temp, n=200))
  curve = cbind(temp$X1, temp$X2, temp$X3)

  # get starting point
  start = curve[1,]

  temp = data.frame(digit.curves(start, curve, numpts - 2, closed = FALSE)) #98 semi landmarks = start + end

  # reintroduce metadata
  temp$species = species
  temp$subject = subject
  temp$hemisphere = hemi
  temp$time = time
  temp$surface = surface
  temp$sulcus = sulcus

  return(temp)
}

#Handles all resampling for a given surface (White or Gray matter)
resampleCurves = function(df, surface, sulcalLengths=NULL) {
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

  # if a reference hashmap was provided, take note. Otherwise all sulci are resampled to 100 points
  REFERENCE = is.NULL(sulcalLengths)
  numpoints = 100

  for (species in list('macaque', 'human')) {
    for (subject in unique(df[df$species == species, 'subject'])) {
      for (time in unique(df[df$species == species & df$subject == subject, 'time'])) {
        for (hemi in list('lh', 'rh')) {
          for (sulcus in sulci) {
            #handle intitial vs second resampling
            if(REFERENCE) {
              numpoints = sulcalLengths[[hemi]][[sulcus]]
            }

            #resample
            tryCatch({resampled = rbind(resample(sulcus, subject, species, hemi, time, df, numpoints), resampled)},
                     error = function(e) {
                       cat(conditionMessage(e), '\n')
                     })
          }
        }

        if (surface == 'pi') {
          for (sulcus in list('pl', 'ol', 'fl')) {

            #handle intitial vs second resampling
            if(REFERENCE) {
              numpoints = sulcalLengths[[hemi]][[sulcus]]
            }

            #resample
            tryCatch({resampled = rbind(resample(sulcus, subject, species, 'central', time, df, numpoints), resampled)},
                     error = function(e) {
                       cat(conditionMessage(e), '\n')
                     })
          }
        }
      }
    }
  }

  # renamed resampled points
  colnames(resampled)[colnames(resampled) == "X1"] ="x"
  colnames(resampled)[colnames(resampled) == "X2"] ="y"
  colnames(resampled)[colnames(resampled) == "X3"] ="z"

  return(resampled)
}

# Excludes a given list of outliers from the analysis
excludeOutliers= function(outliers){
  filteredDf = resampled
  filteredDf$specimen = paste(filteredDf$subject, filteredDf$time, filteredDf$species, sep='_')

  #remove outliers as defined by plotOutliers
  for (outlier in outliers) {
    filteredDf = subset(filteredDf, filteredDf$specimen != outlier)
  }

  return(filteredDf)
}

# basic euclidean distance function for calculating sulcal lengths
euclidDist = function(x1, y1, z1, x2, y2, z2) {
  return(sqrt((x2-x1)**2 + (y2- y1)**2 + (z2-z1)**2))
}

# calculates sulcal length by summing the euclidean distance between every point along the curve
curveDist = function(i1, i2, referenceMatrix) {
  sum = 0
  for(i in seq(i1, i2 - 1)) {
    sum = sum + euclidDist(referenceMatrix[i,][1], referenceMatrix[i,][2], referenceMatrix[i,][3],
                           referenceMatrix[i+1,][1], referenceMatrix[i+1,][2], referenceMatrix[i+1,][3])
  }
  return(sum)
}

# Create a reference hashmap of sulcal lengths
findSulcalLengths = function (surface, sulci, referenceMatrix) {
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

# Transforms the dataframe into a p x k x n array compatible with geomorph functions
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

# Performs Generalized Procrustes Superimposition, minimizing bending energy and utilizing sliding semi landmarks
procrustes = function(num, data, sulcalLengths=NULL) {
  #handle semilandmarks (num = number of landmarks)
  semiLm = c() # landmark to be slid
  slideStart = c() # previous lm
  slideEnd = c()  # next lm

  # If given no reference for sulcus lengths, there is a "real" lm at the start and end of every 100 pt curve
  if (is.NULL(sulcalLengths)) {
    numSemi = num - (2 * num / 100)

    for (i in 1:num) {
      if (i %% 100 != 0 & i %% 100 != 1) {
        semiLm = c(semiLm, i)
        slideStart = c(slideStart, i - 1)
        slideEnd = c(slideEnd, i + 1)
      }
    }
  } else { # OTHERWISE use sulcalLengths as reference
    hemi = sort(keys())
    sulci = sort(keys(sulcalLengths[['lh']]))





  }

  sliders <- cbind(slideStart, semiLm, slideEnd)

  gpa = gpagen(data, curves = sliders, surfaces = NULL, PrinAxes = TRUE,
   max.iter = NULL, ProcD = FALSE, Proj = TRUE, print.progress = FALSE)

  return(gpa)
}

