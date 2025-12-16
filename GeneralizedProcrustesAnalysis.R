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

# Transforms the dataframe into a p x k x n array compatible with geomorph functions
transformToArray = function(df, num) {
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
procrustes = function() {}

