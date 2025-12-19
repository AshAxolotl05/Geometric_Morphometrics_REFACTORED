library('rgl')
library('geomorph')
library('Morpho')
library('abind')
library('magick')
library('roxygen2')

# Contains functions based around the rgl, Morpho and geomorph packages to perform multiple statistical tests,
# plot visual representations of shape variation, and save results as GIFs

#' Saves last rgl window as a GIF int he current directory.
#'
#' @param name A base name to use for the gif, not including any surface or species information.
#' @param species The species that was plotted.
#' @param surface The surface (white matter or gray matter) that was plotted
makeGif = function(name, species, surface) {

  play3d(spin3d(axis = c(0, 0, 1), rpm = 5), duration = 2)
  movie3d(
    movie=paste(name, species, surface, sep='_'),
    spin3d(axis = c(0, 0, 1), rpm = 5), duration = 10,
    type = "gif", clean = T, fps=10, webshot=F, dir='.')
}

#' Add species and age group (based on developmental milestones for humans and macaques up to Late Childhood) as factors to a gpagen object.
#'
#' @param gpa The gpagen object to modify.
#' @param specimens A list of specimen ids in the format [subject]_[time]_[species] (1_12months_macaque). These must be in the same order as specimens in the gpa.
addFactors = function (gpa, specimens) {
  # get subject number
  subject = c()

  # get species
  species = c()
  # get age in months

  months = c()
  age = c()
  for (specimen in specimens) {
    split = unlist(strsplit(specimen, '_'))

    #handle subject and species
    subject = c(subject, split[1])
    species = c(species, split[3])

    #handle months vs weeks
    if(grepl('month', split[2])) {
      month = unlist(strsplit(split[2], 'month'))
      time = as.numeric(month[1])
      months = c(months, time)
    } else if(grepl('week', split[2])) {
      month = unlist(strsplit(split[2], 'week'))
      time = as.numeric(month[1]) / 4
      months = c(months, time)
    }

    if(split[3] == 'macaque') {
      #  sort into age group
      if(as.integer(time) < 4) {
      group = '(Infant)'
      } else if (as.integer(time) < 12) {
      group = '(Toddler)'
      } else if (as.integer(time) < 24) {
      group = '(Childhood)'
      } else if (as.integer(time) <= 36) {
      group = '(Late Childhood)'
      }
    } else if(split[3] == 'human') {
      #  sort into age group
      if(as.integer(time) < 12) {
      group = '(Infant)'
      } else if (as.integer(time) < 36) {
      group = '(Toddler)'
      } else if (as.integer(time) < 72) {
      group = '(Childhood)'
      } else if (as.integer(time) <= 144) {
      group = '(Late Childhood)'
      }
    }

    age = c(age, group)
  }

  #add factors to R analysis
  gpa$species = as.factor(species)
  gpa$ids = specimens
  gpa$age_group = as.factor(age)
  gpa$subject = as.factor(subject)

  return(gpa)
}

#' Performs Pearson's correlation coefficent tests for the first four principal components and centroid size.
#'
#' @param PCA An object of class gm.prcomp to obtain relevant pc scores from.
#' @param gpa An object of class gpagen to obtain centroid sizes from.
pcToSize = function(PCA, gpa) {

  #look for correlation to centroid size
  pc1=PCA$x[, 1]  #for PC1
  print(cor.test(pc1,gpa$Csize))

  pc2=PCA$x[, 2] #for PC2
  print(cor.test(pc2,gpa$Csize))

  pc3=PCA$x[, 3] #for PC2
  print(cor.test(pc3,gpa$Csize))

  pc4=PCA$x[, 4] #for PC2
  print(cor.test(pc4,gpa$Csize))
}

#' Performs analysis of variance tests to investgate the effects of centroid size, age, and potentially species on overall shape.
#'
#' @param gdf A geomorph data.frame containing GPA results as well as age (and species if species=TRUE) as a factor.
#' @param gpa A gpagen object containing results of GPA with age (and species if species=TRUE) as a factor.
#' @species A boolean describing whether species should be included as a covariate(only necesssary for cross-species analysis).
multivariateAnalysis = function(gdf, gpa, species=FALSE) {
  #effects of age group, species and centroid size on shape

  if(species) {
    fit = procD.lm(gpa$coords~gpa$age_group+gpa$species+log(gpa$Csize), data=gdf)
    print(anova(fit))

    #look for interaction between species and allometry
    modelInteraction<-procD.lm(gpa$coords~gpa$species*log(gpa$Csize), data = gdf)
    print(anova(modelInteraction))
  } else {
    fit = procD.lm(gpa$coords~gpa$age_group+log(gpa$Csize), data=gdf)
    print(anova(fit))

    #look for interaction between age and allometry
    modelInteraction<-procD.lm(gpa$coords~gpa$age_group*log(gpa$Csize), data = gdf)
    print(anova(modelInteraction))
  }

  return(modelInteraction)
}

#' Finds the separate consensus configurations for macaques and humans after a combined GPA.
#' The result is a 3D matrix with the macaque configuration as the first 'slice' (accessible by [,,1] and the human configuration second (accessible by [,,2])
#'
#' @param gpa A gpagen object containing the results of cross-species GPA
getBySpeciesConsensus = function(gpa) {
  humanMatrix =  gpa$coords[,,grepl('human', dimnames(gpa$coords)[[3]])]
  macaqueMatrix =  gpa$coords[,,grepl('macaque', dimnames(gpa$coords)[[3]])]

  # get consensus (mean configuration)
  humanConsensus = apply(humanMatrix, 1:2, function(coord) { mean(coord) })

  macaqueConsensus = apply(macaqueMatrix, 1:2, function(coord) { mean(coord) })

  return(abind(macaqueConsensus, humanConsensus, along=3))
}

#' Creates a GIF visually transforming a configuration of landmarks into another and back.
#'
#' @param config1 The first configuration as a 2D matrix of coordinates, formatted so x, y and optionally z are columns.
#' @param config2 The first configuration as a 2D matrix of coordinates, formatted so x, y and optionally z are columns..
#' @param path The path to save the resulting GIF to. Defaults to the current directory.
transformConfigurations = function(config1, config2, path='.') {

  #transform macaque to human
  warpmovie3d(config1, config2, 15, col = 'green', palindrome = TRUE,
              folder = '.', movie = paste('combined', surface, 'consensusPoints', sep='_'))

  #merge pngs to gif
  list.files(path=path,
             pattern = '*.png', full.names = TRUE) %>%
          image_read() %>% # reads each path file
          image_join() %>% # joins image
          image_animate(fps=4) %>% # animates, can opt for number of loops
          image_write(format='gif', path=paste(paste('combined', surface, 'consensusPoints', sep='_'), '.gif', sep='')) # write to current dir

  #teardown pngs
  list.files(path=path,
             pattern = '*.png', full.names = TRUE) %>% unlink()
}

#' Plots deviation associated with the extreme ends of a single principal compoenent. The minimum value is plotted in blue, the maximum in green.
#'
#' @param pca A gm.prcomp object containing PCA information.
#' @param pc The number (1, 2, 3 etc) of the principal component to look at. Default is 1 for pc1.
deviationAlongPC = function(pca, pc=1) {
  comp = paste0('shapes.comp', pc)

  min = get('min', get(comp, pca$shapes))
  max = get('max', get(comp, pca$shapes))

  plot3d(min, col='blue')
  points3d(max, col='green', add=TRUE)

}
