library('rgl')
library('geomorph')
library('Morpho')
library('abind')
library('magick')

# save last rgl window as a gif
makeGif = function(name, species, surface) {
  r3dDefaults$windowRect = c(100, 100, 1000, 1000)

  play3d(spin3d(axis = c(0, 0, 1), rpm = 5), duration = 2)
  movie3d(
    movie=paste(name, species, surface, sep='_'),
    spin3d(axis = c(0, 0, 1), rpm = 5), duration = 10,
    type = "gif", clean = T, fps=10, webshot=F, dir='.')
}

# add species and age as factors to gpagen object
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


multivariateAnalysis = function(gdf, gpa, species=FALSE) {
  #effects of age group, species and centroid size on shape

  if(species) {
    fit = procD.lm(gpa$coords~gpa$age_group+gpa$species+log(gpa_aligned$Csize), data=gdf)
    print(anova(fit))

    #look for interaction between species and allometry
    modelInteraction<-procD.lm(gpa$coords~gpa$species*log(gpa$Csize), data = gdf)
    print(anova(modelInteraction))
  } else {
    fit = procD.lm(gpa$coords~gpa$age_group+log(gpa_aligned$Csize), data=gdf)
    print(anova(fit))

    #look for interaction between age and allometry
    modelInteraction<-procD.lm(gpa$coords~gpa$age_group*log(gpa$Csize), data = gdf)
    print(anova(modelInteraction))
  }

  return(modelInteraction)
}

getBySpeciesConsensus = function(gpa) {
  humanMatrix =  gpa$coords[,,grepl('human', dimnames(gpa$coords)[[3]])]
  macaqueMatrix =  gpa$coords[,,grepl('macaque', dimnames(gpa$coords)[[3]])]

  # get consensus (mean configuration)
  humanConsensus = apply(humanMatrix, 1:2, function(coord) { mean(coord) })

  macaqueConsensus = apply(macaqueMatrix, 1:2, function(coord) { mean(coord) })

  return(abind(macaqueConsensus, humanConsensus, along=3))
}

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

deviationAlongPC = function(pca, pc=1) {
  comp = paste0('shapes.comp', pc)

  min = get('min', get(comp, pca$shapes))
  max = get('max', get(comp, pca$shapes))

  plot3d(min, col='blue')
  points3d(max, col='green', add=TRUE)

}
