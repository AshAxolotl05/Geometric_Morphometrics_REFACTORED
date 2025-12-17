library('rgl')
library('geomorph')
library('Morpho')

# save last rgl window as a gif
makeGif = function(name, surface) {
  r3dDefaults$windowRect = c(100, 100, 1000, 1000)

  play3d(spin3d(axis = c(0, 0, 1), rpm = 5), duration = 2)
  movie3d(
    movie=paste(name, surface, sep='_'),
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


pcToSize = function(PCA) {

  #look for correlation to centroid size
  pc1=PCA$x[, 1]  #for PC1
  print(cor.test(pc1,gpa_aligned$Csize))

  pc2=PCA$x[, 2] #for PC2
  print(cor.test(pc2,gpa_aligned$Csize))

  pc3=PCA$x[, 3] #for PC2
  print(cor.test(pc3,gpa_aligned$Csize))

  pc4=PCA$x[, 4] #for PC2
  print(cor.test(pc4,gpa_aligned$Csize))
}


multivariateAnalysis = function(gpa, species=FALSE) {
  #effects of age group, species and centroid size on shape
  gdf <- geomorph.data.frame(gpa, age=gpa$age_group, species=gpa$species, subject=gpa$subject)

  if(species) {
    fit = procD.lm(gpa$coords~gpa$age_group+gpa$species+log(gpa_aligned$Csize), data=gdf)
    anova(fit)

    #look for interaction between species and allometry
    modelInteraction<-procD.lm(gpa$coords~gpa$species*log(gpa$Csize), data = gdf)
    anova(modelInteraction)
  } else {
    fit = procD.lm(gpa$coords~gpa$age_group+log(gpa_aligned$Csize), data=gdf)
    anova(fit)

    #look for interaction between age and allometry
    modelInteraction<-procD.lm(gpa$coords~gpa$age_group*log(gpa$Csize), data = gdf)
    anova(modelInteraction)
  }

  return(modelInteraction)
}