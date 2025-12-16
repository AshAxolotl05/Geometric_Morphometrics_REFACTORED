source('GeneralizedProcrustesAnalysis.R')

df = read.csv("combined_tracings.csv")

species = readline(prompt="Will this be a human, macaque, or combined analysis?")
surface = readline(prompt="Chose gray matter (pi) or white matter (wh)")



  # 25 1 months had a very incomplete segmentation. 14_20 months is missing sulci
  # 25_1month_macaque, 14_20months_human
  OUTLIERS = list(readline(prompt="List all outlier specimens you want excluded, separated by commas"))