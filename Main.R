source('GeneralizedProcrustesAnalysis.R')

  # 25 1 months had a very incomplete segmentation. 14_20 months is missing sulci
  # 25_1month_macaque, 14_20months_human
  OUTLIERS = list(readline(prompt="List all outlier specimens you wanted excluded, separated by commas"))