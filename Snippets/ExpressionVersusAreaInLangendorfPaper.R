# Download supplementary file 6 from this paper: https://www-nature-com.proxy.library.uu.nl/articles/s41467-018-06639-7#additional-information
# Store it somewhere on your computer 
# Rename the file 'area_vs_expression.xlsx'.
# Get the pathname to the file in the below variable
xlsxFileName = '/Users/j.eding/Work/Projects/SCS in HCM/Langendorf Paper/area_vs_expression.xlsx'

# Make a list of genes that your are interested in: (Mouse notation, so first letter capitalized)
# Example: listOfGenesOfInterest = c('Acta1','Actn2','Atp5d','Chchd3','Clgn','Cryab','Etfb','Fabp3','Fundc2','Hrc','Hspb1','Hspb6','Huwe1','Iscu','Kif5b','Mtrnr2l','Myh7','Myl2','Myl7','Ndufa1','Ndufa4','Pdhb','Tax1bp1','Tbc1d4','Tnnt2')
listOfGenesOfInterest = c('Myh6','Nppa')
                         
# Now run the rest of the script. Graphs will be saved to your desktop

# Then run rest of script.
library(ggplot2)
cellArea <- read.xlsx(
  xlsxFile = xlsxFileName,
  rows = c(1,2),
  rowNames = T
)
cellExpression <- read.xlsx(
  xlsxFile = xlsxFileName,
  startRow = 4,
  colNames=T,
  rowNames=T
)

for(geneOfInterest in listOfGenesOfInterest) {
  # Collect data
  dataToPlot = data.frame(
    Area = as.numeric(cellArea),
    Expression = as.numeric(cellExpression[geneOfInterest,])
  )
  
  # Analyze correlation
  correlation = cor.test(dataToPlot$Area, dataToPlot$Expression)
  
  # Plot
  ggplot(
    data=dataToPlot,
    aes(
      x=Area,
      y=Expression
    )
  ) + geom_point(
    
  ) + geom_smooth(
      method = 'lm',
      se=T,
      color='black',
      na.rm=T
  ) + labs(
    title = paste0(geneOfInterest," vs Area"),
    subtitle = paste0("r = ",round(as.numeric(correlation$estimate),digits=3),"               p = ",round(correlation$p.value,digits=5))
  )
  ggsave(paste0(geneOfInterest,'_vs_area.pdf'),path ="~/Desktop/",useDingbats=F)
}
