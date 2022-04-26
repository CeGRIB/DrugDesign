#Title: Molecular docking, simulation and binding free energy analysis
#of small molecules as PfHT1 inhibitors

#Authors: Afolabi J. Owoloye, Funmilayo C. Ligali, Ojochenemi A. Enejoh, Adesola Z. Musa, Oluwagbemiga O. Aina,
          #Emmanuel T. Idowu, and Kolapo M. Oyebola,

#Quality Control Script


library(ggplot2)
library(ggrepel)
#Import dataset
#adjust to your working directory (Unhash below)
#setwd()
#-----starts here -----
RMSF_DF_ligand <- read.table('L_RMSF.csv', header = T, sep = ',')
RMSF_DF_Protein <- read.table('P_RMSF.txt', header = T, sep = '')
water_bridges <- read.table('PL-Contacts_WaterBridge.dat', header = T, sep = '')
pipi <- read.table('PL-Contacts_Pi-Pi.dat', header = T, sep = '')
picat <- read.table('PL-Contacts_Pi-Cation.dat', header = T, sep = '')
hbond <- read.table('PL-Contacts_HBond.dat', header = T, sep = '')
hphobic <- read.table('PL-Contacts_Hydrophobic.dat', header = T, sep = '')

#write to pdf
pdf(file = '4B_DATA_QC.pdf',useDingbats = F, paper = 'a4r', width = 40, height = 40) #Optional point it to a pdf file

#Root mean sqquare functuation_LIGAND
with(RMSF_DF_ligand, plot(x = PDBResName, y = wrt_Protein, pch = 19, col = 'chocolate', cex = 1, ylim= c(0,7)))
with(RMSF_DF_ligand, lines(x = PDBResName, y = wrt_Protein, pch = 19, col = 'darkgoldenrod', lty='dashed', lwd='0.3'))
with(RMSF_DF_ligand, points(x = PDBResName, y = wrt_Ligand, pch = 19, col = 'darkgreen', cex = 1))
with(RMSF_DF_ligand, lines(x = PDBResName, y = wrt_Ligand, pch = 19, col = 'cyan4', lty='dashed', lwd='0.3'))
with(RMSF_DF_ligand, plot(x = wrt_Protein, y = wrt_Ligand, pch = 19, col = 'cyan4', lty='dashed', lwd='0.3'))

#Root mean sqquare functuation_PROTEIN
ggplot(data = RMSF_DF_Protein, aes(label = ResName, y= Backbone, x=Residue))+
  geom_line(color='brown2', lwd=0.1)+
  geom_point(data = subset(RMSF_DF_Protein, LigandContact == 'Yes'), size=0.2, color='grey')+
  geom_text_repel(data = subset(RMSF_DF_Protein, LigandContact == 'Yes'), aes(label=ResName), size=2, color='grey')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  lims(y=c(0.4,3))


#COMPARE
ggplot(data = RMSF_DF_Protein, aes(label = ResName, y= Backbone, x=CA))+
  geom_line(color='deepskyblue4', lwd=0.1)+
  geom_point(data = subset(RMSF_DF_Protein, LigandContact == 'Yes'), size=0.2, color='grey')+
  geom_text_repel(data = subset(RMSF_DF_Protein, LigandContact == 'Yes'), aes(label=ResName), size=2, color='grey')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#merge contacts with protein flunctuation
bridge_residues <- unique(water_bridges$Residue)
pipi_residues <- unique(pipi$Residue)
picat_residues <- unique(picat$Residue)
hbond_residues <- unique(hbond$Residue)
hphobic_residues <- unique(hphobic$Residue)

#merge function
mergeFix <- function(residueList, colID) { #super specific function for this dataset
  outList <- c()
  for (residue in c(RMSF_DF_Protein$Residue)) {
    if (residue %in% residueList) {
      outList <- c(outList, 'Yes')
    }else {
      outList <- c(outList, 'No')
    }
  }
  
  RMSF_DF_Protein[colID] <- outList
  
  return(RMSF_DF_Protein)
}

RMSF_DF_Protein <- mergeFix(bridge_residues, 'WaterBridge')
RMSF_DF_Protein <- mergeFix(pipi_residues, 'PI_PI')
RMSF_DF_Protein <- mergeFix(picat_residues, 'PI_CAT')
RMSF_DF_Protein <- mergeFix(hbond_residues, 'HBOND')
RMSF_DF_Protein <- mergeFix(hphobic_residues, 'HPHOBIC')


#COMPARE with all contacts
ggplot(data = RMSF_DF_Protein, aes(label = ResName, y= B.factor, x=CA))+
  geom_line(color='deepskyblue4', lwd=0.1)+
  geom_point(data = subset(RMSF_DF_Protein, LigandContact == 'Yes'), size=0.2, color='grey')+
  geom_point(data = subset(RMSF_DF_Protein, WaterBridge == 'Yes'), size=0.2, color='chocolate')+
  geom_point(data = subset(RMSF_DF_Protein, PI_PI == 'Yes'), size=0.2, color='darkorange')+
  geom_point(data = subset(RMSF_DF_Protein, PI_CAT == 'Yes'), size=0.2, color='forestgreen')+
  geom_point(data = subset(RMSF_DF_Protein, HBOND == 'Yes'), size=0.2, color='goldenrod1')+
  geom_point(data = subset(RMSF_DF_Protein, HPHOBIC == 'Yes'), size=0.2, color='red')+
  #geom_text_repel(data = subset(RMSF_DF_Protein, LigandContact == 'Yes'), aes(label=ResName), size=2, color='grey')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(data = RMSF_DF_Protein, aes(label = ResName, y= B.factor, x=CA))+
  geom_point(color='deepskyblue4', lwd=0.1)+
  #geom_point(data = subset(RMSF_DF_Protein, LigandContact == 'Yes'), size=0.2, color='grey')+
  geom_point(data = subset(RMSF_DF_Protein, WaterBridge == 'Yes'), size=0.4, color='chocolate')+
  geom_point(data = subset(RMSF_DF_Protein, PI_PI == 'Yes'), size=0.4, color='darkorange')+
  geom_point(data = subset(RMSF_DF_Protein, PI_CAT == 'Yes'), size=0.4, color='forestgreen')+
  geom_point(data = subset(RMSF_DF_Protein, HBOND == 'Yes'), size=0.4, color='yellow')+
  geom_point(data = subset(RMSF_DF_Protein, HPHOBIC == 'Yes'), size=0.4, color='red')+
  #geom_text_repel(data = subset(RMSF_DF_Protein, LigandContact == 'Yes'), aes(label=ResName), size=2, color='grey')+  theme_bw()+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#plot Histograms
hist(RMSF_DF_Protein$B.factor)

dev.off()

#write merged data into dataframe.txt
write.table(RMSF_DF_Protein, 'RMSF_DF_PROTEIN_ALL_INTERACTION.TSV', sep = '\t', row.names = F, quote = T)
