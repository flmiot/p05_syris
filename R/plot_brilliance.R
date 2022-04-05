# plot emission of cdwo4
library(ggplot2)
library(extrafont)
setwd('/home/ubuntu/Dropbox/MA/TEIL02/Spectra/brilliance_curve')
output = '/home/ubuntu/ownCloud/Latex_Masterarbeit/Latex_Masterarbeit/Abbildungen'
bc.01 = read.table('bc_86p5m_maxK_2p12_1-9_128.d01', skip = 10, 
                   col.names = c("Energy",
                                 "K_Value",
                                 "F.Density",
                                 "GA.Brill.",
                                 "Flux",
                                 "Tot.Power",
                                 "P.Density"))
bc.03 = read.table('bc_86p5m_maxK_2p12_1-9_128.d03', skip = 10, 
                   col.names = c("Energy",
                                 "K_Value",
                                 "F.Density",
                                 "GA.Brill.",
                                 "Flux",
                                 "Tot.Power",
                                 "P.Density"))
bc.05 = read.table('bc_86p5m_maxK_2p12_1-9_128.d05', skip = 10, 
                   col.names = c("Energy",
                                 "K_Value",
                                 "F.Density",
                                 "GA.Brill.",
                                 "Flux",
                                 "Tot.Power",
                                 "P.Density"))
bc.07 = read.table('bc_86p5m_maxK_2p12_1-9_128.d07', skip = 10, 
                   col.names = c("Energy",
                                 "K_Value",
                                 "F.Density",
                                 "GA.Brill.",
                                 "Flux",
                                 "Tot.Power",
                                 "P.Density"))
bc.09 = read.table('bc_86p5m_maxK_2p12_1-9_128.d09', skip = 10, 
                   col.names = c("Energy",
                                 "K_Value",
                                 "F.Density",
                                 "GA.Brill.",
                                 "Flux",
                                 "Tot.Power",
                                 "P.Density"))
bc = rbind(bc.01, bc.03, bc.05, bc.07, bc.09)
bc = data.frame(bc, harm = c(rep("1. Harmonische", 128), rep("3. Harmonische", 128), rep("5. Harmonische", 128), rep("7. Harmonische", 128), rep("9. Harmonische", 128)))
bc = bc[bc$Energy < 50000,]

anX = c(min(bc.01$Energy), min(bc.03$Energy), min(bc.05$Energy), min(bc.07$Energy), min(bc.09$Energy)) 
anY = c(bc.01$GA.Brill.[bc.01$Energy == anX[1]],
        bc.03$GA.Brill.[bc.03$Energy == anX[2]],
        bc.05$GA.Brill.[bc.05$Energy == anX[3]],
        bc.07$GA.Brill.[bc.07$Energy == anX[4]],
        bc.09$GA.Brill.[bc.09$Energy == anX[5]]) 
anY = anY + anY * log(2 * anY / anY) + (anY-anY[1]) * 0.2
anX = anX / 1000 + 3
anX = anX + anX * 0.1
ggplot(data=bc, aes(x = Energy / 1000, y = GA.Brill., group = harm))+
  geom_line()+
  theme_classic()+
  theme(axis.line.x = element_line(arrow=arrow(length = unit(3, 'mm'))), 
        axis.line.y=element_line(colour = "black", arrow = arrow(length = unit(3, 'mm'), type = "open")),
        text=element_text(size=20, family="charter"),
        #axis.ticks.y =  element_blank(),
        #axis.ticks.x =  element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = rel(1)),
        axis.text.y = element_text(colour = "black", size = rel(1))
  )+
  annotate('text', label=c('1.', '3.', '5.', '7.', '9.'), x=anX, y=anY, size=5.5,
            family="charter", fontface =2)+
  annotate('text', label='bold("Harmonische")', x=15, y=1e+19, size=5.5,
           family="charter", parse = TRUE) +
  scale_y_log10()+
  xlab("Energie in keV")+
  ylab("Brillianz in ph/s/mrad²/mm²/0.1% ")
setwd(output)
ggsave('brilliance_petra_u2m.pdf', plot = last_plot(), width = 10, height = 7, scale = 1, device=cairo_pdf)



qecmos = read.table('qe_cmos.dat')
qecmos= data.frame(qecmos, Chip = "CMV20000")
qeccd = read.table('qe_ccd.dat')
qeccd = data.frame(qeccd, Chip = "KAF-09000")
qe = rbind(qecmos, qeccd)

ggplot(data=qe, aes(x = V1, y = V2, linetype = Chip))+
  geom_line()+
  theme_classic()+
  theme(axis.line.x = element_line(arrow=arrow(length = unit(3, 'mm'))), 
        axis.line.y=element_line(colour = "black", arrow = arrow(length = unit(3, 'mm'), type = "open")),
        text=element_text(size=20, family="charter"),
        legend.position = c(0.8, 0.8),
        #axis.ticks.y =  element_blank(),
        #axis.ticks.x =  element_blank(),
        #axis.text.y = element_blank(),
        axis.text = element_text(colour = "black", size = rel(0.8))
  )+
  xlab("Wellenlänge in nm")+
  ylab("Quanteneffizienz in %")



ggsave('filename.pdf', plot = last_plot(), width = 10, height = 7, scale = 1, device=cairo_pdf)




plot(em, type = "l")

font_import('/home/ubuntu/Desktop/ttf')


font_import(pattern = 'Charter')
fonts()
