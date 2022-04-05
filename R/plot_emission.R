# plot emission of cdwo4
library(ggplot2)
setwd('/home/ubuntu/Dropbox/MA/TEIL02/R')
output = '/home/ubuntu/ownCloud/Latex_Masterarbeit/Latex_Masterarbeit/Abbildungen'
em = read.table('cdwo4_emission.dat')


library(extrafont)

ggplot(data=em, aes(x = V1, y = V2))+
  geom_line()+
  theme_classic()+
  theme(axis.line.x = element_line(arrow=arrow(length = unit(3, 'mm'))), 
        axis.line.y=element_line(colour = "black", arrow = arrow(length = unit(3, 'mm'), type = "open")),
        text=element_text(size=20, family="charter"),
        axis.ticks.y =  element_blank(),
        axis.ticks.x =  element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = rel(1))
  )+
  xlab("Wellenlänge in nm")+
  ylab("Intensität a.u.")

ggsave('filename.pdf', plot = last_plot(), width = 10, height = 7, scale = 1, device=cairo_pdf)



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
