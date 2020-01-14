rm(list=ls())
par(mar=c(3,5,4,0.5))

#####Fréquences Attendues VS Observées d'acides aminés
data_Ob = read.table('aa_freq.txt', header = T)
data_Ex = read.table('expectedAaFreq.txt', header = T)
dat = rbind(data_Ob$Frequency_..., data_Ex$Frequency_Expected...)
barplot(dat, beside=T, col = c('salmon', 'lightblue'), ylim = c(0,10),
        main = "Fréquences Observées et Attendues d'acides aminés",
        names.arg =  data_Ex$Amino_Acids, border=NA, ylab = 'Fréquences (%)',
        cex.lab = 2, cex.axis = 2, cex.names = 2, cex.main = 2)
legend('topright', fill = c("salmon", "lightblue"),
       legend = c('Observées', 'Attendues'), cex = 1.5)
 
#####Fréquences Attendues VS Observées de structures secondaires
data_Ob = read.table('struct_freq.txt', header = T)
data_Ex = read.table('expectedStructFreq.txt', header = T)
dat = rbind(data_Ob$Frequency_..., data_Ex$Frequency_Expected...)
barplot(dat, beside=T, col = c('salmon', 'lightblue'),
        main = "Fréquences Observées et Attendues de structures secondaires",
        names.arg =  data_Ex$Secondary_Structure, border=NA, ylab = 'Fréquences (%)',
        cex.lab = 2, cex.axis = 2, cex.names = 2, cex.main = 2)
legend('topright', fill = c("salmon", "lightblue"),
       legend = c('Observées', 'Attendues'), cex = 1.5)

#####Fréquences Attendues VS Observées de type de coudes Beta
data_Ob = read.table('type_freq.txt', header = T)
data_Ex = read.table('expectedTypeFreq.txt', header = T)
dat = rbind(data_Ob$Frequency_..., data_Ex$Frequency_Expected...)
barplot(dat, beside=T, col = c('salmon', 'lightblue'),
        main = "Fréquences Observées et Attendues de type de coudes Beta",
        names.arg =  data_Ex$Type, border=NA, ylab = 'Fréquences (%)',
        cex.lab = 2, cex.axis = 2, cex.names = 1.6, cex.main = 2)
legend('topright', fill = c("salmon", "lightblue"),
       legend = c('Observées', 'Attendues'), cex = 1.5)
