# =========
# Now also create a histogram with total reads per cell
totalReads_noMito = apply(groupedSCS$patientAllMod@expdata,2,sum)
log10TotalReads = log10(totalReads_noMito+1)
theCutoff = log10(1000+1)
binwidth = max(log10TotalReads)/200
breaks = seq(0,max(log10TotalReads)+binwidth,binwidth)
totalReads_A=log10TotalReads[log10TotalReads<theCutoff]
totalReads_B=log10TotalReads[log10TotalReads>=theCutoff]

p=ggplot(rbind(data.frame(reads=totalReads_A,Status='Excluded'),data.frame(reads=totalReads_B,Status='Included')))+
    geom_histogram(aes(x=reads, fill=Status), breaks=breaks)+theme_bw()+shorthand_tsne_joeptheme()+
    xlab('Log10(Total reads+1)')+ylab('# cells')+give_better_textsize_plot(8)
p
shorthand_save(p, thedir = 'analysis_2020-05-04-18h_17m/FigX/', filename = 'TotalReads_Cutoff_Histogram.pdf', mywidth = 7.5, myheight = 7.5)
