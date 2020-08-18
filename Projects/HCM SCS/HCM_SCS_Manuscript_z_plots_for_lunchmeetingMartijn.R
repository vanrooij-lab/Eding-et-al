

# tSNE with clusters 
p = shorthand_XY_plot_clusters(groupedSCS$patientAllMod,'Original')
pp = p + shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
shorthand_save(pp, thedir = 'MW_custom_plots/', filename = 'tSNE_f_all_clusters.png', mywidth = 7.5, myheight = 7.5)

# Show were patients are   
lookup=c(patient1=1, patient2=2, patient3=3, patient4=4, patient5=5)
p=plot_sth_in_XY(comp_1 = groupedSCS$patientAllMod@tsne[,1],
             comp_2 = groupedSCS$patientAllMod@tsne[,2],
             color_by = as.factor(lookup[shorthand_give_p_source(config, groupedSCS$patientAllMod)]),
             name_for_color = 'Patient', print_yes = F,mypointsize=.5)+
             theme_bw()+give_better_textsize_plot(15)
pp = p + shorthand_pretty_tsne(mylegendsize=.3) + xlab(element_blank())+ylab(element_blank())+ggtitle(element_blank())
shorthand_save(pp, thedir = 'MW_custom_plots/', filename = 'tSNE_f_all_patients.png', mywidth = 7.5, myheight = 7.5)



# Some very simple illustrations
N=40

fake_expression1 = rnorm(N)-.5
fake_expression2 = .5*(fake_expression1) + .5*(rnorm(N)-.5)
fake_expression3 = rnorm(N)-.5

p=ggplot(data.frame(expression_1=fake_expression1, expression_2 = fake_expression2),
    aes(x=expression_1, y=expression_2))+
    geom_point()+
    shorthand_pretty_tsne(mylegendsize=.3) + xlab('Gene A')+ylab('Gene B')+ggtitle(element_blank())+
    give_better_textsize_plot(15)+
    geom_smooth(method='lm',color='black')
p
shorthand_save(p, thedir = 'MW_custom_plots/', filename = 'simple_example_corrs1.png', mywidth = 5, myheight = 5)
p=ggplot(data.frame(expression_1=fake_expression1, expression_2 = fake_expression3),
    aes(x=expression_1, y=expression_2))+
    geom_point()+
    shorthand_pretty_tsne(mylegendsize=.3) + xlab('Gene A')+ylab('Gene D')+ggtitle(element_blank())+
    give_better_textsize_plot(15)+
    geom_smooth(method='lm',color='black')
p
shorthand_save(p, thedir = 'MW_custom_plots/', filename = 'simple_example_corrs2.png', mywidth = 5, myheight = 5)

# Fake correlation with cell size
fake_expression4 = rnorm(N)-.5
fake_fsca = .5*(fake_expression4) + .5*(rnorm(N)-.5)

p=ggplot(data.frame(expression=fake_expression4, fsca = fake_fsca),
        aes(x=expression, y=fsca))+
    geom_point()+
    shorthand_pretty_tsne(mylegendsize=.3) + xlab('Gene X')+ylab('FSC-A (cell size)')+ggtitle(element_blank())+
    give_better_textsize_plot(15)+
    geom_smooth(method='lm',color='black')
p
shorthand_save(p, thedir = 'MW_custom_plots/', filename = 'simple_example_corrs_fsca.png', mywidth = 5, myheight = 5)


