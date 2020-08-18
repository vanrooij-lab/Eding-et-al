########################################################################
# Some color palettes

col_Dark2 <- RColorBrewer::brewer.pal(8,name='Dark2')
    # brewer.pal is a function that generates an array of colors

MANY_COLOR_SCHEME <- c('#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabebe', '#469990', '#e6beff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#000075', '#a9a9a9')
    # manual array of colors, obtained from
    # https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/

n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector_60 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


col_YlOrBr<-RColorBrewer::brewer.pal(9,name='YlOrBr')
col_Purples<-RColorBrewer::brewer.pal(9,name='Purples')
col_Spectral<-RColorBrewer::brewer.pal(9,name='Spectral')
col_viridius = viridis_pal(option = "D")(10)
col_viridius_inferno = viridis_pal(option = 'inferno')(100)   
col_viridius_inferno_inv=col_viridius_inferno[seq(length(col_viridius_inferno),1,-1)]
col_blueblack = c('lightblue','black')

TWO_COLOR_SCHEME  <- c('#E5007E','#476AEF') # personal pick

pie(rep(1,length(col_Spectral)),col =    col_Spectral)
pie(rep(1,length(MANY_COLOR_SCHEME)),col =    MANY_COLOR_SCHEME)