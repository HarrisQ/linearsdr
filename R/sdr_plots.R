#######################################################################
#           ggplot Functions for OPCG and MADE
#######################################################################

ggplot_fsdr <- function(y_datta, x_datta, y_on_axis=F, ytype="multinomial",
                        size=1, h_lim=NULL, v_lim=NULL, 
                        h_lab=NULL, v_lab=NULL, main_lab=NULL, 
                        show_legend=T, 
                        y_colors=NULL, y_symbols=NULL,
                        ellipse=F
                        ) {
  
  # y_datta=Y; x_datta=t( B_hat_opcg )%*%(X);
  # y_on_axis=F; ytype="continuous";#"multinomial";
  # size=1; h_lab=NULL; v_lab=NULL; main_lab=NULL;
  # show_legend=T; h_lim=NULL; v_lin=NULL;
  # y_colors=NULL; y_symbols=NULL
  # ellipse=T;
  
  datta_frame0 <- data.frame( t( rbind(y_datta, x_datta) ));
  colnames(datta_frame0) <- c('y', 
                              sapply(1:(dim(datta_frame0)[2]-1), 
                                        function(k) paste0('x',k) )); 
    
  
  p_base = ggplot2::ggplot(datta_frame0);
  
  if(!is.null(h_lim) ) p_base = p_base + xlim(h_lim[1], h_lim[2]); 
  if(!is.null(v_lim) ) p_base = p_base + ylim(v_lim[1], v_lim[2]);  
  
  if( ytype=="multinomial" ) {
    # Colours for discrete Y
    if(!is.null(y_colors)) {
      p_base = p_base + ggplot2::scale_colour_manual(values = y_colors)
    }
    # Shapes for Discrete Y
    if(!is.null(y_symbols)) {
      p_base = p_base + ggplot2::scale_shape_manual(values = y_symbols)
    }
    
    # Add points
    if(!y_on_axis) { # Don't plot Y
      
      # Draw Ellipse
      if(ellipse) {
        pplot = p_base + 
          ggplot2::stat_ellipse(aes(x=x1, y=x2, color = factor(y), group=factor(y)),
                                type="norm", lwd=2, lty=2) #
      } else {
        pplot = p_base + 
          ggplot2::geom_point(aes(x=x1, y=x2, color = factor(y), shape=factor(y)),
                              size=size, show.legend=show_legend );
      }
      
    } else {
      # Plot Y
      pplot = p_base + 
        ggplot2::geom_point(aes(x=x1, y=y, color = factor(y), shape=factor(y)),
                            size=size, show.legend=show_legend  )
    }
    
    
    
  } else if (ytype=="continuous") {
    
    # Colours for continuous Y
    if(!is.null(y_colors)) {
      p_base = p_base + ggplot2::scale_colour_manual(values = y_colors)
    }
    # Shapes for continuous Y
    if(!is.null(y_symbols)) {
      p_base = p_base + ggplot2::scale_shape_manual(values = y_symbols)
    }
    
    
    
    # Add points
    if(!y_on_axis) { # Don't plot Y
      
      # Draw Ellipse
      if(ellipse) {
        pplot = p_base + 
          ggplot2::stat_ellipse(aes(x=x1, y=x2, color = factor(y), group=factor(y)),
                                type="norm")
      } else {
        pplot = p_base +
          ggplot2::geom_point(aes(x=x1, y=x2, color = y),
                              size=size, show.legend=show_legend  )
      }
    } else { # Plot Y
      
      # Draw Ellipse
      if(ellipse) {
        pplot = p_base + 
          ggplot2::stat_ellipse(aes(x=x1, y=y, color = factor(y), group=factor(y)),
                                type="norm")
      } else {
        pplot = p_base +
          ggplot2::geom_point(aes(x=x1, y=y, color = y),
                              size=size, show.legend=show_legend  )
      }
    }
    
    
    
  }

  
  pplot + ggplot2::labs(title = main_lab, x = h_lab, y = v_lab ) +
    ggplot2::theme(legend.position="none",
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black")) 
}

### Tuning plots

tuning_plot_km = function(h_list, tuning_list, x_lab=NULL,y_lab=NULL) {
  ggplot(data.frame(x=h_list, y=as.numeric(tuning_list))  ) +
    geom_line( aes(x , y ), size=1) +
    labs(x = x_lab, y=y_lab ) +
    theme(legend.position="none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
   
### 3D Plot


ggplot3D_fsdr = function(alpha, x_datta, y_datta, ytype="multinomial", 
                         y_color, y_symbol,
                         sdr_method='', size=2, show.legend=F, label_size=5,
                         image=NULL,image_size=0.05) {

    
  # alpha=125; x_datta=(nlopcg_fit$pred_test[1:3,]);
  # y_datta=1:n; ytype="multinomial"; 
  # y_color=clas_col; y_symbol=clas_symb;
  # sdr_method='OPCG'; size=2; show.legend=F; label_size=5;
  # image=image_vec
  
  # x_datta=t(B_hat_opcg)%*%(X); y_datta=Y;
  # y_color=clas_col; y_symbol=clas_symb; sdr_method = 'OPCG';
  # size=3; label_size = 8; show.legend=F;
  
  
  # Internal Functions
  
  # Cabinet Projection, see Wikipedia on 3D projection and Oblique projection
  proj_3to2d <- function(alpha) {
    P <- as.matrix( 
      rbind(c(1,0,.5*cos(alpha)),
            c(0,1,.5*sin(alpha)),
            c(0,0,0)) 
    )
    return(P)
  }
  
  map3to2d <- function(df,alpha) {
    # Takes in a dattaframe with 3 columns
    df_2d=as.data.frame( t(proj_3to2d(alpha)%*%t(df)))[,1:2]
    colnames(df_2d) <- c('x','y')
    return(df_2d)
  }
  
  # Creating Blank Cube plot ----------
  
  # Define the corners of the cube for perspective alpha
  vertices <- function(alpha) {
    corners <- expand.grid(c(-1,1), c(-1,1), c(-1,1) )
    V <- as.data.frame( t(proj_3to2d(alpha)%*%t(corners)))[,1:2]
    # Z <- c(1,2,3,4,1,2,3,4) # Z <- rowSums(sign(V)) but with no diag 
    Z <- c(1,2,3,4,1,2,3,4) 
    cube <- data.frame(V,Z)
    colnames(cube) <- c('x','y','group')
    return(cube)
  }
  
  # vertices(alpha)
  
  p_blank=ggplot2::ggplot(aes(x=x, y=y), data = data.frame( vertices(alpha) ) )+
    ggplot2::geom_segment(aes(x = x[6], y = y[6], xend = x[2], yend = y[2] ), 
                          color='black', linetype=3  )+
    ggplot2::geom_segment(aes(x = x[6], y = y[6], xend = x[5], yend = y[5] ), 
                          color='black'  )+
    ggplot2::geom_segment(aes(x = x[6], y = y[6], xend = x[8], yend = y[8] ), 
                          color='black'  )+
    #
    ggplot2::geom_segment(aes(x = x[1], y = y[1], xend = x[3], yend = y[3] ), 
                          color='black'  )+
    ggplot2::geom_segment(aes(x = x[1], y = y[1], xend = x[5], yend = y[5] ), 
                          color='black'  )+
    ggplot2::geom_segment(aes(x = x[1], y = y[1], xend = x[2], yend = y[2] ), 
                          color='black', linetype=3  )+
    #
    ggplot2::geom_segment(aes(x = x[4], y = y[4], xend = x[3], yend = y[3] ), 
                          color='black'  )+
    ggplot2::geom_segment(aes(x = x[4], y = y[4], xend = x[8], yend = y[8] ), 
                          color='black'  )+
    ggplot2::geom_segment(aes(x = x[4], y = y[4], xend = x[2], yend = y[2] ), 
                          color='black', linetype=3  )+
    # geom_point( color='white' ) + # )+#
    # geom_line( aes(group=y ), color='black' ) +
    # geom_line( aes(group=x ), color='black') +
    # geom_line( aes(group=group ), color='black')+
    ggplot2::geom_text(aes(x = x[1], y = y[1], label = paste(sdr_method,'2') ), size=label_size,
              # data=surface.2d.labels, 
              nudge_x = -.25, nudge_y = 1 )  +
    ggplot2::geom_text(aes(x = x[5], y = y[5], label = paste(sdr_method,'3') ), size=label_size,
              # data=surface.2d.labels, 
              nudge_x = +.15, nudge_y = .25 )  +
    ggplot2::geom_text(aes(x = x[6], y = y[6], label = paste(sdr_method,'1') ), size=label_size,
              # data=surface.2d.labels, 
              nudge_x = -1, nudge_y = -.1 )  +
    theme_void() +
    # labs(x='x', y='y')+
    ggplot2::theme(legend.position="none",
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   legend.background = element_rect(fill = 'transparent')) 
  # p_blank
  # 
  
  
  
  if(ytype=="multinomial"){
    tmp_datta=data.frame( apply( map3to2d(t( x_datta ), alpha), 2, 
                                 function(vec) (vec-mean(vec))/(.55*max(abs(vec))) ), 
                          label=c(y_datta) )
    tmp_plot =
      p_blank +  
      ggplot2::geom_point(data = tmp_datta , aes(y = y,x = x, 
                                                 color = factor(label), shape=factor(label)), 
                          size=size, show.legend=show.legend) +
      ggplot2::scale_colour_manual(values = y_color)+
      ggplot2::scale_shape_manual(values = y_symbol)+
      theme_void() +
      # labs(x='x', y='y')+
      ggplot2::theme(legend.position="none",
                     plot.background = element_blank(),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_blank(),#element_line(colour = "black"),
                     legend.background = element_rect(fill = 'transparent')) 
  
  } else if (ytype=="continuous") {
    tmp_datta=data.frame( apply( map3to2d(t( x_datta ), alpha), 2, 
                                 function(vec) (vec-mean(vec))/(.55*max(abs(vec))) ), 
                          label=c(y_datta), 
                          image=image)
   
    if(!is.null(image)){
      tmp_plot =
        p_blank +  
        ggplot2::geom_point(data = tmp_datta , aes(y = y,x = x), 
                            size=size, show.legend=show.legend)  +
        ggplot2::scale_colour_manual(values = y_color)+
        ggplot2::scale_shape_manual(values = y_symbol)+
        theme_void() +
        # labs(x='x', y='y')+
        ggplot2::theme(legend.position="none",
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(), 
                       axis.line = element_blank(),#element_line(colour = "black"),
                       legend.background = element_rect(fill = 'transparent')) +
        ggimage::geom_image(data=tmp_datta, aes(image=image), size=image_size)
      
    } else {
      tmp_plot = 
        p_blank +  
        ggplot2::geom_point(data = data.frame( apply( map3to2d(t( x_datta ), alpha), 2, 
                                                      function(vec) (vec-mean(vec))/(.55*max(abs(vec))) ), 
                                               label=t(y_datta) ), 
                            aes(y = y,x = x), 
                            size=size, show.legend=show.legend) +
        ggplot2::scale_colour_manual(values = y_color)+
        ggplot2::scale_shape_manual(values = y_symbol)+
        theme_void() +
        # labs(x='x', y='y')+
        ggplot2::theme(legend.position="none",
                       plot.background = element_blank(),
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(), 
                       axis.line = element_blank(),#element_line(colour = "black"),
                       legend.background = element_rect(fill = 'transparent')) 
    }
    
  }
  
  
  return(tmp_plot)
  
  
}

### Plotly function

#### Plots Params ====
# clas.symb <- c("diamond-open" , "cross","square-open", "circle")#,"square" ) # c("square" ,"square-open" , "diamond" , "diamond-open" , "cross" , "x" )
# clas.col <-   c("#83bf93",'#f7b500','#7ab6ff', "#cfa0ff" ) # "#9F3CF4") # "#cfa0ff" "#f7b500")# #"#DFC217" )#, "#7ab6ff")# ,'#efe700',"#28902E") #c('#efe700','#39a079', '#522270') #'Viridis' #
# # purple, yellow, green, blue
# # bg.col <- 'rgb(250, 250, 250)' #'rgb(147, 147, 143)'
# cam.norm <- sqrt( sum( c(0, 0, -10)^2) )/2#1.85
# symb.size <- 6
# ax <- list(
#   zeroline = F,
#   showline = TRUE,
#   showticklabels = FALSE,
#   showgrid = FALSE,
#   mirror = F, #T, # 'all', #T, #"ticks",
#   gridcolor = toRGB("white"),
#   gridwidth = 0,
#   zerolinecolor = toRGB("black"),
#   zerolinewidth = 2,
#   linecolor = toRGB("black"),
#   linewidth = 2,
#   titlefont = list(size=30)
# )

### Plotly Function ====
plotly_plot <- function(y_datta, preds,
                        clas_col=NULL, clas_symb=NULL, symb_size=2,
                        method=NULL, legend=T) {
  
  dattaframe = data.frame(preds ,y_datta);
  colnames(dattaframe) = c('x','y','z','resp');
  
  plot_ly(dattaframe, x=dattaframe$x, y=dattaframe$y , z=dattaframe$z,
          color = as.factor(dattaframe$resp), colors = clas_col,
          symbol = as.factor(dattaframe$resp), symbols = clas_symb,
          marker=list(size=symb_size, sizemode="diameter")
  ) %>%
    add_markers(    ) %>% layout(showlegend = legend)
  
  # %>%
  #   layout(scene = list(xaxis = c(title=paste0(axis_labels_method, as.character(1))),
  #                       yaxis = c(title=paste0(axis_labels_method, as.character(2))),
  #                       zaxis = c(title=paste0(axis_labels_method, as.character(3))) ) )
  # 
}

### Saving plot
save_sdr_plot=function(sdr_plot, filename,width=400,height=400, units="px", 
                       pointsize=12, bg = "white", res = 100) {
  png(filename = filename, width = width, height = height, units = units, 
      pointsize = pointsize, bg = bg,  res = res)
  
  print(sdr_plot)
  dev.off()

}
  

