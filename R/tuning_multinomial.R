#######################################################################
#           Tuning functions for classification
#######################################################################

### Unsupervised km tuning

tuning_km=function(x_datta, n_hlist, edr_list, centers, iter.max=100, nstart = 5){
  
  tune_km = lapply(1:n_hlist, function(k) #k=5; x_datta=x; centers=n_centers
    (kmeans(t( t(edr_list[[k]])%*%x_datta ), centers, iter.max, nstart)$tot.withinss)/
      (kmeans(t( t(edr_list[[k]])%*%x_datta ), centers, iter.max, nstart)$betweenss)
  )
  return(Reduce('c', tune_km))
}

#' This is an internal function called by opcg 
#' 
#'
#'
#' @keywords internal
#' @noRd
#' 
tuning_skm <- function(x, y, d, class_labels, n_cpc, n_hlist, edr_list,
                       iter.max = 100, nstart = 100) {
  # x=X_tune; y=Y_tune; d=km_d; n_classes=m; n_cpc=1; n_hlist=length(h_list);
  # edr_list; nstart = 200; iter.max = 100
  n_classes=length(class_labels);
  if (length(n_cpc)==1) n_cpc_vec=rep(n_cpc, n_classes) else n_cpc_vec=n_cpc;
  
  sapply(1:n_hlist, function(k) {
    # k=2
    km_ss=lapply(1:n_classes, function(l) {
      # l=3
      # SDR Predictors, "n x p" format
      sdr_pred=t( t(edr_list[[k]][,1:d])%*%x[,which(y==class_labels[l])] );
      class_grand_mean=colMeans(sdr_pred);
      
      km=kmeans(sdr_pred, centers=n_cpc_vec[l], iter.max = iter.max, nstart = nstart )
      
      # Sum of WSS for the clusters of a given class l;
      per_class_wss=km$tot.withinss
      # the estimated centers of the given class l;
      per_class_cen=km$centers
      # The BSS of clusters (i.e. centers) by class
      if (n_cpc_vec[l]>1) {
        per_class_cen_wss = sum(diag( var(per_class_cen)*(n_cpc_vec[l]-1) ))
      } else {
        per_class_cen_wss = 0
      }
      
      return( list(per_class_wss, per_class_cen,
                   per_class_cen_wss,
                   class_grand_mean) )
    })
    
    # Want a measure of distance from centers of other classes
    # want a measure of distance from within class centers
    
    # list of all the estimated centers
    class_cen <- lapply(1:n_classes, FUN=function(j) km_ss[[j]][[2]]  )
    # Adding up the WSS by classes
    class_wss <- sum( sapply(1:n_classes, FUN=function(j) km_ss[[j]][[1]]  ) )
    # Adding up the Tot.SS of the centers by class; effectively the supervised BSS
    class_grand_means<- lapply(1:n_classes, FUN=function(j) km_ss[[j]][[4]]  )
    
    # The idea is that if the class has just one cluster, then the variance between
    # intra class centers is 0, because it is just one center/point.
    # The idea behind the class_bss is to add up all the sums of squares for
    # intra-class centers and this is a sort of measure of between center sums
    # of squares.
    
    # We can use class_bss1 for the sum of intra-class between SS
    # We can use class_bss2 for the total variance of class centers
    # We can use class_bss3 for the between class SS (calculate the mean of each class)
    
    class_bss1= sum( sapply(1:n_classes, FUN=function(j) km_ss[[j]][[3]]));
    class_bss2= sum(diag(var( Reduce('rbind', class_cen) )))*(n_classes-1);
    # class_bss3= sum(diag(var( Reduce('rbind', class_grand_means) )))*(n_classes-1);
    
    if ( sum(n_cpc_vec)>n_classes) {
      return( class_wss/ #max(class_bss1, class_bss2) ) #
    } else if(sum(n_cpc_vec)==n_classes) {
      return( class_wss/class_bss2 )#max(class_bss1, class_bss2) )#
    }
    
    
  } )
}

# ### Weighted Supervised km tuning
#' This is an internal function called by opcg 
#' 
#'
#'
#' @keywords internal
#' @noRd
#' 
tuning_skmkm <- function(x, y, d, class_labels, n_cpc, n_hlist, edr_list,
                         iter.max = 100, nstart = 100) {
  # x=X_tune; y=Y_tune; d=km_d; n_classes=m; n_cpc=1; n_hlist=length(h_list);
  # edr_list; nstart = 200; iter.max = 100
  
  # x=X_test; y=Y_test; d=10; class_labels=m_classes; n_cpc=1; n_hlist=length(bw_list);
  # edr_list; nstart = 100; iter.max = 100
  
  n_classes=length(class_labels);

  #Supervised Tuning
  skm_tune=tuning_skm(x, y, d, class_labels, n_cpc, n_hlist, edr_list,
                 iter.max, nstart)

  #Unsupervised Tuning
  # number of centers for unsupervised is n_classes*n_cpc if n_cpc is a
  # number, i.e. each class has same number, otherwise, it should be sum(n_cpc)
  n_centers = if (length(n_cpc)==1) n_classes*n_cpc else sum(n_cpc);
  
  km_tune=tuning_km(x, n_hlist, edr_list,
               centers=n_centers, iter.max=iter.max, nstart = nstart)

  # Standardizing the ratios
  skm_tune_std=(skm_tune-mean(skm_tune))/sd(skm_tune)
  km_tune_std=(km_tune-mean(km_tune))/sd(km_tune)

  return(list(wskm= (.5*skm_tune_std + .5*km_tune_std),
              km_std=km_tune_std,
              skm_std=skm_tune_std ))

}

# 
# tuning_skmkm <- function(x, y, d, n_classes, n_cpc, n_hlist, edr_list,
#                          iter.max = 100, nstart = 100) {
#   # x=X_tune; y=Y_tune; d=km_d; n_classes=m; n_cpc=1; n_hlist=length(h_list); 
#   # edr_list; nstart = 200; iter.max = 100
#   
#   skm_tune=tuning_skm(x, y, d, n_classes, n_cpc, n_hlist, edr_list,
#                       iter.max, nstart)
#   skm_tune_std=(skm_tune-mean(skm_tune))/sd(skm_tune)
#   
#   km_tune=tuning_km(x, n_hlist, edr_list, 
#                     centers=n_classes*n_cpc, iter.max=iter.max, nstart = nstart)
#   km_tune_std=(km_tune-mean(km_tune))/sd(km_tune)
#   
#   return(.5*skm_tune_std + .5*km_tune_std)
#   
# }
############################################################################
###                       K-fold K-Means Tuning
############################################################################
#' K-fold Tuning with K-means
#'
#' This implements the tuning procedure for SDR and classification problems 
#' in the forth coming paper Quach and Li (2021).
#' 
#' The kernel for the local linear regression is fixed at a gaussian kernel.
#' 
#' For large 'p', we strongly recommend using the Conjugate Gradients implement, 
#' by setting method="cg".
#' For method="cg", the hybrid conjugate gradient of Dai and Yuan is implemented,
#' but only the armijo rule is implemented through backtracking, like in Bertsekas'
#' "Convex Optimization Algorithms".
#' A weak Wolfe condition can also be enforced by adding setting c_wolfe > 0 
#' in the control_list, but since c_wolfe is usually set to 0.1 (Wikipedia)
#' and this drastically slows down the algorithm relative to newton for small to 
#' moderate p, we leave the default as not enforcing a Wolfe condition, since we 
#' assume that our link function gives us a close enough initial point that local
#' convergence is satisfactory. Should the initial values be suspect, then maybe
#' enforcing the Wolfe condition is a reasonable trade-off.  
#' 
#' @param d specified the reduced dimension 
#' @param method "newton" or "cg" methods; for carrying out the optimization using
#' the standard newton-raphson (i.e. Fisher Scoring) or using Congugate Gradients 
#' @param parallelize Default is False; to run in parallel, you will need to have
#' foreach and some parallel backend loaded; parallelization is strongly recommended
#' and encouraged.
#' @param control_list a list of control parameters for the Newton-Raphson 
#' or Conjugate Gradient methods
#' @param ytype specify the response as 'continuous', 'multinomial', or 'ordinal' 
#' @param h_list 
#' @param k 
#' @param x_datta 
#' @param y_datta 
#' @param class_labels 
#' @param n_cpc 
#' @param iter.max 
#' @param nstart 
#'
#' @return A list containing both the estimate and candidate matrix.
#' \itemize{
#'   \item opcg - A 'pxd' matrix that estimates a basis for the central subspace.
#'   \item opcg_wls - A 'pxd' matrix that estimates a basis for the central subspace based 
#'   on the initial value of the optimization problem; useful for examining bad starting 
#'   values.
#'   \item cand_mat - A list that contains both the candidate matrix for OPCG and for
#'   the initial value; this is used in other functions for order determination
#'   \item gradients - The estimated local gradients; used in regularization of OPCG
#'   \item weights - The kernel weights in the local-linear GLM. 
#' }
#'  
#' @export
#' 
kfold_km_tuning=function(h_list, k, x_datta, y_datta, d, ytype,
                         class_labels, n_cpc, method="newton", parallelize = F,
                         control_list=list(),
                         iter.max = 100, nstart = 100) {

  n=dim(x_datta)[2]; n_hlist = length(h_list);
  # k = 3;
  # x_datta=X; y_datta=Y; ytype="multinomial"; parallelize=T
  # class_labels=m_classes; iter.max = 100; nstart = 100
  # method="cg"; control_list=list();

  # partition1=1:250#sample(1:n, round(n/2,0), replace = F)
  # partitions=list(sort(partition1), (1:n)[-partition1])
  partitions=sample(1:(k),size=n,replace=T,prob=rep(1/(k),(k) ) );
  # n==length(partitions); table(partitions)

  X_folds=lapply(1:k, function(fold) x_datta[,which(partitions==fold)])
  Y_folds=lapply(1:k, function(fold) matrix(y_datta[,which(partitions==fold)], 1,
                 table(partitions)[fold]) )

  # X_folds=lapply(1:k, function(fold) x_datta[,partitions[[fold]] ])
  # Y_folds=lapply(1:k, function(fold) matrix(y_datta[, partitions[[fold]] ], 1,
  #                length(partitions[[fold]]) ) )

  # edr_whole=list();
  km_whole=list();
  km_std=list();
  skm_std=list();
  for(fold in 1:k) {
    # fold=2;
    X_tune=X_folds[[fold]]; Y_tune=Y_folds[[fold]];

    X_train=x_datta[,-which(partitions==fold) ];
    Y_train=matrix(y_datta[,-which(partitions==fold) ] , 1, ncol=(n-dim(X_tune)[2]) ) ;

    # dim(X_tune)[2] + dim(X_train)[2]
    # table(partitions)

    edr_list=list()
    for (iter in 1:n_hlist) {
      edr_iter=opcg(x_matrix=X_train, y_matrix=Y_train, d=d, bw=h_list[iter],
                    ytype=ytype, method=method, parallelize=parallelize,
                    control_list=control_list)
      edr_list[[iter]]=edr_iter;
    }

    km_tune=tuning_skmkm(X_tune, Y_tune, d, class_labels, n_cpc, n_hlist, edr_list,
                 iter.max, nstart)

    # edr_whole[[fold]]=edr_list;
    km_whole[[fold]]=km_tune$wskm;
    km_std[[fold]]=km_tune$km_std;
    skm_std[[fold]]=km_tune$skm_std;

  }

  # return(list(edr=edr_whole, km=km_whole))
  return(list(wskmf_ratio=Reduce('+', km_whole)/k,
              skm_f_ratio=Reduce('+', skm_std)/k,
              km_f_ratio=Reduce('+', km_std)/k));

}

