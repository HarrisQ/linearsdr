################################################################
#                   Inverse Linear SDR methods
#                        Modified from 
#           Sufficient Dimension Reduction with R (Li, 2018)
################################################################

################################################################
#                       discretize
################################################################
discretize=function(y,h){
  n=length(y);m=round(n/h)
  y=y+.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt=numeric();for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1])
  y1=rep(0,n);y1[y<divpt[1]]=1;y1[y>=divpt[h-1]]=h
  for(i in 2:(h-1)) y1[(y>=divpt[i-1])&(y<divpt[i])]=i
  return(y1)
}


################################################################
#                           SIR
################################################################

sir=function(x,y,nslices,d,ytype, std = T, lambda=0){
  p=ncol(x);n=nrow(x)
  if (std) {
    signrt=matpower_cpp(var(x)+diag(lambda,p,p),-1/2)
    xc=t(t(x)-colMeans(x))
    xst=xc%*%signrt
  } else {
    xst=x
  }
  if(ytype=="continuous") ydis=discretize(y,nslices)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric();exy=numeric()
  for(i in 1:nslices) prob=c(prob,length(ydis[ydis==ylabel[i]])/n) 
  for(i in 1:nslices) exy=rbind(exy,colMeans(xst[ydis==ylabel[i],]))
  sirmat=t(exy)%*%diag(prob)%*%exy;
  
  if (std) {
    beta=signrt%*%eigen_cpp(sirmat)$vectors[,1:d];
  } else {
    beta=eigen_cpp(sirmat)$vectors[,1:d];
  }
  
  return( list(beta=beta, cand_mat=sirmat) )
}

################################################################
#                          save
################################################################
save_sdr=function(x,y,nslices,d,ytype,std=T,lambda=0){
  p=ncol(x);n=nrow(x)
  if (std) {
    signrt=matpower_cpp(var(x)+diag(lambda,p,p),-1/2)
    xc=t(t(x)-colMeans(x))
    xst=xc%*%signrt
  } else {
    xst=x
  }
  if(ytype=="continuous") ydis=discretize(y,nslices)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  prob=numeric() 
  for(i in 1:nslices) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  vxy = array(0,c(p,p,h))
  for(i in 1:nslices) vxy[,,i] = var(xst[ydis==ylabel[i],]) 
  savemat=0
  for(i in 1:nslices){
    savemat=savemat+prob[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p))};
  
  
  if (std) {
    beta=signrt%*%eigen_cpp(savemat)$vectors[,1:d];
  } else {
    beta=eigen_cpp(savemat)$vectors[,1:d];
  }
  return( list(beta=beta, cand_mat=savemat) )
}

################################################################
#                                 dr 
################################################################
# x=Reduce('rbind', y.R.0); y=Reduce('c', phi.R.0)
# x=t(X_std);y=Y;nslices=m;d=d ;ytype="categorical"; std=F; lambda = 7
dr=function(x,y,nslices,d,ytype,std=T,lambda=0){
  p=ncol(x);n=nrow(x)
  if (std) {
    signrt=matpower_cpp(var(x)+diag(lambda,p,p),-1/2)
    xc=t(t(x)-colMeans(x))
    xst=xc%*%signrt
  } else {
    xst=x
  }
  if(ytype=="continuous") ydis=discretize(y,nslices)
  if(ytype=="categorical") ydis=y
  ylabel=unique(ydis)
  # prob=numeric() 
  # for(i in 1:nslices) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)
  # vxy = array(0,c(p,p,nslices)); 
  #* vxy Cannot be done for very large data
  #* eg. p=10000, n=1350
  #* rather than making a cube, maybe 
  #* make a list instead (still gets too large)
  #* or create and overwrite vxy as needed
  # vxy=list()
  # exy=numeric()
  # for(i in 1:nslices) { # at i=19, vxy is already too large.
  #   # vxy[[i]] = var(xst[ydis==ylabel[i],])
  #   # vxy[,,i]=var(xst[ydis==ylabel[i],])
  #   exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))}  # This is fine; just a matrix
  mat1 = matrix(0,p,p);mat2 = matrix(0,p,p); 
  for(i in 1:nslices){ #var(xst[ydis==ylabel[i],])
    prb_i = length(ydis[ydis==ylabel[i]])/n;
    vxy_i = var(xst[ydis==ylabel[i],]); 
    exy_i = apply(xst[ydis==ylabel[i],],2,mean) ;
    
    mat1 = mat1+prb_i*( vxy_i+exy_i%*%t( exy_i ))%*% ( vxy_i+ exy_i %*%t( exy_i ))
    mat2 = mat2+prb_i*exy_i%*%t( exy_i)
    
    # mat1 = mat1+prob[i]*(vxy[,,i]+exy[i,]%*%t(exy[i,]))%*%
    #   (vxy[,,i]+exy[i,]%*%t(exy[i,]))
    # mat2 = mat2+prob[i]*exy[i,]%*%t(exy[i,])
  }
  drmat = 2*mat1+2*mat2%*%mat2+2*sum(diag(mat2))*mat2-2*diag(p)
  
  if (std) {
    beta=signrt%*%eigen_cpp(drmat)$vectors[,1:d];
  } else {
    beta=eigen_cpp(drmat)$vectors[,1:d];
  }
  
  return( list(beta=beta, cand_mat=drmat) )
}



# ################################################################
# #                          pir - monomial basis
# ################################################################
# pir=function(x,y,m,r){
#   xc=t(t(x)-apply(x,2,mean))
#   signrt=matpower_cpp(var(x),-1/2)
#   xstand = xc%*%signrt
#   f=numeric();ystand=(y-mean(y))/sd(y)
#   for(i in 1:m) f = cbind(f, ystand^i) 
#   sigxf=cov(xstand,f);sigff=var(f) 
#   cand=sigxf%*%solve(sigff)%*%t(sigxf)
#   return(signrt%*%eigen_cpp(symmetry(cand))$vectors[,1:r])
# }
# ################################################################
# #                          pir - sin basis
# ################################################################
# pir_sin=function(x,y,m,r){
#   xc=t(t(x)-apply(x,2,mean))
#   signrt=matpower(var(x),-1/2)
#   xstand = xc%*%signrt
#   f=numeric();ystand=(y-mean(y))/sd(y)
#   for(i in 1:m) f = cbind(f, sin(2*pi*ystand^i)) 
#   sigxf=cov(xstand,f);sigff=var(f) 
#   cand=sigxf%*%solve(sigff)%*%t(sigxf)
#   return(signrt%*%eigen(symmetry(cand))$vectors[,1:r])
# }
################################################################
#                        sir ii 
################################################################

# sirii=function(x,y,h,r,ytype="continuous"){                                         
#   p=ncol(x);n=nrow(x)                                              
#   signrt=matpower(var(x),-1/2)                                     
#   xst=t(t(x)-apply(x,2,mean))%*%signrt
#   if(ytype=="continuous") ydis=discretize(y,h)  
#   if(ytype=="categorical") ydis=y 
#   ylabel=unique(ydis)
#   prob=numeric();exy=numeric()                                     
#   for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n)               
#   for(i in 1:h) exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean))         
#   sirmat=t(exy)%*%diag(prob)%*%exy                                 
#   vxy = array(0,c(p,p,h))                                          
#   for(i in 1:h) vxy[,,i] = var(xst[ydis==ylabel[i],])                      
#   savemat=0                                                        
#   for(i in 1:h){                                                   
#     savemat=savemat+prob[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p))} 
#   siriimat=savemat-sirmat%*%t(sirmat)                              
#   return(signrt%*%eigen(siriimat)$vectors[,1:r])}


################################################################
#                               cr 
################################################################
# cr=function(x,y,percent,r){                                  
#   tradeindex12 = function(k,n){                                
#     j = ceiling(k/n)                                             
#     i = k - (j-1)*n                                              
#     return(c(i,j))}                                                                                    
#   mu=apply(x,2,mean);signrt=matpower(var(x),-1/2)              
#   z=t(t(x)-mu)%*%signrt                                        
#   n=dim(x)[1];p = dim(x)[2]                                    
#   ymat=matrix(y,n,n)                                           
#   deltay=c(abs(ymat - t(ymat)))                                
#   singleindex=(1:n^2)[deltay < percent*mean(deltay)]           
#   contourmat=matrix(0,p,p)                                     
#   for(k in singleindex){                                       
#     doubleindex=tradeindex12(k,n)                                
#     deltaz=z[doubleindex[1],]-z[doubleindex[2],]                 
#     contourmat=contourmat+deltaz %*% t(deltaz)}                  
#   signrt=matpower(var(x),-1/2)                                 
#   return(signrt%*%eigen(contourmat)$vectors[,p:(p-r+1)])                   
# }

################################################################
# IHT & IST
################################################################
#x=t(y.R);y=phi.R;h=4;r=2;nboot=1000;method='dr';ytype='continuous'

# 
# standmat=function(x){
#   mu=apply(x,2,mean);sig=var(x);signrt=matpower(sig,-1/2) 
#   return(t(t(x)-mu)%*%signrt)}
# 
# iht=function(x,y,r, can.mat=F){ 
#   signrt=matpower(var(x),-1/2); z=standmat(x);szy=cov(z,y);szz=var(z);
#   szzy= t(z*(y-mean(y)))%*%z/n
#   p=dim(z)[2];
#   imat=szy 
#   for(i in 1:(p-1)) imat=cbind(imat,matpower(szzy,i)%*%szy)
#   ihtmat = imat%*%t(imat)
#   if(can.mat) {
#     return( list(beta=signrt%*%eigen(ihtmat)$vectors[,1:r], 
#                  ihtmat=ihtmat) )
#   } else {
#     return(signrt%*%eigen(ihtmat)$vectors[,1:r])
#   }
# }
# 
# ist=function(x,y,r,h,g,ytype){
#   
#   z=standmat(x); p=dim(z)[2];
#   if(ytype=="continuous") { ydis.g=discretize(y,g); ydis.h=discretize(y,h) } 
#   if(ytype=="categorical") { ydis.g=y; ydis.h=y }
#   ylabel.g=unique(ydis.g); ylabel.h=unique(ydis.h)
#   
#   prob.g=numeric();prob.h=numeric()
#   for(i in 1:g) prob.g=c(prob.g,length(ydis.g[ydis.g==ylabel.g[i]])/n)
#   for(i in 1:h) prob.h=c(prob.h,length(ydis.h[ydis.h==ylabel.h[i]])/n)
#   
#   ezy=numeric()
#   for(i in 1:g) ezy=rbind(ezy,apply(z[ydis.g==ylabel.g[i],],2,mean))
#   sirmat=t(ezy)%*%diag(prob.g)%*%ezy
#   
#   vxy = array(0,c(p,p,h))
#   for(i in 1:h) vxy[,,i] = var(z[ydis.h==ylabel.h[i],]) 
#   savemat=0
#   for(i in 1:h){ savemat=savemat+prob.h[i]*(vxy[,,i]-diag(p))%*%(vxy[,,i]-diag(p)) }
#   
#   imat=sirmat
#   for(i in 1:(p-1)) imat=cbind(imat, matpower(savemat,i)%*%sirmat )
#   
#   return(eigen(imat%*%t(imat))$vectors[,1:r])}


