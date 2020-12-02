# # ###############################################################
# #           Linear SDR Order Determination Methods
# #                        Modified from 
# #           Sufficient Dimension Reduction with R (Li, 2018)
# # ###############################################################
# 
# 
# #################################################################
# # Internal Functions
# #################################################################
# 
# cand_mat=function(x,y,method,nslices,ytype){
#   
#   out_m = switch(method, 
#                  sir=sir(x,y,nslices,d,ytype,lambda=0)$cand_mat,
#                  save=save_sdr(x,y,nslices,r,ytype,lambda=0)$cand_mat,
#                  dr=dr(x,y,nslices,d,ytype,lambda=0)$cand_mat,
#                  opcg=opcg(x,y,d)$cand_mat
#                  ) 
# 
#   return(out_m)}
# 
# # # # x=Reduce('rbind', y.R.0); y=Reduce('c', phi.R.0); r=6; h; method='dr';ytype; k=R.2;
# # mean.candmat <- function(x,y,h,r,ytype,method,k,m=NULL, ensemble=NULL) {
# #   if( is.list(y) ) {
# #     k0 <- length(y)
# #     out <- Reduce('+', lapply(1:k0,
# #                               function(k) candmat(x[[k]],y[[k]],h,r,ytype,method, m, ensemble) ) )/k0
# #   } else {
# #     n=dim(x)[1]/k ;p=dim(x)[2]
# #     # for what we do for now, n/k will always be R.2
# #     # x <- t(Reduce('cbind', y.R.0)); y <- Reduce('c', phi.R.0); method='sir'
# #     x.list <- lapply(0:(k-1), function(k) x[(1+n*k):(n*(k+1)),])
# #     y.list <- lapply(0:(k-1), function(k) y[(1+n*k):(n*(k+1))])
# #     out <- Reduce('+', lapply(1:k, #k<-1
# #                               function(k) candmat(x.list[[k]],y.list[[k]],h,r,ytype,method, m, ensemble) ) )/k
# #   }
# #   return(out)
# # }
# 
# # Non-linear Candidate Matrices
# # nl.candmat=function(x,y,r,ytype,atype,ex,ey,complex.x,complex.y, method){
# #   if(method=="gsir") mat=gsir(x,y,r,ytype,atype,ex,ey,complex.x,complex.y,
# #                               can.mat=T,x.new=NULL)$gsirmat
# #   return(mat)}
# 
# # mean.nl.candmat <- function(x,y,r,ytype,atype,ex,ey,complex.x,complex.y, method, k) {
# #   if( is.list(y) ) {
# #     k0 <- length(y)
# #     out <- Reduce('+', lapply(1:k0,
# #                               function(k) nl.candmat(x.list[[k]],y.list[[k]],
# #                                                      r,ytype,atype,ex,ey,complex.x,complex.y, method) ) )/k0
# #   } else {
# #     n=dim(x)[1]/k ;p=dim(x)[2]
# #     # for what we do for now, n/k will always be R.2
# #     # x <- t(Reduce('cbind', y.R.0)); y <- Reduce('c', phi.R.0); method='sir'
# #     x.list <- lapply(0:(k-1), function(k) x[(1+n*k):(n*(k+1)),])
# #     y.list <- lapply(0:(k-1), function(k) y[(1+n*k):(n*(k+1))])
# #     out <- Reduce('+', lapply(1:k,
# #                               function(k) nl.candmat(x.list[[k]],y.list[[k]],
# #                                                      r,ytype,atype,ex,ey,complex.x,complex.y, method) ) )/k
# #   }
# #   return(out)
# # }
# 
# 
# # Variation in Eigenvalues component
# phi_n=function(kmax,eval){
#   den=1+sum(eval[1:(kmax+1)]);
#   return(eval[1:(kmax+1)]/den)
# }
# 
# # Variation in Eigenvectors sub-component 
# # Measures variation/distance between two Eigen subspaces
# # corresponind the the Eigen matrices given
# prefn0=function(kmax,evec1,evec2){
#   out=numeric();
#   for(k in 0:kmax){
#     if(k==0) out=c(out,0)
#     if(k>0) out=c(out,1-abs(t(evec1[,1:k])%*%evec2[,1:k]))
#     # if(k!=0&k!=1) out=c(out,1-abs(det(t(evec1[,1:k])%*%evec2[,1:k])))}
#   return(out)
# }
# 
# 
# # minimizer=function(a,b) return(a[order(b)][1])
# 
# # x=t(y.R);y=phi.R;h=4;r=2;nboot=1000;method='dr';ytype='continuous'
# 
# ################################################################
# # Ladle estimator
# ################################################################
# 
# # # The number of bootstrap samples
# # set.seed(1991)
# # n_boot=500
# # start_time3=Sys.time();
# # opcg_boot_list=lapply(1:n_boot, function(rep) {
# #   boot_set=sample(1:n,n, replace=T);
# #   eigen_opcg=eigen_cpp(opcg_DD(X[,boot_set], matrix(Y[,boot_set], nrow=1, ncol=n), 
# #                                h=1, ytype="multinomial", method="newton",
# #                                parallelize=T)$opcg_mat); 
# #   return(eigen_opcg);
# # }  ) 
# # end_time3=Sys.time();  end_time3 - start_time3 
# # # about 30 mins on 15 cores for 200
# # # 3 mins with cpp code for 200; # 1,3 mins with package for 200,500
# # 
# # 
# # ### The whole OPG estimate
# # eigen_opcg=eigen_cpp(opcg_DD(X, Y, h=1, ytype='multinomial', parallel = T)$opcg_mat);
# # eigenval_opcg=eigen_opcg$values;
# # eigenvec_opcg=eigen_opcg$vectors;
# # 
# # if(p<=10) kmax=p-2 else if(p>10) kmax=round(p/log(p))
# 
# 
# 
# 
# ladle=function(x,y,method,nslices=NULL,bw=NULL,ytype,nboot,opcg_args=list()){
#   
#   #* for opcg, the extra arguments are
#   #* method="newton", lambda2a=0, lambda2b=0, 
#   #* parallelize=F, control_list=list()
#   #*
#   #*
#   #*
#   #*
#   #*
#   #*
#   #* 
#   #*     
#   #*   
#   
#   # the dimension arg for SDR methods for ladle doesn't matter
#   # just need some specified d to run the command
#   d=2; n=dim(x)[2];p=dim(x)[1]
#   
#   if(p<=10) kmax=p-2 else if(p>10) kmax=round(p/log(p))
#   
#   # Computing the Candidate Matrix for the SDR method
#   out=cand_mat(x,y,nslices,d,ytype,method)
#   
#   eval_full=eigen(out)$values; # calculate full eigenvectors
#   evec_full=eigen(out)$vectors; # calculate full eigenvalues
#   
#   pn=phi_n(kmax,eval_full); #variation of eigenvalues component
# 
#   fn0=0;
#   for(iboot in 1:nboot){
#     # Draws Bootstrap Sample 
#     bootindex=round(runif(n,min=-0.5,max=n+0.5));
#     xs=x[bootindex,]; ys=y[bootindex]
#     
#     #Computes the Candidate Matrix and Eigenvectors
#     mat=candmat(xs,ys,h,r,ytype,method)
#     evec=eigen(mat)$vectors
#     
#     # Adds to the variation in Eigenvectors component
#     fn0=fn0+prefn0(kmax,evec_full,evec)/nboot};
#   # Final Variation in Eigenvectors component
#   fn=fn0/(1+sum(fn0));
#   
#   # Ladle plot values and minimizer
#   gn=pn+fn;
#   rhat=minimizer(0:kmax,gn)
#   return(list(kset=(0:kmax),gn=gn,rhat=rhat))
# }
# 
# 
# 
# ################################################################
# # K-fold Ladle estimator
# ################################################################
# 
# # What was the point of this again? Forgot why I did this.
# 
# # reladle=function(x,y,h,nboot,method,ytype, k,m=NULL,ensemble=NULL){ 
# #   # x=t(y.R.0[[1]]); y=phi.R.0[[1]]; h=10;nboot=200;ytype='continuous'; k=1;method='iht'
# #   # x=t(Reduce('cbind', y.R.0[1:2])); y=Reduce('c', phi.R.0[1:2]); h=10;nboot=200;ytype='continuous'; k=2;method='iht'
# #   # x=t(Reduce('cbind', y.R.0)); y=Reduce('c', phi.R.0); h=10;nboot=200;ytype='continuous'; k=R.2;method='iht'
# #   # x=t(Reduce('cbind', y.R.0)); y=Reduce('c', phi.R.0); h=10;nboot=200;ytype='continuous'; k=R.2;method='sir'
# #   r=2;
# #   n=dim(x)[1];p=dim(x)[2] 
# #   if(p<=10) kmax=p-2 else if(p>10) kmax=round(p/log(p))
# #   # The initial estimate
# #   out=mean.candmat(x,y,h,r,ytype,method,k,m,ensemble)
# #   eval.full=eigen(out)$values;evec.full=eigen(out)$vectors
# #   pn=phin(kmax,eval.full)
# #   
# #   fn0=0;
# #   for(iboot in 1:nboot){
# #     print(iboot)
# #     bootindex = sample(1:n, n, replace=T) # round(runif(n,min=-0.5,max=n+0.5)); # sample runs faster
# #     xs = x[bootindex,]; ys=y[bootindex] 
# #     mat = mean.candmat(xs,ys,h,r,ytype,method,k,m,ensemble) 
# #     eval = eigen(mat)$values;evec=eigen(mat)$vectors 
# #     fn0 = fn0 + prefn0(kmax,evec.full,evec)/nboot} 
# #   fn = fn0/(1+sum(fn0)); 
# #   gn = pn+fn;
# #   rhat = minimizer(0:kmax,gn) 
# #   return(list(kset=(0:kmax),gn=gn,rhat=rhat))
# # }
# 
# 
# 
# 
# 
# # ################################################################
# # # Sequential Test for SIR (Li, 2018)
# # ################################################################
# # #x=t(y.R);y=phi.R;h=4;r=2;nboot=1000;method='dr';ytype='continuous'
# # seqtestsir=function(x,y,h,r,ytype){ 
# #   p=ncol(x);n=nrow(x)
# #   signrt=matpower(var(x),-1/2) 
# #   xc=t(t(x)-apply(x,2,mean))
# #   xst=xc%*%signrt
# #   if(ytype=="continuous") ydis=discretize(y,h) 
# #   if(ytype=="categorical") ydis=y 
# #   yless=ydis;ylabel=numeric() 
# #   for(i in 1:n) {
# #     if(var(yless)!=0) {
# #       ylabel=c(ylabel,yless[1]); 
# #       yless=yless[yless!=yless[1]]
# #     }
# #   }
# #   ylabel=c(ylabel,yless[1])
# #   prob=numeric();exy=numeric() 
# #   for(i in 1:h) prob=c(prob,length(ydis[ydis==ylabel[i]])/n) 
# #   for(i in 1:h) exy=rbind(exy,apply(xst[ydis==ylabel[i],],2,mean)) 
# #   sirmat=t(exy)%*%diag(prob)%*%exy 
# #   test=n*sum(eigen(sirmat)$values[(r+1):p]) 
# #   pval=1-pchisq(test,(p-r)*(h-1-r))
# # return(pval)}
# 
# 
# # ################################################################
# # # Sequential Test for DR (Li, 2018)
# # ################################################################
# # #x=t(y.R);y=phi.R;h=4;r=2;nboot=1000;method='dr';ytype='continuous'
# # seqtestdr=function(x,y,h,r,ytype){ 
# #   p=ncol(x);n=nrow(x) 
# #   xc=t(t(x)-apply(x,2,mean)) 
# #   if(ytype=="continuous") ydis=discretize(y,h) 
# #   if(ytype=="categorical") ydis=y 
# #   yless=ydis;ylabel=numeric() 
# #   for(i in 1:n) { 
# #     if(var(yless)!=0){
# #       ylabel=c(ylabel,yless[1]);
# #       yless=yless[yless!=yless[1]]}} 
# #   ylabel=c(ylabel,yless[1]) 
# #   sk=matrix(0,n,h); 
# #   for(i in 1:n) 
# #     for(k in 1:h) 
# #       if(ydis[i]==ylabel[k]){ sk[i,k]=1
# #       }
# #   pk=apply(sk,2,mean) 
# #   sig=var(xc) 
# #   uk=numeric() 
# #   for(i in 1:h) uk=rbind(uk,apply(xc[ydis==ylabel[i],],2,mean)) 
# #   vk=array(0,c(p,p,h)); 
# #   for(i in 1:h) vk[,,i]=var(xc[ydis==ylabel[i],])-sig
# #   m1=numeric(); 
# #   for(k in 1:h) m1=cbind(m1,pk[k]^(1/2)*vk[,,k]) 
# #   m2=t(uk)%*%diag(pk)%*%uk 
# #   tr=function(a) return(sum(diag(a))) 
# #   m3=numeric(); 
# #   for(k in 1:h) m3=cbind(m3,pk[k]^(1/2)*(tr(m2))^(1/2)*uk[k,])
# #   m=cbind(m1,m2,m3) 
# #   pkst=numeric(); 
# #   for(i in 1:n) pkst=rbind(pkst,sk[i,]-pk) 
# #   ukst=array(0,c(p,h,n)); 
# #   for(i in 1:n) for(k in 1:h) { ukst[,k,i]=(xc[i,]-uk[k,])*sk[i,k]/pk[k]-xc[i,]} 
# #   vkst=array(0,c(p,p,h,n)); 
# #   for(k in 1:h) for(i in 1:n) { 
# #     vkst[,,k,i]=(xc[i,]%*%t(xc[i,])-sig-vk[,,k])*sk[i,k]/pk[k]-uk[k,]%*%t(xc[i,])-xc[i,]%*%t(uk[k,])-xc[i,]%*%t(xc[i,])+sig
# #     } 
# #   m1kst=array(0,c(p,p,h,n)); 
# #   for(k in 1:h) for(i in 1:n) { m1kst[,,k,i]=
# #     (1/2)*pk[k]^(-1/2)*pkst[i,k]*vk[,,k]+pk[k]^(1/2)*vkst[,,k,i]} 
# #   m2st=array(0,c(p,p,n)); 
# #   for(i in 1:n){ for(k in 1:h) {
# #       m2st[,,i]=m2st[,,i]+pkst[i,k]*uk[k,]%*%t(uk[k,])+ pk[k]*ukst[,k,i]%*%t(uk[k,])+pk[k]*uk[k,]%*%t(ukst[,k,i])}
# #   } 
# #   m3kst=array(0,c(p,h,n)); 
# #   for(k in 1:h) for(i in 1:n){ 
# #     m3kst[,k,i]=(1/2)*pk[k]^(-1/2)*pkst[i,k]*(tr(m2))^(1/2)*uk[k,]+ (1/2)*pk[k]^(1/2)*(tr(m2))^(-1/2)*tr(m2st[,,i])*uk[k,]+ pk[k]^(1/2)*(tr(m2))^(1/2)*ukst[,k,i]
# #     }
# #   mst=array(0,c(p,p*h+p+h,n)); 
# #   for(i in 1:n){ 
# #     mst1=numeric(); 
# #     for(k in 1:h) mst1=cbind(mst1,m1kst[,,k,i]) 
# #     mst3=numeric(); 
# #     for(k in 1:h) mst3=cbind(mst3,m3kst[,k,i]) 
# #     mst[,,i]=cbind(mst1,m2st[,,i],mst3)} 
# #   if(r==0) qlef=diag(p); qrig=diag(p*h+p+h) 
# #   svdm=svd(m);u1=svdm$u[,1:r]; v1=svdm$v[,1:r]
# #   qlef=diag(p)-u1%*%t(u1); qrig=diag(p*h+p+h)-v1%*%t(v1)
# #   qmstq=array(0,c(p,p*h+p+h,n)); 
# #   for(i in 1:n) qmstq[,,i]= qlef%*%mst[,,i]%*%qrig
# #   vqmstq=matrix(0,n,p*(p*h+p+h)); 
# #   for(i in 1:n) vqmstq[i,]=c(qmstq[,,i]) 
# #   s=var(vqmstq);omega=eigen(s)$values[1:((p-r)*(p*h+p+h-r))] 
# #   teststat=n*sum(eigen(m%*%t(m))$values[(r+1):p]) 
# #   obs=numeric();nsim=20000 
# #   for(isim in 1:nsim) {
# #     ki=rnorm(length(omega))^2;obs=c(obs,sum(omega*ki))} 
# #   pval=length((1:nsim)[obs>teststat])/nsim;
# #   return(pval)
# # }
