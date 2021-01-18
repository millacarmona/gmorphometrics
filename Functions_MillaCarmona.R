
# 1. BASICS AND MISC ############################################################################

## consensus ------------------------------------------------------------------------

consensus<-function(gpa, index=1:dim(gpa)[3]){
  
  p<-nrow(gpa)
  data2d<-two.d.array(gpa)
  
  if(is.numeric(index)) {
    colmeans<-colMeans(data2d[index,])
    consensus<-matrix(colmeans, nrow=p, byrow=TRUE)
  }
  
  if(is.factor(index)){
    colmeans<-apply(X=data2d, MARGIN=2, FUN=tapply, index, mean)
    consensus<-arrayspecs(A=colmeans, p = p, k = ncol(gpa))
    dimnames(consensus)[[3]]<-levels(index)
  }
  
  return(consensus)
  
  
  #Description
  # This function computes the mean shape of the entire sample, the mean shape of a subset 
  # of the sample, or the mean shape of the levels of a factor (for landmark configurations).
  
  #Arguments:
  # gpa= a k x p x n array of superimposed landmarks
  # index= either a numeric vector indicating the configurations to be averaged or
  #        a factor whose levels are used to average groups of configurations
}


## stack -------------------------------------------------------------------------------

stack<-function(gpa, links=NULL, mshape=TRUE, ...){
  
  longlist<-c() ; for(i in 1:dim(gpa)[3]) longlist<-rbind(longlist, gpa[,,i])
  plot(longlist, col="#708095", axes=FALSE, xlab="", ylab="", ...)
  if(mshape==TRUE) {
    points(consensus(gpa), pch=16, ...)
    for(l in 1:length(links)) lines(consensus(gpa)[links[[l]],])
  }
  
  
  #Description
  # This function plot all the forms/shapes (landmark configurations) in the dataset in a single plot.
  
  #Arguments:
  # gpa= a k x p x n array of superimposed landmarks
  # links= optional. A list with the indices of the coordinates defining the wireframe (Morpho format)
  # mshape= whether to plot the mean configuration
  # ...= additional arguments passed to plot
  
}


## exp_var ---------------------------------------------------------------------------

exp_var<-function(pca, vars=NULL, digits=3, axes=10) {
  
  if(is.null(vars)==TRUE) vars<-pca$sdev^2
  exp<-round(((pca$sdev^2)/sum(vars))*100, digits=digits)
  cumexp<-round((cumsum(pca$sdev^2)/sum(vars))*100, digits=digits)
  return(data.frame(explained=exp[1:min(c(length(exp), axes))], cummulative=cumexp[1:min(c(length(cumexp), axes))]))
  
  
  #Description
  # A function that returns the explained variance and cummulative variance accounted
  # by each axis resulting from a PCA (prcomp format).
  
  #Arguments:
  # pca= the result of a PCA, in prcomp format
  # vars= a vector of variances representing the original variables. If NULL, pca$sdev^2 is used.
  # digits= number of decimal values for the output
  # axes=number of axes to deploy
  
}


## hulls_by_group_2D ---------------------------------------------------------------------

hulls_by_group_2D<-function(xy, fac, col=1:nlevels(fac)) {
  if(length(col)==1) col<-rep(col, nlevels(fac))
  
  for(i in 1:nlevels(fac)){
    x<-xy[fac==levels(fac)[i],1] ; y<- xy[fac==levels(fac)[i],2]
    hullp<-chull(x=x, y=y)
    polygon(x[hullp], y[hullp], border = col[i])}
  
  
  #Description
  # a function to plot convex hulls for different groups in 2D scatterplots created using plot.
  
  #Arguments:
  # xy= coordinates for the scatterplot
  # fac= a factor grouping data points
  # col= colors for each group 
  
}


## hulls_by_group_3D ---------------------------------------------------------------------

hulls_by_group_3D<-function(xyz, fac) {
  
  require(geometry)
  
  for(i in 1:nlevels(fac)){
    matsp<-xyz[fac==levels(fac)[i],1:3]
    surf<-t(convhulln(matsp))
    convex<-rgl.triangles(matsp[surf,1], matsp[surf,2], matsp[surf,3], col=i,alpha=.3)}
  
  
  #Description
  # a function to plot convex hulls for different groups in 3D scatterplots created using rgl's plot3d.
  
  #Arguments:
  # xyz= coordinates for the scatterplot
  # fac= a factor grouping data points
  # col= colors for each group 
  
}


# 2. ADAPTATIONS, WRAPPERS AND PLAGIARISM ##########################################################

## efourier_i2 -----------------------------------------------------------------------

efourier_i2<-function(coe, nb.p=120){
  nb.h<-length(coe)/4
  m_coo<-matrix(0, nrow=nb.p, ncol=2)
  
  a<-c(coe[1:nb.h])
  b<-c(coe[(nb.h+1):(nb.h*2)])
  c<-c(coe[((nb.h*2)+1):(nb.h*3)])
  d<-c(coe[((nb.h*3)+1):(nb.h*4)])
  coefs<-setNames(list(a,b,c,d),c("an","bn","cn","dn"))
  
  m_coo<-efourier_i(coefs, nb.pts=nb.p)
  return(m_coo)
  
  
  #Description
  # The function efourier_i2 performs the inverse Fourier function to transform Fourier 
  # coefficients into (x,y) coordinates. This function does the same as the efourier_i 
  # from Momocs, but deals with a vector of coefficients instead of a list.
  # Modified from Momocs.
  
  #Arguments:
  # coe= a vector with Fourier coefficients
  # nb.p= number of coordinates of the resulting (x,y) matrix of coordinates
}


## fourier2coo ------------------------------------------------------------------------

fourier2coo<-function(ef, n=120){
  cords<-array(0, c(n, 2, nrow(ef$coe)))
  for(i in 1:nrow(ef$coe)) cords[,,i]<-efourier_i2(ef$coe[i,], nb.p = n)
  return(Out(cords))
  
  
  #Description
  # function to transform a sample of closed outlines stored as Fourier coefficients in a
  # OutCoe object (from Momocs) to an array of (x,y) coordinates.
  
  #Arguments:
  # ef= an OutCoe object resulting from efourier
  # n= number of coordinates of the resulting (x,y) matrix of coordinates
}


## panel_ef ------------------------------------------------------------------------

panel_ef<-function(ef, ...){
  require(Momocs)
  cords<-array(0, c(120, 2, nrow(ef$coe)))
  for(i in 1:nrow(ef$coe)) cords[,,i]<-efourier_i2(ef$coe[i,])
  panel(Out(cords), ...)
  
  
  #Description
  # panel function adapted for OutCoe objects. Modified from Momocs.
  
  #Arguments:
  # ef= an OutCoe object
  # ...= additional arguments passed to panel
  
}


## edCurve_sample -------------------------------------------------------------------

edCurve_sample<-function(coord, n){
  newcoord<-array(NA, c(n, dim(coord)[2:3]))
  for(i in 1:dim(coord)[3]) newcoord[,,i]<-equidistantCurve(coord[,,i], n=n)
  return(newcoord)
  
  
  #Description
  # This function applies the Morpho function equidistantCurve to a sample of specimens.
  
  #Arguments:
  # coord= an array of coordinates defining a curve to be resampled
  # n= the number of desired coordinates
}


# edSurface_sample -------------------------------------------------------------------

edSurface_sample<-function(coord, n){
  newcoord<-array(NA, c(n, dim(coord)[2:3]))
  for(i in 1:dim(coord)[3]) newcoord[,,i]<-fastKmeans(coord[,,i], k=n)$centers
  return(newcoord)
  
  
  #Description
  # This function applies the Morpho function fastKmeans to a sample of specimens.
  
  #Arguments:
  # coord= an array of coordinates defining a curve to be resampled
  # n= the number of desired coordinates
}


## readland.tps_list -------------------------------------------------------------------

readland.tps_list<-function(tps.file, n, specID="imageID"){
  
  require(stringr)
  
  file<-scan(tps.file, what="character")
  
  starts<-which(str_detect(file, pattern="LM"))
  stops<-which(str_detect(file, pattern="ID"))
  
  if(specID=="imageID"){
    nameslist<-strsplit(file[which(str_detect(file, pattern="IMAGE"))], split="=")
    specname<-c() ; for(i in 1:length(nameslist)) specname[i]<-nameslist[[i]][2]}
  
  if(specID=="ID"){
    nameslist<-strsplit(file[which(str_detect(file, pattern="ID"))], split="=")
    specname<-c() ; for(i in 1:length(nameslist)) specname[i]<-nameslist[[i]][2]}
  
  confs<-list()  
  for(i in 1:n){
    spec_v<-as.vector(na.omit(as.numeric(file[(starts[i]):(stops[i])])))
    spec_m<-matrix(as.numeric(spec_v), ncol=2, nrow=(length(spec_v)/2), byrow=TRUE)
    confs[[i]]<-spec_m 
    names(confs)[i]<-specname[i]
  }
  
  nlands<-c()
  for(i in 1:length(confs)){
    nlands[i]<-nrow(confs[[i]])
  }
  nlands<-setNames(nlands, names(confs))
  return(list(tps=confs, nlands=nlands))
  
  
  #Description
  # Read tps file as a list. Useful when there are specimens that were mistakenly given a 
  # different number of landmarks.
  
  #Arguments:
  # tps.file= the path to the tps file
  # n= number of specimens
  # specID= ?readland.tps
}


# writeland.tps2 -----------------------------------------------------------------------

writeland.tps2<-function (A, file, curves, scale = NULL, specID = TRUE) 
{
  
  curvec<-c()
  for(i in 1:length(curves)) curvec<-c(curvec, curves[[i]])
  
  A2<-A
  A<-A[-curvec,,]
  
  n <- dim(A)[3]
  k <- dim(A)[2]
  p <- dim(A)[1]
  lmline <- ifelse(k == 2, paste("LM=", p, sep = ""), 
                   paste("LM3=", p, sep = ""))
  file.create(file, showWarnings = TRUE)
  if (!is.null(scale)) {
    scaleline <- paste("SCALE", "=", scale, sep = "")
  }
  
  if(!is.null(curves)){
    
    ncurveslines<-paste("CURVES=", length(curves), sep="")
    curvelines<-c()
    for(j in 1:length(curves)) curvelines<-c(curvelines, paste("POINTS=", length(curves[[j]]), sep=""))
    
  }
  
  for (i in 1:n) {
    write(lmline, file, append = TRUE)
    write.table(A[, , i], file, col.names = FALSE, row.names = FALSE, 
                append = TRUE)
    if(!is.null(curves)){
      write(ncurveslines, file, append = TRUE)
      
      for(j in 1:length(curves)){
        write(curvelines[[j]], file, append = TRUE)
        write.table(A2[curves[[j]], , i], file, col.names = FALSE, row.names = FALSE, 
                    append = TRUE)
      }
    }
    
    if (!is.null(scale)) {
      if (length(scaleline) == 1) {
        write(scaleline, file, append = TRUE)
      }
      if (length(scaleline) > 1) {
        write(scaleline[i], file, append = TRUE)
      }
    }
    
    if (specID == TRUE) {
      if (is.null(dimnames(A)[[3]])) 
        dimnames(A)[[3]] <- c(1:dim(A)[3])
      idline <- paste("ID=", dimnames(A)[[3]][i], 
                      sep = "")
      write(idline, file, append = TRUE)
    }
    write("", file, append = TRUE)
  }
  
  
  #Description
  # The writeland.tps from geomorph function but modified to allow for inclusion of curves in the
  # right tps format.
  
  #Arguments
  # A, file, curves, scale, specID= see ?writeland.tps
  # curves= a list of vectors indicating the coordinates to be trated as cuvres (one or many separate ones)
}


# read.mpp_all ----------------------------------------------------------------------------

read.mpp_all<-function(dir, list=FALSE){
  require(abind)
  require(Morpho)
  
  files<-list.files(dir)
  pps<-list() ; for(i in 1:length(files)) pps[[i]]<-read.mpp(paste(dir, files[i], sep="/"))
  names(pps)<-files
  arr<-abind(pps, along=3)
  if(list==FALSE) return(invisible(arr)) else invisible(return(pps))
  
  
  #Description
  # A function relying on read.mpp function from Morpho for importing a set of 3D landmarks stored
  # as meshlab .pp format. Returns an array or a list with landmark configurations.
  
  #Arguments: 
  # dir= a directory in which ONLY .pp meshes are stored
  # list= logical. Whether return a list (default FALSE imports landmarks as an array)
}


## tps2d -----------------------------------------------------------------------------

tps2d<-function(M, matr, matt)
{p<-dim(matr)[1]; q<-dim(M)[1]; n1<-p+3
P<-matrix(NA, p, p)
for (i in 1:p)
{for (j in 1:p){
  r2<-sum((matr[i,]-matr[j,])^2)
  P[i,j]<- r2*log(r2)}}
P[which(is.na(P))]<-0
Q<-cbind(1, matr)
L<-rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
m2<-rbind(matt, matrix(0, 3, 2))
coefx<-solve(L)%*%m2[,1]
coefy<-solve(L)%*%m2[,2]
fx<-function(matr, M, coef)
{Xn<-numeric(q)
for (i in 1:q)
{Z<-apply((matr-matrix(M[i,],p,2,byrow=T))^2,1,sum)
Xn[i]<-coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+sum(coef[1:p]*(Z*log(Z)))}
Xn}
matg<-matrix(NA, q, 2)
matg[,1]<-fx(matr, M, coefx)
matg[,2]<-fx(matr, M, coefy)
matg


#Description
# A function to perform two-dimensional TPS warping interpolation (for curves) from the
# transformation between a reference and a target shape. Taken directly from Claude (2008).

#Arguments:
# M= matrix with coordinates defining a curve or set of curves which will be warped
# matr= reference configuration
# matt= target configuration

}


## deformGrid2d_amp ---------------------------------------------------------------

deformGrid2d_amp<-function(matrix, tarmatrix, amp=1, ...){
  
  require(geomorph)
  require(Morpho)
  require(abind)
  
  mat2d<-two.d.array(abind:::abind(matrix, tarmatrix, along=3))
  pca<-prcomp(mat2d)
  
  scores<-range(pca$x[,1])*amp
  v<-pca$rotation[,1]
  c<-pca$center
  X<-arrayspecs(t(t(scores%*%t(v))+c), p = nrow(matrix), k = 2)
  
  deformGrid2d(matrix=X[,,1], tarmatrix=X[,,2], ...)
  
  
  #Description
  # This function calls deformGrid2d from Morpho to compare two shapes, but allows
  # amplification of the implied transformation using PCA.
  
  #Arguments
  # matrix= reference matrix with shape coordinates
  # tarmatrix= target matrix with shape coordinates
  # amp= magnification factor
  # ... = further arguments passed to deformGrid2d
  
}


## correct.fourier ------------------------------------------------------------------

correct.fourier<-function(ef){
  print("Click Finish to finish selection")
  options(warn=-1)
  pos<-coo_listpanel(coo_template(fourier2coo(ef))$coo)
  panel_ef(ef)
  
  nb.h<-ncol(ef$coe)/4
  n<-length(ef)
  sel <- rep(FALSE, length(ef))
  while(sum(sel) < n) {
    ans <- identify(pos[!sel,], n = 1, plot = FALSE)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    sel[ans] <- TRUE
    
    a<-c(ef$coe[ans,1:nb.h])
    b<-c(ef$coe[ans,(nb.h+1):(nb.h*2)])
    c<-c(ef$coe[ans,((nb.h*2)+1):(nb.h*3)])
    d<-c(ef$coe[ans,((nb.h*3)+1):(nb.h*4)])
    
    a[1:nb.h %% 2 == 0]<-a[1:nb.h %% 2 == 0]*-1
    b[1:nb.h %% 2 == 0]<-b[1:nb.h %% 2 == 0]*-1
    c[1:nb.h %% 2 == 0]<-c[1:nb.h %% 2 == 0]*-1
    d[1:nb.h %% 2 == 0]<-d[1:nb.h %% 2 == 0]*-1
    
    ef$coe[ans,]<-c(a,b,c,d)
    
    panel_ef(ef)
  }
  options(warn=0)
  return(ef)
  
  
  #Description
  # This function allows to interactively correct the 180° rotation of shapes after EFA. Ispired in SHAPE
  # (Iwata and Ukai 2002)
  
  #Arguments:
  # ef= an OutCoe object resulting from efourier
}


# 3. SHAPE VECTORS ######################################################################

## shapes_vector -------------------------------------------------------------------------

shapes_vector<-function(M, template=NULL, ax=1, amp=2, links=NULL, landmarks=TRUE, pointcols=c(1,2), linecols=c(1,2), lwd=2)
{
  require(geomorph)
  
  p<-ncol(M)/2
  pca<-prcomp(M); 
  X<-pca$x
  v<-pca$rotation
  c<-pca$center
  
  PC.sample<-matrix(0, ncol=min(nrow(M),ncol(M)), nrow=2)
  PC.sample[1,ax]<-0+sd(pca$x[,ax])
  PC.sample[2,ax]<-0-sd(pca$x[,ax])
  PC.amp<-PC.sample*amp
  
  shp.sample<-array(0, c(nrow=p, ncol=2, 2))
  for(i in 1:2) shp.sample[,,i]<-arrayspecs(t(t(PC.amp[i,]%*%t(v))+c), p, 2)
  
  plot(rbind(shp.sample[,,1],shp.sample[,,2]), type="n", axes=FALSE, xlab="", ylab="")
  legend(x="topleft", legend=paste(c(paste("-",amp,"sd", sep=c("", " ")), paste("+",amp,"sd", sep=c("", " "))), paste("PC",ax, sep="")), fill=pointcols, bty="n")
  
  
  if(is.null(template)==TRUE&is.null(links)==FALSE) {
    for(i in 1:nrow(links)) lines(rbind(shp.sample[links[i,1],,1], shp.sample[links[i,2],,1]), col=linecols[1], lwd=lwd)
    for(i in 1:nrow(links)) lines(rbind(shp.sample[links[i,1],,2], shp.sample[links[i,2],,2]), col=linecols[2], lwd=lwd)
  }
  
  if(landmarks==TRUE) for(i in 1:2) points(shp.sample[,,i], pch=21, bg=pointcols[i], cex=1.2)
  
  if(is.null(template)==FALSE) {
    consensus<-arrayspecs(t(t(rep(0, ncol(PC.amp))%*%t(v))+c), p, 2)[,,1]
    warpedconsensus<-rbind(consensus, tps2d(M=template[-(1:p),], template[(1:p),], consensus))
    warped.shp1<-tps2d(M=warpedconsensus[-(1:p),], warpedconsensus[(1:p),], shp.sample[,,1])
    warped.shp2<-tps2d(M=warpedconsensus[-(1:p),], warpedconsensus[(1:p),], shp.sample[,,2])
    
    lines(warped.shp1, lwd=lwd, col=linecols[1])
    lines(warped.shp2, lwd=lwd, col=linecols[2])
  }
  
  
  #Description
  # This function plots the extreme shapes from an axis resulting from a PCA of normalized 
  # landmark data, including the interpolated warped curves (if fed with).
  
  #Arguments:
  # M= a two-dimensional array (i.e., a n x (p x 2) matrix) with normalizaed configurations
  # ax= the PC axes to be represented
  # template= matrix with coordinates defining a curve or set of curves which will be warped,
  #           as well as the landmarks what will be used to interpolate
  # links= wireframe
  # amp=  amplification factor for shape transformation, in standard deviations
  # linecols= wireframe/warped curves color
  # pointcols= landmark colors
  # landmarks= whether to plot landmarks
  # lwd= line width
}


# shapes_vector_3D ----------------------------------------------------------------------------

shapes_vector_3D<-function(pcscores, rotation, center, ax=1, amp=1, sample=2, indep=FALSE, 
                           plot=TRUE, col="gray", points=FALSE, surface=TRUE, ...){
  
  require(akima)
  require(Morpho)
  
  s<-pcscores
  u<-rotation
  c<-center[1:nrow(u)]
  
  sd<-sd(s[,ax])
  lims<-c(-sd,sd)*amp
  range<-seq(from=lims[1], to=lims[2], length.out=sample)
  p<-dim(u)[1]/3
  
  shapes<-array(0, c(p,3,sample))
  for(i in 1:sample){
    PCsample<-rep(0, ncol(u)) ; PCsample[ax]<-range[i]
    shapes[,,i]<-arrayspecs(t(t(PCsample%*%t(u))+c), p, 3)
  }
  
  if(plot==TRUE){
    open3d()
    if(indep==TRUE) mfrow3d(nr=1,nc=sample)
    if(indep==FALSE) layout3d(matrix(1:sample,nrow=1,ncol=sample),sharedMouse = TRUE)
    
    
    if(surface==TRUE){
      for(i in 1:sample){
        next3d()
        model<-shapes[,,i]
        x<-model[,1] ; y<-model[,2] ; z<-model[,3]
        surf<-interp(x,y,z, nx=p, ny=p)
        surface3d(surf$x,surf$y,surf$z, col=col)}
    }
    
    if(points==TRUE){
      for(i in 1:sample){
        next3d()
        plot3d(shapes[,,i], add=TRUE, col=col, ...)}  
      
    }
  }
  
  invisible(shapes)
  
  
  #Description
  # A function to plot the theoretical tridimensional shapes sampling the shape variation of 
  # a PCA. It will also return the configurations for those shapes.
  
  #Arguments:
  # pcscores= scores resulting from the PCA
  # rotation= rotation matrix from the PCA
  # center= the means from the original variables
  # ax= the desired axis
  # amp= amplification factor for shapes
  # sample= number of shapes to draw from the axis at regular intervals
  # indep= should the windows be manipulated independently
  # plot= whether to plot the results
  # points= if TRUE, will plot the semilandmarks
  # surface= if TRUE, will plot the surface
  # ...= further arguments passed to plot3d
}


## plot_allometry_ef ------------------------------------------------------------------------

plot_allometry_ef<-function(ef, size, amp.shp=1, amp.sz=2, plotPLS=FALSE){
  require(Morpho)
  
  pls<-pls2B(y=ef$coe, x=size)
  
  if(plotPLS==TRUE) {
    plot(pls$Xscores,pls$Yscores, cex=size/max(size)*amp.sz, ylab="Y-PLS", xlab="Size")
    abline(lm(pls$Yscores~pls$Xscores), col="red")  
  }
  
  if(sign(cor(size, pls$Yscores))==-1){
    small<-max(pls$Yscores)*amp.shp
    big<-min(pls$Yscores)*amp.shp
  } else {
    small<-min(pls$Yscores)*amp.shp
    big<-max(pls$Yscores)*amp.shp
  }
  
  big.shp<-efourier_i2(big%*%t(pls$svd$v)+apply(X=ef$coe, FUN=mean, MARGIN=2))
  small.shp<-efourier_i2(small%*%t(pls$svd$v)+apply(X=ef$coe, FUN=mean, MARGIN=2))
  
  plot(small.shp, type="l", lwd=2, axes=F, xlab="", ylab="", col="gray")
  points(big.shp, type="l", lwd=2)
  legend(x="bottomright", legend=c("grande","pequeño"), fill=c("black", "gray"), bty="n")
  
  invisible(abind(small.shp, big.shp, along=3))
  if(plotPLS==TRUE) cat("First plot: PLS\nSecond plot: shapes") else print("Plot: Shapes") 
  
  
  #Description
  # The plot_allometry_ef function plots the extreme shapes of the pls-y-scores resulting from 
  # the PLS between a matrix with the normalized fourier coefficients and size.
  
  #Arguments:
  # ef= OutCoe object resulting from the efourier function (from Momocs)
  # size= a numeric vector of sizes of length equal to nrow(ef$coe)
  # amp.shp= amplification factor for the shape transformation
  # amp.sz= amplification factor for plotting size (if plotPLS=TRUE)
  # plotPLS= should the function plot the PLS scatterplot?
  
}


## plot_allometry_ldm -----------------------------------------------------------------------

plot_allometry_ldm<-function(gpa, size, amp.shp=1, amp.sz=2, links=NULL, plotSHPs=TRUE, plotPLS=FALSE, psize=0.1, separated=TRUE, indep=FALSE){
  
  require(Morpho)
  
  consensus<-consensus(gpa)
  pls<-pls2B(y=two.d.array(gpa), x=size)
  
  if(plotPLS==TRUE) {
    plot(pls$Xscores,pls$Yscores, cex=size/max(size)*amp.sz, ylab="Y-PLS", xlab="Size")
    abline(lm(pls$Yscores~pls$Xscores), col="red")  
  }
  
  if(sign(cor(size, pls$Yscores))==-1){
    small<-max(pls$Yscores)*amp.shp
    big<-min(pls$Yscores)*amp.shp
  } else {
    small<-min(pls$Yscores)*amp.shp
    big<-max(pls$Yscores)*amp.shp
  }
  
  big.shp<-arrayspecs(big*t(pls$svd$v)+apply(X=two.d.array(gpa), FUN=mean, MARGIN=2), p=dim(gpa)[[1]], k=dim(gpa)[[2]])[,,1]
  small.shp<-arrayspecs(small*t(pls$svd$v)+apply(X=two.d.array(gpa), FUN=mean, MARGIN=2), p=dim(gpa)[[1]], k=dim(gpa)[[2]])[,,1]
  
  if(dim(gpa)[[2]]==2 & plotSHPs==TRUE){
    plot(small.shp, type="p", lwd=2, axes=F, xlab="", ylab="", col="gray")
    if(is.null(links)==FALSE) for(i in 1:nrow(links)) lines(rbind(small.shp[links[i,1],], small.shp[links[i,2],]), col="gray")
    points(big.shp, type="p", lwd=2)
    if(is.null(links)==FALSE) for(i in 1:nrow(links)) lines(rbind(big.shp[links[i,1],], big.shp[links[i,2],]), col="black")
    legend(x="bottomright", legend=c("grande","pequeño"), fill=c("black", "gray"), bty="n")
  }
  
  if(dim(gpa)[[2]]==3 & plotSHPs==TRUE){
    
    if(separated==TRUE){
      open3d()
      if(indep==TRUE) mfrow3d(nr=1,nc=2)
      if(indep==FALSE) layout3d(matrix(1:2,nrow=1,ncol=2),sharedMouse = TRUE)
      
      propx <- max(consensus[,1])-min(consensus[,1])
      propy <- max(consensus[,2])-min(consensus[,2])
      propz <- max(consensus[,3])-min(consensus[,3])
      
      small.shp<-pcAlign(x=small.shp, y=consensus)
      big.shp<-pcAlign(x=big.shp, y=consensus)
      
      plot3d(small.shp, type="s", aspect=c(propx,propy,propz), size=psize, axes=F, xlab="", ylab="", col="darkgray")
      if(is.null(links)==FALSE) for(i in 1:nrow(links)) lines3d(rbind(small.shp[links[i,1],], small.shp[links[i,2],]), col="gray")
      plot3d(big.shp, type="s", aspect=c(propx,propy,propz), size=psize, axes=F, xlab="", ylab="", col="black")
      if(is.null(links)==FALSE) for(i in 1:nrow(links)) lines3d(rbind(big.shp[links[i,1],], big.shp[links[i,2],]), col="black")
      legend3d(x="bottomright", legend=c("grande","pequeño"), fill=c("black", "gray"), bty="n")
      
      legend3d(x="bottomright", legend=c("grande","pequeño"), fill=c("black", "gray"), bty="n")}
    
    
    if(separated==FALSE){
      
      propx <- max(consensus[,1])-min(consensus[,1])
      propy <- max(consensus[,2])-min(consensus[,2])
      propz <- max(consensus[,3])-min(consensus[,3])
      
      small.shp<-pcAlign(x=small.shp, y=consensus)
      big.shp<-pcAlign(x=big.shp, y=consensus)
      
      plot3d(small.shp, type="s", aspect=c(propx,propy,propz), size=psize, axes=F, xlab="", ylab="", col="darkgray")
      if(is.null(links)==FALSE) for(i in 1:nrow(links)) lines3d(rbind(small.shp[links[i,1],], small.shp[links[i,2],]), col="gray")
      plot3d(big.shp, add=TRUE, type="s", aspect=c(propx,propy,propz), size=psize, axes=F, xlab="", ylab="", col="black")
      if(is.null(links)==FALSE) for(i in 1:nrow(links)) lines3d(rbind(big.shp[links[i,1],], big.shp[links[i,2],]), col="black")
      legend3d(x="bottomright", legend=c("grande","pequeño"), fill=c("black", "gray"), bty="n")}
  }
  
  
  if(plotPLS==TRUE) cat("First plot: PLS\nSecond plot: shapes") else print("Plot: Shapes") 
  return(invisible(abind(small.shp, big.shp, along=3)))
  
  
  #Description
  # The plot_allometry_ldm function plots and returns the extreme shapes of the pls-y-scores 
  # resulting from the PLS between a matrix with the superimposed landmark configurations 
  # and size. Works for 2D and 3D landmark data.
  
  #Arguments:
  # gpa= an array of normalized landmark configurations
  # size= a numeric vector of sizes of length equal to dim(gpa)[[3]]
  # amp.shp= amplification factor for the shape transformation
  # amp.sz= amplification factor for plotting size (if plotPLS=TRUE)
  # plotPLS= should the function plot the PLS scatterplot?
  # plotSHPs= should the function plot extreme shapes?
  
  
}


# 5. OPERATIONS ###########################################################################

## optim_curve -----------------------------------------------------------------------------

optim_curve<-function(m, acc=0.99){
  
  require(Morpho)
  
  totdist<-sum(sapply(2:nrow(m), function(i) {dist(rbind(m[i-1,], m[i,]))}))
  tol<-totdist*(1-acc)
  newdist<-0
  n<-2
  
  while((totdist-newdist)>tol){
    subm<-equidistantCurve(m, n=n)
    newdist<-sum(sapply(2:nrow(subm), function(i) {dist(rbind(subm[i-1,], subm[i,]))}))
    n<-n+1
  }
  
  return(nrow(subm))
  
  
  #Description
  # A function to compute the optimal length of a curve defined by a matrix of coordinates, 
  # i.e., the number of coordinates retaining a certain percentage of the total original distance.
  
  #Arguments:
  # m= a matrix defining a curve (belonging to a single specimen) to be assessed
  # acc= desired accuracy which has to be preserved by the optimal number of coordinates
}


# optim_surf -----------------------------------------------------------------------------

optim_surf<-function(m, acc=0.99, rate=10, n=5){
  
  require(Morpho)
  require(geometry)
  
  totarea<-convhulln(m, output.options=TRUE)$area
  tol<-totarea*(1-acc)
  newarea<-0
  
  while((totarea-newarea)>tol){
    subm<-fastKmeans(m, k=n)$centers
    newarea<-convhulln(subm, output.options=TRUE)$area
    n<-n+rate
  }
  
  return(nrow(subm))
  
  
  #Description
  # A function to compute the optimal number of coordinates to cover a surface defined by a matrix of 
  # coordinates, i.e., the number of coordinates retaining a certain percentage of the total original area.
  
  #Arguments:
  # m= a matrix defining a curve (belonging to a single specimen) to be assessed
  # acc= desired accuracy which has to be preserved by the optimal number of coordinates
  # rate= by how many coordinates should the function compute the area? 
  # n= initial number of coordinates to sample
}


## detrend_shapes ---------------------------------------------------------------------

detrend_shapes<-function(model, xvalue=NULL){
  
  n<-nrow(model$residuals)
  coefs<-model$coefficients
  resids<-model$resid
  
  x<-model$model[,ncol(model$model)]
  
  if(is.null(xvalue)) {
    predicted<-colMeans(model$residuals+model$fitted.values)
    
    predicted_vec<-rep(1, n)%*%t(predicted)
    predicted_mat<-resids+predicted_vec
  } 
  
  if(!is.null(xvalue)) {
    
    if(is.numeric(x)==TRUE){ 
      designmat<-cbind(1, xvalue)
    }
    
    if(is.factor(x)==TRUE){
      designmat<-rep(0, nlevels(x))
      designmat[which(levels(x)==xvalue)]<-1
      if(isFALSE(designmat[1]==1)) designmat[1]<-1
    }
    
    predicted<-as.numeric(designmat%*%coefs)
    
    predicted_vec<-rep(1, n)%*%t(predicted)
    predicted_mat<-resids+predicted_vec
  }
  
  return(predicted_mat)
  
  
  #Description
  # A function to detrendize (i.e., standardize) shape data fitted with a linear model to
  # some external variable (works with both factors and numerics).
  
  #Arguments:
  # model= an object fitted using lm
  # xvalue= value (numeric) or level (character) at which shape data is to be standardized; 
  #         If NULL, the mean of the complete sample is used instead.
  
}


## predict_shapes ---------------------------------------------------------------------

predict_shapes<-function(model, xvalue) {
  
  coefs<-model$coefficients
  designmat<-cbind(1, xvalue)
  predicted_vec<-as.numeric(designmat%*%coefs)
  return(predicted_vec)
  
  
  #Description
  # A function to predict shapes using the linear relationship between shape data and some
  # external variable, fitted using a linear model.
  
  #Arguments:
  # model= an object fitted using lm
  # xvalue= value (numeric) or level (character) at which shape data is to be standardized 
  
}


# 6. MORPHOSPACES ##################################################################################

## mspace1 ---------------------------------------------------------------------------

mspace1<-function(M, axes=c(1,2), links, amp=1, n=5, rotate=0, s=1, adj=0.9, PCA="regular", invertax=NULL, aspX=1, 
                  pch=16, col="black", model.col="#708095", model.lwd=1.4, ldm.cex=0.2, main="", points=T, xlim=NULL, ylim=NULL, ...){
  
  require(geomorph)
  require(spdep)
  
  p<-ncol(M)/2
  
  if(PCA=="regular") {
    pca<-prcomp(M)
  }
  
  if(PCA=="bgPCA") {
    require(Morpho)
    bgpca<-groupPCA(dataarray = M, ...) ; pca<-list()
    pca$rotation<-bgpca$groupPCs
    pca$x<-bgpca$Scores
    pca$center<-colMeans(M)
  }
  
  if(PCA=="phyloPCA") {
    require(phytools)
    ppca<-phyl.pca(Y = M, ...) ; pca<-list()
    pca$rotation<-ppca$Evec
    pca$x<-ppca$S
    pca$center<-colMeans(M)
  }
  
  X<-pca$x ; X[,axes[invertax]]<-X[,axes[invertax]]*-1
  v<-pca$rotation ; v[,axes[invertax]]<-v[,axes[invertax]]*-1
  c<-pca$center
  
  if(is.null(xlim)) xlim<-range(X[,axes[1]])
  if(is.null(ylim)) ylim<-range(X[,axes[2]])
  
  plot(X[,axes[1]], X[,axes[2]], xlim=xlim, ylim=ylim,  main=main,
       xlab=paste("PC", axes[1], sep=""), ylab=paste("PC", axes[2], sep=""), col="white")
  
  m<-matrix(0, nrow=n*n, ncol=length(v[1,]))
  m[,axes[1]]<-rep(seq(from=min(xlim), to=max(xlim), length.out=n), times=c(rep(n, n)))
  m[,axes[2]]<-rep(seq(from=min(ylim), to=max(ylim), length.out=n), n)
  
  PC.sample<-matrix(0, ncol=nrow(M), nrow=length(m[,1]))
  for(i in 1:length(m[,1])){
    PC.sample[i,axes[1]]<-m[i,axes[1]]
    PC.sample[i,axes[2]]<-m[i,axes[2]]}
  
  PC.amp<-PC.sample*amp
  
  if(ncol(PC.amp)>=nrow(v)){PC.amp<-PC.amp[,1:nrow(v)]}
  if(ncol(PC.amp)>=ncol(v)){PC.amp<-PC.amp[,1:ncol(v)]}
  
  shp.sample<-array(0, c(nrow=p, ncol=2, length(m[,1])))
  for(i in 1:length(m[,1])){
    shp.sample[,,i]<-arrayspecs(t(t(PC.amp[i,]%*%t(v))+c), p, 2)
    shp.sample[,,i]<-Rotation(shp.sample[,,i], rotate*0.0174532925199433)
  }
  
  scores.ext<-matrix(0, nrow=4, ncol=ncol(PC.amp))
  for(i in 1:2) scores.ext[((i*2)-1):(i*2),axes[i]]<-range(PC.amp[,axes[i]])
  
  shapes.ext<-array(0, c(nrow=p, ncol=2, 4))
  for(i in 1:nrow(scores.ext)) shapes.ext[,,i]<-arrayspecs(t(t(scores.ext[i,]%*%t(v))+c), p, 2)
  
  dimnames(shapes.ext)[[3]]<-c(paste("PC", axes[1], "min", sep="_"), paste("PC", axes[1], "max", sep="_"),
                               paste("PC", axes[2], "min", sep="_"), paste("PC", axes[2], "max", sep="_"))
  
  for(i in 1:length(m[,1])){
    shp.sample[,1,i]<-shp.sample[,1,i]*aspX
    shp<-shp.sample[,,i]*0.1*s+matrix(rep(PC.sample[i,axes],p),p,2,byrow=T)*adj
    points(shp, pch=16, cex=ldm.cex)
    for(l in 1:length(links)) lines(shp[links[[l]],], col=model.col, lwd=model.lwd)
  }
  
  if(points==T){points(X[,axes], pch=pch, col=col)}
  
  return(invisible(list(PCA=pca, ext.shapes.PC1=shapes.ext[,,1:2], ext.shapes.PC2=shapes.ext[,,3:4])))
  
  
  #Description
  # The mspace1 function plots the first two axes resulting from a principal component analysis
  # (regular, between groups and phylogenetic) and the theoretical shapes (landmark configurations 
  # and associated wireframes/links) representing shape variation in this space. Returns a list with
  # the PCA in prcomp format (scores, eigenvectors, center), and the negative and positive extreme
  # shapes from each PC.
  
  #Arguments:
  # M= a two-dimensional array (i.e., a n x (p x 2) matrix) with normalizaed configurations
  # axes= a numeric vector indicating the PC axes to be represented
  # links= a list with a vector or vectors defining a wireframe
  # amp= numeric; amplifying factor for shape models
  # n= numeric; density of shape models in the plot
  # rotate= numeric; angle (in degrees) to rotate the shape models
  # s= numeric; scaling factor for shape models
  # adj= numeric; ajusting factor for shape models position
  # PCA= character; the function to be used for PCA. if "bgPCA", the groupPCA function from Morpho 
  #      is used; if "phyloPCA", the phyl.pca function from phytools is used instead
  # invertax= numeric; which axis (of those defined by axes) should be inverted?
  # aspX= numeric; a factor for compressing/depressing shape models along the X axis
  # pch= symbol of the scatter points
  # col= color of scatter points
  # model.col= color of the lines for background wireframe models
  # model.lwd= width of the lines for background wireframe models
  # ldm.cex= numeric; size of the points representing the landmark configurations
  # points= logical; wheter to plot the PC scores
  # xlim, ylim, main= standard arguments for the generic plot function
  # ...= further arguments passed to groupPCA (groups) or phyl.pca (tree)
  
}


## mspace2 ---------------------------------------------------------------------------

mspace2<-function(M, axes=c(1,2), amp=1, n=5, rotate=0, s=1, adj=0.9, PCA="regular", invertax=NULL, aspX=1, 
                  pch=16, col="black", model.col="#708095", model.lwd=1, main="", points=T, xlim=NULL, ylim=NULL, ...){
  
  require(geomorph)
  require(spdep)
  
  p<-ncol(M)/2
  
  if(PCA=="regular") {
    pca<-prcomp(M)
  }
  
  if(PCA=="bgPCA") {
    require(Morpho)
    bgpca<-groupPCA(dataarray = M, ...) ; pca<-list()
    pca$rotation<-bgpca$groupPCs
    pca$x<-bgpca$Scores
    pca$center<-colMeans(M)
  }
  
  if(PCA=="phyloPCA") {
    require(phytools)
    ppca<-phyl.pca(Y = M, ...) ; pca<-list()
    pca$rotation<-ppca$Evec
    pca$x<-ppca$S
    pca$center<-colMeans(M)
  }
  
  X<-pca$x ; X[,axes[invertax]]<-X[,axes[invertax]]*-1
  v<-pca$rotation ; v[,axes[invertax]]<-v[,axes[invertax]]*-1
  c<-pca$center
  
  if(is.null(xlim)) xlim<-range(X[,axes[1]])
  if(is.null(ylim)) ylim<-range(X[,axes[2]])
  
  plot(X[,axes[1]], X[,axes[2]], xlim=xlim, ylim=ylim,  main=main,
       xlab=paste("PC", axes[1], sep=""), ylab=paste("PC", axes[2], sep=""), col="white")
  
  m<-matrix(0, nrow=n*n, ncol=length(v[1,]))
  m[,axes[1]]<-rep(seq(from=min(xlim), to=max(xlim), length.out=n), times=c(rep(n, n)))
  m[,axes[2]]<-rep(seq(from=min(ylim), to=max(ylim), length.out=n), n)
  
  PC.sample<-matrix(0, ncol=nrow(M), nrow=length(m[,1]))
  for(i in 1:length(m[,1])){
    PC.sample[i,axes[1]]<-m[i,axes[1]]
    PC.sample[i,axes[2]]<-m[i,axes[2]]}
  
  PC.amp<-PC.sample*amp
  
  if(ncol(PC.amp)>=nrow(v)){PC.amp<-PC.amp[,1:nrow(v)]}
  if(ncol(PC.amp)>=ncol(v)){PC.amp<-PC.amp[,1:ncol(v)]}
  
  shp.sample<-array(0, c(nrow=p, ncol=2, length(m[,1])))
  for(i in 1:length(m[,1])){
    shp.sample[,,i]<-arrayspecs(t(t(PC.amp[i,]%*%t(v))+c), p, 2)
    shp.sample[,,i]<-Rotation(shp.sample[,,i], rotate*0.0174532925199433)
  }
  
  scores.ext<-matrix(0, nrow=4, ncol=ncol(PC.amp))
  for(i in 1:2) scores.ext[((i*2)-1):(i*2),axes[i]]<-range(PC.amp[,axes[i]])
  
  shapes.ext<-array(0, c(nrow=p, ncol=2, 4))
  for(i in 1:nrow(scores.ext)) shapes.ext[,,i]<-arrayspecs(t(t(scores.ext[i,]%*%t(v))+c), p, 2)
  
  dimnames(shapes.ext)[[3]]<-c(paste("PC", axes[1], "min", sep="_"), paste("PC", axes[1], "max", sep="_"),
                               paste("PC", axes[2], "min", sep="_"), paste("PC", axes[2], "max", sep="_"))
  
  for(i in 1:length(m[,1])){
    shp.sample[,1,i]<-shp.sample[,1,i]*aspX
    lines(shp.sample[,,i]*0.1*s+matrix(rep(PC.sample[i,axes],p),p,2,byrow=T)*adj, col=model.col, lwd=model.lwd)
  }
  
  if(points==T){points(X[,axes], pch=pch, col=col)}
  
  return(invisible(list(PCA=pca, ext.shapes.PC1=shapes.ext[,,1:2], ext.shapes.PC2=shapes.ext[,,3:4])))
  
  
  #Description
  # The mspace2 function plots the first two axes resulting from a principal component analysis
  # (regular, between groups and phylogenetic) and the theoretical shapes (open curves resulting from 
  # a semilandmark analysis) and associated wireframes/links) representing shape variation in this space. 
  # Returns a list with the PCA in prcomp format (scores, eigenvectors, center), and the negative and 
  # positive extreme shapes from each PC.
  
  #Arguments:
  # M= a two-dimensional array (i.e., a n x (p x 2) matrix) with normalizaed configurations
  # axes= a numeric vector indicating the PC axes to be represented
  # amp= numeric; amplifying factor for shape models
  # n= numeric; density of shape models in the plot
  # rotate= numeric; angle (in degrees) to rotate the shape models
  # s= numeric; scaling factor for shape models
  # adj= numeric; ajusting factor for shape models position
  # PCA= character; the function to be used for PCA. if "bgPCA", the groupPCA function from Morpho 
  #      is used; if "phyloPCA", the phyl.pca function from phytools is used instead
  # invertax= numeric; which axis (of those defined by axes) should be inverted?
  # aspX= numeric; a factor for compressing/depressing shape models along the X axis
  # pch= symbol of the scatter points
  # col= color of scatter points
  # model.col= color of the lines for background outlines models
  # model.lwd= width of the lines for background outlines models
  # points= logical; wheter to plot the PC scores
  # xlim, ylim, main= standard arguments for the generic plot function
  # ...= further arguments passed to groupPCA (groups) or phyl.pca (tree)
  
}


## mspace3 ---------------------------------------------------------------------------

mspace3<-function(ef, axes=c(1,2), amp=1, n=5, rotate=0, s=0.1, adj=0.9, PCA="regular", invertax=NULL, aspX=1, 
                  pch=16, col="black", model.col="#708095", model.lwd=1, model.p=300, main="", xlim=NULL, ylim=NULL, points=T, ...){
  
  require(geomorph)
  require(spdep)
  
  p<-ncol(M)/2
  
  if(PCA=="regular") {
    pca<-prcomp(M)
  }
  
  if(PCA=="bgPCA") {
    require(Morpho)
    bgpca<-groupPCA(dataarray = M, ...) ; pca<-list()
    pca$rotation<-bgpca$groupPCs
    pca$x<-bgpca$Scores
    pca$center<-colMeans(M)
  }
  
  if(PCA=="phyloPCA") {
    require(phytools)
    ppca<-phyl.pca(Y = M, ...) ; pca<-list()
    pca$rotation<-ppca$Evec
    pca$x<-ppca$S
    pca$center<-colMeans(M)
  }
  
  X<-pca$x ; X[,axes[invertax]]<-X[,axes[invertax]]*-1
  v<-pca$rotation ; v[,axes[invertax]]<-v[,axes[invertax]]*-1
  c<-pca$center
  
  if(is.null(xlim)) xlim<-range(X[,axes[1]])
  if(is.null(ylim)) ylim<-range(X[,axes[2]])
  
  plot(X[,axes[1]], X[,axes[2]], xlim=xlim, ylim=ylim,  main=main,
       xlab=paste("PC", axes[1], sep=""), ylab=paste("PC", axes[2], sep=""), col="white")
  
  m<-matrix(0, nrow=n*n, ncol=length(v[1,]))
  m[,axes[1]]<-rep(seq(from=min(xlim), to=max(xlim), length.out=n), times=c(rep(n, n)))
  m[,axes[2]]<-rep(seq(from=min(ylim), to=max(ylim), length.out=n), n)
  
  PC.sample<-matrix(0, ncol=nrow(M), nrow=length(m[,1]))
  for(i in 1:length(m[,1])){
    PC.sample[i,axes[1]]<-m[i,axes[1]]
    PC.sample[i,axes[2]]<-m[i,axes[2]]}
  
  PC.amp<-PC.sample*amp
  
  if(ncol(PC.amp)>=nrow(v)){PC.amp<-PC.amp[,1:nrow(v)]}
  if(ncol(PC.amp)>=ncol(v)){PC.amp<-PC.amp[,1:ncol(v)]}  
  
  shp.sample<-array(0, c(nrow=p, ncol=2, length(m[,1])))
  for(i in 1:length(m[,1])){
    shp.sample[,,i]<-arrayspecs(efourier_i2(t(t(PC.amp[i,]%*%t(v))+c), model.p), model.p, 2)
    shp.sample[,,i]<-Rotation(shp.sample[,,i], rotate*0.0174532925199433)
  }
  
  scores.ext<-matrix(0, nrow=4, ncol=ncol(PC.amp))
  for(i in 1:2) scores.ext[((i*2)-1):(i*2),axes[i]]<-range(PC.amp[,axes[i]])
  
  shapes.ext<-array(0, c(nrow=p, ncol=2, 4))
  for(i in 1:nrow(scores.ext)) shapes.ext[,,i]<-arrayspecs(efourier_i2(t(t(scores.ext[i,]%*%t(v))+c), model.p), model.p, 2)
  
  dimnames(shapes.ext)[[3]]<-c(paste("PC", axes[1], "min", sep="_"), paste("PC", axes[1], "max", sep="_"),
                               paste("PC", axes[2], "min", sep="_"), paste("PC", axes[2], "max", sep="_"))
  
  for(i in 1:length(m[,1])){
    shp.sample[,1,i]<-shp.sample[,1,i]*aspX
    lines(shp.sample[,,i]*0.1*s+matrix(rep(PC.sample[i,axes],model.p),model.p,2,byrow=T)*adj, col=model.col, lwd=model.lwd)
  }
  
  if(points==T){points(X[,axes], pch=pch, col=col)}
  
  return(invisible(list(PCA=pca, ext.shapes.PC1=shapes.ext[,,1:2], ext.shapes.PC2=shapes.ext[,,3:4])))
  
  
  #Description
  # The mspace3 function plots the first two axes resulting from a principal component analysis
  # (regular, between groups and phylogenetic) and the theoretical shapes (EFA outlines) and 
  # representing shape variation in this space. Returns a list with the PCA in prcomp format 
  # (scores, eigenvectors, center), and the negative and positive extreme shapes from each PC 
  # (as x,y coordinates, not Fourier coefficients!).
  
  #Arguments:
  # ef= an OutCoe object resulting from efourier (from Momocs) and containing a matrix of Fourier
  #     coefficients
  # axes= a numeric vector indicating the PC axes to be represented
  # amp= numeric; amplifying factor for shape models
  # n= numeric; density of shape models in the plot
  # rotate= numeric; angle (in degrees) to rotate the shape models
  # s= numeric; scaling factor for shape models
  # adj= numeric; ajusting factor for shape models position
  # PCA= character; the function to be used for PCA. if "bgPCA", the groupPCA function from Morpho 
  #      is used; if "phyloPCA", the phyl.pca function from phytools is used instead
  # invertax= numeric; which axis (of those defined by axes) should be inverted?
  # aspX= numeric; a factor for compressing/depressing shape models along the X axis
  # pch= symbol of the scatter points
  # col= color of scatter points
  # model.col= color of the lines for background outlines models
  # model.lwd= width of the lines for background outlines models
  # model.p= number of points for background outlines models representation
  # points= logical; wheter to plot the PC scores
  # xlim, ylim, main= standard arguments for the generic plot function
  # ...= further arguments passed to groupPCA (groups) or phyl.pca (tree)
  
}


## mspace4 --------------------------------------------------------------------------

mspace4<-function(M, axes=c(1,2), template, amp=1, n=5, rotate=0, s=1, adj=0.9, PCA="regular", invertax=NULL, aspX=1,
                  pch=1, col="black", model.col="#708095", model.lwd=1, main="", points=TRUE, xlim=NULL, ylim=NULL,...){
  
  require(geomorph)
  require(spdep)
  
  p<-ncol(M)/2
  
  if(PCA=="regular") {
    pca<-prcomp(M)
  }
  
  if(PCA=="bgPCA") {
    require(Morpho)
    bgpca<-groupPCA(dataarray = M, ...) ; pca<-list()
    pca$rotation<-bgpca$groupPCs
    pca$x<-bgpca$Scores
    pca$center<-colMeans(M)
  }
  
  if(PCA=="phyloPCA") {
    require(phytools)
    ppca<-phyl.pca(Y = M, ...) ; pca<-list()
    pca$rotation<-ppca$Evec
    pca$x<-ppca$S
    pca$center<-colMeans(M)
  }
  
  X<-pca$x ; X[,axes[invertax]]<-X[,axes[invertax]]*-1
  v<-pca$rotation ; v[,axes[invertax]]<-v[,axes[invertax]]*-1
  c<-pca$center
  
  if(is.null(xlim)) xlim<-range(X[,axes[1]])
  if(is.null(ylim)) ylim<-range(X[,axes[2]])
  
  plot(X[,axes[1]], X[,axes[2]], xlim=xlim, ylim=ylim,  main=main,
       xlab=paste("PC", axes[1], sep=""), ylab=paste("PC", axes[2], sep=""), col="white")
  
  m<-matrix(0, nrow=n*n, ncol=length(v[1,]))
  m[,axes[1]]<-rep(seq(from=min(xlim), to=max(xlim), length.out=n), times=c(rep(n, n)))
  m[,axes[2]]<-rep(seq(from=min(ylim), to=max(ylim), length.out=n), n)
  
  PC.sample<-matrix(0, ncol=nrow(M), nrow=length(m[,1]))
  for(i in 1:length(m[,1])){
    PC.sample[i,axes[1]]<-m[i,axes[1]]
    PC.sample[i,axes[2]]<-m[i,axes[2]]}
  
  PC.amp<-PC.sample*amp
  
  if(ncol(PC.amp)>=nrow(v)){PC.amp<-PC.amp[,1:nrow(v)]}
  if(ncol(PC.amp)>=ncol(v)){PC.amp<-PC.amp[,1:ncol(v)]}  
  
  
  shp.sample<-array(0, c(nrow=p, ncol=2, length(m[,1])))
  for(i in 1:length(m[,1])){
    shp.sample[,,i]<-arrayspecs(t(t(PC.amp[i,]%*%t(v))+c), p, 2)
    shp.sample[,,i]<-Rotation(shp.sample[,,i], rotate*0.0174532925199433)
  }
  
  scores.ext<-matrix(0, nrow=4, ncol=ncol(PC.amp))
  for(i in 1:2) scores.ext[((i*2)-1):(i*2),axes[i]]<-range(PC.amp[,axes[i]])
  
  shapes.ext<-array(0, c(nrow=p, ncol=2, 4))
  for(i in 1:nrow(scores.ext)) shapes.ext[,,i]<-arrayspecs(t(t(scores.ext[i,]%*%t(v))+c), p, 2)
  
  dimnames(shapes.ext)[[3]]<-c(paste("PC", axes[1], "min", sep="_"), paste("PC", axes[1], "max", sep="_"),
                               paste("PC", axes[2], "min", sep="_"), paste("PC", axes[2], "max", sep="_"))
  
  
  centroid<-arrayspecs(t(t(rep(0, ncol(PC.amp))%*%t(v))+c), p, 2)[,,1]
  warpedcentroid<-rbind(centroid, tps2d(template[-(1:p),], template[1:p,], centroid))
  
  warped.sample<-array(0, c((nrow(warpedcentroid)-p), 2, dim(shp.sample)[3]))
  for(i in 1:dim(shp.sample)[3]){
    warped.sample[,,i]<-tps2d(warpedcentroid[-(1:p),], warpedcentroid[1:p,], shp.sample[,,i])}
  
  p2<-dim(warped.sample)[1]
  for(i in 1:length(m[,1])){
    warped.sample[,1,i]<-warped.sample[,1,i]*aspX
    lines(warped.sample[,,i]*0.1*s+matrix(rep(PC.sample[i,axes],p2),p2,2,byrow=T)*adj, col=model.col, lwd=model.lwd)
  }
  
  if(points==T){points(X[,axes], pch=pch, col=col)}
  
  return(invisible(list(PCA=pca, ext.shapes.PC1=shapes.ext[,,1:2], ext.shapes.PC2=shapes.ext[,,3:4])))
  
  
  #Description
  # The mspace4 function plots the first two axes resulting from a principal component analysis
  # (regular, between groups and phylogenetic) and the theoretical shapes (warped coordinates defining 
  # outlines) and representing shape variation in this space. Returns a list with the PCA in prcomp 
  # format (scores, eigenvectors, center), and the negative and positive extreme shapes from each PC .
  
  
  #Arguments:
  # M= a two-dimensional array (i.e., a n x (p x 2) matrix) with normalizaed configurations
  # axes= a numeric vector indicating the PC axes to be represented
  # amp= numeric; amplifying factor for shape models
  # n= numeric; density of shape models in the plot
  # rotate= numeric; angle (in degrees) to rotate the shape models
  # s= numeric; scaling factor for shape models
  # adj= numeric; ajusting factor for shape models position
  # PCA= character; the function to be used for PCA. if "bgPCA", the groupPCA function from Morpho 
  #      is used; if "phyloPCA", the phyl.pca function from phytools is used instead
  # invertax= numeric; which axis (of those defined by axes) should be inverted?
  # aspX= numeric; a factor for compressing/depressing shape models along the X axis
  # pch= symbol of the scatter points
  # col= color of scatter points
  # model.col= color of the lines for background outlines models
  # model.lwd= width of the lines for background outlines models
  # points= logical; wheter to plot the PC scores
  # xlim, ylim, main= standard arguments for the generic plot function
  # ...= further arguments passed to groupPCA (groups) or phyl.pca (tree)
  
}

