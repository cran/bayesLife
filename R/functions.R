


g.dl6<-function(x,l, p1, p2){
  d1<-x[1]; d2<-x[2]; d3<-x[3]; d4<-x[4]; k1<-x[5]; z<-x[6]

  m1<-sum(d1,.5*d2)
  m2<-sum(d1,d2,d3,.5*d4)
  k2<-z-k1
  k1/(1+exp(-log(p1^2)*(l-m1)/d2))+ k2/(1+exp(-log(p2^2)*(l-m2)/d4))
}

loess.lookup<-function(look){
   #-- spline function
     a<- 2.1828214 
     b<- -0.0335983
     c<- 0.0463050 
     d<- -0.0495180 
   
     x1<- 51.2058
     x2<-65.801 
   
     val<-a+b*look+c*(ifelse(look>x1,look-x1,0))+d*(ifelse(look>x2,look-x2,0))
     
     if(look>77.2) val<-a+b*77.2+c*(ifelse(77.2>x1,77.2-x1,0))+d*(ifelse(77.2>x2,77.2-x2,0))
      
   return(val)
}


dnorm.trunc<-function(x,mean,sd,low,high){
  out<-dnorm(x,mean=mean,sd=sd)/(pnorm(high,mean=mean,sd=sd)-pnorm(low,mean=mean,sd=sd))
  out[x<low]<-0
  out[x>high]<-0
  return(out)
}

rnorm.trunc<-function(mean,sd,low,high){
  temp<--999
  maxit <- 10
  i <- 1
  while((temp<low || temp>high) && i <= maxit) {
     temp<-rnorm(1,mean=mean,sd=sd)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- if(temp<low) low else high
  	warning(paste('Maximum iterations reached in rnorm.trunc(', 
  				mean, ',', sd, '). Value truncated to ', temp, '.', sep=''), immediate.=TRUE)
  }
  return(temp)
}

rgamma.ltrunc<-function(shape,rate,low){
  temp<--999
  maxit <- 10
  i <- 1
  while(temp<low && i <= maxit) {
     temp<-rgamma(1,shape=shape,rate=rate)
     i <- i+1
  }
  if (i > maxit) {
  	temp <- low
  	warning(paste('Maximum iterations reached in rgamma.ltrunc(', shape, ',', rate, ').', sep=''), immediate.=TRUE)
  }
  return(temp)
}


