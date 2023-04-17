library("statmod")
library("mvtnorm")
library("gtools")
#S is the number of stages for each experimantial arm
#L is the lower boundary for each arm
#U is the upper boundary for each arm
# t is the unknown which we are intergrating out
# r is the ratio martix
# r0 is the ratio for control
# n is the number of patients ratio matrix
# n0 is the number of patients on control 
#delta1 is the effect of the clinically relevant treatment
#delta0 is the uninteresting effect
#f is the clinically relevant treatment
#J is the point where treatment f is taken forward 
#sh say when treatment h is added
#sd is standard derivation

#boundary functions

po=function(rf)
{
  rfNA=rf[!is.na(rf)]
  LrfNA=length(rfNA)
  UfNA=rep(1,LrfNA)  
  LfNA=rep(-1,LrfNA)
  LfNA[LrfNA]=UfNA[LrfNA]
  
  NNA=length(rf)-LrfNA #number of NA needed
  Uf=c(UfNA,rep(NA,NNA))
  Lf=c(LfNA,rep(NA,NNA))
  
  return(list("Uf"=Uf,"Lf"=Lf))
}

po0=function(rf)
{
  rfNA=rf[!is.na(rf)]
  LrfNA=length(rfNA)
  UfNA=rep(1,LrfNA)  
  LfNA=rep(0,LrfNA)
  LfNA[LrfNA]=UfNA[LrfNA]
  
  NNA=length(rf)-LrfNA #number of NA needed
  Uf=c(UfNA,rep(NA,NNA))
  Lf=c(LfNA,rep(NA,NNA))
  
  return(list("Uf"=Uf,"Lf"=Lf))
}

obf=function(rf)
{
  rfNA=rf[!is.na(rf)]
  MrfNA=max(rfNA)
  LrfNA=length(rfNA)
  rfNA1=rfNA[1]
  UfNA=1/sqrt(rfNA/rfNA1)
  LfNA=-1/sqrt(rfNA/rfNA1)
  LfNA[LrfNA]=UfNA[LrfNA]
  
  NNA=length(rf)-LrfNA #number of NA needed
  Uf=c(UfNA,rep(NA,NNA))
  Lf=c(LfNA,rep(NA,NNA))
  
  return(list("Uf"=Uf,"Lf"=Lf))
}

obf0=function(rf)
{
  rfNA=rf[!is.na(rf)]
  MrfNA=max(rfNA)
  LrfNA=length(rfNA)
  rfNA1=rfNA[1]
  UfNA=1/sqrt(rfNA/rfNA1)
  LfNA=rep(0,LrfNA)
  LfNA[LrfNA]=UfNA[LrfNA]
  
  NNA=length(rf)-LrfNA #number of NA needed
  Uf=c(UfNA,rep(NA,NNA))
  Lf=c(LfNA,rep(NA,NNA))
  
  return(list("Uf"=Uf,"Lf"=Lf))
}

tri=function(rf) 
{
  rfNA=rf[!is.na(rf)]
  MrfNA=max(rfNA)
  LrfNA=length(rfNA)
  rfNA1=rfNA[1]
  rfNAa1=rfNA/rfNA1
  MrfNAa1=MrfNA/rfNA1
  UfNA=1*(1+rfNAa1/MrfNAa1)/sqrt(rfNAa1)
  LfNA=-1*(1-3*rfNAa1/MrfNAa1)/sqrt(rfNAa1)
  
  NNA=length(rf)-LrfNA #number of NA needed
  Uf=c(UfNA,rep(NA,NNA))
  Lf=c(LfNA,rep(NA,NNA))
  
  return(list("Uf"=Uf,"Lf"=Lf))  
}




#
gausspoints=function(N,a,b) 
{
  gq = gauss.quad(N)
  w = gq$weights*(b-a)/2
  n = a + (b-a)*(gq$nodes + 1)/2
  return(list("w"=w,"n"=n))
}



# Alpha code
editedbound=function(L,U,sh,S,r,r0,t,lsh)
{
  
  editL=matrix(nrow=lsh,ncol = S) #new lower boundaries
  editU=matrix(nrow=lsh,ncol = S) #new lower boundaries
  for(h in 1:lsh)
  {
    r0sh=0    #finds out if it is 0 or not as R0[0] gives error 
    if(sh[h]>0)  
    {
      r0sh=r0[(sh[h])]  
    }
    
    for(j in 1:S) 
    {
      partl=L[h,j]*sqrt(1+(r[h,j]/(r0[(sh[h]+j)]-r0sh))) #bit that using lower boundary
      partu=U[h,j]*sqrt(1+(r[h,j]/(r0[(sh[h]+j)]-r0sh))) #bit that using upper boundary
      partt1=sqrt(r[h,j])/(r0[(sh[h]+j)]-r0sh) 
      partt2=0 #the t part in the sum
      for(i in 1:j)
      {
        tpart=t[(sh[h]+i)]
        frpart=r0[(sh[h]+i)] #first r part 
        srpart=0 #second r part
        if((sh[h]+i-1)>0)
        {
          srpart=r0[(sh[h]+i-1)]  
        }
        partt2i=tpart*sqrt(frpart-srpart) #the t part for this i
        partt2=partt2+partt2i
      } 
      editL[h,j]=partl+partt1*partt2
      editU[h,j]=partu+partt1*partt2
    }
  }
  return(list(editL,editU))
}



ALPhaforgivent=function(L,U,sh,S,r,r0,t,acm,lsh)
{
  eb= editedbound(L,U,sh,S,r,r0,t,lsh)#edited boundaries
  leb= eb[[1]]
  ueb= eb[[2]] 
  A=length(sh)   #number of arms 
  prodpart=rep(NA,A) 
  for(h in 1:A)
  {
    cm=acm[[h]]  
    sumpart=rep(NA,S) #sum part
    sumpart[1]=pmvnorm(lower=-Inf, upper = leb[h,1],mean = 0, sigma = 1)[1]
    for(j in 2:S)
    {
      sumpart[j]=pmvnorm(lower=c(leb[h,1:(j-1)],-Inf) , upper = c(ueb[h,1:(j-1)],leb[h,j])  ,mean = rep(0,j), corr = cm[1:j,1:j])[1]
    }
    prodpart[h]=sum(sumpart)
  }
  ans=prod(prodpart) #1 minus FWER given t 
  return(ans)
}  



FWER=function(L,U,sh,S,r,r0,N,acm,lsh,lt)
{ 
  
  quad=gausspoints(N=N, -6, 6)
  perm=permutations(n=N, r=lt,repeats.allowed = TRUE)
  permv=as.vector(t(perm)) #perm vector
  tt=quad$n[permv]
  tw=quad$w[permv] #test weights 
  #testing stuff
  counter=0
  
  for(i in 1:(N^lt))
  {
    fv=ALPhaforgivent(L,U,sh,S,r,r0,t=tt[(lt*(i-1)+1):(lt*i)],acm,lsh)
    wei=prod(tw[(lt*(i-1)+1):(lt*i)]) #weight part
    normpart=prod(dnorm(tt[(lt*(i-1)+1):(lt*i)]))
    counter=counter+fv*wei*normpart
  }
  Fwercal=1-counter
  return(Fwercal)
}




#power code 
editedboundpower=function(L,U,sh,S,r,r0,n,n0,delta0,f,delta1,sd,t,lsh)
{
  
  editL=matrix(nrow=lsh,ncol = S) #new lower boundaries
  editU=matrix(nrow=lsh,ncol = S)
  delta=rep(delta0,lsh) #so have delta 0 for all apart from f 
  delta[f]=delta1
  for(h in 1:lsh)
  {
    r0sh=0    #finds out if it is 0 or not as R0[0] gives error 
    if(sh[h]>0)  
    {
      r0sh=r0[(sh[h])]  
    }
    
    for(j in 1:S) 
    {
      partl=L[h,j]*sqrt(1+(r[h,j]/(r0[(sh[h]+j)]-r0sh))) #bit that using lower boundary
      partu=U[h,j]*sqrt(1+(r[h,j]/(r0[(sh[h]+j)]-r0sh))) #bit that using upper boundary
      partt1=sqrt(r[h,j])/(r0[(sh[h]+j)]-r0sh)
      partt2=0 #the t part in the sum
      for(i in 1:j)
      {
        tpart=t[(sh[h]+i)]
        frpart=r0[(sh[h]+i)] #first r part 
        srpart=0 #second r part
        if((sh[h]+i-1)>0)
        {
          srpart=r0[(sh[h]+i-1)]  
        }
        partt2i=tpart*sqrt(frpart-srpart) #the t part for this i
        partt2=partt2+partt2i
      } 
      partdelta=delta[h]*sqrt(n[h,j])/sd
      
      editL[h,j]=partl+partt1*partt2-partdelta
      editU[h,j]=partu+partt1*partt2-partdelta
    }
  }
  return(list(editL,editU))
}

# edited bounds for treatment f
#eb is edited bounds found above
editedboundf=function(f,eb,v,J,r)
{
  leb= eb[[1]] 
  ueb= eb[[2]] 
  Lf=rep(NA,J-1) #lower for f 
  Uf=rep(NA,J-1) #upper for f 
  for(j in 1:J-1)  
  {
    Lf[j]=sqrt(r[f,J]/(r[f,J]-r[f,j]))*(leb[f,j]-v*sqrt(r[f,j]/r[f,J])) 
    Uf[j]=sqrt(r[f,J]/(r[f,J]-r[f,j]))*(ueb[f,j]-v*sqrt(r[f,j]/r[f,J])) 
  }
  return(list(Lf,Uf))
}



corralationmatix=function(r) #given h
{
  A=dim(r)[1]
  S=dim(r)[2]
  AllCormat=vector(mode = "list", length = A)
  for(h in 1:A)
  {  
    Cormat=diag(1,nrow = S)
    for(j1 in 1:S) #j
    {
      for(j2 in 1:S) #j*
      {
        Cormat[j1,j2]=min(sqrt(r[h,j1]/r[h,j2]),sqrt(r[h,j2]/r[h,j1]))  
      }
    }
    AllCormat[[h]]=Cormat  
  }
  return(AllCormat)  
}

corralationmatixf=function(r,f,J) 
{
  Cormat=diag(1,nrow = (J-1))
  for(j1 in 1:(J-1)) #j
  {
    for(j2 in 1:(J-1)) #j*
    {
      inside1=(r[f,j1]*(r[f,J]-r[f,j2]))/(r[f,j2]*(r[f,J]-r[f,j1]))  #inside the sqrt part
      inside2=(r[f,j2]*(r[f,J]-r[f,j1]))/(r[f,j1]*(r[f,J]-r[f,j2]))  #inside the sqrt part
      Cormat[j1,j2]=min(sqrt(inside1),sqrt(inside2))  
    }
  }
  return(Cormat)  
}


treatmentfatJ=function(r,f,J,eb,v,fcm,lsh) #find probability of treatment f gets to stage J 
{
  startpart=1
  if(J>1)
  {
    cm=fcm[[((J-2)*lsh+f)]]
    ebf=editedboundf(f=f,eb=eb,r=r,v=v,J=J) #edited boundaries for f
    lebf= ebf[[1]] #edited lower boundaries for f 
    uebf= ebf[[2]] # edited upper boundaries for f
    if(J == 2)
    {
      startpart=pmvnorm(lower=lebf, upper = uebf,mean = 0, sigma = 1)[1]
    }
    else
    {
      startpart=pmvnorm(lower=lebf, upper = uebf,mean = rep(0,(J-1)) , corr = cm)[1]
    }
    
  }
  return(startpart)  
}



Betaforgivent=function(L,U,sh,S,r,r0,n,n0,delta0,delta1,f,sd,t,lsh,v,J,acm,fcm)
{ 
  eb= editedboundpower(L,U,sh,S,r,r0,n,n0,delta0,f,delta1,sd,t,lsh)#edited boundaries
  leb= eb[[1]]
  ueb= eb[[2]] #upper edited boundaries 
  A=length(sh)   #number of arms 
  prodpart=rep(NA,A) #product part
  for(h in 1:A)
  {
    ip=sh[f]+J-sh[h]
    if(h==f | ip <= 0)
    {
      sumpart=1  
    }
    else
    {
      #S is also the same as Y
      mYsf=min(ip,S) #min of Y,sf+J-s(h)
      cm=acm[[h]]  #cor matrix for treatment h 
      sumpart=rep(NA,mYsf) #sum part
      if(mYsf==1)
      {
        odotUB=sqrt(n[h,ip]/n[f,J])*v+(sqrt(n[h,ip])*(delta1-delta0))/sd   
        dotUB=max(ueb[h,1],odotUB) 
        sumpart[1]=pmvnorm(lower=-Inf, upper = dotUB,mean = 0, sigma = 1)[1]
      }
      else
      {
        sumpart[1]=pmvnorm(lower=-Inf, upper = leb[h,1],mean = 0, sigma = 1)[1]
        for(j in 2:mYsf)
        {
          if(ip==j)
          {
            odotUB=sqrt(n[h,ip]/n[f,J])*v+(sqrt(n[h,ip])*(delta1-delta0))/sd   
            dotUB=max(ueb[h,j],odotUB)  
            sumpart[j]=pmvnorm(lower=c(leb[h,1:(j-1)],-Inf) , upper = c(ueb[h,1:(j-1)],dotUB)  ,mean = rep(0,j), corr = cm[1:j,1:j])[1]
          }
          else
          {  
            sumpart[j]=pmvnorm(lower=c(leb[h,1:(j-1)],-Inf) , upper = c(ueb[h,1:(j-1)],leb[h,j])  ,mean = rep(0,j), corr = cm[1:j,1:j])[1]
          }        
        }
      }
    }
    prodpart[h]=sum(sumpart)
    
    
  }
  OT=prod(prodpart)   

  BJ=treatmentfatJ(r=r,f=f,J=J,eb = eb,v = v,fcm=fcm,lsh=lsh) #before treatment is taken forward at J
  
  ans=BJ*OT
  return(ans)
  
  
}  


vstuff=function(n,n0,delta0 ,delta1,f,sd,lsh,J,L,U,sh,S,r,r0,ts)
{
  

  lhsI1=(delta1/sd) #Left hand side of indicator function part 1
  lhsI2v=rep(NA,J) #left hand side of vector indicator function part 2
  
  if(sh[f]==0)
  {
    n0sf=0  # n0 at the point s(f)
  }  
  else
  {
    n0sf=n0[sh[f]]
  }
  
  for(i in 1:J)
  {
    if((sh[f]+i-1)==0)
    {
      n0sfi=0  
    }
    else
    {
      n0sfi=n0[(sh[f]+i-1)]  
    }
    
    lhsI2v[i]=ts[sh[f]+i]*sqrt(n0[sh[f]+i]-n0sfi)  
  }
  
  vmin=sqrt(n[f,J])*(U[f,J]*sqrt((n[f,J]^(-1))+(n0[sh[f]+J]-n0sf)^(-1))-lhsI1+sum(lhsI2v)/(n0[sh[f]+J]-n0sf))  #Left hand side of indicator function
  return(vmin)
}


###### Find power for J

PowerforJ=function(L,U,sh,S,r,r0,n,n0,delta0,delta1,f,sd,J,N,acm,lsh,fcm)
{ 
  
  lt=sh[f]+J 
 
  
  quad=gausspoints(N=N, -6, 6)
  perm=permutations(n=N, r=lt+1,repeats.allowed = TRUE)
  permv=as.vector(t(perm)) #perm vector
  allpoints=quad$n[permv]
  tt=allpoints[c(rep(TRUE,lt),FALSE)]
  
  wquad=quad$w[permv] #test weights
  ttw=wquad[c(rep(TRUE,lt),FALSE)]
  
  
  counter=0
  for(i in 1:(N^(lt+1)))
  {
    
    if(i %% N ==1)
    {
      minforv=vstuff(n = n,n0 = n0,delta0 = delta0,delta1 = delta1,f = f,sd = sd,lsh = lsh,J = J,L=L,U=U,sh=sh,S=S,r=r,r0=r0,ts=tt[(lt*(i-1)+1):(lt*i)])
      quadforv=gausspoints(N=N, minforv, 6)
      cim1=i-1 #current i minus 1
      ctv=quadforv$n #current test v
      ctvw=quadforv$w  #current test v weight
    }
    
    fv=Betaforgivent(n = n,n0 = n0,delta0 = delta0,delta1 = delta1,f = f,sd = sd,lsh = lsh,J = J,L=L,U=U,sh=sh,S=S,r=r,r0=r0,v=ctv[i-cim1],t=tt[(lt*(i-1)+1):(lt*i)],acm=acm,fcm=fcm)#function given value # in this we have egnored the dphi(t) as this cancels with the normal distribution
   
    wei=prod(ttw[(lt*(i-1)+1):(lt*i)])*ctvw[i-cim1] #weight part
    normpart=prod(dnorm(tt[(lt*(i-1)+1):(lt*i)]))*dnorm(ctv[i-cim1])
    counter=counter+wei*normpart*fv
  }
  powerJ=counter
  return(powerJ)
}



LFCPowerf=function(L,U,sh,S,r,r0,n,n0,delta0,delta1,sd,N,acm,lsh,fcm,f)
{
  
  Jpower=rep(NA,S)
  for(j in 1:S)
  {
    Jpower[j]=PowerforJ(f = f,delta1 = delta1,r=r,r0 = r0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd,J = j,N=N, acm=acm, lsh=lsh,fcm=fcm)
  }
  Fpower=sum(Jpower)
  print(Fpower)
  return(Fpower)
}





corralationmatixfet=function(r,r0,sh) #given h # correlation matrix for each treatment
{
  A=dim(r)[1]
  S=dim(r)[2]
  AllCormat=vector(mode = "list", length = A)
  for(h in 1:A)
  {  
    Cormat=diag(1,nrow = (S))
    for(j1 in 1:(S)) #j
    {
      for(j2 in 1:(S)) #j*
      {
        if(sh[h]==0)
        {
          r0sh=0  
        }
        else
        {
          r0sh=r0[(sh[h])]
        }
        part1=sqrt( (r[h,j1])^(-1) + (r0[(sh[h]+j1)]-r0sh)^(-1))
        part2=sqrt( (r[h,j2])^(-1) + (r0[(sh[h]+j2)]-r0sh)^(-1))
        part3a=(1/r[h,j1])+(1/(r0[(sh[h]+j1)]-r0sh))
        part3b=(1/r[h,j2])+(1/(r0[(sh[h]+j2)]-r0sh))
        Cormat[j1,j2]=min(part3a,part3b)/(part1*part2) 
      }
    }
    AllCormat[[h]]=Cormat  
  }
  print(AllCormat)
  return(AllCormat)
}


currentagivenhUL=function(acm,sh,L,U,f,S,seed.number)
{
  set.seed(seed.number)
  cm=acm[[f]]  #cor matrix for treatment h 
  sumpart=rep(NA,S) #sum part
  sumpart[1]=pmvnorm(lower=-Inf, upper = L[1],mean = 0, sigma = 1)[1]
  if((S)>1)
  {
    for(j in 2:(S))
    {
      sumpart[j]=pmvnorm(lower=c(L[1:(j-1)],-Inf) , upper = c(U[1:(j-1)],L[j])  ,mean = rep(0,j), corr = cm[1:j,1:j])[1]
    }
  }
  onemw=sum(sumpart) #one minus w 
  return(onemw)
}


giveaapproximation=function(alpha,r,r0,sh,S,tolerance,a1,BshU,BshL,seed.number,amax)
{
  acm=corralationmatixfet(r=r,r0=r0,sh=sh)
  lsh=length(sh)
  U1=a1*BshU[1,]    
  L1=a1*BshL[1,]
  onemw=currentagivenhUL(U=U1,L=L1,acm=acm,sh=sh,f=1,S=S,seed.number=seed.number)
  #test for the same boundary shape
  a2=a1
  U2=a2*BshU[2,]    
  L2=a2*BshL[2,]
  testa2=currentagivenhUL(U=U2,L=L2,acm=acm,sh=sh,f=2,S=S,seed.number=seed.number)
  
  if(abs(testa2-onemw)<tolerance)
  {
    a2=a1
    print("sameshape")
  }
  else
  {
  
    avector=rep(NA,lsh)
    ua2=amax
    la2=0
    a2=(ua2+la2)/2
    U2=a2*BshU[2,]    
    L2=a2*BshL[2,]
    chpv=currentagivenhUL(U=U2,L=L2,acm=acm,sh=sh,f=2,S=S,seed.number=seed.number) #current halfway point value 
  
    while((ua2-la2)>tolerance/1000)
    {
     
      if(onemw < chpv)
      {
        la2=la2
        ua2=a2
        a2=(ua2+la2)/2
        U2=a2*BshU[2,]    
        L2=a2*BshL[2,]
        chpv=currentagivenhUL(U=U2,L=L2,acm=acm,sh=sh,f=2,S=S,seed.number=seed.number) #current halfway point value 
        
      }
      
      else if(onemw > chpv)
      {
        la2=a2
        ua2=ua2
        a2=(ua2+la2)/2
        U2=a2*BshU[2,]    
        L2=a2*BshL[2,]
        chpv=currentagivenhUL(U=U2,L=L2,acm=acm,sh=sh,f=2,S=S,seed.number=seed.number) #current halfway point value 
      }
      
    }
  }
  print("a2")
  print(a2)
  return(a2)  
}



boundfinder=function(f,sh,S,r,r0,N,alpha,tolerance,power,sd,aratio,amax,BshL,BshU)
{
  set.seed(1)
  lsh=length(sh) #how many treatments there are
  acm=corralationmatix(r) #all correlation matrices
  A=length(sh)
  ua=rep(amax,A)  #upper value for a 
  lt=max(sh)+S
  la=rep(0,A)    #lower value for a          
  a=aratio*(ua[1]+la[1])/(2*aratio[1])  
  L=matrix(NA,nrow = lsh,ncol = S , byrow = T) 
  U=matrix(NA,nrow = lsh,ncol = S , byrow = T)
  wmaxar=which(aratio==max(aratio))[1] 
  
  
  for(h in 1:lsh)
  {
    L[h,]=a[h]*BshL[h,]
    U[h,]=a[h]*BshU[h,]
  }
  print("a")
  print(a)
  print("BshU")
  print(BshU)
  
  halfpointvalue=FWER(L,U,sh,S,r,r0,N=N,acm=acm,lsh=lsh,lt=lt)
  
  print("halfpointvalue")
  print(halfpointvalue)
  
  counter=1
  while((ua[wmaxar]-la[wmaxar])>tolerance)
  {  
    counter=counter+1
    if(alpha > halfpointvalue)
    {
      la=la
      ua=a
      a=aratio*(ua[1]+la[1])/(2*aratio[1])  
      L=matrix(NA,nrow = lsh,ncol = S , byrow = T) 
      U=matrix(NA,nrow = lsh,ncol = S , byrow = T)
      for(h in 1:lsh)
      {
        L[h,]=a[h]*BshL[h,]
        U[h,]=a[h]*BshU[h,]
      }  
      halfpointvalue=FWER(L,U,sh,S,r,r0,N=N,acm=acm,lsh=lsh,lt=lt)
    }
    
    else if(alpha < halfpointvalue)
    {
      la=a
      ua=ua
      a=aratio*(ua[1]+la[1])/(2*aratio[1]) 
      L=matrix(NA,nrow = lsh,ncol = S , byrow = T) 
      U=matrix(NA,nrow = lsh,ncol = S , byrow = T)
      for(h in 1:lsh)
      {
        L[h,]=a[h]*BshL[h,]
        U[h,]=a[h]*BshU[h,]
      }  
      halfpointvalue=FWER(L,U,sh,S,r,r0,N=N,acm=acm,lsh=lsh,lt=lt)
    }
    print("halfpointvalue")
    print(halfpointvalue)
  }  
  
  
  L=matrix(NA,nrow = lsh,ncol = S , byrow = T) 
  U=matrix(NA,nrow = lsh,ncol = S , byrow = T)
  for(h in 1:lsh)
  {
    L[h,]=a[h]*BshL[h,]
    U[h,]=a[h]*BshU[h,]
  }  

  
  delta=rep(0,lsh)
  
  return(list(a,U,L))
}





exactafinder=function(sh,S,r,r0,N,alpha,tolerance,sd,amax)
{
  lsh=length(sh) #how many treatments there are
  caratio=rep(1,lsh) # current aratio
  oaratio=rep(Inf,lsh)
  
  BshU=matrix(NA,nrow=nrow(r),ncol = ncol(r)) #boundary shape
  BshU[1,]=tri(rf=r[1,])$Uf #change this code here for different shapes!
  BshU[2,]=tri(rf=r[2,])$Uf #change this code here for different shapes!
  
  
  BshL=matrix(NA,nrow=nrow(r),ncol = ncol(r)) #boundary shape
  BshL[1,]=tri(rf=r[1,])$Lf #change this code here for different shapes!
  BshL[2,]=tri(rf=r[2,])$Lf #change this code here for different shapes!
  
  
  boundpart=boundfinder(sh=sh,S=S,r=r,r0=r0,N=N,alpha=alpha,tolerance=tolerance,aratio = caratio,f=f,sd=sd,amax=amax,BshL=BshL,BshU=BshU)
  
  ca=boundpart[[1]]
  print("ca")
  print(ca)
  oa=rep(Inf,length(ca))
  counter=1
  while(sum(abs(oa-ca)>tolerance)>0)
  {
    print(counter)
    counter=counter+1
    oa=ca
    a2=giveaapproximation(alpha=alpha,r=r,r0=r0,sh=sh,S=S,tolerance=tolerance,a1=ca[1],BshU=BshU,BshL=BshL,seed.number=1,amax=amax)
    caratio=c(1,a2/ca[1])
    boundpart=boundfinder(sh=sh,S=S,r=r,r0=r0,N=N,alpha=alpha,tolerance=tolerance/100,aratio = caratio,f=f,sd=sd,amax=amax,BshL=BshL,BshU=BshU)
    ca=boundpart[[1]]
  }
  
  print("U")
  print(boundpart[[2]])
  
  print("caratio")
  print(caratio)
  return(list("U"=boundpart[[2]],"L"=boundpart[[3]],"aratio"=caratio))
}


findboundMamsandn=function(sh,S,N,alpha,tolerance,power,delta0,delta1,sd,nmax,amax)
{
  set.seed(1)
  dubhap=0
  lsh=length(sh) #how many treatments there are
  n1=1
  n2=1
  r=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
  r0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
  
  acm=corralationmatix(r) #all correlation matrices
  alist=c() #to stop it going into an invite loop if by changing the desegn it keeps jumping (it will stick with the worse U) 
  
  fcm=vector(mode = "list",length = lsh*(S-1))
  for(f in 1:lsh)
  {
    for(j in 2:S)
    {
      fcm[[((j-2)*lsh+f)]]=corralationmatixf(r=r,f=f,J=j)
    }
  }
  
  startbound=exactafinder(sh=sh,S=S,r=r,r0=r0,N=N,alpha=alpha,tolerance=tolerance,sd=sd,amax=amax)
  print("startbound")
  print(startbound)
  L=startbound$L
  U=startbound$U
  alist=c(alist,U[2])
    
  nl=0 
  nu=nmax
  i1=(nu+nl)/2
  n1=i1
  n2=i1
  n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
  n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
  cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=1)
  print(cp)
  print(i1)
  
  
  while((nu-nl)>(tolerance))
  {  
    if(cp<power)
    {
      nl=i1 
      nu=nu
      i1=(nu+nl)/2
      n1=i1
      n2=i1
      n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
      n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
      cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=1)
      print(cp)
      print(i1)
    }
    
    if(cp>power)
    {
      nl=nl 
      nu=i1
      i1=(nu+nl)/2
      n1=i1
      n2=i1
      n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
      n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
      cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=1)
      print(cp)
      print(i1)
    }
    
  } 
  
  nl=1 
  nu=nmax
  i2=(nu+nl)/2
  n1=n1
  n2=i2
  n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
  n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
  cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=2)
  print(cp)
  print(i2)

  while((nu-nl)>(tolerance))
  {  
    if(cp<power)
    {
      nl=i2 
      nu=nu
      i2=(nu+nl)/2
      n1=n1
      n2=i2
      n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
      n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
      cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=2)
      print(cp)
      print(i2)
    }
    
    if(cp>power)
    {
      nl=nl 
      nu=i2
      i2=(nu+nl)/2
      n1=n1
      n2=i2
      n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
      n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
      cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=2)
      print(cp)
      print(i2)
    }
    
  } 
  
  n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
  n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
  
  Ci1=n1
  Ci2=n2
  Oi1=Inf
  Oi2=Inf
  
  while(round(Oi1) != round(Ci1) | round(Oi2) != round(Ci2))
  {
    Oi1=Ci1
    print("Oi1")
    print(Oi1)
    Oi2=Ci2
    print("Oi2")
    print(Oi2)
    r=n
    r0=n0
    
    startbound=exactafinder(sh=sh,S=S,r=r,r0=r0,N=N,alpha=alpha,tolerance=tolerance,sd=sd,amax=amax)
    L=startbound$L
    U=startbound$U  
    
    
    
    nl=0 
    nu=nmax
    i1=(nu+nl)/2
    n1=i1
    n2=n2
    n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
    n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
    cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=1)
    print(cp)
    print(i1)
    
    
    while((nu-nl)>(tolerance))
    {  
      if(cp<power)
      {
        nl=i1 
        nu=nu
        i1=(nu+nl)/2
        n1=i1
        n2=n2
        n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
        n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
        cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=1)
        print(cp)
        print(i1)
      }
      
      if(cp>power)
      {
        nl=nl 
        nu=i1
        i1=(nu+nl)/2
        n1=i1
        n2=n2
        n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T) #need to change this code here for different number of stages or arms
        n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
        cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=1)
        print(cp)
        print(i1)
      }
      
    } 
    
    print("n1 now finished first part")
    print(n1)
    
    n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
    n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
    
    print("start of new run for t3")
    nl=1 
    nu=nmax
    i2=(nu+nl)/2
    n1=n1
    n2=i2
    n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
    n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
    cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=2)
    print(cp)
    print(i2)
    
    while((nu-nl)>(tolerance))
    {  
      if(cp<power)
      {
        nl=i2 
        nu=nu
        i2=(nu+nl)/2
        n1=n1
        n2=i2
        n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
        n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)#need to change this code here for different number of stages or arms
        cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=2)
        print(cp)
        print(i2)
      }
      
      if(cp>power)
      {
        nl=nl 
        nu=i2
        i2=(nu+nl)/2
        n1=n1
        n2=i2
        n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)#need to change this code here for different number of stages or arms
        n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2) #need to change this code here for different number of stages or arms
        cp=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=2)
        print(cp)
        print(i2)
      }
      
    } 
    
    n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T) #need to change this code here for different number of stages or arms
    n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2) #need to change this code here for different number of stages or arms
    
    Ci1=n1
    Ci2=n2
  }
  
  print("found expact n1 and n2")
  print(c(n1,n2))
  
  n1=ceiling(n1)
  n2=ceiling(n2)
  n=matrix(c(n1,2*n1,n2,2*n2),nrow = lsh,ncol = S, byrow = T)
  n0=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)
  
  r=n
  r0=n0
  
  startbound=exactafinder(sh=sh,S=S,r=r,r0=r0,N=N,alpha=alpha,tolerance=tolerance,sd=sd,amax=amax)
  L=startbound$L
  U=startbound$U  

  
  t1power=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=1)
  print("t1power")
  print(t1power)
  t3power=LFCPowerf(delta1 = delta1,r=n,r0 = n0,n=n,n0=n0,sh=sh,S=S,U=U,L=L,delta0 = delta0,sd = sd, N=N ,acm=acm,lsh=lsh,fcm=fcm,f=2)
  
  print("t3power")
  print(t3power)
  
  return(list(n1,n2,U,L))  
}




gcEBS=function(sh,S,N,alpha,tolerance,power,delta0,delta1,sd,nmax,amax)
{
boundsandn=findboundMamsandn(sh=sh,S=S,N=N,alpha=alpha,tolerance=tolerance,power=power,delta0=delta0,delta1=delta1,sd=sd,nmax=nmax,amax=amax)
return(list(boundsandn))
}

ptm=proc.time()
set.seed(1)
each2stagestritri=gcEBS(sh=c(0,1),S=2,N=20,alpha=0.025,tolerance=0.0001,power=0.8,delta0=-log(0.99),delta1=-log(0.69),sd=1,nmax=400,amax=40)
each2stagestritri
proc.time() - ptm

save(each2stagestritri,file = "each2stagestritri.RData")
