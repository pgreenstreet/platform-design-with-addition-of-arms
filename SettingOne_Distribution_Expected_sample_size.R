library("mvtnorm") #pvnorm
library("statmod")
library("mvtnorm")
library("gtools")
#S is the number of stages for each experiential arm
#L is the lower boundary for each arm
#U is the upper boundary for each arm
# t is the unknown which we are integrating out
# r is the ratio matrix
# r0 is the ratio for control
# n is the number of patients ratio matrix
# n0 is the number of patients on control 
#delta1 is the effect of the clinically relevant treatment
#delta0 is the uninteresting effect
#f is the clinically relevant treatment
#J is the point where treatment f is taken forward 
#sh say when treatment h is added
#sd is standard derivation

gausspoints<-function(N,a,b) 
{
  gq = gauss.quad(N)
  w = gq$weights*(b-a)/2
  n = a + (b-a)*(gq$nodes + 1)/2
  return(list("w"=w,"n"=n))
}

corralationmatixE=function(r) #given h
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

editedboundE=function(L,U,sh,S,r,r0,delta,n,n0,t,lsh,sd)
{
  
  editL=matrix(nrow=lsh,ncol = S) #new lower boundaries
  editU=matrix(nrow=lsh,ncol = S) #new lower boundaries
  for(h in 1:lsh)
  {
    r0sh=0    
    if(sh[h]>0)  
    {
      r0sh=r0[(sh[h])]  
    }
    
    for(j in 1:S) 
    {
      partl=L[h,j]*sqrt(1+(r[h,j]/(r0[(sh[h]+j)]-r0sh)))-delta[h]*sqrt(n[h,j])/sd #bit that using lower boundary
      partu=U[h,j]*sqrt(1+(r[h,j]/(r0[(sh[h]+j)]-r0sh)))-delta[h]*sqrt(n[h,j])/sd #bit that using upper boundary
      partt1=sqrt(r[h,j])/(r0[(sh[h]+j)]-r0sh) # the t part at the begining
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





#inside the integral
in_int=function(L,U,sh,S,r,r0,delta,n,n0,t,acm,lsh,sd,k,q)
{
  eb= editedboundE(L=L,U=U,sh=sh,S=S,r=r,r0=r0,delta=delta,n=n,n0=n0,t=t,lsh=lsh,sd=sd)#edited boundaries
  leb= eb[[1]]
  ueb= eb[[2]] #upper edited boundaries 
  A=length(sh)   #number of arms 
  prodpart=rep(NA,A) #product part
  for(h in 1:A)
  {
    cm=acm[[h]]  #cor matrix for treatment h 
    
    vkh=k[h] #value of k[h]
    if(q[h]==1)
    {
      if(vkh==1)
      {
        endpart=pmvnorm(lower=ueb[h,1], upper = Inf,mean = 0, sigma = 1)[1]
      }
      
      if(vkh>1)
      {
        endpart=pmvnorm(lower=c(leb[h,1:(vkh-1)],ueb[h,vkh]) , upper = c(ueb[h,1:(vkh-1)],Inf)  ,mean = rep(0,vkh), corr = cm[1:vkh,1:vkh])
      }
      
    }
    
    if(q[h]==Inf)
    {
      if(vkh==1)
      {
        endpart=pmvnorm(lower= -Inf, upper = leb[h,1],mean = 0, sigma = 1)[1]
      }
      
      if(vkh>1)
      {
        endpart=pmvnorm(lower=c(leb[h,1:(vkh-1)],-Inf) , upper = c(ueb[h,1:(vkh-1)],leb[h,vkh])  ,mean = rep(0,vkh), corr = cm[1:vkh,1:vkh])
      }
      
    }
    
    prodpart[h]=endpart
  }
  ans=prod(prodpart) #1 minus FWER given t 
  return(ans)
}  

ntpart=function(n,n0,sh,lsh,k,q)
{
  controln=NA
  activen=rep(NA,lsh)
  #active treatments 
  for(f in 1:lsh)
  {
    minstarmsf=min(k[f]+sh[f],(k+sh)*q)-sh[f]
    if(minstarmsf>0)
    {
      activen[f]=n[f,minstarmsf]
    }
    else
    {
      activen[f]=0  
    }
    
  }
  
  controlminstar=min(max((k+sh)),(k+sh)*q)
  controln=n0[controlminstar]
  
  totaln=controln+sum(activen)
  return(c(controln,activen,totaln))
}

#probability expected sample size distribution
PESSD=function(L,U,sh,S,r,r0,N,delta,n,n0,acm,lsh,lt,sd,k,q)
{ 
  shstar=max(sh)
  quad=gausspoints(N=N, -6, 6)
  perm=permutations(n=N, r=lt,repeats.allowed = TRUE)
  permv=as.vector(t(perm)) #perm vector
  tt=quad$n[permv]
  tw=quad$w[permv] #test weights
  
  probbit=rep(NA, (N^lt))
  
  for(i in 1:(N^lt))
  {
      wei=prod(tw[(lt*(i-1)+1):(lt*i)]) #weight part
      normpart=prod(dnorm(tt[(lt*(i-1)+1):(lt*i)]))
      probbit[i]=normpart*wei*in_int(L=L,U=U,sh=sh,S=S,r=r,r0=r0,delta=delta,n=n,n0=n0,t=tt[(lt*(i-1)+1):(lt*i)],acm=acm,lsh=lsh,sd=sd,k = k,q=q)#function given value # in this we have egnored the dphi(t) as this cancels with the normal distribution
  }
  
  ans=sum(probbit)
  return(ans)

}

Expectedsamplesize=function(sh,S,r,r0,N,gn,U,L,delta,sd)
{
  #set.seed(1)
  i=gn #will tell you how much to times r by
  n=r*i
  n0=r0*i
  lsh=length(sh) #how many treatments there are
  acm=corralationmatixE(r) #all correlation matrices
  lt=max(sh)+S #the number of t we will integrate over
  A=length(sh)
  perplist=permutations(n = S*2,r=lsh,v = c(-S:-1,1:S),repeats.allowed = TRUE)
  d1pl=dim(perplist)[1] 
  resultmatrix=matrix(data = NA,nrow = d1pl, ncol=(lsh+3))
  for(i in 1:d1pl)
  {
    print(i)
    cpl=perplist[i,] #current perplist
    k=abs(cpl)
    qm11= sign(cpl) 
    q=replace(qm11,qm11==-1,Inf)
    cntpart=ntpart(n=n,n0=n0,sh=sh,lsh=lsh,k=k,q=q)
    cPESSD=PESSD(L=L,U=U,sh=sh,S=S,r=r,r0=r0,N=N,delta=delta,n=n,n0=n0,acm=acm,lsh=lsh,lt=lt,sd=sd,k=k,q=q)
    resultmatrix[i,]=c(cntpart,cPESSD)
  }
  return(resultmatrix)
}


load("each2stagestritri.RData")
boundsandn=each2stagestritri[[1]]
n1=boundsandn[[1]]
n2=boundsandn[[2]]
nt=matrix(c(n1,2*n1,n2,2*n2),nrow = 2,ncol = 2, byrow = T)
n0t=c(n1,n1+max(n1,n2),n1+max(n1,n2)+n2)
Ut=boundsandn[[3]]
Lt=boundsandn[[4]]
gnt=1

ptm <- proc.time()
set.seed(1)
ansNull=Expectedsamplesize(gn=gnt,sh=c(0,1),S=2,r=nt,r0=n0t,N=20,U=Ut,L=Lt,delta=c(0,0),sd=1)
sum(ansNull[,4]*ansNull[,5]) #expected sample size
set.seed(1)
ansT1=Expectedsamplesize(gn=gnt,sh=c(0,1),S=2,r=nt,r0=n0t,N=20,U=Ut,L=Lt,delta=c(-log(0.69),-log(0.99)),sd=1)
sum(ansT1[,4]*ansT1[,5]) #expected sample size
set.seed(1)
ansT2=Expectedsamplesize(gn=gnt,sh=c(0,1),S=2,r=nt,r0=n0t,N=20,U=Ut,L=Lt,delta=c(-log(0.99),-log(0.69)),sd=1)
sum(ansT2[,4]*ansT2[,5]) #expected sample size
proc.time() - ptm



sampledis2tritri=list(ansNull,ansT1,ansT2)


#the null sample size distribution
ans=sampledis2tritri[[1]]

sum(ans[,5])
dataans=data.frame(ans)
names(dataans)  

#total sample size
orderans=dataans[order() ,]

names(dataans)[1] <- "CT"
names(dataans)[2] <- "T1"
names(dataans)[3] <- "T2"
names(dataans)[4] <- "TN"
names(dataans)[5] <- "PN"
orderans=dataans[order(dataans$TN) ,]

TotalN_dataTN=aggregate(. ~ TN, data=dataans, FUN=sum)
TotalN_dataTN$TN
TotalN_dataTN$PN
cumsum(TotalN_dataTN$PN)
plot(TotalN_dataTN$TN,TotalN_dataTN$PN,type = "s")
plot(TotalN_dataTN$TN,cumsum(TotalN_dataTN$PN),xlim=c(90,550),xlab="",ylab="",ylim=c(0,1),type = "s",lwd = 2)
points(TotalN_dataTN$TN,cumsum(TotalN_dataTN$PN),lwd=2)

#control treatment sample size
TotalN_dataCT=aggregate(. ~ CT, data=dataans, FUN=sum)
plot(TotalN_dataCT$CT,TotalN_dataCT$PN,type = "s")
plot(TotalN_dataCT$CT,cumsum(TotalN_dataCT$PN),xlim=c(0,250),ylim=c(0,1),xlab="",ylab="",type = "s",col="Blue",lwd = 2)
points(TotalN_dataCT$CT,cumsum(TotalN_dataCT$PN),col="Blue",lwd=2)

# active treatment 2 sample size
TotalN_dataT2=aggregate(. ~ T2, data=dataans, FUN=sum)
points(TotalN_dataT2$T2,cumsum(TotalN_dataT2$PN),xlab="",ylab="",type = "s",col="Dark Green",lwd = 2)
points(TotalN_dataT2$T2,cumsum(TotalN_dataT2$PN),col="Dark Green",lwd=2)

# active treatment 1 sample size
TotalN_dataT1=aggregate(. ~ T1, data=dataans, FUN=sum)
points(TotalN_dataT1$T1,cumsum(TotalN_dataT1$PN),xlab="",ylab="",type = "s",col="purple",lwd = 2)
points(TotalN_dataT1$T1,cumsum(TotalN_dataT1$PN),col="purple",lwd=2)


points(TotalN_dataCT$CT,cumsum(TotalN_dataCT$PN),xlab="",ylab="",type = "s",col="Blue",lwd = 2)
points(TotalN_dataCT$CT,cumsum(TotalN_dataCT$PN),col="Blue",lwd=2)
