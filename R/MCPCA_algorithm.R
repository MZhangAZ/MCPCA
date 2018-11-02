#################################
# run MCPCA
#################################

n= dim(G)[1]
p=dim(G)[2]
q= 10 # KyFan norm number
Phi_zero= G

eig_zero= prcomp(Phi_zero, center = F, scale. = F)

v_zero= eig_zero$rotation[,1:q]

rho_zero=sum((eig_zero$sdev^2)[1:q])

rm(eig_zero)
######################################################

P_old= Phi_zero
v_old= v_zero
rho= rho_zero

rm(Phi_zero)
rm(v_zero)
rm(rho_zero)

######################################################

for (j in 1:1000){
  
  if (j==1){
    v= v_old
    P= P_old
  }
  
  for (k in 1:p){
    w= numeric(n)
    
    for (r in 1:q){
      
      w.temp= v[k,r]*( t(t(P[,-k]) * v[-k,r])  )
      w.temp= apply(w.temp, 1, sum)
      w= w+ w.temp 
    }
    
    if (sum(w)!=0){
      w=(w-mean(w))/sd(w)
      P[,k]=w
    }
  }  
  
  #  K= (1/n)*t(P)%*%P
  
  eig= prcomp(P, center = F, scale. = F)
  v= eig$rotation[,1:q]
  rho_update=sum((eig$sdev^2)[1:q])
  
  if (abs(rho-rho_update)<=10e-5) {
    break
  }else{
    rho=rho_update
  }
}


mcpc= P%*%eig$rotation[,1:10] # top 10 PCs
