
















delta = dt*thin/(m+1)
lik3 <- function(par){
  #log-likelihood
  l = 0
  for (j in 1:1) {
    #defining transition covariance matrix
    Q = diag(delta*par[4],2,2)
    
    #control matrix
    B = diag(delta*par[4]/2,2,2)
    
    #initial covariance guess
    P = diag(10*Q, diag(100, 3, 3))
    
    #initial state
    z = c(X[1, ], 0, 0 , 0)
    
    
    for (i in 2:nrow(X)) {
      #control vector
      u = bilinearGrad(z, covlist) 
      
      #transition matrix
      F_k = cbind(rbind(diag(1,2,2) + (delta*par[4]/2)*hessian(z, covlist, par), matrix(rep(0, 6), nrow = 3, ncol = 2)),
                  rbind(u, diag(1,3,3)))
      
      #predicted state estimate
      z_p = c(z[1:2] + B %*% u %*% par[1:3], z[3:5])
      
      #predicted estimate covariance 
      P = F_k %*% P %*% t(F_k) + diag(Q, diag(0, 3, 3))
      
      
      #updated state estimate
      z =  c(X[i, ], P[3:5, 1:2] %*% solve(P[1:2, 1:2] %*% (X[i, ] - z_p[1:2])))
      

      #updated estimate covariance
      P = diag(diag(0,2,2), P[1:3, 1:3] - P[1:3, 1:2] %*% solve(P[1:2, 1:2]) %*% P[1:2, 1:3])
      
      #adding likelihood contribution of i-th state
      l = l - dmvnorm(c(X[i, ] - z_p), mean = c(0,0), sigma = P[1:2, 1:2], log = T)
      
      
      for (k in 1:m) {
        #control vector
        u = bilinearGrad(z, covlist) 
        
        #transition matrix
        F_k = cbind(rbind(diag(1,2,2) + (delta*par[4]/2)*hessian(z, covlist, par), matrix(rep(0, 6), nrow = 3, ncol = 2)),
                    rbind(u, diag(1,3,3)))
        
        #predicted state estimate
        z_p = c(z[1:2] + B %*% u %*% par[1:3], z[3:5])
        
        #predicted estimate covariance 
        P = F_k %*% P %*% t(F_k) + diag(Q, diag(0, 3, 3))
        
        #updated state estimate
        z = z_p 
        
        #updated estimate covariance
        P = P
        
      }
      
    }
  }
  return(l)
  
}






































