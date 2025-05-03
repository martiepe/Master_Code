for (ik in 1:100) {
  for (jk in 1:5) {
    beta <- c(4,2,-0.1)
    thin = 100
    dt = 0.01
    delta = dt*thin
    N = 49
    M = c(5, 10, 50, 100, 200)[jk]
    n_obs = 5000
    
    Tmax = n_obs*thin*dt
    
    
    
    time <- seq(0, Tmax, by = dt)
    alltimes <- list()
    for(i in 1:ntrack)
      alltimes[[i]] <- time
    
    # Generate tracks 
    alldat <- lapply(alltimes, function(times) {
      simLangevinMM(beta = beta, gamma2 = speed, times = times,
                    loc0 = c(0, 0), cov_list = covlist, silent = TRUE)
    })
    
    # Add ID column
    for(zoo in 1:ntrack)
      alldat[[zoo]] <- cbind(ID = rep(zoo, length(time)), alldat[[zoo]])
    
    
    #thinning tracks
    X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
    n = nrow(X)
    X = X[(0:(n%/%thin -1))*thin +1, ]
    
    
    gradArray = bilinearGradArray(X, covlist)
    
    times = alldat[[1]]$t[(0:(n%/%thin -1))*thin +1]
    ID = alldat[[1]]$ID[(0:(n%/%thin -1))*thin +1]
    
    gradArray = bilinearGradArray(X, covlist)
    gammasq = langevinUD(locs=X, times=times, ID=ID, grad_array=gradArray)$gamma2Hat
    
    
    
    #simulate brownian bridges
    sigma = matrix(nrow = N, ncol = N)
    for (k in 1:N) {
      for (m in 1:k) {
        sigma[k,m] = gammasq*delta*(1 - k/(N+1))*(m/(N+1))
        sigma[m,k] = gammasq*delta*(1 - k/(N+1))*(m/(N+1))
      }
    }
    
    
    ##### Constructing Brownian bridges ############
    mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
    mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
    
    B <- array(data = NA, c(2, nrow(X)-1, M, N))
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1, M))
    for (i in 1:(nrow(X)-1)) {
      mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
      mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
      
      # Generate all M sample tracks at once
      B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, mu_x, sigma = sigma)
      B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, mu_y, sigma = sigma)
      
      
      P[i, 1:M] = 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], mu_x, sigma = sigma) * 
                       mvnfast::dmvn(B[2, i, 1:M, 1:N], mu_y, sigma = sigma))
    }
    
    
    
    #find the gradient at the bridge nodes
    Grad <- array(data = NA, c(ncov +1,nrow(X)-1, M, N, 2))
    for (i in 1:(nrow(X)-1)) {
      for(j in 1:M){
        #grad = bilinearGradArray(t(B[1:2, i, j, 1:N]), covlist)
        
        grads <- bilinearGradVec(cbind(B[1, i, j, ], B[2, i, j, ]), covlist)
        
        for (k in 1:(ncov+1)) {
          Grad[k, i, j, 1:N, 1] = grads[k , ,1]
          Grad[k, i, j, 1:N, 2] = grads[k , ,2]
        }
      }
    }
    
    
    #vectorized and paralellized likelihood and gradient function
    lik_grad <- function(par, cl){
      #log-likelihood
      l = 0
      lik_grad = c(0,0,0,0)
      #number of simulations
      
      compute <- function(i){
        L_k = 1
        lik_grad_k = 0
        # Add endpoints to all samples (M x (N+2) matrices)
        # calc initial gradient
        x_samples = B[1, i, , ]
        y_samples = B[2, i, , ]
        
        grad_0 <- array(data = t(gradArray[i, , ]), c(3,1,2))
        
        u_0 <- (delta*par[4]/((N+1)*2)) * 
          (par[1] * grad_0[1,,] + par[2] * grad_0[2,,] + par[3] * grad_0[3,,])
        
        full_x <- cbind(X[i,1], x_samples, X[i+1,1])
        full_y <- cbind(X[i,2], y_samples, X[i+1,2])
        # likelihood of all locations
        L_k <- sapply(seq(M), function(j) {
          grads <- Grad[ , i, j, , ]
          us <- (delta*par[4]/((N+1)*2)) * 
            (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
          us <- rbind(u_0, us)
          prod(dmvn((cbind(full_x[j,0:N+2], full_y[j,0:N+2]) - 
                       cbind(full_x[j,0:N+1], full_y[j,0:N+1])) - us, 
                    matrix(c(0,0)),
                    diag(delta*par[4]/(N+1), 2, 2)))
        })*P[i, ]
        
        
        lik_grad_k <- sapply(seq(M), function(j){
          grads <- Grad[ , i, j, , ]
          us <- (delta*par[4]/((N+1)*2)) * 
            (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
          us <- rbind(u_0, us)
          
          g = cbind(grad_0[,1,],array(aperm(Grad[ , i, j, , ], c(1, 3, 2)), dim = c(3, 2 * N)))
          
          D = matrix(t((cbind(full_x[j,0:N+2], full_y[j,0:N+2]) - 
                          cbind(full_x[j,0:N+1], full_y[j,0:N+1])) - us), ncol = 1)
          
          rbind(g%*%D, -(N+1)/par[4] + (N+1)/(2*delta*par[4]^2)*t(D)%*%D + (1/(2*par[4]))* t(t(g)%*%par[1:3]) %*% D)
          
        })
        
        return(c(log(sum(L_k/M)), -sum(lik_grad_k[1, ]*L_k)/(2*sum(L_k)), -sum(lik_grad_k[2, ]*L_k)/(2*sum(L_k)), -sum(lik_grad_k[3, ]*L_k)/(2*sum(L_k)), -sum(lik_grad_k[4, ]*L_k)/(sum(L_k))))
      }
      
      results <- parLapply(cl, 1:(nrow(X)-1), compute)
      
      l = sum(unlist(results)[(1:(nrow(X)-1))*5 -4])
      lik_grad[1] = sum(unlist(results)[(1:(nrow(X)-1))*5 -3])
      lik_grad[2] = sum(unlist(results)[(1:(nrow(X)-1))*5 -2])
      lik_grad[3] = sum(unlist(results)[(1:(nrow(X)-1))*5 -1])
      lik_grad[4] = sum(unlist(results)[(1:(nrow(X)-1))*5])
      
      
      return(list(l = -l, g = lik_grad))
      
    }
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(16)
    clusterExport(cl, c("X", "M", "N", "delta", "P", "B", "Grad", "gradArray", "dmvn",  "lik_grad", "bilinearGradArray"))
    
    o = optim(par = c(0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B")
    
    stopCluster(cl)
    
    print(o$par)
    params[ik*5+jk-5, 1:4] = o$par
    params[ik*5+jk-5, 5] = N
  }
  
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], M = params[5])
  save(df,file="varying_M_estimates.Rda")
  
  
  print(ik)
}
