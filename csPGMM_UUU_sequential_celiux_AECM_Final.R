require(ggplot2)
require(reshape2)
library(Matrix)
library(LaplacesDemon)
library(ggplot2)
library(caTools)
library(mclust)
library(pgmm)
library(nnet)
library(fossil)
library(MASS)
library(sys)
library(MixSim)
library(IMIFA)
library(parallel)
library(doParallel)

#Creating a Parsimonious Gaussian Mixture Model(UUU) Alternating Expectation - Conditional Maximization Algorithm:
pgmmAECM <- function(x,G = 1,Q = 2,epoch=100, convergence_value = 0.1,zstart=1,zlist=c(),seed=123456){
  #Input:
  #x : Dataframe(N,P). Where n is the number of observations for creating the clusters and p is the features of each observation.
  #epoch : int. Number of iterations.
  #G: int. Number of clusters.
  #Q: int. Number of factors.
  #convergence_value: convergence score. Aitken accelaration is used to assess convergence.
  #zstart: int. kmeans starting values (1) or user defined starting values(2)
  #zlist: vector. User defined starting values. Applicable only when zstart = 2
  #seed: int. seed for kmeans clustering. Applicable only when zstart = 1
  
  #Output:
  #predicted_cluster: final cluster membership of observations.
  #initial_cluster: initial cluster membership of observations.
  #predicted_score: matrix of posterior probability of cluster membership.
  #log_likelihood_list: vector of log-likelihood.
  #tau: prior probability of mixture components.
  #mu: mean of mixture components.
  #Sigma: Parsimonious covariance matrix of mixture components.
  #lambda: loadings matrix.
  #psi: diagonal covariance matrix of error/specific factor.
  #converged: whether the model is converged or not.
  #bic: Bayesian information criterion.
  #model_time: time taken to fit the model.
  
  #computing the total number of training pixels
  start_time = Sys.time()
  N=nrow(x)
  #computing the number of features for training pixels
  P=ncol(x)
  #Creating a initial cluster with kmeans clustering
  if (zstart == 1){
    set.seed(seed)
    init_cluster = kmeans(x,G,nstart = 10)$cluster
  }
  else if(zstart==2){
    init_cluster = zlist
  }
  #Creating the initial data with initial clusters for maximum log-likelihood estimation
  initial_data=cbind(x,init_cluster)
  colnames(initial_data)=c(colnames(x),"init_cluster")
  #Computing the initial parameters
  tau = matrix(0,G,1)
  mu = matrix(0,G,P)
  lambda = array(0,c(P,Q,G))
  psi = array(0,c(P,P,G))
  sigma = array(0,c(P,P,G))
  beta = matrix(0,Q,P)
  theta = matrix(0,Q,Q)
  s = array(0,c(P,P,G))
  Log_Likelihood = c()
  predicted_cluster=c()
  z=matrix(0,N,G)
  z_hat=matrix(0,N,G)
  z_hat=unmap(initial_data[,"init_cluster"])
  covw_computed = covw(x,z_hat)
  #Initialization
  for (g in 1:G){
    #Computing the probability of belonging to each class
    tau[g,]=nrow(initial_data[initial_data[,'init_cluster']==g,])/N
    #Computing the mean of each class
    mu[g,]=covw_computed$mean[,g]
    #Computing the covariance matrix of each class
    s[,,g]=covw_computed$S[,,g]
    eigen_g=eigen(s[,,g])
    for (q in 1:Q){
      for (p in 1:P){
        if (eigen_g$vectors[p,1] < 0){
          lambda[p,q,g] = -exp((0.5*log(eigen_g$values[q])+log(-eigen_g$vectors[p,1])))
        }
        else{
          lambda[p,q,g] = exp((0.5*log(eigen_g$values[q])+log(eigen_g$vectors[p,1])))
        }
      }
    }
    diag(psi[,,g]) = abs(diag(s[,,g]-lambda[,,g]%*%t(lambda[,,g])))
    # Updating sigma_g based on latent factors
    sigma[,,g] = lambda[,,g]%*%t(lambda[,,g]) + psi[,,g]
  }
  # AECM Algorithm
  for (e in 0:epoch){
    if(e%%100 == 0){
      print(e)      
    }
    #Step 1
    #Conditional Maximization
    covw_computed = covw(x,z_hat)
    for (g in 1:G){
      #Computing the probability of belonging to each class
      tau[g,]=(sum(z_hat[,g]))/N
      #Computing the mean of each class
      mu[g,]=covw_computed$mean[,g]
    }
    #Creating a temporary variable
    z_old = z_hat
    #Expectation
    if(e!=0){
      #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
      for (g in 1:G){
        z[,g]=log(tau[g,])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log=TRUE)
      }
      #Computing Log-Sum-Exp
      z_max = apply(z,1,max)
      z_sum=log(apply(exp(z-z_max),1,sum))+z_max
      z_hat=exp(z-z_sum)
    }
    # Step 2:
    # Update s,beta,theta and Conditional Maximization
    covw_computed = cov_weighted(x,G,mu,z_hat)
    for (g in 1:G){
      s[,,g]=covw_computed[,,g]
      # Using woodbury identity for computing the Inverse of lambda%*%t(lambda) + psi
      psi_inv = solve(psi[,,g])
      beta = t(lambda[,,g]) %*% (psi_inv - (psi_inv %*% lambda[,,g] %*% solve(diag(Q)+t(lambda[,,g]) %*% psi_inv %*% lambda[,,g]) %*% t(lambda[,,g]) %*% psi_inv))
      theta = diag(Q) - (beta %*% lambda[,,g]) + (beta %*% s[,,g] %*% t(beta))
      # Estimating new lambda
      lambda[,,g] = s[,,g] %*% t(beta) %*% solve(theta)
      # Estimating new psi
      diag(psi[,,g]) = (diag(s[,,g]-lambda[,,g]%*%beta%*%s[,,g]))
      # Estimating new sigma based on latent factors
      sigma[,,g] = lambda[,,g]%*%t(lambda[,,g]) + psi[,,g]
    }
    # Expectation
    #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
    for (g in 1:G){
      z[,g]=log(tau[g,])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log = TRUE)
    }
    #Computing Log-Sum-Exp
    z_max = apply(z,1,max)
    z_sum=log(apply(exp(z-z_max),1,sum))+z_max
    z_hat=exp(z-z_sum)
    print(sum(z_sum))
    Log_Likelihood =c(Log_Likelihood,sum(z_sum))
    write_time = Sys.time()
    if(e >= 3){
      log_lik_change = (Log_Likelihood[length(Log_Likelihood)]-Log_Likelihood[length(Log_Likelihood)-1])/(Log_Likelihood[length(Log_Likelihood)-1])
      if(abs(log_lik_change)< convergence_value){
        converged = 'Converged'
        print('Cluster solution converged')
        break
      }
    }
    if (e == epoch){
      converged = 'not converged'
      print('Cluster Solution not Converged')
      break
    }
  }
  #bic_model = 2*tail(Log_Likelihood,1) - PGMM_dfree(Q,P,G,'UUU') * log(N)
  bic_model = 'not_computed'
  final_prediction=data.frame(z_hat)
  colnames(final_prediction) = c(1:G)
  predicted_cluster = colnames(final_prediction)[max.col(final_prediction, ties.method = "first")]
  end_time = as.numeric(Sys.time() - start_time,units = "mins")
  return(list(predicted_cluster=predicted_cluster,mu=mu,sigma=sigma,lambda=lambda,psi=psi))
}

cov_weighted = function(x,G=1,mu,z_hat){
  s = array(0,c(ncol(x),ncol(x),G))
  for(g in 1:G){
    x_temp = sweep(x,2,mu[g,],'-')
    x_temp_2 = sweep(x_temp,1,z_hat[,g],'*')
    s[,,g] = (t(x_temp_2) %*% as.matrix(x_temp))/sum(z_hat[,g])
  }
  return(s)
}

aitken_acceleration = function(l1,l2,l3){
  a_k = (l3 - l2)/(l2-l1)
  if(a_k<1.0){
    l_inf = l2 + ((l3-l2)/(1-a_k))
  }
  else{
    l_inf = "invalid"
  }
  return(l_inf)
}


neighbour_observations_1layer <-function(n_y,n_x,circle,start_value){
  W<-Matrix(0,nrow = (n_y*n_x),ncol = c(n_y*n_x))
  n = 1
  for(x in 1:n_x){
    print(x)
    for(y in 1:n_y){
      temp = c()
      if(x-1>0){
        if(y-1>0){
          temp = c(temp,(y-1)+((x-2)*n_y) + start_value)
        }
        temp = c(temp,(y)+((x-2)*n_y) + start_value)
        if(y+1<=n_y){
          temp = c(temp,(y+1)+((x-2)*n_y) + start_value)
        }
      }
      if(x+1<=n_x){
        if(y-1>0){
          temp = c(temp,(y-1)+((x)*n_y) + start_value)
        }
        temp = c(temp,(y)+((x)*n_y) + start_value)
        if(y+1<=n_y){
          temp = c(temp,(y+1)+((x)*n_y) + start_value)
        }        
      }
      if(y-1>0){
        temp = c(temp,(y-1)+((x-1)*n_y) + start_value)
      }
      if(y+1<=n_y){
        temp = c(temp,(y+1)+((x-1)*n_y) + start_value)
      } 
      W[n,sort(temp)] = 1
      n = n+1
    }
  }
  return(W)
}
spatial_tau_estimation <- function(reg_coeff,neighbour_pixels_samedifferentclass,N,G){
  log_odds = matrix(0,nrow = N,ncol = G)
  for(g in 1:G){
    log_odds[,g] = (cbind(1,neighbour_pixels_samedifferentclass[,g])%*%reg_coeff[g,])[,1]
  }
  log_odds = exp(log_odds)
  #Handling +Inf values
  log_odds[!is.finite(log_odds)] = exp(709.78)
  tau <- log_odds/rowSums(log_odds)
  return(tau)
}

spatial_multinomial_logistic_log_likelihood<-function(params,neighbour_pixels_samedifferentclass,N,G,z_hat){
  reg_coeff = matrix(0,nrow = G,ncol = 2)
  new_gamma = matrix(params,nrow=G-1,ncol = 2,byrow = FALSE)
  reg_coeff[c(2:G),]<- new_gamma
  probabilities = spatial_tau_estimation(reg_coeff,neighbour_pixels_samedifferentclass,N,G)
  log_likelihood = 0
  for(g in 1:G){
    temp = log(probabilities[,g])
    temp[!is.finite(temp)] = log(2.5e-324)
    log_likelihood = log_likelihood + (t(z_hat[,g]) %*% temp)[,1]
  }
  print(log_likelihood)
  return(-log_likelihood)
}

# Define the gradient of the logistic regression log likelihood function
spatial_multinomial_gradient <- function(params,neighbour_pixels_samedifferentclass,N,G,z_hat){
  gradient = matrix(0,nrow = G-1,ncol = 2)
  reg_coeff = matrix(0,nrow = G,ncol = 2)
  new_gamma = matrix(params,nrow=G-1,ncol = 2,byrow = FALSE)
  reg_coeff[c(2:G),]<- new_gamma
  probabilities = spatial_tau_estimation(reg_coeff,neighbour_pixels_samedifferentclass,N,G)
  for(g in c(2:G)){
    gradient[g-1,] = (t(cbind(1,neighbour_pixels_samedifferentclass[,g]))%*%(z_hat[,g]-probabilities[,g]))[,1]
  }
  return(-as.matrix(gradient))
}

log_likelihood_calculation <- function(reg_coeff,mu,sigma,x,neighbour_pixels_samedifferentclass,N,G){
  z = matrix(0,nrow = N,ncol = G)
  tau = spatial_tau_estimation(reg_coeff,neighbour_pixels_samedifferentclass,N,G)
  #Expectation
  #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
  for (g in 1:G){
      z[,g]=log(tau[,g])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log=TRUE)
  }
  z_max = apply(z,1,max)
  z_sum=log(apply(exp(z-z_max),1,sum))+z_max
  z_hat=exp(z-z_sum)
  #Computing the log likelihood of the prediction
  return(sum(z_sum))
}

mode_approx_per_image <- function(pixel_order,W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff){
  density = z_map
  z_hat = z_map
  z = z_map
  tau = z_map
  for(g in c(1:G)){
    density[,g] = dmvn(as.matrix(x),mu[g,],sigma[,,g],log=TRUE)
  }
  for(n in pixel_order){
    print(n)
    neighbour_pixels_sameclass[n,] = c(W[n,]%*%z_map)
    neighbour_pixels_differentclass[n,] = (sum(neighbour_pixels_sameclass[n,])-neighbour_pixels_sameclass[n,])
    neighbour_pixels_samedifferentclass[n,] = neighbour_pixels_sameclass[n,] - neighbour_pixels_differentclass[n,]
    for(g in c(1:G)){
      tau[n,g] = exp((cbind(1,neighbour_pixels_samedifferentclass[n,g])%*%reg_coeff[g,])[,1])
      if(!is.finite(tau[n,g])){
        tau[n,g] = exp(709.78)
      }
    }
    tau[n,] = tau[n,]/sum(tau[n,])
    for(g in c(1:G)){
      z[n,g] = c(log(tau[n,g])) + density[n,g]
    }
    #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
    z_max = max(z[n,])
    z_sum = log(sum(exp(z[n,]-z_max))) + z_max
    z_hat[n,] = exp(z[n,]-z_sum)
    z_map[n,] = 0
    z_map[n,which.max(z_hat[n,])] = 1
  }
  return(z_map[c(min(pixel_order):max(pixel_order)),])
}

#Creating a Parsimonious Gaussian Mixture Model(UUU) Alternating Expectation - Conditional Maximization Algorithm:
spatial_pgmmAECM_UUU <- function(x,colnum,W,n_images,image_axs_list,image_obs_list,G = 1,Q = 2,epoch=100, convergence_value = 1e-8,zstart=1,zlist=c(),seed=123456){
  #Input:
  #x : Dataframe(N,P). Where n is the number of observations for creating the clusters and p is the features of each observation.
  #epoch : int. Number of iterations.
  #G: int. Number of clusters.
  #Q: int. Number of factors.
  #convergence_value: convergence score. Aitken accelaration is used to assess convergence.
  #zstart: int. kmeans starting values (1) or user defined starting values(2)
  #zlist: vector. User defined starting values. Applicable only when zstart = 2
  #seed: int. seed for kmeans clustering. Applicable only when zstart = 1
  
  #Output:
  #predicted_cluster: final cluster membership of observations.
  #initial_cluster: initial cluster membership of observations.
  #predicted_score: matrix of posterior probability of cluster membership.
  #log_likelihood_list: vector of log-likelihood.
  #tau: prior probability of mixture components.
  #mu: mean of mixture components.
  #Sigma: Parsimonious covariance matrix of mixture components.
  #lambda: loadings matrix.
  #psi: diagonal covariance matrix of error/specific factor.
  #converged: whether the model is converged or not.
  #bic: Bayesian information criterion.
  #model_time: time taken to fit the model.
  
  #computing the total number of training pixels
  start_time = Sys.time()
  N=nrow(x)
  #computing the number of features for training pixels
  P=ncol(x)
  result = list()
  reg_ord_list = list()
  rev_ord_list = list()
  left_right_list = list()
  right_left_list = list()
  for(i in c(1:n_images)){
    reg_ord_list = append(reg_ord_list,list(list(c(image_obs_list[[i]][1]:image_obs_list[[i]][2]))))
    rev_ord_list = append(rev_ord_list,list(list(c(image_obs_list[[i]][2]:image_obs_list[[i]][1]))))
    left_right_list = append(left_right_list,list(list(as.vector(matrix(c(image_obs_list[[i]][1]:image_obs_list[[i]][2]),ncol = image_axs_list[[i]][1],nrow = image_axs_list[[i]][2],byrow = TRUE)))))
    right_left_list = append(right_left_list,list(list(as.vector(matrix(c(image_obs_list[[i]][2]:image_obs_list[[i]][1]),ncol = image_axs_list[[i]][1],nrow = image_axs_list[[i]][2],byrow = TRUE)))))
  }
  #print(left_right_list)
  #Creating a initial cluster with kmeans clustering
  if (zstart == 1){
    set.seed(seed)
    init_cluster = kmeans(x,G,nstart = 10)$cluster
  }
  else if(zstart==2){
    init_cluster = zlist
  }
  else if(zstart==3){ 
    #pgmm_result = pgmmAECM(x,G = G,Q = Q,epoch=2000, convergence_value = 1e-8,zstart=1,seed=123456)
    #init_cluster = pgmm_result$predicted_cluster
    loc = paste('originalpuffedcereal_cPGMM_r25_d20_',as.character(colnum),'.RData',sep = '')
    load(loc)
    init_cluster = as.numeric(s[[3]])
    pgmm_result = s
  }
  #Creating the initial data with initial clusters for maximum log-likelihood estimation
  initial_data=cbind(x,init_cluster)
  colnames(initial_data)=c(colnames(x),"init_cluster")
  #Computing the initial parameters
  tau = matrix(0,N,G)
  mu = matrix(0,G,P)
  lambda = array(0,c(P,Q,G))
  psi = array(0,c(P,P,G))
  sigma = array(0,c(P,P,G))
  beta = matrix(0,Q,P)
  theta = matrix(0,Q,Q)
  s = array(0,c(P,P,G))
  reg_coeff = matrix(0,ncol = 2,nrow = G)
  neigbour_pixels_sameclass = matrix(0,ncol = G,nrow = N)
  neighbour_pixels_differentclass = matrix(0,ncol = G,nrow = N)
  Log_Likelihood = c()
  predicted_cluster=c()
  z=matrix(0,N,G)
  z_hat=matrix(0,N,G)
  z_hat=unmap(initial_data[,"init_cluster"])
  z_map = z_hat
  #Initialization
  neighbour_pixels_sameclass = W%*%z_map
  neighbour_pixels_differentclass = (rowSums(neighbour_pixels_sameclass)-neighbour_pixels_sameclass)
  neighbour_pixels_samedifferentclass = neighbour_pixels_sameclass - neighbour_pixels_differentclass
  same_differentclass_diff <- mapply(function(row, col) neighbour_pixels_sameclass[row, col]-neighbour_pixels_differentclass[row, col], 
                                     row = seq_len(nrow(W)), 
                                     col = as.numeric(init_cluster))    
  reg_coeff[c(2:G),] = coef(multinom(init_cluster~same_differentclass_diff),matrix = TRUE)
  out = optim(par = reg_coeff[c(2:G),],fn = spatial_multinomial_logistic_log_likelihood,gr = spatial_multinomial_gradient,
              neighbour_pixels_samedifferentclass = neighbour_pixels_samedifferentclass,z_hat = z_hat, G = G,N = N,lower = c(rep(-Inf,2*(G-1))), upper = c(rep(+Inf,2*(G-1))),method = 'L-BFGS-B',control = list(factr = 1e7,maxit = 2000))
  reg_coeff[c(2:G),]<- out$par
  covw_computed = covw(x,z_hat)
  if(zstart==3){
     for (g in 1:G){
        #Computing the mean of each class
        mu[g,] = pgmm_result[[10]][g,]
        sigma[,,g] = pgmm_result[[13]][,,g]
        psi[,,g] = pgmm_result[[12]][,,g]
        lambda[,,g] = pgmm_result[[11]][,,g]
     }
  }
  else{
     for (g in 1:G){
        #Computing the mean of each class
        mu[g,]=covw_computed$mean[,g]
        #Computing the covariance matrix of each class
        s[,,g]=covw_computed$S[,,g]
        eigen_g=eigen(s[,,g])
        for (q in 1:Q){
          for (p in 1:P){
            if (eigen_g$vectors[p,1] < 0){
               lambda[p,q,g] = -exp((0.5*log(eigen_g$values[q])+log(-eigen_g$vectors[p,1])))
            }
            else{
               lambda[p,q,g] = exp((0.5*log(eigen_g$values[q])+log(eigen_g$vectors[p,1])))
            }
          }
        }
        diag(psi[,,g]) = abs(diag(s[,,g]-lambda[,,g]%*%t(lambda[,,g])))
        # Updating sigma_g based on latent factors
        sigma[,,g] = lambda[,,g]%*%t(lambda[,,g]) + psi[,,g]
     }
  }
  save_num = 0
  interim_out = 0
  interim_out_2 = 0
  set.seed(100)
  # AECM Algorithm
  for (e in 1:epoch){
    print(e)
    if(e%%4==1){
        mode_approx = foreach(m = reg_ord_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    }
    else if(e%%4==2){
        mode_approx = foreach(m = rev_ord_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    }
    #else if(e%%4==3){
    #    mode_approx = foreach(m = reg_ord_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(sample(m[[1]]),W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    #}
    else if(e%%5==3){
        mode_approx = foreach(m = left_right_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    }
    else if(e%%4==0){
        mode_approx = foreach(m = right_left_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    }
    if(e%%100 == 0){
      print(e)      
    }
    #Expectation
    z_map = mode_approx
    neighbour_pixels_sameclass = W%*%z_map
    neighbour_pixels_differentclass = (rowSums(neighbour_pixels_sameclass)-neighbour_pixels_sameclass)
    neighbour_pixels_samedifferentclass = neighbour_pixels_sameclass - neighbour_pixels_differentclass
    tau = spatial_tau_estimation(reg_coeff,neighbour_pixels_samedifferentclass,N,G)
    #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
    for (g in 1:G){
      z[,g]=log(tau[,g])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log=TRUE)
    }
    #Computing Log-Sum-Exp
    z_max = apply(z,1,max)
    z_sum=log(apply(exp(z-z_max),1,sum))+z_max
    z_hat=exp(z-z_sum)
    #Step 1
    #Conditional Maximization
    covw_computed = covw(x,z_hat)
    for (g in 1:G){
      #Computing the mean of each class
      mu[g,]=covw_computed$mean[,g]
    }
    out = optim(par = reg_coeff[c(2:G),],fn = spatial_multinomial_logistic_log_likelihood,gr = spatial_multinomial_gradient,
                neighbour_pixels_samedifferentclass = neighbour_pixels_samedifferentclass,z_hat = z_hat, G = G,N = N,lower = c(rep(-Inf,2*(G-1))), upper = c(rep(+Inf,2*(G-1))),method = 'L-BFGS-B',control = list(factr = 1e7,maxit = 2000))
    reg_coeff[c(2:G),]<- out$par
    #Mode Field Approximation
    if(e%%4==1){
        mode_approx = foreach(m = reg_ord_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    }
    else if(e%%4==2){
        mode_approx = foreach(m = rev_ord_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    }
    #else if(e%%4==3){
    #    mode_approx = foreach(m = reg_ord_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(sample(m[[1]]),W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    #}
    else if(e%%5==3){
        mode_approx = foreach(m = left_right_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    }
    else if(e%%4==0){
        mode_approx = foreach(m = right_left_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,x,mu,sigma,G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,reg_coeff))
    }
    #Expectation
    z_map = mode_approx
    neighbour_pixels_sameclass = W%*%z_map
    neighbour_pixels_differentclass = (rowSums(neighbour_pixels_sameclass)-neighbour_pixels_sameclass)
    neighbour_pixels_samedifferentclass = neighbour_pixels_sameclass - neighbour_pixels_differentclass
    tau = spatial_tau_estimation(reg_coeff,neighbour_pixels_samedifferentclass,N,G)
    #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
    for (g in 1:G){
      z[,g]=log(tau[,g])+dmvn(as.matrix(x),mu[g,],sigma[,,g],log=TRUE)
    }
    #Computing Log-Sum-Exp
    z_max = apply(z,1,max)
    z_sum=log(apply(exp(z-z_max),1,sum))+z_max
    z_hat=exp(z-z_sum)
    # Step 2 CM:
    # Update s,beta,theta and Conditional Maximization
    covw_computed = cov_weighted(x,G,mu,z_hat)
    for (g in 1:G){
      s[,,g]=covw_computed[,,g]
      # Using woodbury identity for computing the Inverse of lambda%*%t(lambda) + psi
      psi_inv = solve(psi[,,g])
      beta = t(lambda[,,g]) %*% (psi_inv - (psi_inv %*% lambda[,,g] %*% solve(diag(Q)+t(lambda[,,g]) %*% psi_inv %*% lambda[,,g]) %*% t(lambda[,,g]) %*% psi_inv))
      theta = diag(Q) - (beta %*% lambda[,,g]) + (beta %*% s[,,g] %*% t(beta))
      # Estimating new lambda
      lambda[,,g] = s[,,g] %*% t(beta) %*% solve(theta)
      # Estimating new psi
      diag(psi[,,g]) = (diag(s[,,g]-lambda[,,g]%*%beta%*%s[,,g]))
      # Estimating new sigma based on latent factors
      sigma[,,g] = lambda[,,g]%*%t(lambda[,,g]) + psi[,,g]
      #print(sigma[,,g])
    }
    log_lik = log_likelihood_calculation(reg_coeff,mu,sigma,x,neighbour_pixels_samedifferentclass,N,G)
    Log_Likelihood = c(Log_Likelihood,log_lik)
    if(e%%10==0){
	result = list()
        save_num = save_num + 1
    }
    result[[e]] = list(e,z_hat)
    file_name = paste('originalpuffedcerealcs_csPGMM_UUU_1e_8_',as.character(colnum),'_G4_epoch_',as.character(save_num),'_final_v3.RData',sep = '')
    save(result,file = file_name)
    end_time = as.numeric(Sys.time() - start_time,units = "mins")
    file_name_2 = paste('originalpuffedcerealcs_csPGMM_UUU_1e_8_',as.character(colnum),'_G4_last_params_final_v3.RData',sep = '')
    last_params = list(z_map,z_hat,reg_coeff,mu,sigma,tau,Log_Likelihood,lambda,psi,beta,theta,end_time)
    save(last_params,file = file_name_2)
    if(e >= 3){
      log_lik_change = (Log_Likelihood[length(Log_Likelihood)]-Log_Likelihood[length(Log_Likelihood)-1])/(Log_Likelihood[length(Log_Likelihood)-1])
      if(abs(log_lik_change)<1e-6 & interim_out == 0){
	interim_result = list(e,z_hat)
        file_name3 = paste('originalpuffedcerealcs_csPGMM_UUU_1e_8_',as.character(colnum),'_G4_1e_6_interim_result_final_v3.RData',sep = '')
    	save(interim_result,file = file_name3)
        interim_out = interim_out + 1
      }
      if(abs(log_lik_change)<1e-7 & interim_out_2 == 0){
	interim_result = list(e,z_hat)
        file_name3 = paste('originalpuffedcerealcs_csPGMM_UUU_1e_8_',as.character(colnum),'_G4_1e_7_interim_result_final_v3.RData',sep = '')
    	save(interim_result,file = file_name3)
        interim_out_2 = interim_out_2 + 1
      }
      if(abs(log_lik_change)< convergence_value){
        converged = 'Converged'
        print('Cluster solution converged')
        break
      }
    }
    if (e == epoch){
      converged = 'not converged'
      print('Cluster Solution not Converged')
      break
    }
  }
  #bic_model = 2*tail(Log_Likelihood,1) - PGMM_dfree(Q,P,G,'UUU') * log(N)
  bic_model = 'not_computed'
  final_prediction=data.frame(z_hat)
  colnames(final_prediction) = c(1:G)
  predicted_cluster = colnames(final_prediction)[max.col(final_prediction, ties.method = "first")]
  end_time = as.numeric(Sys.time() - start_time,units = "mins")
  return(list(predicted_cluster=predicted_cluster,initial_cluster = init_cluster, probability_score=z_hat,log_likelihood_list = Log_Likelihood,tau=tau,reg_coeff = reg_coeff,mu=mu,sigma=sigma,lambda=lambda,psi=psi,converged = converged,bic = bic_model,time = end_time))
}

# Data
Data = read.csv('Data.csv')
Data$X<-NULL

load('cols_list.RData')
n_cols_list = 1:length(cols_list)

load('puffed_cereal_neighbourmatrix_W.RData')
Image_x_y = list(c(67,53),c(62,64),c(65,55),c(61,34),c(39,69),c(56,54),c(56,49),c(78,54),c(50,44))
obs_list = list(c(1,3551),c(3552,7519),c(7520,11094),c(11095,13168),c(13169,15859),c(15860,18883),c(18884,21627),c(21628,25839),c(25840,28039))


G = 4
N = nrow(Data)
P = ncol(Data)
Q = 2
e = 10000

exp_funcs = c('pgmmAECM','cov_weighted','aitken_acceleration','neighbour_observations_1layer','spatial_tau_estimation','spatial_multinomial_logistic_log_likelihood','spatial_multinomial_gradient','log_likelihood_calculation','mode_approx_per_image')

no_cores <- detectCores()
cl <- makeCluster(120)  
print(length(cl))
registerDoParallel(cl)  
clusterEvalQ(cl, {
  require(ggplot2)
  require(reshape2)
  library(Matrix)
  library(LaplacesDemon)
  library(ggplot2)
  library(caTools)
  library(mclust)
  library(pgmm)
  library(nnet)
  library(fossil)
  library(MASS)
  library(sys)
  library(MixSim)
  library(IMIFA)
  library(parallel)
  library(doParallel)
})


print('simulation begins')
r_index = 1
for(r in c(25)){
      start_time = Sys.time()
      #spgmmAECM
      result_final <- foreach(m = n_cols_list,.export = exp_funcs) %dopar% spatial_pgmmAECM_UUU(Data[,cols_list[[m]][[1]]],m,W,9,Image_x_y,obs_list,G = G,Q = Q,epoch=e, convergence_value = 1e-8,zstart=3,seed=123456)
      print('Models computed')
      end_time = as.numeric(Sys.time() - start_time,units = "mins")
      print(end_time)
      parallel_list = list(no_models = r,parallel_time = end_time)
      n_result = 1
      for(s in result_final){
        location = paste('originalpuffedcerealcs_csPGMM_UUU_1e_8_r',as.character(r),'_d',as.character(20),'_',as.character(n_result),'_final_v3.RData',sep = '')
        save(s,file = location)
        n_result=n_result+1
      }
      location = paste('originalpuffedcerealcs_csPGMM_UUU_1e_8_r',as.character(r),'_d',as.character(20),'_creation_time_final_v3.RData',sep = '')
      save(parallel_list,file=location)
}

stopCluster(cl)


data = read.csv('Mixture_Image.csv')
data$X<-NULL
load('mixture_image_cereal_neighbourmatrix_W.RData')
N = nrow(data)
G = 4
z=matrix(0,N,G)
z_hat=matrix(0,N,G)
pred_cluster_list = list()
for(ss in c(1:25)){
  loc = paste('originalpuffedcereal_csPGMM_UUU_1e_8_',as.character(ss),'_G4_last_params_final.RData',sep = '')
  load(loc)
  for (g in 1:G){
    z[,g]=log(0.25)+dmvn(as.matrix(data[,cols_list[[ss]][[1]]]),last_params[[4]][g,],last_params[[5]][,,g],log=TRUE)
  }
  z_max = apply(z,1,max)
  z_sum=log(apply(exp(z-z_max),1,sum))+z_max
  z_hat=exp(z-z_sum)
  final_prediction=data.frame(z_hat)
  colnames(final_prediction) = c(1:G)
  init_cluster = colnames(final_prediction)[max.col(final_prediction, ties.method = "first")]
  table(init_cluster)
  neigbour_pixels_sameclass = matrix(0,ncol = G,nrow = N)
  neighbour_pixels_differentclass = matrix(0,ncol = G,nrow = N)
  z_hat=unmap(init_cluster)
  z_map = z_hat
  #Initialization
  neighbour_pixels_sameclass = W%*%z_map
  neighbour_pixels_differentclass = (rowSums(neighbour_pixels_sameclass)-neighbour_pixels_sameclass)
  neighbour_pixels_samedifferentclass = neighbour_pixels_sameclass - neighbour_pixels_differentclass
  same_differentclass_diff <- mapply(function(row, col) neighbour_pixels_sameclass[row, col]-neighbour_pixels_differentclass[row, col], 
                                     row = seq_len(nrow(W)), 
                                     col = as.numeric(init_cluster)) 
  
  reg_ord_list = list()
  reg_ord_list = append(reg_ord_list,list(list(c(1:43621))))
  mode_approx = foreach(m = reg_ord_list,.combine = rbind,.export='mode_approx_per_image')%dopar%(mode_approx_per_image(m[[1]],W,z_map,data[,cols_list[[ss]][[1]]],last_params[[4]],last_params[[5]],G,neighbour_pixels_sameclass,neighbour_pixels_differentclass,neighbour_pixels_samedifferentclass,last_params[[3]]))
  
  z_map = mode_approx
  neighbour_pixels_sameclass = W%*%z_map
  neighbour_pixels_differentclass = (rowSums(neighbour_pixels_sameclass)-neighbour_pixels_sameclass)
  neighbour_pixels_samedifferentclass = neighbour_pixels_sameclass - neighbour_pixels_differentclass
  tau = spatial_tau_estimation(last_params[[3]],neighbour_pixels_samedifferentclass,N,G)
  #Estimating the probability of each training pixels belonging to each class based on the estimated parameters
  for (g in 1:G){
    z[,g]=log(tau[,g])+dmvn(as.matrix(data[,cols_list[[ss]][[1]]]),last_params[[4]][g,],last_params[[5]][,,g],log=TRUE)
  }
  #Computing Log-Sum-Exp
  z_max = apply(z,1,max)
  z_sum=log(apply(exp(z-z_max),1,sum))+z_max
  z_hat=exp(z-z_sum)
  final_prediction=data.frame(z_hat)
  colnames(final_prediction) = c(1:G)
  pred_cluster = colnames(final_prediction)[max.col(final_prediction, ties.method = "first")]
  table(pred_cluster)
  pred_cluster_list[[ss]] = pred_cluster
}
save(pred_cluster_list,file = 'pred_cluster_list.RData')

table(pred_cluster_list[[1]])
table(pred_cluster_list[[25]])
temp = rep(0,43621)
temp[which(pred_cluster_list[[25]]==1)] = 2
table(temp)
table(pred_cluster_list[[1]],pred_cluster_list[[25]])
table(pred_cluster_list[[1]],temp)

pred_cluster_list[[25]] = temp

save(pred_cluster_list,file = 'pred_cluster_list.RData')

pred_cluster_dataframe = matrix(0,nrow = 43621,ncol = 25)
for(i in c(1:25)){
  pred_cluster_dataframe[,i] = as.numeric(pred_cluster_list[[i]])
}

# Function to find most frequent value(s) and how many times it occurs
max_freq_value <- function(x) {
  freq_table <- table(x)
  max_count <- max(freq_table)
  max_values <- names(freq_table)[freq_table == max_count]
  
  # If there's more than one value with same max frequency, list all
  list(values = max_values, count = max_count)
}

# Apply to each row and get results
results <- apply(pred_cluster_dataframe, 1, max_freq_value)

# Transform into a readable data frame
output <- data.frame(
  row = paste0("Row", seq_along(results)),
  value = sapply(results, function(x) paste(x$values, collapse = ", ")),
  count = sapply(results, function(x) x$count)
)

print(output)
write.csv(output,file = 'mixture_image_csPGMM_predicted_cluster_count.csv')
write.csv(pred_cluster_dataframe,file = 'mixture_image_csPGMM_all_predicted_cluster.csv')
