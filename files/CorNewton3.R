CorNewton3 = function(G, b, I, J, tau){
  
  ######### This code is designed to solve #############
  ##       min    0.5*<X-G, X-G>
  ##         s.t. X_ij = b_k     for (i,j) in (I,J)
  ##              X>=tau*I       X is PD
  ##
  ###########    Based on the algorithm  in  ##############
  ####  "A Quadratically Convergent Newton Method for  ####
  ####   Computing the Nearest Correlation Matrix"     ####
  ##########    By Houduo Qi and Defeng Sun    ############
  ###   SIAM J. Matrix Anal. Appl. 28 (2006) 360--385.  ###
  
  #################### Parameters 
  #   Input
  #   G       the given symmetric correlation matrix
  #   b       the right hand side of equality constraints
  #   I       row indices of the fixed elements
  #   J       column indices of the fixed elements
  #
  #   Output
  #   X         the optimal primal solution
  #   y         the optimal dual solution 
  ###################
  
  t0 <- proc.time()
  n  = dim(G)[1] 
  k  = length(b)
  
  print(sprintf(' ******************************************************** '), quote = F)
  print(sprintf( '          The Semismooth Newton-CG Method                '), quote = F)
  print(sprintf(' ******************************************************** '), quote = F)
  print(sprintf(' The information of this problem is as follows: '), quote = F)
  print(sprintf(' Dim. of    sdp      constr  = %d ',n), quote = F)
  print(sprintf(' Num. of equality    constr  = %d ',k), quote = F)
  
  G = (G + t(G))/2   # make G symmetric
  b0 = b
  
  G = G - tau*diag(1,n)      # reset G
  Ind = which(I == J)     
  b0[Ind] = b[Ind] - tau  # reset the diagonal part of b0
  
  # set parameters  
  Iter_Whole = 200
  Iter_inner = 20       # maximum num of Line Search in Newton method
  maxit      = 200      # maximum num of iterations in PCG
  
  error_tol = 1.0e-6    # termination tolerance
  tol       = 1.0e-2    # relative accuracy for CGs
  sigma     = 1.0e-4    # tolerance in the line search of the Newton method
  
  k1      = 0
  f_eval  = 0
  iter_cg = 0
  
  prec_time = 0
  pcg_time  = 0
  eig_time  = 0
  
  num = 5
  f_hist = array(0, dim = c(num, 1))
  
  # initial point
  y  = array(0, dim = c(k, 1))
  x0 = y
  
  X = matrix(0, ncol=n, nrow=n)
  
  for (i in 1:k){
    X[I[i], J[i]] = y[i]
  }
  
  X = 0.5*(X + t(X))
  X = G + X;
  X = (X + t(X))/2;
  
  t1 <- proc.time()
  Eigen_X <- Myeigen(X)
  eig_etime = proc.time() - t1
  eig_time = eig_time + eig_etime[['elapsed']]*1e0
  P <- Eigen_X$vectors
  D <- Eigen_X$values
  
  lambda = Re(D)
  P = Re(P)
  
  gradient <- gradient(y,I, J, lambda, P, X, b0)
  f0 = gradient$f
  Fy = gradient$Fy
  
  f_eval = f_eval + 1
  f = f0
  b = b0 - Fy
  norm_b = norm(b, "2")
  
  f_hist[1] = f
  
  Omega12 = omega_mat(lambda)
  
  val_G = 0.5*sum(colSums(G*G))
  Initial_f = val_G - f0
  print(sprintf("Initial Dual Objective Function value  = %.3f", Initial_f), quote = F)
  
  etime = proc.time() - t0
  etime = etime[['elapsed']]*1e0
  
  print(sprintf('   Iter.   No. of CGs     Step length      Norm of gradient     func. value      time_used '), quote = F)
 print(sprintf('    0          -              -              %3.2e         %10.9e      %3.3f ', norm_b,f, etime), quote = F)
  
  while(norm_b > error_tol && k1 < Iter_Whole){
    
    t2 <- proc.time()
    c = precond_matrix(I, J, Omega12, P)
    prec_etime = proc.time() - t2
    prec_time = prec_time + prec_etime[['elapsed']]*1e0
  
    t3 <- proc.time()
    pre = pre_cg(b, I, J, tol, maxit, Omega12, P, c)
    
    d = pre$p
    flag = pre$flag
    relres = pre$relres
    iterk = pre$iterk
    
    pcg_etime = proc.time() - t3
    pcg_time = pcg_time + pcg_etime[['elapsed']]*1e0
    iter_cg = iter_cg + iterk
    
    slope = sum((Fy-b0)*d)
    
    y = x0 + d
    
    X = matrix(0, n, n)
    for(i in 1:k){
      X[I[i], J[i]] = y[i]
    }
    
    X = 0.5*(X + t(X))
    X = G+X
    X = 0.5*(X + t(X))
    
    t1 <- proc.time()
    Eigen_X <- Myeigen(X)
    eig_etime = proc.time() - t1
    eig_time = eig_time + eig_etime[['elapsed']]*1e0
    
    P <- Eigen_X$vectors
    D <- Eigen_X$values
    
    lambda = Re(D)
    P = Re(P)
    
    gradient <- gradient(y,I, J, lambda, P, X, b0)
    f = gradient$f                                    
    Fy = gradient$Fy
    
    k_inner = 0
    while(k_inner <= Iter_inner && f > f0 + sigma*0.5^k_inner*slope + 1e-6){
      
      # backtracking
      y = x0 + 0.5^k_inner*d
      
      X = matrix(0, n, n)
      for(i in 1:k){
        X[I[i], J[i]] = y[i]
      }
      
      X = 0.5*(X + t(X))
      X = G+X
      X = 0.5*(X + t(X))
      
      
      t1 <- proc.time()
      Eigen_X <- Myeigen(X)
      eig_etime = proc.time() - t1
      eig_time = eig_time + eig_etime[['elapsed']]*1e0
      P <- Eigen_X$vectors
      D <- Eigen_X$values
      
      lambda = Re(D)
      P = Re(P)
      
      gradient <- gradient(y,I, J, lambda, P, X, b0)
      f = gradient$f                                   
      Fy = gradient$Fy
      
      k_inner = k_inner + 1
      
    } # end of line search
    
    k1 = k1 + 1
    f_eval = f_eval + k_inner + 1
    
    x0     = y
    f0     = f
    b      = b0 - Fy
    norm_b = norm(b, "2")
    
    Omega12 = omega_mat(lambda)
    
    etime = proc.time() - t0
    etime = etime[['elapsed']]*1e0
    
    print(sprintf('   %2.0d         %2.0d           %3.2e          %3.2e         %10.9e      %3.3f ',k1,iterk,0.5^k_inner,norm_b,f, etime), quote = F)
    
    # slow convergence test
    if( k1 < num ){
      
      f_hist[k1+1] = f
      
    } else{
      
      for(i in 1:(num-1)){
        
        f_hist[i] = f_hist[i+1]
        
      }
      
      f_hist[num] = f
    }
    
    if(k1 >= (num - 1) && f_hist[1] - f_hist[num] < 1e-7){
      
      print(' Progress is too slow! :( ', quote = F)
      break;
    }
    
  } # End of while loop
  
  # Optimal solution X*
  Ip = which(lambda > 0)
  r = length(Ip)
  
  if(r == 0){
    X = matrix(0, n, n)
  }else if(r==n){
    X = X
  }else if(r <= n/2){
    lambda1 = sqrt(lambda[Ip])
    
    P1 = P[, Ip]
    if(r > 1){
      P1 = P1 %*% diag(lambda1, nrow=r)
      X = P1 %*% t(P1)   # Optimal solution X*
    }else{
      X = lambda1^2 %*% P1 %*% t(P1)
    }
  }else{
    lambda2 = -lambda[(r+1):n]
    lambda2 = sqrt(lambda2)
    
    P2 = P[, (r+1):n]
    
    P2 = P2 %*% diag(lambda2, nrow=n-r)
    X = X + P2 %*% t(P2)   # Optimal solution X*
  }
  
  X = 0.5*(X + t(X))
  
  # optimal primal and dual objective values
  Final_f = val_G - f
  val_obj = sum(((X-G)*(X-G)))/2
  
  # convert to original X
  X = X + tau*diag(1,n)
  
  time_used = proc.time() - t0
  time_used = time_used[['elapsed']]*1e0
  
  sprintf('')
  print(sprintf(' ================ Final Information ================= '), quote = F)
  print(sprintf(' Total number of iterations      = %2.0f ',k1), quote = F)
  print(sprintf(' Number of func. evaluations     = %2.0f ',f_eval), quote = F)
  print(sprintf(' Number of CG Iterations         = %2.0f ',iter_cg), quote = F)
  print(sprintf(' Primal objective value          = %f ', val_obj), quote = F)
  print(sprintf(' Dual objective value            = %f ', Final_f), quote = F)
  print(sprintf(' Norm of Gradient                = %3.2e ', norm_b), quote = F)
  print(sprintf(' Rank of X-tau*I             ====== %8.0f ', r), quote = F)
  print(sprintf(' Computing time for preconditioner     = %3.1f ',prec_time), quote = F)
  print(sprintf(' Computing time for CG iterations      = %3.1f ',pcg_time), quote = F)
  print(sprintf(' Computing time for eigen-decom        = %3.1f ',eig_time), quote = F)
  print(sprintf(' Total Computing time (secs)           = %3.1f ',time_used), quote = F)
  print(sprintf(' ====================================================== '), quote = F)
  
  ######################################
  ####  End of the main program   ######
  ######################################
  
  return(list(
    X = X,
    y = y
  ))
}

##  **************************************
##  ******** All Sub-routines  ***********
##  **************************************

# To generate F(y)
gradient = function(y, I, J, lambda, P, X, b0)
{
  
  n = dim(P)[1]
  k = length(y)
  
  
  const_sparse = 2 # min(5,n/50)
  
  f = 0
  Fy = array(0, dim = c(k, 1))
  
  
  I1 = which(lambda > 1.0e-18)
  r = length(I1)
  
  if(r > 0){
    
    if(r == n){
      
      f = sum(lambda*lambda)
      
      i = 1
      
      while( i <= k){
        
        Fy[i] = X[I[i], J[i]]
        
        i = i+1
        
      }
      
    }else if(r <= n/2){
      
      lambda1 = lambda[I1]
      
      f = sum(lambda1*lambda1)
      
      lambda1 = sqrt(lambda1)
      
      P1 = P[, I1, drop=F]
      
      
      
      if (r > 1){
        
        P1 = P1 %*% diag(lambda1, nrow=r)
        
      }else{
        
        P1 = lambda1*P1
        
      }
      
      P1T = t(P1)
      
      if(k <= const_sparse*n){   # sparse form
        i = 1
        while(i <= k){
          Fy[i] = P1[I[i], ] %*% P1T[, J[i]]
          i = i+1
        }
      }else{   # dense form
        P = P1%*%P1T
        i = 1
        while(i <= k){
          Fy[i] = P[I[i], J[i]]
          i = i+1
        }
      }
      
    }else{   # n/2<r<n
      
      lambda2 = -lambda[(r+1):n]
      f = sum(lambda*lambda) - sum(lambda2*lambda2)
      
      lambda2 = sqrt(lambda2)
      
      P2 = P[, (r+1):n, drop=F]
      P2 = P2 %*% diag(lambda2, nrow = n-r)
      
      P2T =t(P2)
      
      if(k <= const_sparse*n){   # sparse form
        i = 1
        
        while(i <= k){
          Fy[i] = X[I[i], J[i]] + P2[I[i], ]%*%P2T[, J[i]]
          i = i+1
        }
        
      }else{   # dense form
        P = P2%*%P2T
        i = 1
        while(i <= k){
          Fy[i] = X[I[i], J[i]] + P[I[i], J[i]]
          i = i+1
        }
      }
      
      
    }
  }
  
  f = 0.5*f - sum(b0*y)
  
  return(list(f = f, Fy = Fy))
} # End of gradient.R

# To generate the essential part of the first-order difference of d
omega_mat <- function(lambda){
  
  # We compute omega only for 1<=|idx|<=n-1
  n = length(lambda)
  idx.idp = which(lambda > 0)
  idx.idm = setdiff(1:n, idx.idp)
  
  r = length(idx.idp)
  
  if ( sum(idx.idp >0) )
  {
    if (r==n)
    {
      Omega12 <- matrix(1,n,n)
    }
    else 
    {
      s <- n-r
      
      Omega12 <- matrix(1,nrow = r,ncol = s) 
      for (i in 1:r)
      {
        for (j in 1:s)
        {
          Omega12[i,j] <- lambda[i]/(abs(lambda[i]) + abs(lambda[r+j]));
        }
      }
    }
  }
  else
  {
    Omega12 <- matrix(nrow=0, ncol=0)
  }
  return(Omega12)
} # End of omega_mat.R

####################### PCG method  ##################################
### This is exactly the algorithm given by Hestenes and Stiefel (1952)
### An iterative method to solve A(x)=b  
### The symmetric positive definite matrix M is a preconditioner for A.
### See Pages 527 and 534 of Golub and va Loan (1996)
pre_cg <- function(b, I, J, tol, maxit, Omega12, P, c){
  
  k1 = length(b)
  dim_n = dim(P)[1]
  dim_m = dim(P)[2]
  flag <- 1
  iterk <- 0
  relres <- 1000   # give a big value on relres
  
  
  r <- b  # initial x0=0 
  n2b <- norm(b,'2')   # norm of b
  tolb <- max(tol, min(0.1, n2b))*n2b   # relative tolerance   
  
  if(n2b > 1e2){
    maxit = min(1, maxit)
  }
  
  p <- array(0,c(k1,1))
  
  # preconditioning
  z <- r/c   # z = M\r; here M =diag(c); if M is not the identity matrix
  rz1 <- sum(r*z)  
  rz2 <- 1 
  d <- z
  
  # CG iteration
  for (k in 1:maxit)
  {
    if (k > 1)
    {
      beta <- rz1/rz2
      d <- z + beta*d
    }
    
    w <- Jacobian_matrix(d, I, J, Omega12, P)   # W =A(d)
    
    if(k1 > dim_n){   # if there are more constraints than n
      w = w + 1e-2 * min(1, 0.1*n2b)*d   # perturb it to avoid numerical singularity
    }
    
    denom <- sum(d*w)
    iterk <- k
    relres <- norm(r,'2')/n2b   # relative residue=norm(r)/norm(b)
    
    if (denom <= 0 )
    {
      p <- d/norm(d)   # d is not a descent direction
      break   # exit
    }
    else
    {
      alpha <- rz1/denom
      p <- p + alpha*d
      r <- r - alpha*w
    }
    
    z <- r/c   # z = M\r; here M =diag(c); if M is not the identity matrix   
    
    if (norm(r,'2') <= tolb)   # Exit if Hp=b solved within the relative tolerance
    { 
      iterk <- k
      relres <- norm(r,'2')/n2b   # relative residue =norm(r)/norm(b)         
      flag <- 0
      break
    }
    rz2 <- rz1
    rz1 <- sum(r*z)
  }   # End of pre_cg.R
  
  list(p = p, flag = flag, relres = relres, iterk = iterk)  
  
}

# To generate the Jacobian product with x: F'(y)(x)
Jacobian_matrix = function(x, I, J, Omega12, P){
  n = dim(P)[1]
  k = length(x)
  Dim <-  dim(Omega12) 
  r <- Dim[1] 
  s <- Dim[2]
  
  
  if(r == 0){
    Ax = 1e-10*x
  }else if(r == n){
    Ax = (1 + 1e-10)*x
  }else{
    Ax = array(0,c(k,1))
    
    P1 = P[, 1:r]
    P2 = P[, (r+1):n]
    
    Z = matrix(0, ncol=n, nrow=n)
    
    for(i in 1:k){
      Z[I[i], J[i]] = x[i]
    }
    
    Z = 0.5*(Z + t(Z))
    
    const_sparse = 2 # min(5, n/50)
    
    if(k <= const_sparse*n){   # sparse form
      if(r < n/2){
        H1 = t(P1) %*% Z
        
        Omega12 = Omega12 * (H1 %*% P2)
        
        el1 = (H1%*%P1)%*%t(P1) + Omega12%*%t(P2)
        el2 = t(Omega12)%*%t(P1)
        
        H = rbind(el1, el2)
        
        i = 1
        while(i <=k){
          Ax[i] = P[I[i], ]%*% H[, J[i]]
          Ax[i] = Ax[i] + 1e-10 * x[i] # add a small perturbation
          
          i = i+1
        }
        
      }else{ # if r>=n/2, use a complementary formula. ##### giusto
        H2 = t(P2)%*%Z
        
        Omega12 = matrix(1, nrow=r, ncol=s) - Omega12
        
        Omega12 = Omega12 * t(H2%*%P1)
        
        
        el1 = Omega12 %*% t(P2)
        el2 = t(Omega12) %*% t(P1) + (H2 %*% P2 ) %*% t(P2)
        
        H = rbind(el1, el2)
        
        i = 1
        
        while(i <= k){   # AA^* is not the identity matrix
          
          if( I[i] == J[i]){
            Ax[i] = x[i] - P[I[i], ] %*% H[, J[i]]
          }else{
            Ax[i] = 0.5 * x[i] - P[I[i], ] %*% H[, J[i]]
          }
          
          Ax[i] = Ax[i] + 1e-10 * x[i] 
          
          i = i + 1
        }
        
      }
    }else{ # dense form
      #Z = full(Z); to use the full form
      # dense form
      if( r < n/2){
        H1 = t(P1)%*%Z
        Omega12 = Omega12 * (H1%*%P2)
        H = P1 %*% ((H1%*%P1) %*% t(P1) + 2*Omega12%*%t(P2))
        H = (H + t(H))/2
        
        i = 1
        while( i <= k){
          Ax[i] = H[I[i], J[i]]
          Ax[i] = Ax[i] + 1e-10 * x[i] 
          
          i = i + 1
        }
      }else{ # if r>=n/2, use a complementary formula.
        H2 = t(P2) %*% Z
        Omega12 = matrix(1, nrow=r, ncol=s) - Omega12
        Omega12 = Omega12 * t(H2%*%P1)
        H = P2 %*% ( 2*(t(Omega12) %*% t(P1)) + (H2 %*% P2) %*% t(P2) )
        H = (H + t(H))/2
        H = Z - H
        
        i = 1
        while( i <= k){   #  AA^* is not the identity matrix
          Ax[i] = H[I[i], J[i]]
          Ax[i] = Ax[i] + 1e-10 * x[i] 
          
          i = i + 1
        }
        
      } 
    } 
    
  } 
  
  
  return(Ax)
  
}   # End of Jacobian_matrix.R

# To solve the eigenvalue problem and sort eigenvalues and eigenvectors in descending order
Myeigen <- function(X,n)
{
    Eigen_X <- eigen(X,symmetric = TRUE)
    P <- Eigen_X$vectors
    lambda <- Eigen_X$values
    
    if (is.unsorted(-lambda))   #### arrange the eigenvalues in non-increasing order
    {
        lambda_sort <- sort(lambda,decreasing = TRUE,index.return = TRUE)
        lambda <- lambda_sort$x
        idx <- lambda_new$ix
        P <- P[,idx]
    }
    
    list(values = lambda, vectors = P)
    
}


# To generate the (approximate) diagonal preconditioner
precond_matrix <- function(I, J, Omega12,P){
  n = dim(P)[1]
  k = length(I)
  Dim <-  dim(Omega12) 
  r <- Dim[1] 
  s <- Dim[2]
  
  c <- array(1,c(k,1))
  
  H = t(P)
  H = H*H
  const_prec = 1
  
  if(r < n){
    if(r > 0){
      if(k <= const_prec*n){   # compute the exact diagonal preconditioner
        Ind = which(I != J)
        k1 = length(Ind)
        if(k1 > 0){
          H1 = matrix(0, nrow=n, ncol=k1)
          for(i in 1:k1){
            H1[, i] = P[I[Ind[i]], ]*t(P[J[Ind[i]], ])
          }
        }
        
        if(r < n/2){
          H12 = t(H[1:r, ,drop=F]) %*% Omega12
          if(k1 > 0){
            H12_1 = t(H1[1:r, ,drop=F]) %*% Omega12
          }
          
          d <- array(1,c(r,1))
          
          j = 0
          for(i in 1:k){
            if(I[i] == J[i]){
              c[i] = sum(H[1:r, I[i]]) %*% (t(d) %*% H[1:r, J[i]])
              c[i] = c[i] + 2*(H12[I[i], ]%*%H[(r+1):n, J[i]])
            }else{
              j = j+1
              c[i] = sum(H[1:r, I[i]]) %*% (t(d) %*% H[1:r, J[i]])
              c[i] = c[i] + 2*(H12[I[i], ] %*% H[(r+1):n, J[i]])
              c[i] = sum(H[1:r, j]) %*% (t(d) %*% H[1:r, j])
              c[i] = c[i] + 2*(H12[j, ] %*% H[(r+1):n, j]) 
              c[i] = c[i]/2
            }
            if(c[i] < 1e-8){
              c[i] = 1e-8
            }
          }
          
        }else{ # if r>=n/2, use a complementary formula
          Omega12 = matrix(1, nrow=r, ncol=s) - Omega12
          H12 = Omega12 %*% H[(r+1):n, ]
          if(k1 > 0){
            H12_1 = Omega12 %*% H1[(r+1):n, ]
          }
          
          d <- array(1,c(s,1))
          dd <- array(1,c(n,1))
          
          j = 0
          for(i in 1:k){
            if(I[i] == J[i]){
              c[i] = sum(H[(r+1):n, I[i]]) %*% (t(d) %*% H[(r+1):n, J[i]])
              c[i] = c[i] + 2*(t(H[1:r, I[i]]) %*% H12[, J[i]])
              alpha = sum(H[, I[i]])
              c[i] = alpha * (t(H[, J[i]]) %*% dd) - c[i]
            }else{
              j = j + 1
              c[i] = sum(H[(r+1):n, I[i]]) %*% (t(d) %*% H[(r+1):n, J[i]])
              c[i] = c[i] + 2*(t(H[1:r, I[i]]) %*% H12[, J[i]])
              alpha = sum(H[, I[i]])
              c[i] = alpha * (t(H[, J[i]]) %*% dd) - c[i]
              
              tmp = sum(H1[(r+1):n, j])%*%(t(d)%*%H1[(r+1):n, j])
              tmp = tmp + 2*(t(H1[1:r, j]) %*% H12_1[, j])
              alpha = sum(H1[, j])
              tmp = alpha * (t(H[, j]) %*% dd) - tmp
              
              
              c[i] = (tmp + c[i])/2
            }
            if(c[i] < 1e-8){
              c[i] = 1e-8
            }
          }
        }
        
        
      }else{   # approximate the diagonal preconditioner
        HH1 = H[1:r, , drop=F]
        HH2 = H[(r+1):n, , drop=F]
        
        if(r == 1){
          
          H0 = t(HH1) %*% (Omega12 %*% HH2)
          tmp = sum(HH1)
          H0 = H0 + t(H0) + tmp^2
          
        }else if(r < n/2){
          
          H0 = t(HH1) %*% (Omega12 %*% HH2)
          tmp = colSums(HH1)
          H0 = H0 + t(H0) + tcrossprod(tmp)
          
        }else if(n - r == 1){
          
          Omega12 = matrix(1, nrow=r, ncol=s) - Omega12
          H0 = t(HH2) %*% (t(Omega12) %*% HH1)
          tmp = sum(HH2)
          H0 = H0 + t(H0) + tmp^2
          tmp = colSums(H)
          H0 = tcrossprod(tmp) - H0 
        }else{
          
          Omega12 = matrix(1, nrow=r, ncol=s) - Omega12
          H0 = t(HH2) %*% (t(Omega12) %*% HH1)
          tmp = colSums(HH2) 
          H0 = H0 + t(H0) + tcrossprod(tmp) 
          tmp = colSums(H) 
          H0 = tcrossprod(tmp)  - H0
          
        }
        
        i = 1
        while(i <= k){
          if(I[i] == J[i]){
            c[i] = H0[I[i], J[i]]
          }else{
            c[i] = 0.5*H0[I[i], J[i]]
          }
          if(c[i] < 1e-8){
            c[i] = 1e-8
          }
          i = i+1
        } 
      }
    }
    
  }else{
    
    tmp = colSums(H)
    H0 = tmp%*%t(tmp)
    
    i = 1
    while(i <= k){
      if(I[i] == J[i]){
        c[i] = H0[I[i], J[i]]
      }else{
        c[i] = 0.5 * H0[I[i], J[i]]
      }
      if(c[i] < 1e-8){
        c[i] = 1e-8
      }
      i = i + 1
    } 
  } 
  
  return(c)
} #   End of precond_matrix.R
