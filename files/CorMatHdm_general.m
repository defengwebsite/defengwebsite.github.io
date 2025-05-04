function   [X0,z_e0,z_l0,z_u0,Gamma0,val_obj] = CorMatHdm_general(G,H,e,I_e,J_e,l,I_l,J_l,u,I_u,J_u,tau,TOL1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This code is designed to solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              min  0.5 ||H o ( X - G )||^2        % "o" is the Hadamard product symbol
%%              s.t.  X_ij   = e_ij,   for (i,j) in (I_e,J_e)
%%                    X_ij  >= l_ij,   for (i,j) in (I_l,J_l)
%%                    X_ij  <= u_ij,   for (i,j) in (I_u,J_u)
%%                      X   >= tau*I   (symmetric and positive semi-definite; tau can be zero)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% based on the algorithm  in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% "An augmented Lagrangian dual approach for the H-weighted nearest correlation matrix problem" %%%%
%%%%%%%%%%%%%%%%%%%%%%%  By Houduo Qi and Defeng Sun,  March 2008  %%%%%%%%%%%%%
%   Parameters:
%
%   Input
%   G       the given symmetric correlation matrix
%   H       the weight matrix for G
%   e       the right hand side of equality constraint
%   I_e     row indices of the fixed elements
%   J_e     column indices of the fixed elements
%   l       the right hand side of lower bounded constraint
%   I_l     row indices of the lower bounded elements
%   J_l     column indices of the lower bounded elements
%   u       the right hand side of upper bounded constraint
%   I_u     row indices of the upper bounded elements
%   J_u     column indices of the upper bounded elements
%   tau     the lower bound for the smallest eigevalue of X0 (can be zero)
%   TOL1    stopping crierion for the KKT system: 1.0e-4 ~ 1.0e-6
%
%   Output
%   X0:       the calibrated correlation matrix
%   z_e0:     Lagrangian dual variable corresponds to X_ii  = e_ij
%   z_l0:     Lagrangian dual variable corresponds to X_ij >= l_ij
%   z_u0:     Lagrangian dual variable corresponds to X_ij <= u_ij
%   Gamma0:   Lagrangian dual variable corresponds to X >= 0 (psd)
%  val_obj:   final objective function value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Last updated on September 12, 2009 %%%%%%%%%%%%%%%%%




t0 = clock;
n = length(G);
k_e = length(e);       % the number of the fixed elements
k_l = length(l);       % the number of the elements with lower bounds
k_u = length(u);       % the number of the elements with upper bounds
k = k_e+k_l+k_u;

fprintf('\n *****************************************************  \n')
fprintf( '  The  Augmented Lagrangian Dual Method for the H-norm Case       ')
fprintf('\n ******************************************************  \n')

fprintf('\n The information of this problem is as follows: \n')
fprintf(' Dim. of    sdp      constr  = %d \n',n)
fprintf(' Num. of equality    constr  = %d \n',k_e)
fprintf(' Num. of lower bound constr  = %d \n',k_l)
fprintf(' Num. of upper bound constr  = %d \n',k_u)
fprintf(' The lower bounds: [ %2.1e, %2.1e ] \n',min(l),max(l))
fprintf(' The upper bounds: [ %2.1e, %2.1e ] \n',min(u),max(u))

G = (G + G')/2;      % make sure that G is symmetric
H = (H + H')/2;      % make sure that H is symmetric

H_ave = sum(sum(H))/n^2;
H = H/H_ave;                 
H2H = H.*H;

%% set parameters
max_lambda = 1.0e5;                     % maximum penalty parameter
lambda = 10;                            % 0.5; or 1.0e2;   
lambda = min(lambda, max_lambda);       % the initial penalty parameter
rho = 1.4;                              % 10; the ration to increase the penalty parameter

mu = 1.0e-12;                           % parameter used in the line search (it can be other small positive numbers)
rho1 = 0.0;                             % rho1>=0 controls the regularized quadratic term
tol =  1.0e-1;                          % 1.0e-2; tolerance for the CG method 1.0e-1 ~5.0e-3 (relative error)

%TOLRel1 = TOL1;
TOL2 =  5.0e-1*TOL1;                     % tolerance of || DL ||=0, slightly smaller than outside tolerance
TOLrel2 = 1.0e1;                        % initial tolerance of || DL ||=0
 
maxit =  200;         % 200; maximum step the CG method
maxit1 = 200;         % maximum steps for k: maximal number of outer iterations
maxit2 = 50;          % maximum steps for j: maximal number of inner iterations 
maxit3 = 20;          % maximum steps for t: maximal number of line searches

eig_time = 0;
pcg_time = 0;
eye_n = eye(n);

% initial value
X0 = G;
z_e0 = zeros(k_e,1);
z_l0 = zeros(k_l,1);
z_u0 = zeros(k_u,1);
Gamma0 = zeros(n,n);

x0z_e = z_e0;
x0z_l = z_l0;
x0z_u = z_u0;

    Z_e0 = zeros(n,n);
    for i = 1:k_e
         Z_e0(I_e(i),J_e(i)) = z_e0(i);
    end
    Z_e0 = 0.5*(Z_e0 + Z_e0');
   
    Z_l0 = zeros(n,n);
    for i = 1:k_l
         Z_l0(I_l(i),J_l(i)) = z_l0(i);
    end
    Z_l0 = 0.5*(Z_l0 + Z_l0');

    Z_u0 = zeros(n,n);
    for i = 1:k_u
         Z_u0(I_u(i),J_u(i)) = z_l0(i);
    end
    Z_u0 = 0.5*(Z_u0 + Z_u0');

 
  
    
Gamma = Gamma0 - lambda*(X0 - tau*eye_n);
Gamma = (Gamma + Gamma')/2;

eig_t0 = clock;
[P,D] =  eig(Gamma);
eig_time =  eig_time + etime(clock,eig_t0);
d = real(diag(D));
P = real(P);
if d(1) < d(n)
    d = d(n:-1:1);
    P = P(:,n:-1:1);
end


        Ip = find(d>0);
         r = length(Ip);

         if (r==0)
              Gamma = zeros(n,n);
         elseif (r==n)
              Gamma =  Gamma;
         elseif (r<=n/2)
             d1 = d(Ip);

             d1 = d1.^0.5;
             P1 = P(:,Ip);
             if r >1
                 P1 = P1*sparse(diag(d1));
                  Gamma = P1*P1'; %
             else
                  Gamma = d1^2*P1*P1';
             end
         else
             d2 = -d(r+1:n);
             d2 = d2.^0.5;
             P2 = P(:,r+1:n);
             P2 = P2*sparse(diag(d2));
             Gamma = Gamma + P2*P2'; % 
         end

          Gamma = (Gamma + Gamma')/2;





 

Err_est0 = mytrace((Gamma0 - Gamma), (Gamma0 - Gamma));
Err_est0 = Err_est0^0.5/(min(100, lambda^0.5));      

Gamma = H2H.*(X0 - G);                                 % a temporary matrix 
val_obj = mytrace(Gamma, X0 - G)/2;                    % objective function value at initial point

Gamma = H2H.*(X0 - G) - Z_e0 - Z_l0 + Z_u0 - Gamma0;   % the rhs side of the first KKT equation
Err_est1 = mytrace(Gamma, Gamma);
Err_est = max(Err_est1^0.5, Err_est0);

i = 1;
while (i <= k_e)
    x0z_e(i) = X0(I_e(i), J_e(i));
    i =i+1;
end 
Err_est2 = (e - x0z_e)'*(e - x0z_e);
Err_est = max(Err_est2^0.5, Err_est); 

i=1;
while (i<=k_l)
    x0z_l(i) = X0(I_l(i), J_l(i));
    i=i+1;
end
z_l  = x0z_l - l;
z_l  = z_l0 - max(0, z_l0 - z_l);
Err_est3 = norm(z_l);
Err_est = max(Err_est3, Err_est); 

i=1;
while (i<=k_u)
    x0z_u(i)= X0(I_u(i), J_u(i));
    i=i+1;
end
z_u  = u - x0z_u;
z_u  = z_u0 - max(0, z_u0 - z_u);
Err_est4 = norm(z_u);
Total_err = max(Err_est4, Err_est);     

Total_err0 = Total_err;  

fprintf('\n')
fprintf('The initial absolute error: %d \n', Total_err)
 
 if Total_err > TOL1
     
     fprintf('The initial correlation test is not good enough, continue...... \n');
     
     k = 0;
     count_LS = 0;     % the number of linear systems to be solved
     count_cg = 0;     % the number of CG steps to be needed
     
     [f0,Fx,d_l,d_u] = grad(G,H,X0,z_e0,z_l0,z_u0,Gamma0,Gamma,e,I_e,J_e,l,I_l,J_l,u,I_u,J_u,lambda,d,P,X0,rho1);  % the gradient of the augmented Lagrangian function
     b = -Fx;
     f_eval = 1;
     
     while ( k < maxit1 )
         
         j = 0;
         Tol2_k = min(TOLrel2, 0.5*Total_err);
         Tol2_k = max(TOL2, Tol2_k);

         norm_Lag = mytrace(b, b);
         norm_Lag = norm_Lag^0.5;
         
         fprintf('\n ********** One Level Augmented Lagrangian Function Minimization ***************')
         fprintf('\n  Inner Iter.   Num. of CGs      Step length        Norm of Gradient')
         
         norm_Lag = 1.0e18; % Must do at least one step
         
         while (norm_Lag > Tol2_k) & (j<maxit2)

              
            Omega12 = omega_mat(P, d); % generate the Omeage matrix (the off-diagonal part only)
            M = precond_matrix(Omega12, P); 
            pcg_t0 = clock;
            [dx,flag,relres,iterk] = p_cg(b,tol,maxit,H,M,I_e,J_e,I_l,J_l,I_u,J_u,lambda,Omega12,P,d_l,d_u,rho1);
            %dx=b; %corresponds to the steepest descent direction method
            pcg_time = pcg_time + etime(clock,  pcg_t0);    
            count_cg = count_cg + iterk;
            count_LS = count_LS + 1;         % the total number of linear system solved
             
            tmp =  mu*mytrace(Fx,dx);

             alpha =1.0;
             Ok = false;
             line_num = 0;
             
             while ( Ok == false ) & ( line_num < maxit3 )
                 % update X
                 delta_x = alpha*dx;

                 X = X0 + delta_x;
                 X = (X + X')/2;
                 Gamma = Gamma0 - lambda*( X-tau*eye_n );
                 Gamma = (Gamma + Gamma')/2;
                 eig_t0 = clock;
                 [P,D] =  eig(Gamma);
                 eig_time =  eig_time + etime(clock,eig_t0);
                 P = real(P);
                 d = real(diag(D));
                 
                 if d(1) < d(n)
                     d = d(n:-1:1);
                     P = P(:,n:-1:1);
                 end
                 
                 
                 
                 
                 
                 [f,Fx,d_l,d_u] = grad(G,H,X,z_e0,z_l0,z_u0,Gamma0,Gamma,e,I_e,J_e,l,I_l,J_l,u,I_u,J_u,lambda,d,P,X0,rho1);
  
                 tmp1 = alpha*tmp;
                 dlag = f - f0 - tmp1;

                 if  ( dlag <1e-6 )           % dlag <= 0 (theoretically)
                     Ok = true;
                 else
                     alpha = alpha/2;
                 end
                 
                 line_num = line_num + 1;
             end   % End of line search 
             
             f_eval = f_eval + line_num;
             
             j = j+1;  
             f0 = f;
             X0 = X;
             % compute the first derivative of Lagrangian
             b = -Fx;
             norm_Lag = mytrace(b,b);
             norm_Lag = norm_Lag^0.5;
             
             fprintf('\n   %2.0d              %2.0d            %2.1e            %3.2e',j,iterk,alpha,norm_Lag);
                  
         end   % j while loop: this is the end of one level Augmented Lagrangian function minimization
            
         k=k+1;   
         
         %%%%%% Update the dual variables and the penalty parameter 
         %%% Checking if convergence has reached 
         i = 1;
         while (i<=k_e)
             x0z_e(i) = X(I_e(i),J_e(i));
             i =i+1;
         end
         z_e = z_e0 + lambda*(e - x0z_e);
         
         i=1;
         while (i<=k_l)
             x0z_l(i) = X(I_l(i),J_l(i));
             i=i+1;
         end
         z_l = z_l0 - lambda*(x0z_l-l);
         z_l = max(0,z_l);

         i=1;
         while (i<=k_u)
             x0z_u(i) = X(I_u(i),J_u(i));
             i=i+1;
         end
         z_u = z_u0 - lambda*(u - x0z_u);
         z_u = max(0, z_u);

         
       
         Ip = find(d>0);
         r = length(Ip);

         if (r==0)
              Gamma = zeros(n,n);
         elseif (r==n)
              Gamma =  Gamma;
         elseif (r<=n/2)
             d1 = d(Ip);

             d1 = d1.^0.5;
             P1 = P(:,Ip);
             if r >1
                 P1 = P1*sparse(diag(d1));
                  Gamma = P1*P1'; %
             else
                  Gamma = d1^2*P1*P1';
             end
         else
             d2 = -d(r+1:n);
             d2 = d2.^0.5;
             P2 = P(:,r+1:n);
             P2 = P2*sparse(diag(d2));
             Gamma = Gamma + P2*P2'; % 
         end

          Gamma = (Gamma + Gamma')/2;

         
         
         
         
         

         Err_est1 = mytrace((Gamma0 - Gamma), (Gamma0 - Gamma));
         Err_est  = Err_est1^0.5/(min(100,lambda^0.5));      % the test on the complementary condition
                 
         Err_est2 = (e - x0z_e)'*(e - x0z_e);       
         Err_est = max(Err_est2^0.5, Err_est);
         
         Err_est3 = norm(z_l0 - z_l);
         Err_est3 = Err_est3/(min(100, lambda^0.5)) ;
         Err_est = max(Err_est3, Err_est);

         Err_est4 = norm(z_u0 - z_u);
         Err_est4 = Err_est4/(min(100, lambda^0.5)) ;
         Err_est = max(Err_est4, Err_est);

         Total_err = max(norm_Lag, Err_est);

         if Total_err <= TOL1
             z_e0 = z_e;
             z_l0 = z_l;
             z_u0 = z_u;

             X0 = X0 + (Gamma - Gamma0)/lambda;    % to get a positive semidefinite matrix
             X0 = (X0 + X0')/2;
             Gamma0 = Gamma;
             
             i=1;
             while (i<=k_l)
                 x0z_l(i)= X(I_l(i), J_l(i));
                 i=i+1;
             end

             i=1;
             while (i<=k_u)
                 x0z_u(i)= X(I_u(i), J_u(i));
                 i=i+1;
             end
             
             gap = mytrace(X0 - tau*eye_n, Gamma0);
             
             break;   % successful already
         else

             klambda =0;
             if (lambda < max_lambda)
                 if Total_err >(0.25)*Total_err0      % update lambda only if convergence is not fast
                     lambda = rho*lambda;              % new lambda
                     lambda = min(lambda, max_lambda);
                     klambda = 1;                 
                 end                
             end
             
             Total_err0 = Total_err;

             if klambda == 1  % Only when lambda is updated we need recompute Gamma

                 Gamma = Gamma0 - lambda*(X -tau*eye_n);
                 Gamma = (Gamma + Gamma')/2;

                 eig_t0 = clock;
                 [P,D] =  eig(Gamma);
                 eig_time =  eig_time + etime(clock, eig_t0);
                 P =real(P);
                 d = real(diag(D));
                 if d(1) < d(n)
                     d = d(n:-1:1);
                     P = P(:,n:-1:1);
                 end
                
                 
                 
                 Ip = find(d>0);
                 r = length(Ip);

                 if (r==0)
                     Gamma0 = zeros(n,n);
                 elseif (r==n)
                     Gamma0 =  Gamma;
                 elseif (r<=n/2)
                     d1 = d(Ip);

                     d1 = d1.^0.5;
                     P1 = P(:,Ip);
                     if r >1
                         P1 = P1*sparse(diag(d1));
                         Gamma0 = P1*P1'; %
                     else
                         Gamma0 = d1^2*P1*P1';
                     end
                 else
                     d2 = -d(r+1:n);
                     d2 = d2.^0.5;
                     P2 = P(:,r+1:n);
                     P2 = P2*sparse(diag(d2));
                     Gamma0 = Gamma + P2*P2'; %
                 end

                 Gamma0 = (Gamma0 + Gamma0')/2;

         

                 z_e0 = z_e0 + lambda*(e - x0z_e);   % new z_e0

                 z_l0 = z_l0 - lambda*(x0z_l - l);
                 z_l0 = max(0,z_l0);

                 z_u0 = z_u0 - lambda*(u - x0z_u);
                 z_u0 = max(0, z_u0);
                 
             else
                 Gamma0 = Gamma;
                 z_e0 = z_e;
                 z_l0 = z_l;
                 z_u0 = z_u;
             end
             
             Gamma = Gamma0 - lambda*(X0 - tau*eye_n);
             Gamma = (Gamma + Gamma')/2;
             eig_t0 = clock;
             [P,D] = eig(Gamma);
             eig_time =  eig_time + etime(clock,eig_t0);
             P = real(P);
             d = real(diag(D));
             if d(1) < d(n)
                     d = d(n:-1:1);
                     P = P(:,n:-1:1);
             end
                 
             [f0,Fx,d_l,d_u] = grad(G,H,X0,z_e0,z_l0,z_u0,Gamma0,Gamma,e,I_e,J_e,l,I_l,J_l,u,I_u,J_u,lambda,d,P,X0,rho1);
             b = -Fx;
             
        end % Total_err loop     
        
         fprintf('\n')
         tt = etime(clock, t0);
         [hh,mm,ss] = time(tt);              
         fprintf('\n   Level of AugLag Min.       Total Absolute Error           Time_used \n ')
         fprintf('        %2.0d                         %2.1e                  %d:%d:%d ',k,Total_err,hh,mm,ss)
         fprintf('\n ---------- End of One Level Augmented Lagrangian Function Optimization ---------- \n')
         fprintf('\n')
         
     end  % End of while
     
         tt = etime(clock, t0);
         [hh,mm,ss] = time(tt);              
         fprintf('\n   Level of AugLag Min.       Total Absolute Error           Time_used \n ')
         fprintf('        %2.0d                         %2.1e                  %d:%d:%d ',k,Total_err,hh,mm,ss)
         fprintf('\n ---------- End of One Level Augmented Lagrangian Function Optimization ---------- \n')
     
 else
     
     fprintf(' The ininial correlation test succeeds already: stop! \n');
     
 end   % End of if

 Feas_err1 = (e - x0z_e)'*(e - x0z_e);
 Feas_err1 = Feas_err1^0.5;
 RelFeas_err1 = Feas_err1/(1+ norm(e));
 Feas_err2 = norm(min(0, x0z_l-l));
 RelFeas_err2 = Feas_err2/(1+norm(l));
 Feas_err3 = norm(min(0, u- x0z_u));
 RelFeas_err3 = Feas_err3/(1+norm(u));

    Z_e0 = zeros(n,n);
    for i = 1:k_e
         Z_e0(I_e(i),J_e(i)) = z_e0(i);
    end
    Z_e0 = 0.5*(Z_e0 + Z_e0');
   
    Z_l0 = zeros(n,n);
    for i = 1:k_l
         Z_l0(I_l(i),J_l(i)) = z_l0(i);
    end
    Z_l0 = 0.5*(Z_l0 + Z_l0');

    Z_u0 = zeros(n,n);
    for i = 1:k_u
         Z_u0(I_u(i),J_u(i)) = z_u0(i);
    end
    Z_u0 = 0.5*(Z_u0 + Z_u0');

 Gamma = H.*(X0 - G);
 val_obj = mytrace(Gamma, Gamma)/2;
 val_obj = H_ave^2*val_obj;                       % scale back to the original function value
 
Gamma = H.*Gamma - Z_e0 - Z_l0 + Z_u0 - Gamma0;   % the residue of the first equation in the KKT system
KKT_err = mytrace(Gamma, Gamma);
KKT_err = KKT_err^0.5;                            % the norm of the residue of the first KKT equation 

H2G = H2H.*G;
scale_kkt = mytrace(H2G, H2G);
scale_kkt = scale_kkt^0.5;
RelKKT_err = KKT_err/(1+scale_kkt);    % the norm of the relative residue of the first KKT equation 
gap = gap*H_ave^2;
Relgap = gap/(1+abs(val_obj));         % relative gap

%%% The Lagrangian multipliers to the original problem %%%%%%
z_e0 = H_ave^2*z_e0;
z_l0 = H_ave^2*z_l0;
z_u0 = H_ave^2*z_u0;
Gamma0 = H_ave^2*Gamma0;

%%% compute the rank of X0 and Gamma0
r_G = length( find(d>1.0e-8) );

eig_t0 = clock;
[P,D] =  eig(X0 - tau*eye_n);
eig_time =  eig_time + etime(clock, eig_t0);
 
r_X0 = length( find(real(diag(D))>1.0e-8) );


time_used = etime(clock,t0);

fprintf('\n')
fprintf('\n ================ Final Information ================= \n');
fprintf(' Total number of iterations      = %2.0f \n',k)
fprintf(' Number of linear systems        = %2.0f \n',count_LS)
fprintf(' Total number of CG steps        = %2.0f \n',count_cg)
fprintf(' Number of func. evaluations     = %2.0f \n',f_eval)
fprintf(' Final penalty parameter         = %2.0f \n',lambda)

fprintf(' Primal objective value            = %3.2e \n', val_obj)

fprintf(' Norm of KKT equation                = %3.2e \n', KKT_err)
fprintf(' Relative norm of  KKT equation      = %3.2e \n', RelKKT_err)
fprintf(' Equation feasibility                = %3.2e \n', Feas_err1)
fprintf(' Relative equation feasibility       = %3.2e \n', RelFeas_err1)
fprintf(' Lower bound feasibility             = %3.2e \n', Feas_err2)
fprintf(' Relative lower bound feasibility    = %3.2e \n', RelFeas_err2)
fprintf(' Upper bound feasibility             = %3.2e \n', Feas_err3)
fprintf(' Relative upper bound feasibility    = %3.2e \n', RelFeas_err3)
fprintf(' Gap < X-tau*I, Gamma0 >             = %3.2e \n', abs(gap))
fprintf(' Relative gap                        = %3.2e \n', abs(Relgap))

fprintf(' Rank of X0-(tau0*I)            = %2.0d \n', r_X0)
fprintf(' Rank of Gamma0                 = %2.0d \n', r_G)
fprintf(' Computing time for CG          = %3.1f \n',pcg_time)
fprintf(' Computing time for eigen-decom = %3.1f \n',eig_time)
fprintf(' Total computing time (secs)    = %3.1f \n',time_used)
fprintf(' ====================================================== \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% end of the main program %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%  Subroutines  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Time format function %%%
function [h,m,s] = time(t)
t = round(t); 
h = floor(t/3600);
m = floor(rem(t,3600)/60);
s = rem(rem(t,60),60);
return
%%% End of time.m



%%% Trace function %%%
function val = mytrace(A,B)
if (nargin == 1)
   val = sum(diag(A));
elseif (nargin == 2)
   val = sum(sum(A.*B));
end
return
%%% end of mytrace.m %%%

%% To generate the gradient of the augmented lagrangian function %%%
function  [f,Fx,d_l,d_u] = grad_old(G,H,X,z_e,z_l,z_u,Gamma,e,I_e,J_e,l,I_l,J_l,u,I_u,J_u,lambda,d,P,X0,rho1)
%%%% output
% f:   the augmented Lagrangian function value
% Fx:  gradient of the augmented Lagrangian function 
% d_l: the derivative at (z_l - lambda*(xz_l-l))_+
% d_u: the derivative at (z_u - lambda*(u-xz_u))_+

n = length(X);
k_e = length(z_e);
k_l = length(z_l);
k_u = length(z_u);

f = 0.0;
Fx = zeros(n,n);
d_l = ones(k_l,1);
d_u = ones(k_u,1);

xz_e = z_e;
x0z_e = xz_e;
xz_l = z_l;
xz_u = z_u;

i=1;
while (i<=k_e)
    x0z_e(i) = X(I_e(i),J_e(i));
    i=i+1;
end
xz_e = z_e + lambda*(e-x0z_e);
Z_e = sparse(I_e,J_e,xz_e,n,n);
Z_e = 0.5*(Z_e+Z_e');

i=1;
while (i<=k_l)
    xz_l(i) = X(I_l(i),J_l(i));
    xz_l(i) = z_l(i) - lambda*( xz_l(i)-l(i) );
    if xz_l(i) < 0
        d_l(i) = 0;
        xz_l(i) = 0;
    end
    i=i+1;
end 
Z_l = sparse(I_l,J_l,xz_l,n,n);
Z_l = 0.5*(Z_l+Z_l');
   

i=1;
while (i<=k_u)
    xz_u(i) = X(I_u(i),J_u(i));
    xz_u(i) = z_u(i) - lambda*( u(i)-xz_u(i) );
    if xz_u(i)<0
        xz_u(i) = 0;
        d_u(i) = 0;
    end
    i=i+1;
end  
Z_u = sparse(I_u,J_u,xz_u,n,n);
Z_u = 0.5*(Z_u+Z_u');

W = P';
for i=1:n 
     W(i,:) = max(d(i),0)*W(i,:);             % W = P*diag(max(d,0))*P^T
end
W = P*W; 
 
Fx = H.*(X-G);
f = f + 0.5*mytrace(Fx,Fx);
f = f + z_e'*(e-x0z_e) + ((e-x0z_e)'*(e-x0z_e))*lambda/2;
f = f + (xz_l'*xz_l - z_l'*z_l)/(2*lambda);
f = f + (xz_u'*xz_u - z_u'*z_u)/(2*lambda);
f = f + (mytrace(W,W)-mytrace(Gamma,Gamma))/(2*lambda);
f = f + mytrace(X-X0,X-X0)*rho1/(2*lambda);     % regularized term

Fx = H.*Fx;
Fx = Fx - Z_e - Z_l + Z_u - W;
Fx = Fx + (X-X0)*(rho1/lambda);
return
%%% end of grad.m %%%



%%% To generate the gradient of the augmented lagrangian function %%%
function  [f,Fx,d_l,d_u] = grad(G,H,X,z_e,z_l,z_u,Gamma0,Gamma,e,I_e,J_e,l,I_l,J_l,u,I_u,J_u,lambda,d,P,X0,rho1)
%%%% output
% f:   the augmented Lagrangian function value
% Fx:  gradient of the augmented Lagrangian function 
% d_l: the derivative at (z_l - lambda*(xz_l-l))_+
% d_u: the derivative at (z_u - lambda*(u-xz_u))_+

n = length(X);
k_e = length(z_e);
k_l = length(z_l);
k_u = length(z_u);

f = 0.0;
 Fx = zeros(n,n);
d_l = ones(k_l,1);
d_u = ones(k_u,1);

xz_e = z_e;
x0z_e = xz_e;
xz_l = z_l;
xz_u = z_u;

i=1;
while (i<=k_e)
    x0z_e(i) = X(I_e(i),J_e(i));
    i=i+1;
end
xz_e = z_e + lambda*(e-x0z_e);


    Z_e = zeros(n,n);
    for i = 1:k_e
         Z_e(I_e(i),J_e(i)) = xz_e(i);
    end
    Z_e = 0.5*(Z_e + Z_e');

i=1;
while (i<=k_l)
    xz_l(i) = X(I_l(i),J_l(i));
    xz_l(i) = z_l(i) - lambda*( xz_l(i)-l(i) );
    if xz_l(i) < 0
        d_l(i) = 0;
        xz_l(i) = 0;
    end
    i=i+1;
end 
    Z_l  = zeros(n,n);
    for  i = 1:k_l
         Z_l(I_l(i),J_l(i)) = xz_l(i);
    end
    Z_l  = 0.5*(Z_l + Z_l');

i=1;
while (i<=k_u)
    xz_u(i) = X(I_u(i),J_u(i));
    xz_u(i) = z_u(i) - lambda*( u(i)-xz_u(i) );
    if xz_u(i)<0
        xz_u(i) = 0;
        d_u(i) = 0;
    end
    i=i+1;
end  
    Z_u  = zeros(n,n);
    for i = 1:k_u
         Z_u (I_u(i),J_u(i)) = xz_u (i);
    end
    Z_u = 0.5*(Z_u  + Z_u');


Ip = find(d>0);
    r = length(Ip);
    
    if (r==0)
        W =zeros(n,n);
    elseif (r==n)
        W = Gamma;
    elseif (r<=n/2)
        d1 = d(Ip);
        d1 = d1.^0.5;
        P1 = P(:, 1:r);
        if r >1
            P1 = P1*sparse(diag(d1));
             W = P1*P1'; %
        else
            W = d1^2*P1*P1';
        end
       
    else %% r >n/2
        d2 = -d(r+1:n);
        d2 = d2.^0.5;
        P2 = P(:, r+1:n);
        P2 = P2*sparse(diag(d2));
         W = Gamma + P2*P2'; %
    end
    
    W = (W + W')/2;
 
 
Fx = H.*(X - G);
f = f + 0.5*mytrace(Fx, Fx);
f = f + z_e'*(e - x0z_e) + ((e - x0z_e)'*(e - x0z_e))*lambda/2;
f = f + (xz_l'*xz_l - z_l'*z_l)/(2*lambda);
f = f + (xz_u'*xz_u - z_u'*z_u)/(2*lambda);
f = f + (mytrace(W,W) - mytrace(Gamma0,Gamma0))/(2*lambda);
f = f + mytrace(X-X0,X-X0)*rho1/(2*lambda);     % regularized term

Fx = H.*Fx;
Fx = Fx - Z_e - Z_l + Z_u - W;
Fx = Fx + (X-X0)*(rho1/lambda);
return
%%% end of grad.m %%%

%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the essential part of the first -order difference of d
%%%%%%%
function Omega12 = omega_mat(P,lambda)
%We compute omega only for 1<=|idx|<=n-1
n = length(lambda);
idx.idp = find(lambda>0);
idx.idm = setdiff([1:n],idx.idp);
 
r = length(idx.idp);
 
if ~isempty(idx.idp)
    if (r == n)
        Omega12 = ones(n,n);
    else
        s = n-r;
        dp = lambda(1:r);
        dn = lambda(r+1:n);
        Omega12 = (dp*ones(1,s))./(abs(dp)*ones(1,s) + ones(r,1)*abs(dn'));
        %  Omega12 = max(1e-15,Omega12);

    end
else
    Omega12 =[];
end

    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of omega_mat.m %%%%%%%%%%


%%%%%% To generate the first -order difference of $\Pi_{S^n_+}( )$ at d
function [Omega12,P1,P2] = omega_mat_old(P,d)

idx.idp = find(d>0);
idx.idm = find(d<=0);
n = length(d);
r = length(idx.idp);

if isempty(idx.idp)
    Omega12 = [];
elseif (r == n)
    Omega12 = ones(n,n);
else
    s = n-r;
    if idx.idp(r)< idx.idm(1) % to know largest eigenvalue comes first or not
    dp = d(1:r);
    dn = d(r+1:n); 
    else
    dp = d(s+1:n);
    dn = d(1:s);
    end
    Omega12 = (dp*ones(1,s))./(abs(dp)*ones(1,s) + ones(r,1)*abs(dn'));
end

% ***** perturbation *****
% Omega12 = max(1e-15,Omega12);
if ~isempty(idx.idp)
    P1 = P(:,idx.idp); P2 = P(:,setdiff([1:n],idx.idp));
else
    P1 = []; P2 = P;
end
return
%%%% end of omega_mat.m %%%%%%%%%%



%%%%%% PCG method %%%%%%%
%%%%%%% This is exactly the algorithm given by Hestenes and Stiefel (1952)
% An iterative method to solve A(x) =b  
% The symmetric positive definite matrix M is a preconditioner for A.
%%%%%% See Pages 527 and 534 of Golub and va Loan (1996)
function [p,flag,relres,iterk] = p_cg(b,tol,maxit,H,M,I_e,J_e,I_l,J_l,I_u,J_u,lambda,Omega12,P,d_l,d_u,rho1)                                   
% Initializations
n = length(b);      % the dimension of the unknown 
k_e = length(I_e);
k_l = length(I_l);
k_u = length(I_u);

r = b;              % We take the initial guess x0=0 to save time in calculating A(x0) 
n2b = mytrace(b,b);  
n2b = n2b^0.5;      % norm of b
 tolb = max(tol,min(0.5,n2b))*n2b; 
 %tolb = tol * n2b;  % relative tolerance 

p = zeros(n,n);
 
flag= 1;
iterk = 0;
relres= 1000;   % To give a big value on relres

% Precondition 


 

z_e = ones(k_e,1);
z_l = ones(k_l,1);
z_u = ones(k_u,1);

 
z_l = d_l.*z_l;
z_u = d_u.*z_u;

            Z = zeros(n,n);
            for i=1:k_e
                Z(I_e(i), J_e(i)) = z_e(i);
            end
            for i=1:k_l
                Z(I_l(i), J_l(i)) = z_l(i) + Z(I_l(i), J_l(i));
            end
            for i=1:k_u
                Z(I_u(i), J_u(i)) = z_u(i) + Z(I_u(i), J_u(i));  %%% upper bound
            end
            Z = 0.5*(Z + Z');
   

M =  (H.*H) + lambda* (Z + M)  + (rho1/lambda); 

z = r./M;  %%%%% z = M\r; here M =diag(ww); if M is not the identity matrix 

rz1 = mytrace(r,z); 
rz2 = 1; 
dX = z;
% CG iteration
for k = 1:maxit
   if k > 1
       beta = rz1/rz2;
       dX = z + beta*dX;
   end
   
   %diag_d =diag(dX);
   
   w = Jacobian_mat(dX,H,I_e,J_e,I_l,J_l,I_u,J_u,lambda,Omega12,P,d_l,d_u,rho1);
   w = w + min(1.0e-2, 1.0e1*n2b)*dX; %w = A(d); Perturbed
  
  
   denom = mytrace(dX,w);
   iterk = k;
   norm_r = mytrace(r,r);
   norm_r = norm_r^0.5;
   relres = norm_r/n2b;              %relative residue =norm(z) / norm(b)
   if denom <= 0 
       norm_d = mytrace(dX,dX);
       norm_d =norm_d^0.5;
       p = dX/norm_d; % d is not a descent direction
       break % exit
   else
       alpha = rz1/denom;
       p = p + alpha*dX;
       r = r - alpha*w;
   end
   z = r./M; %  z = M\r if M is not the identity matrix;
   
   norm_r = mytrace(r, r);
   norm_r = norm_r^0.5;
   if norm_r <= tolb % Exit if Hp=b solved within the relative tolerance
       iterk =k;
       relres = norm_r/n2b;          %relative residue =norm(z) / norm(b)
       flag =0;
       break
   end
   rz2 = rz1;
   rz1 = mytrace(r,z);
end
return
%%% end of pre_cg.m %%%%%%%%%%%



%%%%%% To generate the Jacobain product with dX for the Augmented Lagrangian Function %%%%%%%
function JFX = Jacobian_mat(dX,H,I_e,J_e,I_l,J_l,I_u,J_u,lambda,Omega12,P,d_l,d_u,rho1)

n = length(dX);
k_e = length(I_e);
k_l = length(I_l);
k_u = length(I_u);

JFX = zeros(n,n);

z_e = zeros(k_e,1);
z_l = zeros(k_l,1);
z_u = zeros(k_u,1);

i = 1;
while (i<=k_e)
    z_e(i) = dX(I_e(i),J_e(i));
    i =i+1;
end 

i =1;
while (i<=k_l)
    z_l(i) = dX(I_l(i),J_l(i));
    i =i+1;
end 
z_l = d_l.*z_l;


i =1;
while (i<=k_u)
    z_u(i) = dX(I_u(i),J_u(i));
    i =i+1;
end 
z_u = d_u.*z_u;

            Z = zeros(n,n);
            for i=1:k_e
                Z(I_e(i), J_e(i)) = z_e(i);
            end
            for i=1:k_l
                Z(I_l(i), J_l(i)) = z_l(i) + Z(I_l(i), J_l(i));
            end
            for i=1:k_u
                Z(I_u(i), J_u(i)) = z_u(i) + Z(I_u(i), J_u(i));  %%% upper bound
            end
            Z = 0.5*(Z + Z');
   
JFX =  (H.*H).*dX + lambda*Z + (rho1/lambda)*dX; 

JFX = JFX + lambda*Jacobian_mat0(dX,Omega12,P);
return
%%% end of Jacobian_mat.m %%%






%%%%%% To generate the Jacobain product with dX for the Augmented Lagrangian Function %%%%%%%
function JFX = Jacobian_mat_old(dX,H,I_e,J_e,I_l,J_l,I_u,J_u,lambda,Omega12,P1,P2,d_l,d_u,rho1)

n = length(dX);
k_e = length(I_e);
k_l = length(I_l);
k_u = length(I_u);

JFX = zeros(n,n);

z_e = zeros(k_e,1);
z_l = zeros(k_l,1);
z_u = zeros(k_u,1);

i = 1;
while (i<=k_e)
    z_e(i) = dX(I_e(i),J_e(i));
    i =i+1;
end 
Z_e = sparse(I_e,J_e,z_e,n,n);
Z_e = 0.5*(Z_e+Z_e');

i =1;
while (i<=k_l)
    z_l(i) = dX(I_l(i),J_l(i));
    i =i+1;
end 
z_l = d_l.*z_l;
Z_l = sparse(I_l,J_l,z_l,n,n);
Z_l = 0.5*(Z_l+Z_l');
   

i =1;
while (i<=k_u)
    z_u(i) = dX(I_u(i),J_u(i));
    i =i+1;
end 
z_u = d_u.*z_u;
Z_u = sparse(I_u,J_u,z_u,n,n);
Z_u = 0.5*(Z_u+Z_u');
   
JFX =  (H.*H).*dX + lambda*Z_e + lambda*Z_l + lambda*Z_u + (rho1/lambda)*dX; 

JFX = JFX + lambda*Jacobian_mat0(dX,Omega12,P1,P2);
return
%%% end of Jacobian_mat.m %%%


%%%%%% To generate the Jacobain product with x: F'(y)(x) %%%%%%%
%%%%%%%

function JFX = Jacobian_mat0(X,Omega12,P)
n = length(X);
JFX = zeros(n);


[r,s] = size(Omega12);

if (r>0)
     P1 = P(:,1:r);
     H1 = P1; 
   if (r< n/2)
         P2 = P(:,r+1:n);
         H1=X*P1;
         Omega12 = Omega12.*(H1'*P2);


        %H =[(H1'*P1)*P1'+ Omega12*P2';Omega12'*P1']; %%%  H= [Omega o (P^T*diag(x)*P)]*P^T
        %H =[(H1'*P1)*P1'+ Omega12*P2';Omega12'*P1']; %%%  H= [Omega o (P^T*diag(x)*P)]*P^T
        JFX = P1*((H1'*P1)*P1' + 2.0*Omega12*P2');
        JFX = (JFX + JFX')/2;
        
        
    else % r >n/2
        if (r==n)
            JFX = X;
        else
             
            P2 = P(:,r+1:n);
            H2 =  P2;
            
             
            H2 = X*P2;
    
         
            Omega12 = ones(r,s)-Omega12;
            Omega12 = Omega12.*(P1'*H2);

           % H =[Omega12* (P2)';Omega12'*(P1)'+ ( (P2)'*H2)* (P2)']; %%% Assign H*P' to H= [(ones(n,n)-Omega) o (P^T*diag(x)*P)]*P^T

            JFX = P2*(2.0*Omega12'*P1' + (P2'*H2)* (P2)');
            JFX = (JFX + JFX')/2;
            
            JFX = X - JFX;

        end
    end
end
return
 
%%%%%%%%%%%%%%%
%end of Jacobian_mat0.m%%%


%%%%%% To generate the Jacobain product: F'(X)(dX) (generalized  Jacobian) for the projector %%%%%%%
function JFX = Jacobian_mat0_old(dX,Omega12,P1,P2)

n = length(dX);
JFX = zeros(n);
 
if isempty(P2)
    JFX = dX;
else
    if ~isempty(P1)
        tmp0 = P1'*dX;
        tmp1 = (tmp0*P1)*P1';

        tmp2 = Omega12.*(tmp0*P2);
        tmp2 = P1*tmp2*P2';

        JFX = P1*tmp1 + tmp2 + tmp2';
        JFX = (JFX+JFX')/2; %check symmetry
    end
end
return
%%% end of Jacobian_mat0.m %%%

 


%%%%%% To generate the Jacobain product: F'(X)(dX) (generalized  Jacobian) for the projector %%%%%%%
function JFX = Jacobian_mat0_old_old(dX,Omega12,P1,P2)

n = length(dX);
JFX = zeros(n);
 
if isempty(P2)
    JFX = dX;
else
    if ~isempty(P1)
        tmp0 = P1'*dX;
        tmp1 = (tmp0*P1)*P1';

        tmp2 = Omega12.*(tmp0*P2);
        tmp2 = P1*tmp2*P2';

        JFX = P1*tmp1 + tmp2 + tmp2';
        JFX = (JFX+JFX')/2; %check symmetry
    end
end
return
%%% end of Jacobian_mat0.m %%%



%%% To generate the (approximate) diagonal preconditioner
function H0 = precond_matrix(Omega12,P) 
n     = length(P);
[r,s] = size(Omega12);
H0  = zeros(n,n);

H = P';
H = H.*H;
if (r<n)
    if (r>0)  
            HH1 = H(1:r,:);
            HH2 = H(r+1:n,:);

            if (r<n/2)
                H0 = HH1'*(Omega12*HH2);
                tmp = sum(HH1);
                H0 = H0 + H0'+ tmp'*tmp;
            else
                Omega12 = ones(r,s) - Omega12;
                H0 = HH2'*((Omega12)'*HH1);
                tmp  = sum(HH2);
                H0 = H0 + H0' + tmp'*tmp;
                tmp = sum(H);
                H0 = tmp'*tmp - H0;
            end       
    end  %End of second if
    
else % if r=n
   H0 = ones(n,n); % exact diagonal matrix
end  %End of the first if
return
%%% End of precond_matrix.m 






