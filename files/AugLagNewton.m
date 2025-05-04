
function [X0,val_obj] =AugLagNewton(C,m,n,tau,tau0,TOL1)
%%%%%%%%%%%%%%%%%%%                           %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This code is for the Band Correlation Stress Testing of the Augmenetd Lagrangian Dual approach described in
%%%%%%  Houduo Qi and  Defeng Sun, "Correlation Stress Testing for Value-at-Risk: An Unconstrained
%%%%%%  Convex Optimization Approach", Department of Mathematics, National
%%%%%%  University of Singapore, March 2007 
%%%%%%%%%%%%%%%%%%%%%%%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 
% for solving 
%%  min  0.5 || X -  C ||^2 
%%  s.t.  X_ii = 1, i=1, ..., n
%%        X_ij = C_ij, 1<= i < j <=m
%%        X- tau*I >=0 (symmetric and positive semi-definite)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%
% C: stressed correlation matrix (n by n)
% C(1:m,1:m) is fixed
% tau0: the threshold for positive definiteness of  C  
% tau: the lower bound for the smallest eigevalue of X0 
% TOL1: stopping crierion for the KKT system: 1.0e-5 ~ 1.0e-6
% X0: the calibrated correlation matrix
% val_obj: final objective function value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Last updated on March 20, 2007 %%%%%%%%%%%%%%%%%%%%%%%%%%%%


t0 =cputime;
if m >n-1
    disp(' ---m >=n --- ')
end
m =min(m,n-1);
m1 =n-m; % C(1:m,1:m) fixed
if m<1
    tau =tau0;
else
    tau =min(tau, 0.90*tau0);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs =ones(n,1);

if m>0
    rhs1=ones(m,1);
    C1 = C(1:m,1:m);

    [P1,D1] = eig(C1);
    d1 = diag(D1);
    d1 =real(d1);
    if min(d1) <tau0
        fprintf('\n')
        fprintf('The smallest eignevalue of the initial C(1:m,1:m) < %d \n', tau0)
        disp(' ---Start the Pre-processing--- ')
        [C1,y]=CorNewton1(C1,rhs1,tau0);
        C(1:m,1:m) = C1(1:m,1:m); % update C to make sure  C(1:m,1:m) is PD
        disp(' ---Pre-processing finished--- ')
    else
        fprintf('\n')
        fprintf('The smallest eignevalue of the initial matrix C(1:m,1:m) >= %d \n', tau0)
        disp(' ---No pre-processing needed--- ')
    end
end
l =n*m1 - (m1+1)*m1/2;
l =round(l); % the number of stressed entries
c = zeros(l,1); 
b =zeros(l,1);
x =  zeros(l,1);
x0 = zeros(l,1);

y0 = zeros(n,1);
A0 = eye(n);
h = ones(l,1);
%h = h;  % Corresponding to the H-norm case 
ww =  ones(l,1); %The M =diag(ww) diagonal preconditioner
w =  h.*h;


Gamma0 =zeros(n,n);
X0 = zeros(n,n);

%% set parameters
max_lambda = 1.0e5; % maximum penalty parameter
lambda = 1.0e2;     % penalty parameter
rho = 10;           % the ration to increase the penalty parameter
mu = 1.0e-4;        % parameter used in the line search


tol0 = 1e-12;      % tolerance for the nearest correlation approximation 
tol  = 1e-6;       % tolerance for the CG method

TOL2 = 1e-6;       % tolerance of || DL ||=0 
TOLRel1 =TOL1;
TOL2 = 1.0e-7;     % tolerance of || DL ||=0  
TOLrel2 = min(1.0e-3, TOL2*(1+max_lambda));      % relative tolerance of || DL ||=0  depebding on lambda
%TOL3 = 1e-8;      % tolerance of || dx ||


maxit =min(20000, max(10000, (l+1)));    % maximum step the CG method
maxit1 = 200;     % maximum steps for k:  maximal number of outer iterations
maxit2 = 100;      % maximum steps for j: maximal number of innerer iterations 
maxit3 = 40;      % maximum steps for t:  maxima number of line searches




%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial data (X, y, Lambda)
%%%%%%% We first solve the following problem to get an innitial guess %%%%%%%%%%%%%% 
%%%%%%%%  min 0.5*<X-C, X-C>
%%%%%%%   s.t. X_ii =1, i=1,2,...,n
%%%%%%%        X>=0 (symmetric and positive semi-definite) %%%%%%%%%%%%%%%
%%%%%%%%
%%%%%%  based on the algorithm  in %%%%%
%%%%%%  ``A Quadratically Convergent Newton Method for %%%
%%%%%%    Computing the Nearest Correlation Matrix %%%%%
%%%%%%%   By Houduo Qi and Defeng Sun  %%%%%%%%%%%%
%%%%%%%   SIAM J. Matrix Anal. Appl. 28 (2006)360--385.
%%%%%%%  
%%%%%% Last modified date:  December 26, 2006  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Diagonal Preconditioner is added %%%%%%%%%%%%%%%%%%%%%%%
%%%%%% The  input argument is the matrix C   %%%%%%%%
%%%%%% The outputs are the optimal primal and dual solutions %%%%%%%%

disp('Start the ininial correlation test approximation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X0,y0] = CorNewton1(C,rhs,tau0);
val_obj =sum(sum((X0-C).*(X0-C)))/2;

disp('The ininial correlation test approximation finished')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma0 = X0-C-diag(y0);  %% The guessed optimal Lagrangian multiplier matrix
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% use vector c to represent the stressed parts%%%%%%%%%%%%%%%%%%%%%%
for i=m+1:n
    j = m*(i-m-1)+(i-m-1)*(i-m-2)/2;
    j = round(j);
    c(j+1:j+i-1)= C(i,1:i-1)'; % the ith row of the lower triangular part of C
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eig_time =0;
pcg_time =0;

if m>0
%%%%%%%% Generate the consant term in A0 -tau*I+ A(x)%%%%
A0(1:m,1:m) = C(1:m,1:m);
A0 = A0 -tau*eye(n);
A0 = (A0+A0')/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Gamma =eye(n) +  Constraint_fuction(h,m,n);
 
% Gamma = Gamma.*(X0-C);
% val_obj = sum(sum(Gamma.*Gamma))/2;
C1 = X0(1:m,1:m)-C1;
initial_err = sum(sum(C1.*C1))/2;
else
    initial_err =0;
end 
 

 fprintf('\n')
 fprintf('The ininial correlation test approximation ===%d \n', initial_err)
if initial_err > tol0
 fprintf('\n')
 fprintf('The ininial correlation test is not good enough !!!!!!!!!! %d \n');
 fprintf('\n')
 fprintf('\n')
%fprintf('Newton: Norm of Gradient %d \n',norm_b) 
 fprintf('The Augmented Lagrangian method is activated!!!!!!!!!! %d \n');
 fprintf('\n')
 
 
%%%%%%%% Otherwise compute the initial x0 guess  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=m+1:n
    j = m*(i-m-1)+(i-m-1)*(i-m-2)/2;
    j =round(j);
    x0(j+1:j+i-1)= X0(i,1:i-1)'; % the ith row of the lower triangular part of C
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 =c;
% Gamma0 = eye(n);

k=0;
count_LS=0; % number of linear systems solved


% compute the projection of $(Lambda -lambda*A(x0))$ on S^n_+
     
     Gamma = A0 + Constraint_fuction(x0,m,n);
     
     %Gamma = Gamma0-lambda*Gamma;
     Gamma = Gamma0/lambda-Gamma;
     Gamma = (Gamma +Gamma')/2;
     eig_t0 =cputime;
     [P,D] = eig(Gamma);
     eig_time =  eig_time + cputime -eig_t0;
     d = diag(D);
     P =real(P);
     d =real(d);
  
     [f0,Fx] = grad(x0,c,Gamma0,lambda,d,P,m,n,h);  %%% the gradient of the augmented Lagrangian function 
      b = -Fx;
     f_eval =1;
     
while (k < maxit1)    
     j=0;  
     f_eval0 =0;
     Tol2_k = max(TOL2, TOLrel2/(1+lambda));
  
     norm_Lag = norm(b);
     fprintf('AugLagNewton: Norm of the Gradient === %d \n', norm_Lag)
    
    norm_Lag = 1.0e10; %Must do at least one step
    
     while (norm_Lag > Tol2_k) & (j<maxit2) 
        % compute dx   
        Omega = omega_mat(d);
       
        pcg_t0 =cputime;
        [dx,flag,relres,iterk] = p_cg(b,tol,maxit,ww,Omega,P,m,lambda,h);  
        pcg_time = pcg_time + cputime - pcg_t0;
        %one_lstime = cputime -pcg_t0;
        
        j=j+1; % the number of linear system solved at this level
        count_LS = count_LS + 1; % the total number of linear system solved 
        fprintf('AugLagNewton: Number of CG Iterations %d \n', iterk)
        tmp =  mu*Fx'*dx;      
        
        alpha = 1;
        Ok = false;
        
        t=0;
        
        while (Ok==false) & (t<maxit3) 
            % update dy
            delta_x = alpha*dx;
            x = x0 + delta_x;
          
     %         Gamma = Gamma0- lambda*Constraint_fuction(x,m,C);
    
             Gamma = Gamma0/lambda - (A0 + Constraint_fuction(x,m,n));
             Gamma = (Gamma +Gamma')/2;
            
             eig_t0 =cputime;
            [P,D] = eig(Gamma);
             eig_time =  eig_time + cputime -eig_t0;
            d = diag(D);
            P = real(P);
            d =real(d);
          
            
            [f,Fx] = grad(x,c,Gamma0,lambda,d,P,m,n,h);  %%% the function of the augmented Lagrangian function 
            
            f_eval0 = f_eval0 +1;
            
            tmp1 =alpha*tmp;
            
            dlag = f - f0 -tmp1;
               
            if  ( dlag < 1e-8)
                Ok = true;              
            else                     
                alpha = alpha/2;  
            end                  
            t=t+1;
         end  % This is the end of line search for one level Augmented Lagrangian function minimization
        
        % compute the first derivative of Lagrangian
        b = -Fx;
        norm_Lag = norm(b);
        fprintf('AugLagNewton: Number of Line Search Steps =========== %d \n', t)
        fprintf('AugLagNewton: Steplength ================= %d \n', alpha)
        fprintf('AugLagNewton: Norm of the Gradient === %d \n', norm_Lag)
        
        f0 = f;    
        x0=x;
   
    end   %%%% $j$ while loop:::This is the end of one level Augmented Lagrangian function minimization
     
   f_eval = f_eval + f_eval0;
    
 
  fprintf('AugLagNewton: Number of  Steps at this level ================== %d \n', j)
  fprintf('.... End of Another Level Aug Lagrangian Optimiztion  :):):):)...:):)...');  
  fprintf('\n');
   
  k=k+1;
 %%% Update the dual variable and penalty parameter %%% 
   %%%%% Checking if convergence reached %%%%%%   
    Gamma = P';
    for i=1:n
       Gamma(i,:) =max(0,d(i))*Gamma(i,:);
    end
    Gamma = P*Gamma; % (Gamma0/lambda - A(x))_+
    Gamma = (Gamma + Gamma')/2;
    Gamma0 =lambda*Gamma; %% New Gamma0
    
     %Gamma = Gamma0 - lambda*Constraint_fuction(x0,m,C);
     Gamma = Gamma0/lambda - (A0+Constraint_fuction(x0,m,n));
     Gamma = (Gamma + Gamma')/2;
     eig_t0 = cputime;
     [P,D] = eig(Gamma);
      eig_time =  eig_time + cputime -eig_t0;
     d = diag(D);
     P =real(P);
     d =real(d);
     
     f_eval =f_eval + 1;
     Gamma =P';
    i=1;
     while (i<=n)
       Gamma(i,:)=max(d(i),0)*Gamma(i,:);
       i=i+1;
     end
 
     Gamma= P*Gamma; %%%% (Gamma0_new/lambda - A(x))_+
     Gamma = (Gamma + Gamma')/2;

 Err_est = mytrace(Gamma0/lambda-Gamma, Gamma0/lambda-Gamma);
 
 
 Total_err = max(norm_Lag,  Err_est^0.5); 
 
 
 fprintf('\n');
 fprintf('AugLagNewton: Total Relative Error =============== %d \n', Total_err);   
  if Total_err <= TOLRel1
     Gamma = Gamma0 - (A0+Constraint_fuction(x0,m,n));
     %Gamma = Gamma0 - Constraint_fuction(I,J,x0,C);
     eig_t0 = cputime;
     [P1,D1] = eig(Gamma);
     eig_time  = eig_time + (cputime- eig_t0); 
     d1 = diag(D1);
     P1 = real(P1);
     d1 = real(d1);
     
     f_eval =f_eval + 1;
     
     Gamma =P1';
    i=1;
     while (i<=n)
       Gamma(i,:)=max(d1(i),0)*Gamma(i,:);
       i=i+1;
     end
     Gamma = P1*Gamma; %%%% (Gamma0_new - lambda*A(x))_+
     Gamma = (Gamma+Gamma')/2;
     
     Err_est = mytrace(Gamma0-Gamma, Gamma0-Gamma);
     
     Total_err = max(norm_Lag,  Err_est^0.5) % True error
 
     if Total_err <= TOL1
       break;
     else
       klambda =0;
       if (lambda < max_lambda)
         lambda = rho*lambda;  %new lambda
         klambda=1;
       end
    
        Gamma = Gamma0/lambda - (A0 + Constraint_fuction(x0,m,n));
        Gamma = (Gamma + Gamma')/2;
      if klambda == 1  %% Only when lambda is updated we need recompute Gamma
          eig_t0 = cputime;
          [P,D] = eig(Gamma);
          
          eig_t1 = cputime;
          eig_time  = eig_time + (eig_t1- eig_t0); 
            d = diag(D); 
            P =real(P);
            d =real(d);
            
          f_eval =f_eval + 1;
     end

     [f0,Fx] =   grad(x0,c,Gamma0,lambda,d,P,m,n,h);  %%% the gradient of the augmented Lagrangian function 
   
     b = -Fx;  
     end
     
   else
   klambda =0;
    if (lambda < max_lambda)
       lambda = rho*lambda;
       klambda=1;
    end
    
    % Gamma = Gamma0 - lambda*Constraint_fuction(x0,m,C);
     Gamma = Gamma0/lambda - (A0+Constraint_fuction(x0,m,n));
     Gamma = (Gamma + Gamma')/2;
    if klambda == 1  %% Only when lambda is updated we need recompute Gamma
           eig_t0 =cputime;
        [P,D] = eig(Gamma);
           eig_time = eig_time +cputime -eig_t0;
         d = diag(D); 
         
         P =real(P);
         d =real(d);
         
         f_eval =f_eval + 1;
     end
     [f0,Fx] = grad(x0,c,Gamma0,lambda,d,P,m,n,h);  %%% the gradient of the augmented Lagrangian function 
     
   
     b = -Fx;
     %norm_Lag = norm(b);
     
  
  end %convergence checking loop 
   
    

end % k while loop

 fprintf('\n');

 fprintf('\n');
 fprintf('Number of  Aug Lagrangian Optimiztion Levels solved===%d \n', k); 
 fprintf('Number of  linear system solved ===%d \n',  count_LS); 
 fprintf('Number of  function value calculations ===%d \n',  f_eval); 
 fprintf('AugLagNewton: Error ================= %d \n', Total_err);
 
 X0 = A0 + Constraint_fuction(x0,m,n) +tau*eye(n);  % restore the original X
 Gamma =eye(n) +  Constraint_fuction(h,m,n);
 
 Gamma = Gamma.*(X0-C);
 val_obj = sum(sum(Gamma.*Gamma))/2;
 fprintf('\n')

%fprintf('Newton: Norm of Gradient %d \n',norm_b)
fprintf('AugLagNewton: Final Objective Function value==== %d \n', val_obj);
 
else 
     fprintf('The ininial correlation test succeeds already!!!!!!!!!! %d \n');
end %%% corresponds to the first "if"  
eig_t0 =cputime;
[P,D] = eig(X0);
 eig_time = eig_time +cputime -eig_t0;
d = diag(D);
d =real(d);
lambda_min = min(d);


time_used= cputime-t0;
fprintf('\n');
fprintf('AugLagNewton: the least eigenvalue of the correlation matrix = %d \n', lambda_min);
fprintf('AugLagNewton: cputime used================= %d \n', time_used);
fprintf('AugLagNewton: Eig decomposition time used   ================ %d \n', eig_time);
fprintf('AugLagNewton: Conjugate Gradient time used  ================ %d \n', pcg_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                       %%
%%% end of the main program %%
%%%                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%% trace function %%%
function val = mytrace(A,B);

if (nargin == 1)
   val = sum(diag(A));
elseif (nargin == 2)
   val = sum(sum(A.*B));
end
return
%%% end of trace funcion %%%%%%

%%%%%% To generate the linear part of the constraint funcion with x: A(x) %%%%%%%
%%%%%%%

function Ax = Constraint_fuction(x,m,n)

Ax =zeros(n,n);

%%%%%
for i=m+1:n
    j = m*(i-m-1)+(i-m-1)*(i-m-2)/2;
    j =round(j);
    Ax(i,1:i-1)= x(j+1:j+i-1)'; % the ith row of the lower triangular part of Ax
    Ax(1:i-1,i) = x(j+1:j+i-1);  %  the ith column of the upper triangular part of Ax
end
%Ax(1:m,1:m)= C(1:m,1:m);
Ax=(Ax+Ax')/2;
return
 
%%%%%% End of generatig the constraint function with x: A(x) %%%%%%%
%%%%%%%

%%%%%% To generate the adjoint  the linear part of the afine funcion A(x) %%%%%%%
%%%%%%%

function x =  AtX(X,m,n)

m1 =n-m;
l = m1*n-(m1+1)*m1/2;
l =round(l); % the number of stressed entries
x =zeros(l,1);


for i=m+1:n
    j = m*(i-m-1)+(i-m-1)*(i-m-2)/2;
    j =round(j);
    x(j+1:j+i-1)= X(i,1:i-1)'; % the ith row of the lower triangular part of X
end 

x = 2*x; %%% A*X = 2h with h as the vector fomed by the strict upper trianglar part of the stressted entries
return
 
%%%%%% End of generatig the constraint funcion with x: A(x) %%%%%%%
%%%%%%%


%%%%%%
%%%%%% To generate F(x) %%%%%%%
%%%%%%%

function [f,Fx]= grad(x,c,Gamma0,lambda,d,P,m,n,h)


l = length(x);

f=0.0;
Fx =zeros(l,1);

%Im =find(d<0);
%Ip =find(d>=0);
%dp=max(0,d);
%H =diag(dp); %% H =P *diag(dp)* P^T
%  H =H*P'; %%% Assign H*P' to H
 H=P';
 i=1;
 while (i<=n)
     H(i,:)=max(d(i),0)*H(i,:);
     i=i+1;
 end
 
 H = P*H; 
 
 
 Fx = h.*(x-c);
 f = f + 0.5*mytrace(Fx,Fx);
 
 h = h.*h;
 Fx = x-c;
 Fx = Fx.*h;
 
 Fx = Fx -lambda*AtX(H,m,n);

 i=1;
 while (i<=n)
     f =f+lambda*(max(d(i),0))^2/2;
     i=i+1;
 end
 
f =f -mytrace(Gamma0,Gamma0)/(2*lambda);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of grad.m %%%%%%

%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the first -order difference of d
%%%%%%%
function omega = omega_mat(d)
n = length(d); 
omega =ones(n,n);
%Im =find(lambda<0);
%Ip =find(lamba>=0);
i=1;
while (i<=n)
    j=1;
    while (j<=n)
        if abs(d(i)-d(j))>1.0e-10
            omega(i,j) = (max(0,d(i))-max(0,d(j)))/(d(i)-d(j));
             elseif max(d(i),d(j))<=1.0e-15
                  omega(i,j)=0;
         end
      j=j+1;
    end
    i=i+1;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of omega_mat.m %%%%%%%%%%

%%%%%% PCG method %%%%%%%
%%%%%%% This is exactly the algorithm by  Hestenes and Stiefel (1952)
%%%%%An iterative method to solve A(x) =b  
%%%%%The symmetric positive definite matrix M is a
%%%%%%%%% preconditioner for A. 
%%%%%%  See Pages 527 and 534 of Golub and va Loan (1996)

function [p,flag,relres,iterk] = p_cg(b,tol,maxit,ww,Omega,P,m,lambda,h);

% Initializations
n = length(b); % the dimension of the unknown 

h = h.*h;

r = b;  %We take the initial guess x0=0 to save time in calculating A(x0) 
n2b =norm(b);    % norm of b
tolb = tol * n2b;  % relative tolerance 
p = zeros(n,1);
flag=1;
iterk =0;
relres=1000; %%% To give a big value on relres
% Precondition 
z =r./ww;  %%%%% z = M\r; here M =diag(ww); if M is not the identity matrix 
rz1 = r'*z; 
rz2 = 1; 
d = z;
% CG iteration
for k = 1:maxit
   if k > 1
       beta = rz1/rz2;
       d = z + beta*d;
   end
   
   w= h.*d +lambda*Jacobian_mat(d,Omega,P,m); %w = A(d); 
   denom = d'*w;
   iterk =k;
   relres =norm(z)/n2b;              %relative residue =norm(z) / norm(b)
   if denom <= 0 
       sssss=0
       p = d/norm(d); % d is not a descent direction
       break % exit
   else
       alpha = rz1/denom;
       p = p + alpha*d;
       r = r - alpha*w;
   end
   z = r./ww; %  z = M\r if M is not the identity matrix ;
   if norm(z) <= tolb % Exit if Hp=b solved within the relative tolerance
       iterk =k;
       relres =norm(z)/n2b;          %relative residue =norm(z) / norm(b)
       flag =0;
       break
   end
   rz2 = rz1;
   rz1 = r'*z;
end

return

%%%%%%%% %%%%%%%%%%%%%%%
%%% end of pre_cg.m%%%%%%%%%%%




%%%%%% To generate the Jacobain product with y: F'(x)(d) %%%%%%%
%%%%%%%

function JFd= Jacobian_mat(d,Omega,P,m)

[n,n1] =size(Omega); 
l = length(d); 

JFd = zeros(l,1);

P1 = P(1:m,1:m);
P2 = P(1:m,m+1:n);
P3 = P(m+1:n,1:m);
P4 = P(m+1:n,m+1:n);

H = Constraint_fuction(d,m,n);   % Ad  is produced
%H = H - Constraint_fuction(JFd,m,C); % Ad is produced

H2 = H(1:m,m+1:n);
H4 = H(m+1:n,m+1:n);

S2 = H2*P4; 
S3 = H2'*P1;
%S4 = H2'*P2;

W1 = P3'*S3;
W1 = (W1'+W1);
W1 = W1 + P3'*H4*P3;

S4 = (H2'*P2+ H4*P4); 

W2 = P1'*S2 + P3'*S4;

W4 = P2'*S2 + P4'*S4;

W1 = Omega(1:m,1:m).*W1;
W2 = Omega(1:m,m+1:n).*W2;
W4 = Omega(m+1:n,m+1:n).*W4;

S2 = W1*P3' + W2*P4';
S4 = W2'*P3' + W4*P4';

H2 = P1*S2+P2*S4;
H4 = P3*S2 + P4*S4;

H(1:m,m+1:n)   = H2;
H(m+1:n,m+1:n) = H4;
H(m+1:n,1:m)   = H2';

H =(H+H')/2; 

JFd = AtX(H,m,n);  %%% A*{ P*[Omega o (P^T *A(d)*P)]*P^T }
 return
 
%%%%%%%%%%%%%%%
%end of Jacobian_mat.m%%%








