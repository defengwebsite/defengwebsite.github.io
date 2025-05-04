function [X,y] = CorMatWnorm(G,W,b,tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%              This code is for computing 
%%%%%% "the W-weighted version of the nearest correlation matrix problem"
%%%%%%                          based on 
%%%%%%  Houduo Qi and Defeng Sun, "A Quadratically Convergent Newton Method
%%%%%%  for Computing the Nearest Correlation Matrix"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%                        for solving 
%%              min  0.5 ||W^1/2*(X - G)*W^1/2 ||^2    
%%              s.t.  X_ii =1,       i=1,2,..., n   
%%                    X >= tau*I    (symmetric and positive semi-definite; tau can be zero)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input:
% G: the estimated potentially inconsistent correlation matrix (n by n)
% W: the weight matrix for G
% b: = ones(n,1)
% tau: the lower bound for the smallest eigevalue of X0 (can be zero)

% Output:
% X: the calibrated correlation matrix
% y: Lagrangian dual variable corresponding to X_ii =1
% Warning: accuracy may not be guaranteed %%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Last modified on April 5, 2008  %%%%%%%%%%%%%%%%%%%%%

fprintf('\n --- Semimsooth Newton method with diagonal preconditioner--- \n')

t0 = cputime;
n = length(G);
G = (G + G')/2;        % make G symmetric
G0 = G;
W = real(W);
W = (W + W')/2;        % make W symmetric
 

if nargin == 4
    G  =  G - tau*eye(n);        % reset G
end

if nargin == 3
  tau =0;       
end

if nargin == 2
  tau = 0;  
  b   = ones(n)
end

if nargin == 1
  tau = 0;  
  b   = ones(n);
  W = eye(n);
end

b0 =b;


b0 = b0-tau*ones(n,1);


eig_time =0;
t1 = cputime;
[P,D] = eig(W);
eig_time = eig_time + cputime - t1;

lambda = diag(D);
P = real(P);
lambda = real(lambda);

if max(lambda)<=0
    disp('Warning: Warning: Warning: W is not positive definite: Stop Stop Stop')
   
    return;
end

W_half = P';
W_half_inv =P';
i=1;
while (i<n+1)
    W_half(i,:) = lambda(i)^0.5*W_half(i,:);
    W_half_inv(i,:) = lambda(i)^(-0.5)*W_half_inv(i,:);
    i=i+1;
end
W_half = P * W_half;
W_half = (W_half + W_half')/2;
W_half_inv = P * W_half_inv;           

W_half_inv  = (W_half_inv + W_half_inv)/2;

G = W_half*G*W_half;
G = (G+G')/2;

 
Iter_Whole = 200;
Iter_inner = 20;       % Maximum number of Line Search in Newton method
maxit =200;            % Maximum number of iterations in PCG

error_tol = 1.0e-6;    % termination tolerance
sigma = 1.0e-4;        % tolerance in the line search of the Newton method
tol = 1.0e-2;           % relative accuracy for CGs

k = 0;
f_eval = 0;
iterk = 0;
Inner = 0;
y = zeros(n,1);        % initial point
%y = b -diag(G);

Fy = zeros(n,1);
d = zeros(n,1);

CG_num = 0;
CG_time = 0;

prec_time = 0;

% approximate preconditioner
c = ones(n,1);
%%% no need update for the following case
 

val_G = sum(sum(G.*G))/2;

C = W_half_inv;
i=1;
while (i<n+1)
    C(i,:) = y(i)*C(i,:);
    i=i+1;
end


C = G + W_half_inv * C; % C = G + W^{-1/2}*diag(y)*W^{-1/2}
C = (C+C')/2;

t1 = cputime;
[P,D] = eig(C);
eig_time = eig_time + cputime - t1;

lambda = diag(D);
P = real(P);
lambda = real(lambda);

W0P = W_half_inv * P;

[f,Fy] = gradient(y,lambda,W0P,b0,n);
f_eval = f_eval + 1;

f0 = f;
x0 = y;
Omega = omega_mat(P,lambda,n);

b = -(Fy-b0);
norm_b = norm(b);

Initial_f = val_G - f0;

  fprintf('Newton: Initial Dual Objective Function value ==== %d \n', -f0)
  fprintf('Newton: Norm of Gradient for the initial point ==== %d \n',norm_b)

fprintf('\n  Iter.   Num. of CGs      Step length      Norm of Gradient      time_used ')

while (norm_b > error_tol & k< Iter_Whole)
    
    % update approximate preconditioner
    t1 = cputime;
    c = precond_matrix(Omega,W0P);
    prec_time = prec_time + cputime - t1;
    
    
    
    t2 = cputime;
    [d,flag,relres,iterk] = pre_cg(b,tol,maxit,c,Omega,W0P);
    CG_time = CG_time + cputime - t2;
    CG_num = CG_num + iterk;

    if (flag~=0);      % if CG is unsuccessful, use the negative gradient direction
        disp('..... Not a full Newton step......')
    end
    
    slope = (Fy-b0)'*d; 

    y = x0 + d;   
    
        C = W_half_inv;
        i=1;
        while (i<n+1)
          C(i,:) = y(i)*C(i,:);
          i=i+1;
        end

      C = W_half_inv * C + G;
      C = (C+C')/2;
    
    t1 = cputime;
    [P,D] = eig(C);
    eig_time = eig_time + cputime - t1;
    
    lambda = diag(D);
    P = real(P);
    
    lambda = real(lambda);
    
    W0P = W_half_inv * P;
    
    
    [f,Fy] = gradient(y,lambda,W0P,b0,n);

    k_inner=0;
    
    while ( k_inner <= Iter_inner & ( f -f0 -  sigma*0.5^k_inner*slope )/max(1,abs(f0)) > 1.0e-8 )

        k_inner = k_inner+1;
        
        y = x0 + 0.5^k_inner*d;              % backtracking
       
        C = W_half_inv;
        i=1;
        while (i<n+1)
          C(i,:) = y(i)*C(i,:);
          i=i+1;
        end

        C = W_half_inv * C +G;
        C = (C+C')/2;

        
        t1 = cputime;
        [P,D] = eig(C);
        eig_time = eig_time + cputime - t1;
        
        lambda = diag(D);
        P = real(P);
        lambda = real(lambda);
        
        W0P = W_half_inv * P;
        
        [f,Fy] = gradient(y,lambda,W0P,b0,n);
     
       % fprintf('\n Newton: function value =========================== %d \n',f)
        
    end      % loop for while
    
    k = k+1;
    f_eval = f_eval + k_inner+1;
    
    x0 = y;
    f0 = f;
    b = -(Fy-b0);
    norm_b = norm(b);

    Omega = omega_mat(P,lambda,n);
    
    tt = cputime - t0;
    [hh,mm,ss] = time(tt);
    
    fprintf('\n   %2.0d         %2.0d            %2.1e            %3.2e           %d:%d:%d ',k,iterk,0.5^k_inner,norm_b,hh,mm,ss)

end   % end loop for while

 
C = P';
i=1;
while (i<n+1)
    C(i,:) = max(0,lambda(i))*C(i,:);
    i=i+1;
end

X = P*C;  

X = (X+X')/2;


%Final_f = val_G - f; 

% recover X to the original one
X = W_half_inv*X*W_half_inv;
X = X + tau*eye(n);
val_obj = sum(sum((X-G0).*(X-G0)))/2;


time_used = cputime-t0;

fprintf('\n')
fprintf('\n ================ Final Information ================= \n');
fprintf(' Total number of iterations      = %2.0f \n',k);
fprintf(' Number of func. evaluations     = %2.0f \n',f_eval)
fprintf(' Number of CG iterations         = %2.0f \n',CG_num)
fprintf(' Primal objective value          = %d \n',val_obj)
%fprintf(' Dual objective value            = %d \n',Final_f)
fprintf(' Norm of Gradient                = %3.2e \n', norm_b)
fprintf(' CPU time for preconditioner     = %3.1f \n',prec_time)
fprintf(' CPU time for PCG iterations     = %3.1f \n',CG_time)
fprintf(' CPU time for eigen-decom        = %3.1f \n',eig_time)
fprintf(' Total CPU time (secs)           = %3.1f \n',time_used)
fprintf(' ====================================================== \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% end of the main program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%% To change the format of time 
function [h,m,s] = time(t)
t = round(t); 
h = floor(t/3600);
m = floor(rem(t,3600)/60);
s = rem(rem(t,60),60);
%%% End of time.m



%%%%%% To generate F(y) %%%%%%%
function [f,Fy] = gradient(y,lambda,W0P,b0,n)

f = 0.0;
Fy = zeros(n,1);

H = W0P'; % W0P = W^{-1/2}*P
i=1;
while (i<=n)
    H(i,:) = max(lambda(i),0)*H(i,:);
    i=i+1;
end
i=1;
while (i<=n)
    Fy(i) = W0P(i,:)*H(:,i);
    i = i+1;
end  

i = n;
while (i>=1)
    f = f + (max(lambda(i),0))^2;
    i = i-1;
end

f = 0.5*f - b0'*y;

return

%%%%% end of gradient.m %%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the first -order difference of d  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Omega = omega_mat(P,lambda,n)
%We compute omega only for 1<=|idx|<=n-1
idx.idp = find(lambda>0);
idx.idm = setdiff([1:n],idx.idp);
n =length(lambda);
r = length(idx.idp);
Omega = zeros(n);

if ~isempty(idx.idp)
    if (r == n)
        Omega = ones(n,n);
    else
        s = n-r;
        if idx.idp(1)< idx.idm(1)
            dp = lambda(1:r);
            dn = lambda(r+1:n);
            Omega12 = (dp*ones(1,s))./(abs(dp)*ones(1,s) + ones(r,1)*abs(dn'));
            Omega12 = max(1e-15,Omega12);
            Omega =[ones(r) Omega12;Omega12' zeros(s)];
        else
            dp = lambda(s+1:n);
            dn = lambda(1:s);
            Omega12 = (dp*ones(1,s))./(abs(dp)*ones(1,s) + ones(r,1)*abs(dn'));
            Omega12 = max(1e-15,Omega12);
            Omega =[zeros(s) Omega12';Omega12 ones(r)];
        end

    end
end

%%***** perturbation *****
return

%%%% end of omega_mat.m %%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    PCG method     %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% This is exactly the algorithm by  Hestenes and Stiefel (1952)
%%%%%%%  An iterative method to solve A(x) =b  
%%%%%%%  The symmetric positive definite matrix M is a
%%%%%%%  preconditioner for A. 
%%%%%%%   See Pages 527 and 534 of Golub and va Loan (1996)

function [p,flag,relres,iterk] = pre_cg(b,tol,maxit,c,Omega,W0P)

n = length(W0P);

% Initializations
r = b;             % We take the initial guess x0=0 to save time in calculating A(x0) 
n2b =norm(b);      % norm of b
tolb = tol * n2b;  % relative tolerance 
p = zeros(n,1);
flag=1;
iterk =0;
relres=1000;       % To give a big value on relres

% Precondition 
z =r./c;           % z = M\r; here M =diag(c); if M is not the identity matrix 
rz1 = r'*z; 
rz2 = 1; 
d = z;
% CG iteration
for k = 1:maxit
   if k > 1
       beta = rz1/rz2;
       d = z + beta*d;
   end
   w = Jacobian_matrix(d,Omega,W0P); % W =A(d)
   denom = d'*w;
   iterk =k;
   relres =norm(r)/n2b;              %relative residue =norm(r) / norm(b)
   if denom <= 0 
       sssss=0
       p = d/norm(d); % d is not a descent direction
       break % exit
   else
       alpha = rz1/denom;
       p = p + alpha*d;
       r = r - alpha*w;
   end
   z = r./c;                  %  z = M\r; here M =diag(c); if M is not the identity matrix ;
   if norm(r) <= tolb         %  Exit if Hp=b solved within the relative tolerance
       iterk =k;
       relres =norm(r)/n2b;   % relative residue =norm(z) / norm(b)
       flag =0;
       break
   end
   rz2 = rz1;
   rz1 = r'*z;
end

return

%%%%%%%% end of pre_cg.m %%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the Jacobain product with x: F'(y)(x) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ax = Jacobian_matrix(x,Omega,W0P)

n = length(W0P);
Ax = zeros(n,1);


H = W0P;
i=1;
while (i<=n)
    H(i,:) = x(i)*H(i,:);
    i=i+1;
end
H = W0P'*H;
H = Omega.*H;          % H =[Omega o ( P^T*(W^(-1/2)*diag(x)*W^(-1/2))*P )]

H = H*W0P';
i=1;
 while (i<=n)
       Ax(i) = W0P(i,:)*H(:,i);
       Ax(i) = Ax(i) + 1.0e-10*x(i); % add a small perturbation
       i=i+1;
 end

 return

%%%%%% end of Jacobian_matrix.m 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  To generate the diagonal preconditioner  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We chooose e1-en to compute the diagonal elements of Jacobian matrix

function c = precond_matrix(Omega,W0P)

n = length(W0P);
c = ones(n,1);

 
H = W0P';
H = H.*H;
Omega = Omega*H;
H = H';

for i=1:n
    c(i) = H(i,:)*Omega(:,i);
    if c(i) <= 1.0e-8
        c(i) =1.0e-8;
    end
end
return
%%% end of precond_matrix.m

 













