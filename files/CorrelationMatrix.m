%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  This code is designed to solve                                 
%%%%%%%%%    min  0.5*<X-G, X-G>
%%%%%%%%%    s.t. X_ii =b_i, i=1,2,...,n
%%%%%%%%%         X>=tau*I (symmetric and positive semi-definite)          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Based on the algorithm developed in                            
%%%%%%%%%     "A Quadratically Convergent Newton Method for               
%%%%%%%%%       Computing the Nearest Correlation Matrix"                 
%%%%%%%%%              By Houduo Qi and Defeng Sun                         
%%%%%%%%%  SIAM Journal on Matrix Analysis and  Applications 28 (2006) 360--385.
%%%%%%%%%  
%%%%%%%%%  Last modified date: August 30, 2019  by Qian LI (18090081g@connect.polyu.hk)
%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                          
%%%%%%%%%  The  input arguments  G, b>0, tau>=0, and tol (tolerance error)            
%%%%%%%%%                                                                         
%%%%%%%%%     For correlation matrix, set b =ones(n,1)                     
%%%%%%%%%                                                                  
%%%%%%%%%     For a positive definite matrix                               
%%%%%%%%%         set tau = 1.0e-5 for example                             
%%%%%%%%%         set tol = 1.0e-6 or lower if no very high accuracy required     
%%%%%%%%%  The outputs are the optimal primal and dual solutions 
%%%%%%%%%  Diagonal Preconditioner is added         
%%%%%%%%%  Send your comments and suggestions to   
%%%%%%%%%  hdqi@soton.ac.uk  or defeng.sun@polyu.edu.hk 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
%%%%%%%%% Warning: Though the code works extremely well, it is your call to use it or not.  


function [X,y] = CorrelationMatrix(G,b,tau,tol)

 fprintf('*********************************************************************************');
 fprintf('\n*  --- Semismooth Newton-CG method starts  --- ');
 fprintf('\n*  Developed by Houduo Qi and Defeng Sun  ');
 fprintf('\n*  Based on the algorithm in `A Quadratically Convergent Newton Method for ');
 fprintf('\n*  Computing the Nearest Correlation Matrix'' ' );
 fprintf('\n*  SIAM J. Matrix Anal. Appl. 28 (2006) 360--385. ');
 fprintf('\n* --- This version: August 30, 2019 ---- ');
 fprintf('\n*********************************************************************************');
 t0 = clock;
 [n,m] =size(G);

%global  b0

 G =(G + G')/2;  % make G symmetric
     b0 = ones(n,1);     %% default
     c = ones(n,1);     %% diagonal preconditioner
 error_tol = 1.0e-6; %% default for the termination accuracy
 
if nargin==1 
   tau = 0;
end


if nargin==2
   b0 = b; % reset b0   
   tau = 0;
end

if nargin==3
      b0 =  b;
   if tau > 0
      b0 =  b0 - tau*ones(n,1); % reset b0 
      G  =  G - tau*eye(n); % reset G
   end
end

if nargin==4
     b0 = b; 
   if tau > 0
      b0 =  b0 - tau*ones(n,1); % reset b0 
      G  =  G - tau*eye(n); % reset G
   end
   error_tol = max(1.0e-12,tol); %reset error tolerance
end

Res_b = zeros(300,1);
norm_b0 =norm(b0);

y = zeros(n,1);       %Initial point
%y=b0-diag(G);              
 

k=0;
f_eval = 0;

Iter_Whole = 500;
Iter_inner = 20; % Maximum number of Line Search in Newton method
maxit = 200; %Maximum number of iterations in PCG
iterk = 0;
Inner = 0;
tol = 1.0e-2; %relative accuracy for CGs


sigma_1 = 1.0e-4; %tolerance in the line search of the Newton method


prec_time = 0;
pcg_time = 0;
eig_time =0;

 

  val_G = sum(sum(G.*G))/2;

  X = G + diag(y);
 
  eig_time0 = clock;
  [P,lambda] = Mymexeig(X);   %% X= P*diag(D)*P'
  eig_time = eig_time + etime(clock,eig_time0);
  [f0,Fy] = gradient(y,lambda,P,b0,n);
 
  Initial_f = val_G - f0;
  
  
  X = PCA(X,lambda,P,b0,n);  %% generate a feasible primal solution by using PCA
  val_obj = sum(sum((X - G).*(X - G)))/2;
  gap = (val_obj - Initial_f)/(1+ abs(Initial_f) + abs(val_obj));
 
  f = f0;
  f_eval = f_eval + 1;
  b = b0 - Fy;
  norm_b = norm(b);
  time_used = etime(clock,t0); 
  eta = norm_b/(1+norm_b0);
  Omega12 = omega_mat(P,lambda,n);
  x0 = y;

 fprintf('\n matrix dimension n = %d', n);
 fprintf('\n ---------------------------------------------------------');
 fprintf('\n  iter     pobj          dobj            relgap          etaorg        eta      time | cg_its inner_its');
 fprintf('\n   0   %- 5.4e     %-5.4e      %-3.2e        %3.2e     %3.2e     %3.1f', ...
     val_obj, Initial_f, gap, norm_b, eta, time_used);

 %while (abs(gap) > error_tol  & k< Iter_Whole) %% based on the relative duality gap for practitioners
while (eta > error_tol & k< Iter_Whole) %% based on the relative gradient for academic research
 size_plus = length(Omega12); % check if there positive eigenvalues
  prec_time0 = clock;
  if size_plus >0
   c = precond_matrix(Omega12,P,n); % comment this line for  no preconditioning
  end 
  prec_time = prec_time + etime(clock, prec_time0);
  
  size_plus = length(Omega12); % check if there positive eigenvalues
  
  if size_plus >0
     pcg_time0 = clock;
     [d,flag,relres,iterk]  = pre_cg(b,tol,maxit,c,Omega12,P,n);
      pcg_time = pcg_time + etime(clock,pcg_time0);
      %d = b0-Fy; gradient direction
 %fprintf('Newton-CG: Number of CG Iterations == %d \n', iterk)
  
     if (flag~=0); % if CG is unsuccessful, use the negative gradient direction
        % d =b0-Fy;
        fprintf('\n *********************************************************************************\n')
        
         disp('..... Warning: This step is not a completed Newton-CG step .....')
     end
      
  else
         d  =b0-Fy;
         fprintf('\n *********************************************************************************\n')
      disp('..... Warning: This step uses the negative gradient direction .....')
  end
  
 slope = (Fy-b0)'*d; %%% nabla f d
 

     y = x0 + d; %temporary x0+d 

      X = G + diag(y);

      eig_time0 = clock;
      [P,lambda] = Mymexeig(X); % Eig-decomposition: X =P*diag(D)*P^T
      eig_time = eig_time + etime(clock,eig_time0); 
      [f,Fy] = gradient(y,lambda,P,b0,n);
        

     k_inner = 0;
     while(k_inner <=Iter_inner & f> f0 + sigma_1*0.5^k_inner*slope + 1.0e-6)
         k_inner = k_inner+1;
         y = x0 + 0.5^k_inner*d; % backtracking   
         
         X = G + diag(y);

         
         eig_time0 = clock;
             [P,lambda] = Mymexeig(X); % Eig-decomposition: X =P*diag(D)*P^T
             eig_time = eig_time + etime(clock,eig_time0); 
            [f,Fy] = gradient(y,lambda,P,b0,n);
      end % loop for while
      f_eval = f_eval + k_inner+1;
      x0 = y;
      f0 = f;
       val_dual = val_G - f0;
       X = PCA(X,lambda,P,b0,n);
       val_obj = sum(sum((X - G).*(X - G)))/2;
       gap = (val_obj - val_dual)/(1+ abs(val_dual) + abs(val_obj));
      
     k=k+1;  
     b = b0 - Fy;
     norm_b = norm(b);
     eta = norm_b/(1 + norm_b0);
     time_used = etime(clock,t0);
      Res_b(k) = norm_b;
      fprintf('\n   %d   %- 5.4e     %-5.4e      %-3.2e        %3.2e     %3.2e     %3.1f | %d       %d', ...
     k, val_obj, val_dual, gap, norm_b, eta, time_used, iterk,k_inner);
    
     Omega12 = omega_mat(P,lambda,n);

 end %end loop for while i=1;
 

 rank_X = length(find(max(0,lambda)>0));  
 Final_f = val_G - f;
 if tau > 0
    X = X + tau*eye(n); 
 end
 time_used = etime(clock,t0);
fprintf('\n =====================================================')
fprintf('\n eta = %3.2e, etaorg = %3.2e', eta, norm_b)
fprintf('\n per eig time = %5.4f', eig_time/f_eval)
fprintf('\n time per iteration = %5.4f', time_used/k)
fprintf('\n');
fprintf('Newton-CG: Number of Iterations ============ %d \n', k)
fprintf('Newton-CG: Number of Function Evaluations == %d \n', f_eval)
fprintf('Newton-CG: Final Dual Objective Function value ========== %d \n',Final_f)
fprintf('Newton-CG: Final primal Objective Function value ======== %d \n', val_obj)
fprintf('Newton-CG: The final relative duality gap ========================  %d \n',gap)
fprintf('Newton-CG: The rank of the Optimal Solution - tau*I ================= %d \n',rank_X)

fprintf('Newton-CG: computing time for computing preconditioners == %d \n', prec_time)
fprintf('Newton-CG: computing time for linear systems solving (cgs time)            == %d \n', pcg_time)
fprintf('Newton-CG: computing time for  eigenvalue decompositions (calling eig time)== %d \n', eig_time)
fprintf('Newton-CG: computing time used for equal weights calibration ============================== %d \n',time_used)



%%% end of the main program



%%%%%%
%%%%%% To generate F(y) %%%%%%%
%%%%%%%

function [f,Fy]= gradient(y,lambda,P,b0,n)
 
%  
% f = 0.0;
% Fy =zeros(n,1);
 
rankX = sum(lambda > 0);

if rankX == 0
   Fy = zeros(n,1);
   f = 0;
else
   if rankX < n
      P1 = P(:,1:rankX);
      lambdanew = lambda(1:rankX);
   else
      P1 = P;
      lambdanew = lambda;
   end
   tmp = spdiags(lambdanew.^(0.5),0,rankX,rankX);
   Ptmp = P1*tmp;
   Fy = sum(Ptmp.*Ptmp,2);
   f = norm(lambdanew)^2;
end
f =0.5*f -b0'*y;






return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of gradient.m %%%%%%


%%%% use PCA to generate a primal feasible solution %%%
function [X] = PCA(X,lambda,P,b0,n)
 
Ip = find(lambda>0);
r = length(Ip);
 
if (r==0)
    X =zeros(n,n);
elseif (r==n)
    X = X;
elseif (r<=n/2)
    lambda1 = lambda(Ip);
    lambda1 = lambda1.^0.5;
    P1 = P(:, 1:r);
    if r >1
        P1 = P1*spdiags(lambda1,0,r,r); %sparse(diag(lambda1));
        X = P1*P1'; %
    else
        X = lambda1^2*P1*P1';
    end
    
else
    lambda2 = -lambda(r+1:n);
    lambda2 = lambda2.^0.5;
    P2 = P(:, r+1:n);
    P2 = P2*sparse(diag(lambda2));
    X = X + P2*P2';
end
 
  
  %%% To make X positive semidefinite with diagonal elements exactly b0
     d = diag(X);
     d = max(b0, d); 
     for i = 1:n
        X(i,i) = d(i);
     end
     %X = X - diag( diag(X)) + diag(d); %%% make the diagonal vector to be b0, still PSD
     d = d.^(-0.5);
     d = d.*(b0.^0.5);
%      X = X.*(d*d'); 
D = spdiags(d,0,n,n);
X = (D*X)*D;
     
     return
 %%%%%%%%%%%%%% end of PCA %%%
 


%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the first -order difference of lambda
%%%%%%%

%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the essential part of the first -order difference of d
%%%%%%%
function Omega12 = omega_mat(P,lambda,n)
%We compute omega only for 1<=|idx|<=n-1
idx.idp = find(lambda>0);
idx.idm = setdiff([1:n],idx.idp);
n =length(lambda);
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

    %%***** perturbation *****
    return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of omega_mat.m %%%%%%%%%%

%%%%%% PCG method %%%%%%%
%%%%%%% This is exactly the algorithm by  Hestenes and Stiefel (1952)
%%%%%An iterative method to solve A(x) =b  
%%%%%The symmetric positive definite matrix M is a
%%%%%%%%% preconditioner for A. 
%%%%%%  See Pages 527 and 534 of Golub and va Loan (1996)

function [p,flag,relres,iterk] = pre_cg(b,tol,maxit,c,Omega12,P,n);
% Initializations
r = b;  %We take the initial guess x0=0 to save time in calculating A(x0) 
n2b =norm(b);    % norm of b
tolb = tol * n2b;  % relative tolerance 
p = zeros(n,1);
flag=1;
iterk =0;
relres=1000; %%% To give a big value on relres
% Precondition 
z =r./c;  %%%%% z = M\r; here M =diag(c); if M is not the identity matrix 
rz1 = r'*z; 
rz2 = 1; 
d = z;
% CG iteration
for k = 1:maxit
   if k > 1
       beta = rz1/rz2;
       d = z + beta*d;
   end
   %w= Jacobian_matrix(d,Omega,P,n); %w = A(d); 
   w = Jacobian_matrix(d,Omega12,P,n); % W =A(d)
   denom = d'*w;
   iterk =k;
   relres = norm(r)/n2b;              %relative residue = norm(r) / norm(b)
   if denom <= 0 
       sssss = 0;
       p = d/norm(d); % d is not a descent direction
       break % exit
   else
       alpha = rz1/denom;
       p = p + alpha*d;
       r = r - alpha*w;
   end
   z = r./c; %  z = M\r; here M =diag(c); if M is not the identity matrix ;
   if norm(r) <= tolb % Exit if Hp=b solved within the relative tolerance
       iterk =k;
       relres = norm(r)/n2b;          %relative residue =norm(r) / norm(b)
       flag =0;
       break
   end
   rz2 = rz1;
   rz1 = r'*z;
end

return

%%%%%%%% %%%%%%%%%%%%%%%
%%% end of pre_cg.m%%%%%%%%%%%




%%%%%% To generate the Jacobain product with x: F'(y)(x) %%%%%%%
%%%%%%%

function Ax = Jacobian_matrix(x,Omega12,P,n)
Ax =zeros(n,1);
[r,s] = size(Omega12);
if (r>0 && r<n)
    P1=P(:,1:r);
    P2=P(:,(r+1):n);
    O12=Omega12.*(P1'*bsxfun(@times,P2,x));
    PO=(P1*O12).*P2;
    om=ones(n,1);
    hh=2*PO*om(1:s,1);         %hh=2*diag(P1(Omega12.*(P1'*H*P2))P2')
    if (r<n/2)
        PP1=P1*P1';
        PP11=PP1.^2;
        Ax=PP11*x+hh+1e-10*x;  %Ax=P1*P1'*H*P1*P1'+hh+1e-10*h
    else
        PP2=P2*P2';
        PP22 = PP2.^2;
        Ax=x+PP22*x+hh-2*x.*diag(PP2)+1e-10*x;  %Use the property P1*P1'+P2*P2'=E
    end
else
    if (r==n)
        Ax =(1+1.0e-10)*x;
    end
end

return
 
%%%%%%%%%%%%%%%
%end of Jacobian_matrix.m%%%

%%%%%% To generate the diagonal preconditioner%%%%%%%
%%%%%%%

function c = precond_matrix(Omega12,P,n)

[r,s] =size(Omega12);
c = ones(n,1);

if (r>0)

    if (r< n/2)
        H = P';
        H = H.*H;

        H12 = H(1:r,:)'*Omega12;
        d =ones(r,1);
        for i=1:n
            c(i) = sum(H(1:r,i))*(d'*H(1:r,i));
            c(i) = c(i) + 2.0*(H12(i,:)*H(r+1:n,i));
            if c(i) < 1.0e-8
                c(i) =1.0e-8;
            end
        end
    else  % if r>=n/2, use a complementary formula
        if (r < n)
            H = P';
            H = H.*H;
            Omega12 = ones(r,s)-Omega12;
            H12 = Omega12*H(r+1:n,:);
            d =ones(s,1);
            dd = ones(n,1);
            
            for i=1:n
                c(i) = sum(H(r+1:n,i))*(d'*H(r+1:n,i));
                c(i) = c(i) + 2.0*(H(1:r,i)'*H12(:,i));
                alpha = sum(H(:,i));
                c(i) = alpha*(H(:,i)'*dd)-c(i);
                if c(i) < 1.0e-8
                    c(i) =1.0e-8;
                end
            end
        end

    end
end

return

 
%%%%%%%%%%%%%%%
%end of precond_matrix.m%%%

%
function [P,lambda] = Mymexeig(X)

if (exist('mexeig')==3)
    
    [P,D] = mexeig(full(X));
 if size(D,2)>1
     lambda =diag(D);
 else
     lambda =D;
 end
    
else
    
    [P,D]  = eig(full(X));
    
    lambda = diag(D);

end   
D=real(D);
    
P          = real(P);
lambda     = real(lambda);
if issorted(lambda)
    lambda = lambda(end:-1:1);
    P      = P(:,end:-1:1);
elseif issorted(lambda(end:-1:1))
    return;
else
    [lambda, Inx] = sort(lambda,'descend');
    P = P(:,Inx);
end

return

%*********************** End of Mymexeig.m ********************************