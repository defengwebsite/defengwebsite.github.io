clear all
  
tau =.0e-8;


%% initial of the problem 
%% Case I 
% %%%%%% RiskMetrics Data (n =387)
% load x.mat
% G = subtract(x);
% G  = (G +G')/2;
% [n, n_c]=size(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%
 n = 1000;

%% Case II
% % % % %%%%%%%%%%%%%%%%%%%%%
% n = 500;
% x = 10.^[-4:4/(n-1):0];
% G = gallery('randcorr',n*x/sum(x));
% %%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:n
%     G(i,i) =1;
% end




%% Case I
 E = 2.0*rand(n,n)-ones(n,n);
   %E = 2.0*rand(n,n);
% E = randn(n,n);
  E = triu(E) + triu(E,1)';
% E = (E+E')/2;


 
 
  
  


alpha = .1;
 % G = (1-alpha)*G+ alpha*E;
% norm_E = norm(E, 'fro')
 
 G = E;
for i=1:n
    G(i,i) =1;
end


B = rand(n,n);
W = B'*B;
W = (W + W')/2;
%cond_W0 = cond(W)
%W =  n^0.5* W/norm(W,'fro'); % scaled such that its F-norm is n^0.5
[P,D] = eig(W);
 
al = 0;
d = (1-al )*rand(n,1) + al *ones(n,1);
% P = eye(n); %%% diagonal weighting
W = P*sparse(diag(d));
W = W*P';

 
 
%  w = rand(n,1);
%  W = diag(w);
 
 cond_W = cond(W) 
 

%%% number of fixed off diagonal elements in each row
lh =  50;          
lh = min(lh,n-1);
I_h = [];
J_h = [];
for i = 1:n-lh
    r = rand(n-i,1);
    [r,ind] = sort(r);
    I_h = [I_h; i*ones(lh,1)];
    J_h = [J_h; i+ind(n-i-lh+1:n-i)];
end
for i = ((n-lh)+1):(n-1)
    I_h = [I_h; i*ones(n-i,1)];
    J_h = [J_h;[(i+1):n]'];
end
k_h = length(I_h);
%%%% off-diagonal value
h = zeros(k_h,1);

%%% for fixed diagonal entries
I_d = [1:1:n]';
J_d = I_d;
%%%% diagonal value
rhs = ones(n,1);      % diagonal elements
alpha0 = 1;
rhs = alpha0*rhs + (1-alpha0)*rand(n,1);
 

I = [I_d;I_h];
J = [J_d;J_h];
b = [rhs;h];
k = length(b);




 %[X,y] =  CorMatWnorm(G,W,rhs,tau); 
 fprintf('\n')
 disp('=====end========end==========end===========end=============end==========end==========')
  
 [X,y] =  CorNewton3_Wnorm(G,b,I,J,tau,W);
 

  
    

   
 
   

   
 
 
 