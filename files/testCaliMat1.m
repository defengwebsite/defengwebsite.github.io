clear all


tau   = .0e-8;
OPTIONS.tau0 = tau;


%% G  

% % Case I 
% % %%%% RiskMetrics Data (n =387)
%   load x.mat
%   G = subtract(x);
%   G  = (G +G')/2;
%   [n, n_c]=size(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%


% Case II
% % % %%%%%%%%%%%%%%%%%%%%%
%  n = 200;
% x = 10.^[-4:4/(n-1):0];
% G = gallery('randcorr',n*x/sum(x));
% %%%%%%%%%%%%%%%%%%%%%%%%%


  n = 1000;   


%%  
% Case I  (-1,1)
 E = 2.0*rand(n,n)-ones(n,n);
 
% Case II  (0,1)
% E = rand(n,n);

% case III  (0,2)
% E = 2*rand(n,n);


E = triu(E) + triu(E,1)';
for i=1:n
    E(i,i) =1;
end


 alpha =  1;
 %G = (1-alpha) * G + alpha * E;
 
 G = E; 

 


%%                               
lh =   10 ;       % number of fixed off diagonal elements in each row
ll =    50;       % number of off diagonal elements of lower bounds in each row
lu =   100;      % number of off diagonal elements of upper bounds in each row

ll = min(ll,n-1); 
lu = min(lu,n-1);
lh = min(lh,n-1);


%% I_e,J_e
%%%% for fixed  diagonal entries
I_d = [1:1:n]';
J_d = I_d;

%%%% for fixed off-diagonal entries
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


I_b = [I_d;I_h];
J_b = [J_d;J_h];
k_b = length(I_b);




%%  I_l,J_l 
%%%  entries with lower bounds
I_l = [];
J_l = [];
for i = 1:n-ll
    r = rand(n-i,1);
    [r,ind] = sort(r);
    I_l = [I_l; i*ones(ll,1)];
    J_l = [J_l; i+ind(n-i-ll+1:n-i)];
end
for i = ((n-ll)+1):(n-1)
    I_l = [I_l; i*ones(n-i,1)];
    J_l = [J_l;[(i+1):n]'];
end
k_l = length(I_l);




%%  I_u,J_u      
%%%%%  entries with upper bounds
I_u = [];
J_u = [];
for i = 1:n-lu
    r = rand(n-i,1);
    [r,ind] = sort(r);
    I_u = [I_u; i*ones(lu,1)];
    J_u = [J_u; i+ind(n-i-lu+1:n-i)];
end
for i = ((n-lu)+1):(n-1)
    I_u = [I_u; i*ones(n-i,1)];
    J_u = [J_u;[(i+1):n]'];
end
k_u = length(I_u) ;







%% to generate the bound b,l & u
%%%%%%% b
rhs = ones(n,1);  % diagonal elements

alpha0 = .1 ;
rhs = alpha0 * rhs + (1-alpha0) * rand(n,1);

h = zeros(k_h,1);
b = [rhs;h];

%%%%%%% l
l = -0.10*ones( k_l,1);
%l = 0.50 * (2*rand(k_l,1)-ones(k_l,1));
%l = 1.0 * (rand(k_l,1) - ones(k_l,1));

%%%%%%% u
u = 0.10* ones(k_u,1);
%u = 1.0 * (rand(k_l,1) - ones(k_l,1));





max_l = max(l);
min_l = min(l);
max_u = max(u);
min_u = min(u);

 
 








[X,z_b,z_l,z_u] = CaliMat1(G,b,I_b,J_b,l,I_l,J_l,u,I_u,J_u,OPTIONS);


disp('=====end========end==========end===========end=============end==========end==========')




 


    

 
 
 
 
 
 



 

 

 


 






 
