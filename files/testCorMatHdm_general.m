clear all

%% initial of the problem 

 

TOL1 = 5.0e-6;        % tolerance of |KKT|<= TOL1 : 1.0e-5 ~ 1.0e-6
tau =  .0e-8;         % tau0 to make sure X-tau*I is positive semidefinite



% % % % % % %%%%%% RiskMetrics Data (n =387)
% load x.mat
% G =subtract(x);
% G  = (G +G')/2;
% [n, n_c]=size(G);
%%%%%%%   %%%%%%%%%%%%%%%%%%%



 


% % % %%%%%%%%%%%%%%%%%%%%%
n = 1000;
x = 10.^[-4:4/(n-1):0];
G =gallery('randcorr',n*x/sum(x));
% %%%%%%%%%%%%%%%%%%%%%%%%%%


%G0 =2.0*rand(n,n)-ones(n,n); %Case I
%G0 =2.0*rand(n,n); %Case II
%G0 =(G0+G0')/2;

%for i=1:n
%    G(i,i) =1;
%end
% H = rand(n,n);
% for i =1:n-1
%     H(i+1:n,i) = H(i, i+1:n);  % so H is likely to have small numbers
% end
% H = (H + H')/2;
% H = 5.0*H;
% H = 0.1*ones(n)+9.9*H;
% H = 0.01*ones(n)+99.99*H;



H0 =sprand(n,n,0.5);
H0 = triu(H0) + triu(H0,1)'; % W0 is likely to have small numbers
H0 = (H0 + H0')/2;

H1 = rand(n,n);
H1 = triu(H1) + triu(H1,1)'; % W1 is likely to have small numbers
H1 = (H1 + H1')/2;
%H =5.0*H1;
%H = .1*ones(n,n) + 9.9*H1;
%H = .01*ones(n,n) + 99.99*H1;





%H0 =H1;
%H0 = 5.0*H1;


%H0 =max(.1,H1) + 9*H0; %%% H0 is in [0.1, 10]
H0 =0.01*ones(n,n) + 99.99*H0; %%% H0 is in [0.01, 100]
%H =2*H1;
%H = ones(n,n)+ 0.001*H1;
H =0.1*ones(n,n)+9.9*H1; %%% H is in [.1,10]
%H =0.01*ones(n,n)+99.99*H1;


%H =1*ones(n,n)+0.1*H1; %%% H is around ones
%H =ones(n,n);

% % %%%%%%%%%%%%%%%%%%%% Assign weights H0  on partial elemens 
s =sprand(n,1,min(10/n,1));
I = find(s>0);
d =sprand(n,1,min(10/n,1));
J = find(d>0);
if length(I) >0 & length(J)>0
 H(I,J) = H0(I,J);
 H(J,I) = H0(J,I);
end
H = (H + H')/2;
%%%%%%%%%%%%% end of  assignings weights from one only on partial elemens 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%H =ones(n,n);


Hmin =min(min(H))
Hmax =max(max(H))

% d =10000*rand(n,1);  Diagonal weight only
%H =ones(n,n);
% for i=1:n
%    H(i,i) =d(i);
% end

E =2.0*rand(n,n)-ones(n,n); %Case I
%E = 2.0*rand(n,n);
 %E =sprand(n,n,100/n^2);
E = triu(E) + triu(E,1)';
 
   for i=1:n
      E(i,i) =1;
   end

alpha =.1;
 G =(1-alpha)*G+ alpha*E;
 
G =(G+G')/2;
G = min(G,ones(n,n));
G = max(-ones(n,n),G);

%G =2*rand(n,n)-ones(n,n);
%G =(G+G')/2;

for i=1:n
    G(i,i) =1;
end




%%                               
lh =       10 ;       % number of fixed off diagonal elements in each row
ll =       50;       % number of off diagonal elements of lower bounds in each row
lu =     50;      % number of off diagonal elements of upper bounds in each row

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

alpha0 =  1 ;
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

 

 

   
    
 
   [X0,y0,Gamm0,z_l0,z_u0,val_obj] = CorMatHdm_general(G,H,b,I_b,J_b,l,I_l,J_l,u,I_u,J_u,tau,TOL1);
  % [X0,y0,Gamm0,z_l0,z_u0,val_obj] = CorMatHdmMex_general(G,H,b,I_b,J_b,l,I_l,J_l,u,I_u,J_u,tau,TOL1);
 
 
 
   disp('=====end========end==========end===========end=============end==========end==========')
 
 
 
 