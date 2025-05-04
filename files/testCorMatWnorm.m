clear all

%% initial of the problem 

 
%% Case I 
%%%%%% RiskMetrics Data (n =387)
% load x.mat
% G = subtract(x);
% G  = (G +G')/2;
% [n, n_c]=size(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Case II
% % % % %%%%%%%%%%%%%%%%%%%%%
n = 500;
x = 10.^[-4:4/(n-1):0];
G =gallery('randcorr',n*x/sum(x));
% %%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n
    G(i,i) =1;
end




%% Case I
%E =2.0*rand(n,n)-ones(n,n);
 
E = randn(n,n);

 E = triu(E) + triu(E,1)';
 E =(E+E')/2;


 
%G =2.0*rand(n,n); %Case II
 
  


alpha = .2;
G =(1-alpha)*G+ alpha*E;
% norm_E = norm(E, 'fro')
 
% for i=1:n
%     G(i,i) =1;
% end



  
  %[X,y] = CorNewton2(G,rhs,0);
  
tau =1.0e-8;
B = rand(n,n);
W = B'*B;
W = (W+W')/2;
W = n^0.5 * W / norm(W,'fro'); % scaled
W = (W + W')/2 + 1.0e-1*eye(n);

  w =rand(n,1);
  W = diag(w);



b = ones(n,1);
    
 fprintf('\n')
 disp('=====end========end==========end===========end=============end==========end==========')
 
 [X,y] =  CorMatWnorm(G,W,b,tau);
 
 %   [X,y] = CorNewton_Wnorm(G,W,b,tau);
    
%      [X,y] = Correlation_Newton_Tian(G,W);
%   
%    [X,y] = Correlation_Newton_Diag(G,W,b,tau);
       %[X,y] = Correlation_Newton_Chol(G,W,b,tau); 
%    [X1,X,y] = CNewton_W(G,W,tau);
%    
%    [X,y] = CorNewton_Cai(G,W,b,tau);
   % [Xw,y] = CorNewton_Chen(G,W,b,tau);
      %[X,y] =  C_Newton_Dong(G,W,b,tau);
      %[X,y] = CorNewton5_Goh(G,W); 
     % [X,y] = WeightedCorNewton_Jin(G,W);
     %  [X,y] = Correlation_Newton_Miao(G,W);
   
     
    
    

   
   % [X,y] = Correlation_Newton_Luo(G,b,W,tau);
   %[X,y] = CorNewton_Wu(G,W);
   
  %  [X,y] = Correlation_Newton_Shen(G,W,b,tau);
  %  [X,y] = Correlation_Newton_Yang(G,W);
    
 % [X,y] = CorMatrix_Wang(G,tau,W);
   
   
   
 % CorMat_Jiang(G,W,b,tau);
   
 
   %[X,y] = CorNewton_Chan(G,W);
   %[X,y] = CorNewton_Zhang(G,W);
   % [X,y] = CorNewton_Song3(G,W);
   
      
  fprintf('\n')
 
  
  %[X,z_e,z_l,z_u] = CaliMat0(G,b,I_b,J_b,l,I_l,J_l,u,I_u,J_u);
 
 