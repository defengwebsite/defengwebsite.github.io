%% initial of the problem 

 

TOL1 = 5.0e-6;        % tolerance of |KKT|<= TOL1 : 1.0e-5 ~ 1.0e-6
tau = 1.0e-8;         % tau0 to make sure X-tau*I is positive semidefinite


% % % %%%%%% RiskMetrics Data (n =387)
%load x.mat
%G=subtract(x);
%G = (G+G')/2;
%[n, n_c]=size(G);
%%%%%%%   %%%%%%%%%%%%%%%%%%%


n =  20;

%G =2.0*rand(n,n)-ones(n,n); %Case I
%    G =2.0*rand(n,n); %Case II
% G = triu(G) + triu(G,1)';

 % % %%%%% generate a random G  
 x = 10.^[-4:4/(n-1):0];
 G =gallery('randcorr',n*x/sum(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%  G =(G+G')/2;
 
% W = rand(n,n);
% for i =1:n-1
%     W(i+1:n,i) = W(i, i+1:n);  % so W is likely to have small numbers
% end
% W = (W+W')/2;
% W =5.0*W;
% W =0.1*ones(n)+9.9*W;
% W =0.01*ones(n)+99.99*W;



W0 = sprand(n,n,0.5);
W0 = triu(W0) + triu(W0,1)'; % W0 is likely to have small numbers
W0 = (W0+W0')/2;

W1 = rand(n,n);
W1 = triu(W1) + triu(W1,1)'; % W1 is likely to have small numbers
W1 = (W1 + W1')/2;
%W =W1;
%W = .1*ones(n,n) + 9.9*W1;
%W = .01*ones(n,n) + 99.99*W1;





%W0 =W1;
%W0 = 5.0*W1;


%W0 =max(.1,W1) + 9*W0; %%% W0 is in [0.1, 10]
W0 = 0.01*ones(n,n) + 99.99*W0; %%% W0 is in [0.01, 100]
%W =2*W1;
W = 0.1*ones(n,n)+9.9*W1; %%% W is in [.1,10]
%W = ones(n,n)+ 0.001*W1; %%% W is in [.1,10]
%W =0.01*ones(n,n)+99.99*W1;


%W =1*ones(n,n)+0.1*W1; %%% W is around ones
%W =ones(n,n);

% % %%%%%%%%%%%%%%%%%%%% Assign weights W0  on partial elemens 
s =sprand(n,1,min(10/n,1));
I = find(s>0);
d =sprand(n,1,min(10/n,1));
J = find(d>0);
if length(I) >0 & length(J)>0
 W(I,J) = W0(I,J);
 W(J,I) = W0(J,I);
end
W = (W+W')/2;
%%%%%%%%%%%%% end of  assignings weights from one only on partial elemens 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%W =ones(n,n)+.5*E;
Wmin = min(min(W))
Wmax = max(max(W))

% d =10000*rand(n,1);  Diagonal weight only
% W =ones(n,n);
% for i=1:n
%    W(i,i) =d(i);
% end


%E = 2.0*rand(n,n) - ones(n,n); %Case I
 E = 2.0*rand(n,n);
 %E =sprand(n,n,100/n^2);
E = triu(E) + triu(E,1)';
%E =(E + E')/2;
 for i=1:n
    E(i,i) =1;
 end

%
 %E =rand(n,n);
 %E = (E +E')/2;
%E =E/(sum(sum((E.*E))))^0.5;
%for i=1:n
%    E(i,i) =1;
%end

alpha =  .1;
G =(1-alpha)*G+ alpha*E;
G =(G + G')/2;
%G = min(G, ones(n,n));
%G = max(-ones(n,n),G);

% G =2*rand(n,n)-ones(n,n);
% G = -(G+G')/2;

for i=1:n
    G(i,i) =1;
end
% G =2*rand(n,n);
% G =(G+G')/2;
% W =ones(n,n);

%norm_E = norm(E, 'fro')

% norm_G = norm(G, 'fro');

 
rhs =ones(n,1);
[X0,y0,Gamma0, val_obj] = CorMatHdmMex(G,W,rhs,tau,TOL1); 
%[X0,y0,Gamma0, val_obj] = CorMatHdm(G,W,rhs,tau,TOL1); 
% [X0,y0,Gamma0, val_obj]= CorMatHdmMex_Sep9(G,W,rhs,tau,TOL1); 
 %[X0,y0,Gamma0, val_obj]= CorMatHdmMex_March5_08(G,W,rhs,tau,TOL1); 
  [X0,y0,Gamma0, val_obj]= CorMatHdm_March5_08(G,W,rhs,tau,TOL1); 
%b =ones(n,1);
%[X,y] = CorNewton1(G,b,tau);
%[X,y] = CorNewton2(G,b,tau);




