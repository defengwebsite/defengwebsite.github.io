%test
clear all
tau0=1.0e-4;    % to test if the smallest eigenvalue of G(1:m,1:m) is at least tau0
tau = 0.5*tau0; %tau must be smaller than tau0
tol = 1.0e-6; %stopping criterion

% % %%%%%%%%%%%%%%%% Generate a random correlation matrix with given eigenvalues%%%%%%%%%%%%%%%%
    n = 250;
    m = 51;
    n1 = round(n/10);
    d1 = rand(n1,1)/n;  %% To get larger eigenvalues 
    k =min(10, n-m);
    d3 = n*rand(k,1);

    d = rand(n,1);

    d(1:n1) =d1;

    d(n-k+1:n,1)=d3;  %% To get larger eigenvalues 


    %sum_d =sum(d)

    d = n*d/sum(d); % sum(d) =n;

    min_d = min(d)
    max_d = max(d)

    %G =gallery('randcorr',n);
   G =gallery('randcorr',d); %generate a correlation matrix with given eigenvalues

    G =(G+G')/2;

% Correlation matrix used in RiskMetrics

% load x.mat
% G=extract(x);
% [n, n_c]=size(G);
% m=51% G(1:m,1:m) fixed

% Generate a random stress matrix
stress_matrix = 2.0*rand(n,n)-ones(n,n);
stress_matrix = 0.5*(stress_matrix + stress_matrix');
for i=1:n
    stress_matrix(i,i) =1;
end
stress_matrix(1:m,1:m) = G(1:m,1:m);

% %%%%%%%%%%%%%%%% Produce a perturbed correlation matrix with C_1  %%%%%%%%%%%%%%%%
 alpha_purt =.10; 
 G =(1-alpha_purt)*G+ (alpha_purt)*stress_matrix; 
% G(1:m,1:m) =(1-alpha_purt)*G(1:m,1:m)+ (alpha_purt)*stress_matrix(1:m,1:m); 
 G =(G+G')/2;
% %%%%%%%%%%%%%%%%%%%%%%%

%% G(1:m,1:m) is fixed
% %%%%%%%%%%%%%%%%%%%%%%%
% G: stressed correlation matrix (n by n)
% G(1:m,1:m) is fixed
% tau0: the threshold for positive definiteness of  G  
% tau: the lower bound for the smallest eigevalue of X  
% tol: stopping crierion for the KKT system: 1.0e-5 ~ 1.0e-6
% X: the calibrated correlation matrix
% y and z: dual solutions
% obj_value: final objective function value


[X,y,z,obj_value] = LagDualNewton(G,m,n,tau,tau0,tol);


