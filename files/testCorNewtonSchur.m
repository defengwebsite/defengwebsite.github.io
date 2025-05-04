%
clear all
tau0=1.0e-4;     % to test if the smallest eigenvalue of G(1:m,1:m) is at least tau0
tau = 0.5*tau0;  % tau must be smaller than tau0
tol = 1.0e-6;    % stopping criterion

% % %%%%%%%%%%%%%%%% Generate a random correlation matrix with given eigenvalues%%%%%%%%%%%%%%%%
    n = 300;
    m = 100;
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


% % Correlation matrix used in RiskMetrics
% 
% load x.mat
% G=extract(x);
% [n, n_c]=size(G);
% m=10;

C=G;

% Generate a random stress matrix
stress_matrix = 2.0*rand(n,n)-ones(n,n);
stress_matrix = 0.5*(stress_matrix + stress_matrix');
for i=1:n
    stress_matrix(i,i) =1;
end
stress_matrix(1:m,1:n) = C(1:m,1:n);
stress_matrix((m+1):n,1:m) = C((m+1):n,1:m);

% %%%%%%%%%%%%%%%% Produce a perturbed correlation matrix with C_1  %%%%%%%%%%%%%%%%
 alpha_purt =.10; 
 C =(1-alpha_purt)*C+ (alpha_purt)*stress_matrix; 
 C =(C+C')/2;
% %%%%%%%%%%%%%%%%%%%%%%%

%% Both C(1:m,1:n) and C(1:n,1:m) are fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%
% G: the unstresed correlation matrix (n by n)
% C: stressed correlation matrix (n by n)
% C(1:m,1:m) is fixed
% tau0: the threshold for positive definiteness of  C  
% tau: the lower bound for the smallest eigevalue of X 
% tol: stopping crierion for the KKT system: 1.0e-5 ~ 1.0e-6
% X: the calibrated correlation matrix
% val_obj: final objective function value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,val_obj] = CorNewtonSchur(G,C,m,n,tau,tau0,tol);



