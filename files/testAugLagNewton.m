clear all

%% initial of the problem 

tau0 = 1.0e-4;
tau = .50e-4;         % tau0>tau to make sure the existence of interior points
TOL1 = 1.0e-6;        % tolerance of |KKT|<= TOL1

% %%%%%%%%%%%%%%%% Generate a random correlation matrix with given eigenvalues%%%%%%%%%%%%%%%%
    n = 100;
    m = 90;
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

    %C =gallery('randcorr',n);
    C =gallery('randcorr',d); %generate a correlation matrix with given eigenvalues

    C =(C+C')/2;

% Correlation matrix used in RiskMetrics
% 
%  load x.mat
%  C=extract(x);
%  C =(C+C')/2;
%  [n, n_c]=size(C);
%  m =300; % m<=386





% Generate a random stress matrix
stress_matrix = 2.0*rand(n,n)-ones(n,n);
stress_matrix = 0.5*(stress_matrix + stress_matrix');
for i=1:n
    stress_matrix(i,i) =1;
end
% Keep the first m by m principal submatrix  the same as that of C
stress_matrix(1:m, 1:m)=C(1:m, 1:m);
%%%

% %%%%%%%%%%%%%%%% Produce a perturbed correlation matrix with C_1  %%%%%%%%%%%%%%%%
 alpha_purt =.10; 
 C =(1-alpha_purt)*C+ (alpha_purt)*stress_matrix; 
 C =(C+C')/2;
% %%%%%%%%%%%%%%%%%%%%%%%
% C: stressed correlation matrix (n by n)
% C(1:m,1:m) is fixed
% tau0: the threshold for positive definiteness of  C  
% tau: the lower bound for the smallest eigevalue of X0 
% TOL1: stopping crierion for the KKT system: 1.0e-5 ~ 1.0e-6
% X0: the calibrated correlation matrix
% val_obj: final objective function value

[X0, val_obj]=AugLagNewton(C,m,n,tau,tau0,TOL1);





