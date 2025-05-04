clear all
error_tol=1.0e-6;%1.0e-7; % termination tolerance (relative error)
tol =1.0e-6; %1.0e-6% Tolerance for PCG
%%%%%% Weight for Obj. 
c1=1e0;
c2=1e-0;

%%%%%%%%%%%       Generate the prescribed eigendata    %%%%%%%%%%%%%%%

%%%% Generate k random  eigenvalues with k/2 complex-value ones %%%%%%%
n=1000;
k=30;
Lambda=zeros(k,k); 
alpha=-abs(randn(k,1));
s=floor(1/4*k); % Number of complex eigenvalues
beta=randn(s,1);
Lambda=diag(alpha);
i=1;
while i<=s
    Lambda(2*(i-1)+1:2*i,2*(i-1)+1:2*i)=[alpha(i) beta(i);-beta(i) alpha(i)];
    i=i+1;
end

%%%%   Generate the eigenvector matrix randomly  %%%
X=randn(n,k);
[Q,R]=qr(X);% QR factorization of the eigenvectors X
R=R(1:k,:); % R factor
i=1;
while i<=k
    if norm(R(i,i))<1.0e-6
        R(i,i) =1.0e-6;
    end
    i=i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Strictly feasible solution given by Bai, Chu and Sun      %%%%%%
%%%%% See Example 5.1 in our paper 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Produce a strictly feasible solution with the given eigendata  %%%%
Rinv=inv(R);

Mtrue=[Rinv'*Rinv zeros(k,n-k);zeros(n-k,k) eye(n-k)];
Mtrue=c1^0.5*Mtrue; 

Ctrue=-[Rinv'*(Lambda+Lambda')*Rinv zeros(k,n-k);zeros(n-k,n)];
Ctrue=c2^0.5*Ctrue;

Ktrue=Lambda*Rinv; Ktrue=[Ktrue'*Ktrue zeros(k,n-k);zeros(n-k,k) eye(n-k)];    

W1=2*rand(n)-ones(n,n);W1=(W1'+W1)/2;
W2=2*rand(n)-ones(n,n);W2=(W2'+W2)/2;
W3=2*rand(n)-ones(n,n);W3=(W3'+W3)/2;

%%%% Generate a perturbed solution (Ma,Ca,Ka) as the analytic solution %%%%
tau=0.1; 
Ma=Mtrue+tau*W1; % generate the given symmetric Ma, Ca,and Ka
Ca=Ctrue+tau*W2;
Ka=Ktrue+tau*W3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%
% Ma: the analytic mass matrix (n by n)
% Ca: the analytic damping matrix (n by n)
% Ka: the analytic stiffness matrix (n by n)
% R: R factor of the given eigenvector matrix X (n by k)
% Lambda: the given eigenvalue matrix (k by k)
% c1 and c2: weighted parameters
% error_tol: stopping crierion for the method: 1.0e-7
% tol:  Tolerance for PCG
% M,C,K: final solution
% val_obj: final objective function value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; [M,C,K,y,z]= IQEP_Newton(Ma,Ca,Ka,R,Lambda,n,k,c1,c2);toc

 tic; [M1,C1,K1,y1,z1]= IQEP_Newton_old(Ma,Ca,Ka,R,Lambda,n,k,c1,c2);toc; 
%tic; [M,C,K,y,z]= IQEP_Newton_old(Ma,Ca,Ka,R,Lambda,n,k,c1,c2);toc

%tic; [M,C,K,y,z]= IQEP_Newton(Ma,Ca,Ka,R,Lambda,n,k,c1,c2);toc