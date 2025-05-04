%**************************************************************
%%%%%% This code is designed to solve   %%%%%%%%%%%%%%%%%%
%%%%%% min 0.5*||M-Ma||^2 + 0.5*||C-Ca||^2 + 0.5*||K-Ka||^2
%%%%%% s.t 1/sqrt(c1)*M*[R;zeros(n-k,k)]*Lambda^2 + 1/sqrt(c2)*C*[R;zeros(n-k,k)]*Lambda + K*[R;zeros(n-k,k)] = 0
%%%%%%      M>=0 (PSD), C=C^T, K>=0 (PSD) (PSD: Symmetric and positive semi-definite)
%%%%%%
%%%%%% based on the algorithm in %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% "A Dual Optimization Approach to Inverse Quadratic Eigenvalue %%%%%
%%%%%% Problems with Partial Eigenstructure"  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% By Zheng-Jian Bai, Delin Chu, and Defeng Sun    %%%%%%%%%%%%%%%%%%%
%%%%%% SIAM J. Sci. Comput. 29 (2007), pp. 2531-2561.
%%%%%%
%%%%%% Last Modified Date: July 15, 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%         by Ying Cui and Defeng Sun       %%%%%%%%%%%%%%%%%%%%%%
%%%%%% A new block preconditioner  inv(A1'A1) is added to the first part
%%%%%% of the newton system
%%%%%% The  input argument is the given symmetric Ma, Ca,and Ka,  %%%%%%%%
%%%%%% the eigendata (Lambda,R), where R is the R-factor of QR
%%%%%% factoriztion of the matrix X of k given eigenvectors, the
%%%%%% problem dimension n and the number k of given eigenpairs, and the
%%%%%% positive weighting parameters c1 and c2
%%%%%% The output are the optimal primal and dual solutions %%%%%%%%
%%%%%%%
%%%%%%% Send your comments and suggestions to    %%%%%%
%%%%%%% zjbai@xmu.edu.cn, matchudl@nus.edu.sg  or matsundf@nus.edu.sg %%%%
%%%%%
%%%%% Warning: Accuracy may not be guaranteed!!!!! %%%%%%%%
function [M,C,K,y,z]=IQEP_Newton(Ma,Ca,Ka,X,Lambda,n,k,c1,c2,error_tol)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Newton Method Starts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

t0=cputime;

tstart = clock;
m=0;
f_eval =0;
Iter_Whole=500;
Iter_inner =20; % Maximum number of Line Search in Newton method
maxit0 =min(10000,n*k); % Maximum number of iterations in PCG at the 1st step
maxit =min(10000,n*k); % Maximum number of iterations in PCG
%Inner =0;
tol_cg = 5.0e-3;%1.0e-6; % Initial tolerance for PCG
tol_cg_max = 5.0e-6;
%error_tol=1.0e-7; %1.0e-6;%1.0e-7; % termination tolerance (relative error)
sigma_1=1.0e-4; %tolerance in the line search of the Newton method
regterm_const = 1.0e-5;
precond_num =1;

[Q,R]=qr(X);% QR factorization of the eigenvectors X
Q1 = Q(:,1:k);
Q2 = Q(:,(k+1):n);

R=R(1:k,:);

i=1;
while i<=k
    if norm(R(i,i))<1.0e-6
        R(i,i) =1.0e-6;
    end
    i=i+1;
end

Ma = sqrt(c1)*Q'*Ma*Q;  Ca = sqrt(c2)*Q'*Ca*Q; Ka = Q'*Ka*Q;

S=R*Lambda/R; S2=S*S;

%%%%%% Preconditioners

precond2=c1^(-1)*S2'*S2 + c2^(-1)*S'*S + eye(k);
precond2 = precond2/2;   %%% preconditioner for A2PA2'


if k <= 150              %%% preconditioner for A1PA1'
    precond.choice = 1;
    A1A1T = AAT1(S2,S,c1,c2,k);
    A = A1A1T + regterm_const*eye(k^2);
    precond.A = inv(A);
else
    precond.choice = 2;     %%% diagonal preconditioner for k >150
    precond.A = 2*precond2;
end



run_rel_error =0;

if run_rel_error ==0;
    mck0=1;
else
    mck0 = max(sqrt((norm(Ma,'fro')/sqrt(c1))^2 + (norm(Ca,'fro')/sqrt(c2))^2+(norm(Ka,'fro'))^2),1);
end


Res_b=zeros(300,1);

%%%%%% The denominator used in the stopping criterion %%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b01 =-zeros(k,k);
b02=zeros(k,n-k);
y=zeros(k,k);
z=zeros(k,n-k);
F1=zeros(k,k);
F2=zeros(k,n-k);

%%%%%% Initial guess (It should find optimal solutions with equations only)
y02=b01-c1^(-0.5)*Ma(1:k,1:k)*S2-c2^(-0.5)*Ca(1:k,1:k)*S-Ka(1:k,1:k);
if precond.choice == 1
    y02 = reshape(y02,k^2,1);
    
    y = precond.A*y02;
    %y = linsysolve(precond,y02);
    y = reshape(y,k,k);
    
    fprintf('    Initial step:           ');
else
    [y,flag,relres,iterk] = pre_cg0(y02,tol_cg,maxit0,eye(k),precond.A,S2,S,k,c1,c2); % initial y
    fprintf('Initial step: CGs = %5.0d,  ',  iterk)
end
%fprintf('---------Number of CG Iterations for the initial guess === %d \n', iterk)
z02=b02-c1^(-0.5)*S2'*Ma(1:k,k+1:n)-c2^(-0.5)*S'*Ca(1:k,k+1:n)-Ka(1:k,k+1:n);
z=precond2\z02; % initial z


M =Ma;C =Ca; K =Ka;
[A11,A12,A21,A22,A31,A32]=Adjoint(y,z,S2,S,c1,c2); % ajoint of A
M(1:k,:) = M(1:k,:) + [A11 A12];
M(k+1:n,1:k) = M(k+1:n,1:k) + A12';
M  =(M+M')/2; % New M

C(1:k,:) = C(1:k,:) + [A21 A22];
C(k+1:n,1:k) =C(k+1:n,1:k)+  A22';
C=(C+C')/2; % New C

K(1:k,:) = K(1:k,:) + [A31 A32];
K(k+1:n,1:k) = K(k+1:n,1:k) + A32';
K=(K+K')/2; % New K
eig_total = 0;
tstart_eig = clock;
[lam1,P1,lam3,P3]=ED(M,K);% Eigendecomposition of M,K
eig_total = eig_total + etime(clock,tstart_eig);

[f0,F1,F2]=gradient(C,lam1,P1,lam3,P3,n,k,S2,S,c1,c2);
const_matrix =0.5*(sum(sum(Ma.*Ma))+sum(sum(Ca.*Ca))+sum(sum(Ka.*Ka)));
f0=f0-sum(sum(b01.*y))- sum(sum(b02.*z)) -const_matrix; %the objective function value
fprintf(' fval = %- 12.8e,  ',-f0)
f_eval =f_eval + 1;
b1=b01-F1;
b2=b02-F2;
norm_b=sqrt(norm(b1,'fro')^2 + norm(b2,'fro')^2);
rel_b=norm_b/mck0;
fprintf(' Grad = %3.2e, ',norm_b)
if run_rel_error ==1;
    fprintf(' RelGrad = %3.2e, ',rel_b)
end
fprintf(' Eig Time = %5.1f ',eig_total);
fprintf(' Clock Time = %5.1f, CPU Time = %5.1f, \n',etime(clock,tstart), cputime-t0);%[y,flag,relres,iterk] = pre_cg0(y02,tol_cg,maxit0,precond,S2,S,k,c1,c2); % initial y

Omega1 = omega_matrix(lam1,n);
Omega3 = omega_matrix(lam3,n);

m1 =size(Omega1,1); m3 =size(Omega3,1);
QQ1.ll =P1(1:k,1:m1); QQ1.rl=P1(k+1:n,1:m1);
QQ1.lr = P1(1:k,m1+1:n); QQ1.rr=P1(k+1:n,m1+1:n);
QQ3.ll =P3(1:k,1:m3);  QQ3.rl=P3(k+1:n,1:m3);
QQ3.lr = P3(1:k,m3+1:n); QQ3.rr=P3(k+1:n,m3+1:n);
if m1>=n/2 && m1 < n
    Omega1 = ones(m1,n-m1) - Omega1;
end
if m3 >= n/2 && m3 <n
    Omega3 = ones(m3,n-m3) - Omega3;
end

x01=y;
x02=z;
d1=zeros(k,k);
d2=zeros(k,n-k);
while (rel_b>error_tol && m<= Iter_Whole)
    
    regterm = regterm_const*norm_b;
    
    if m>5 && mod(m,5)==1
        tol_cg = max(tol_cg/2,tol_cg_max);
        
    end
    
    
    
    rank_ave = (m1 + m3)/2;
    iter_num= max(30,rank_ave*(n-rank_ave)*k^2/n^2);
    
    if  m >= 1 && precond.choice ==1 &&  iterk > iter_num && norm_b > 1.0e-5
        
        precond_num = precond_num + 1;
        fprintf('Generate  new  preconditioner ...')
        precond_time = clock;
        A = AVAT1(Omega1,Omega3,QQ1,QQ3,S2,S,c1,c2,n,k,m1,m3,regterm);
        precond.A = inv(A);
        fprintf('time used =  %5.1f, \n',etime(clock,precond_time))
    end
    [d1,d2,flag,relres,iterk]=pre_cg(b1,b2,tol_cg,maxit,precond,precond2,Omega1,P1,QQ1,Omega3,P3,QQ3,n,k,S2,S,c1,c2,regterm,m1,m3);
    fprintf(' Iter = %3.0d,  CGs = %5.0d,  ',m+1,iterk);
    
    if (flag~=0) % flag =0 means CG successful, otherwise not successful
        d1 = -(F1-b01);
        d2 = -(F2-b02);
    end % end loop for if
    slope =sum(sum((F1-b01).*d1))+sum(sum((F2-b02).*d2)); %%% nabla f d
    y =x01+d1; %temporary x0+d
    z =x02+d2;
    M =Ma; C= Ca; K= Ka;
    [A11,A12,A21,A22,A31,A32]=Adjoint(y,z,S2,S,c1,c2); % Ajoint of A
    M(1:k,:) = M(1:k,:) + [A11 A12];
    M(k+1:n,1:k) = M(k+1:n,1:k) + A12';
    M  =(M + M')/2; % New M
    
    C(1:k,:) = C(1:k,:) + [A21 A22];
    C(k+1:n,1:k) = C(k+1:n,1:k)+  A22';
    C=(C + C')/2; % New C
    
    K(1:k,:) = K(1:k,:) + [A31 A32];
    K(k+1:n,1:k) = K(k+1:n,1:k) + A32';
    K=(K+K')/2; % New K
    tstart_eig = clock;
    [lam1,P1,lam3,P3]= ED(M,K); % Eigendecomposition of M,K
    
    eig_total = eig_total + etime(clock,tstart_eig);
    
    
    [f,F1,F2] =gradient(C,lam1,P1,lam3,P3,n,k,S2,S,c1,c2);
    f=f-sum(sum(b01.*y))- sum(sum(b02.*z)) -const_matrix; %objective function value
    k_inner=0;
    while(k_inner <=Iter_inner & f>f0+sigma_1*0.5^k_inner*slope+1.0e-8)
        k_inner=k_inner+1;
        y = x01+0.5^k_inner*d1; % backtracking line search
        z = x02+0.5^k_inner*d2;
        
        M =Ma; C= Ca; K=Ka;
        [A11,A12,A21,A22,A31,A32] =Adjoint(y,z,S2,S,c1,c2); % Ajoint of A
        M(1:k,:) = M(1:k,:) + [A11 A12];
        M(k+1:n,1:k) = M(k+1:n,1:k) + A12';
        M  =(M+M')/2; % New M
        
        C(1:k,:) = C(1:k,:) + [A21 A22];
        C(k+1:n,1:k) =C(k+1:n,1:k)+  A22';
        C=(C+C')/2; % New C
        
        K(1:k,:) = K(1:k,:) + [A31 A32];
        K(k+1:n,1:k) = K(k+1:n,1:k) + A32';
        K=(K+K')/2; % New K
        
        
        [lam1,P1,lam3,P3]=ED(M,K); % Eigen-Decomposition of M,K
        [f,F1,F2] =gradient(C,lam1,P1,lam3,P3,n,k,S2,S,c1,c2);
        f=f-sum(sum(b01.*y))- sum(sum(b02.*z)) -const_matrix;  %objective function value
    end % end loop for while
    % fprintf('---------Number of Inner Iterations at the current step::::: %d \n', k_inner)
    f0=f;
    fprintf(' fval = %- 12.8e,  ',-f)
    f_eval =f_eval+k_inner+1;
    x01=y;
    x02=z;
    m=m+1;
    b1=b01-F1;
    b2=b02-F2;
    norm_b=sqrt(norm(b1,'fro')^2+norm(b2,'fro')^2);
    rel_b=norm_b/mck0;
    fprintf(' Grad = %3.2e, ',norm_b)
    if run_rel_error ==1;
        fprintf(' RelGrad = %3.2e, ',rel_b)
    end
    fprintf(' Eig Time = %5.1f ',eig_total);
    fprintf(' Clock Time = %5.1f, CPU Time = %5.1f, \n',etime(clock,tstart), cputime-t0);
    Res_b(m)=norm_b;
    Omega1 = omega_matrix(lam1,n);
    Omega3 = omega_matrix(lam3,n);
    
    
    m1 =size(Omega1,1); m3 =size(Omega3,1);
    QQ1.ll =P1(1:k,1:m1); QQ1.rl=P1(k+1:n,1:m1);
    QQ1.lr = P1(1:k,m1+1:n); QQ1.rr=P1(k+1:n,m1+1:n);
    QQ3.ll =P3(1:k,1:m3);  QQ3.rl=P3(k+1:n,1:m3);
    QQ3.lr = P3(1:k,m3+1:n); QQ3.rr=P3(k+1:n,m3+1:n);
    
    if m1>=n/2 && m1 <n
        Omega1 = ones(m1,n-m1) - Omega1;
    end
    if m3 >= n/2 && m3 <n
        Omega3 = ones(m3,n-m3) - Omega3;
    end
end %end loop for while

fprintf('Newton: Number of Iterations %d \n', m)
fprintf('Newton: Number of Function Evaluations %d \n', f_eval)
%fprintf('---------Norm of Initial Gradient:::::: %d \n',res_initial)

M = 1/sqrt(c1)*Q*M*Q';
M = (M+M')/2;
C = 1/sqrt(c2)*Q*C*Q';
C = (C+C')/2;

K = Q*K*Q';
K = (K+K')/2;

time_used=cputime-t0;

fprintf('Newton: Eig Time = %3.1f, Clock Time = %3.1f, CPU Time = %3.1f , \n',eig_total, etime(clock,tstart),time_used)

%%% end of the main program

%%%%%%
%%%%%% To generate F(Y,Z) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%

function [f,F1,F2]= gradient(C,lam1,P1,lam3,P3,n,k,S2,S,c1,c2)
f=0.0;
F1 =zeros(k,k);
F2=zeros(k,n-k);
% use  eigendecomposition of M
H1=P1';
i=1;
while (i<=n)
    H1(i,:)=max(lam1(i),0)*H1(i,:);
    i=i+1;
end % end loop for while
H1=P1(1:k,:)*H1;

H2=C(1:k,:);

% use eigendecomposition  of K
H3=P3';
i=1;
while (i<=n)
    H3(i,:)=max(lam3(i),0)*H3(i,:);
    i=i+1;
end % end loop for while
H3=P3(1:k,:)*H3;


F1=c1^(-1/2)*H1(:,1:k)*S2+c2^(-1/2)*H2(:,1:k)*S+H3(:,1:k);% F1(Y,Z)
F2=c1^(-1/2)*S2'*H1(:,k+1:n)+c2^(-1/2)*S'*H2(:,k+1:n)+H3(:,k+1:n);% F2(Y,Z)

f=(norm(max(lam1,zeros(n,1))))^2;
f=f+(norm(C,'fro'))^2;
f=f+(norm(max(lam3,zeros(n,1))))^2;
f =0.5*f;  % f is the first part of the objective function

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% end of gradient.m %%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Form the adjoint of A(Y,Z): A^*(Y) %%%%%%%%%
%%%%%%

function [A11,A12,A21,A22,A31,A32]=Adjoint(Y,Z,S2,S,c1,c2)
A11=Y*S2'; A11=A11+A11';
A11 = (0.5*c1^(-0.5))*A11;
A12=S2*Z; A12 = (0.5*c1^(-0.5))* A12;


%A1=0.5*c1^(-0.5)*[H1 G1;G1' zeros(n-k)]; % A1^*(Y,Z)

A21=Y*S';A21=A21+A21';
A21 = (0.5*c2^(-0.5))*A21;
A22=(S*Z); A22 = (0.5*c2^(-0.5))*A22;

%A2=0.5*c2^(-0.5)*[H2 G2;G2' zeros(n-k)]; % A2^*(Y,Z)

A31=Y;A31=A31+A31';
A31  = 0.5*  A31 ;
A32  = 0.5* Z;
%A3=0.5*[H3 Z;Z' zeros(n-k)]; % A3^*(Y,Z)

return
% function [A1,A2,A3]=Adjoint(Y,Z,Lambda_square,Lambda,R,S2,S,c1,c2,n,k)
% A1= zeros(n,n); A2=A1; A3=A1;
% H1=R*(Lambda_square*Y)*R'; H1=H1+H1';
% H1 = (0.5*c1^(-0.5))*H1;
% G1=S2*Z; G1 = (0.5*c1^(-0.5))* G1;
% A1(1:k,:) =  [H1 G1];
% A1(k+1:n,1:k) =  G1';
%
% %A1=0.5*c1^(-0.5)*[H1 G1;G1' zeros(n-k)]; % A1^*(Y,Z)
%
% H2=R*(Lambda*Y)*R';H2=H2+H2';
% H2 = (0.5*c2^(-0.5))*H2;
% G2=(S*Z); G2 = (0.5*c2^(-0.5))* G2;
% A2(1:k,:) =[H2 G2];
% A2(k+1:n,1:k) =  G2';
% %A2=0.5*c2^(-0.5)*[H2 G2;G2' zeros(n-k)]; % A2^*(Y,Z)
%
% H3=R*Y*R';H3=H3+H3';
% A3(1:k,:) = 0.5* [H3 Z];
% A3(k+1:n,1:k) = 0.5* Z';
% %A3=0.5*[H3 Z;Z' zeros(n-k)]; % A3^*(Y,Z)
%
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% end of Adjoint.m %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   Eigenvalue Decomposition (ED) of M and K   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,lambda] = Myeig(X)
if issparse(X); X = full(X); end
[P,D]  = eig(X);
lambda = diag(D);

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





function [lam1,P1,lam3,P3]=ED(M,K)
[P1,lam1]=Myeig(M); % eig of M
[P3,lam3]=Myeig(K); % eig of K
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% end of ED.m




%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the essential part of the first -order difference of d
%%%%%%%
function Omega12 = omega_matrix(lambda,n)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PCG method %%%%%%%
%%%%%% This is exactly the algorithm by  Hestenes and Stiefel (1952)
%%%%%% An iterative method to solve A(x) =b
%%%%%% The symmetric positive definite matrix M is a preconditioner for A.
%%%%%%  See Pages 527 and 534 of Golub and va Loan (1996)

function [p1,p2,flag,relres,iterk]=pre_cg(b1,b2,tol,maxit,precond,precond2,Omega1,P1,Q1,Omega3,P3,Q3,n,k,S2,S,c1,c2,regterm,m1,m3);
% Initializations
r1 = b1;  %We take the initial guess x0=0 to save time in calculating A(x0)
r2=b2;
n2b1 =norm(b1,'fro');    % norm of b1
n2b2 =norm(b2,'fro');    % norm of b2
n2b=sqrt(n2b1^2+n2b2^2); % norm of b
tolb = max(tol *n2b,1.0e-12) ;
p1 = zeros(k,k);
p2=zeros(k,n-k);
flag=1;
iterk =0;
relres=1000; %%% To give a big value on relres
% Precondition
%z1 =(M1\r1)/M3;  %%%%% z = M\r; if M is not the identity matrix
if precond.choice == 1
    r1temp = reshape(r1,k^2,1);
    %z1 = linsysolve(L,r1temp);
    z1 = precond.A*r1temp;
    z1 = reshape(z1,k,k);
else
    z1 = r1/precond.A;
end
%precond2 = eye(k);
z2=precond2\r2;
rz1 = sum(sum(r1.*z1))+sum(sum(r2.*z2));
rz2 = 1;
d1 = z1;
d2=z2;
% CG iteration
m=1;
for m =1:maxit
    if m > 1
        beta = rz1/rz2;
        d1 = z1 + beta*d1;
        d2 = z2 + beta*d2;
    end
    
    [w1,w2]= Jacobian_matrix(d1,d2,Omega1,P1,Q1,Omega3,P3,Q3,n,k,S2,S,c1,c2,regterm,m1,m3); %w = A(d);
    %     %%%%% test %%%%%
    %   d2test = zeros(k,n-k);
    %   [w1test,w2test]= Jacobian_matrix(d1,d2test,Omega1,P1,Omega3,P3,n,k,S2,S,norm_b,Lambda2,Lambda,R,c1,c2); %w = A(d);
    %   A = AVAT1(Omega1,Omega3,P1,P3,S2,S,c1,c2,n,k,1e-4);
    %   d1test = reshape(d1,k^2,1);
    %   Ad1= A*d1test;
    %   Ad1 = reshape(Ad1,k,k);
    %   err = norm(Ad1-w1test,'fro');
    %    fprintf('err = %8.7e\n',err)
    %    keyboard;
    %%%%%%%%%%%%%%%%%%%%%
    denom =  sum(sum(d1.*w1))+sum(sum(d2.*w2));
    iterk =m;
    norm_z1=norm(z1,'fro');
    norm_z2=norm(z2,'fro');
    norm_z=sqrt(norm_z1^2+norm_z2^2);
    relres =norm_z/n2b;              %relative residue =norm(z) / norm(b)
    if denom <= 0
        norm_d1=norm(d1,'fro');
        norm_d2=norm(d2,'fro');
        norm_d=sqrt(norm_d1^2+norm_d2^2);
        p1 = d1/norm_d; % d is not a descent direction
        p2=d2/norm_d;
        break % exit
    else
        alpha = rz1/denom;
        p1 = p1 + alpha*d1;
        p2 = p2 + alpha*d2;
        r1 = r1 - alpha*w1;
        r2 = r2 - alpha*w2;
    end
    %z1 =(M1\r1)/M3; %  z = M\r; if M is not the identity matrix ;
    if precond.choice == 1
        r1temp = reshape(r1,k^2,1);
        %z1 = linsysolve(L,r1temp);
        z1 = precond.A*r1temp;
        z1 = reshape(z1,k,k);
    else
        z1 = r1/precond.A;
    end
    z2=precond2\r2;
    if norm_z <= tolb  % Exit if Hp=b solved within the relative tolerance
        slope = -sum(sum(b1.*p1)) - sum(sum(b2.*p2)); %%% nabla f d
        if slope <0
            iterk =m;
            relres =norm_z/n2b;          %relative residue =norm(z) / norm(b)
            flag =0;
            break
        end
    end
    rz2 = rz1;
    rz1 = sum(sum(r1.*z1))+sum(sum(r2.*z2));
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% end of pre_cg.m %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the Jacobain product with x: F'(y,z)(x) %%%%%%%
%%%%%%%

function [Ax1,Ax2]= Jacobian_matrix(x1,x2,Omega1,P1,Q1,Omega3,P3,Q3,n,k,S2,S,c1,c2,regterm,m1,m3)


%A12_new = zeros(n,n);
%---------------------------
[A_11,A_12,A_21,A_22,A_31,A_32]=Adjoint(x1,x2,S2,S,c1,c2); % Ajoint of A: A^*(Y)
% use ED of M




if m1==0
    A1 =sparse(n,n);
elseif m1 == n
    A1 = [A_11,A_12];
elseif m1 < n/2
    A1temp1 = Q1.rl'*A_12;
    A1temp2 = A1temp1 + Q1.ll'*A_11;
    
    A11 = Q1.ll'*A1temp1';
    A11 = A11 + A1temp2*Q1.ll;   %%%% compute P1'*Aty*P1
    
    A12 = Q1.ll'*A_12'*Q1.lr;
    A12 = A12 + A1temp2 *Q1.lr;    %%%% compute P1'*Aty*P2
    A12 = Omega1.*A12;
    
    A1 = [Q1.ll*A11 + Q1.rl*A12', Q1.ll*A12]*P1';  %%%% compute P1*(Omega.*(P1'*Aty*P1))*P';
    
else
    %Omega1 = ones(m1,n-m1) - Omega1;
    
    A1temp1 = A_12*Q1.rr;
    A1temp2 = A1temp1 + A_11*Q1.lr;
    
    A12 = Q1.ll'*A1temp2;
    A12 = A12 + (A_12*Q1.rl)'*Q1.lr;     %%%% compute P1'*Aty*P2
    % keyboard;
    A12 = Omega1.*A12;
    
    A13 = Q1.lr'* A1temp2;
    A13 = A13 + A1temp1'*Q1.lr;     %%%% compute P2'*Aty*P2
    
    A1 = [Q1.lr*A12', Q1.ll*A12 + Q1.lr*A13']*P1'; %%%% compute P1*((E - Omega).*(P1'*Aty*P1))*P';
    % A1try = Q1.lr*(Q1.ll*A12)' + (Q1.ll*A12 + Q1.lr*A13')*Q1.rl';
    A1 = [A_11,A_12] - A1;
    
    
end

A2=[A_21 A_22];%A2(1:k,:);




if m3==0
    A3 =sparse(n,n);
elseif m3 == n
    A3 = [A_31,A_32];
elseif m3 < n/2
    A3temp1 = Q3.rl'*A_32;
    A3temp2 = A3temp1 + Q3.ll'*A_31;
    
    A31 = Q3.ll'*A3temp1';
    A31 = A31 + A3temp2*Q3.ll;   %%%% compute P1'*Aty*P1
    
    A32 = Q3.ll'*A_32'*Q3.lr;
    A32 = A32 + A3temp2 *Q3.lr;    %%%% compute P1'*Aty*P2
    A32 = Omega3.*A32;
    
    A3 = [Q3.ll*A31 + Q3.rl*A32', Q3.ll*A32]*P3';  %%%% compute P1*(Omega.*(P1'*Aty*P1))*P';
    
else
    % Omega3 = ones(m3,n-m3) - Omega3;
    A3temp1 = A_32*Q3.rr;
    A3temp2  = A3temp1 + A_31*Q3.lr;
    
    A32 = Q3.ll'* A3temp2;
    A32 = A32 + (A_32*Q3.rl)'*Q3.lr;   %%%% compute P1'*Aty*P2
    A32 = Omega3.*A32;
    
    A33 = Q3.lr'* A3temp2;
    A33 = A33 + A3temp1'*Q3.lr;    %%%% compute P2'*Aty*P2
    
    A3 = [Q3.lr*A32', Q3.ll*A32 + Q3.lr*A33']*P3';   %%%% compute P1*((E - Omega).*(P1'*Aty*P1))*P';
    A3 = [A_31,A_32] - A3;
    
    
end


% A3 = Omega3.*A3;
% A3 =(P3_1*A3)*P3';
Ax1=c1^(-0.5)*A1(:,1:k)*S2+c2^(-0.5)*A2(:,1:k)*S+A3(:,1:k); % F1'(Y,Z)
Ax2=c1^(-0.5)*S2'*A1(:,k+1:n) + c2^(-0.5)*S'*A2(:,k+1:n)+A3(:,k+1:n);% F2'(Y,Z)



Ax1 =Ax1 +  regterm*x1;
Ax2 =Ax2 +  regterm*x2;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% end of Jacobian_matrix.m %%%%%%


function A = AAT1(S2,S,c1,c2,k)

A =zeros(k^2,k^2);
S2TS2 = S2'*S2; STS = S'*S;
for j=1:k
    
    A_11 = S2(:,j);
    A_12 = S2TS2(j,:);
    A_21 = S(:,j);
    A_22 = STS(j,:);
    
    for i = 1:k
        
        A1 = A_11*S2(i,:);
        A1(i,:) = A1(i,:) + A_12;
        A1 = A1/(2*c1);
        
        A2 = A_21*S(i,:);
        A2(i,:) = A2(i,:) + A_22;
        A2 = A2/(2*c2);
        
        
        A3 = zeros(k,k);
        A3(i,j) = 0.5;
        A3(j,i) = A3(j,i)+ 0.5;
        xx= A1 + A2 + A3;
        
        
        A(:,k*(j-1)+i) = reshape(xx,k^2,1);
        
    end
end
%A = A + sigma*eye(k^2,k^2);

function [p1,flag,relres,iterk] = pre_cg0(b1,tol,maxit,M1,M2,S2,S,k,c1,c2);
% Initializations
r1 = b1;  % we take the initial guess x0=0 to save time in calculating A(x0)
n2b1 =norm(b1,'fro');    % norm of b1
n2b=n2b1; % norm of b
tolb = tol;  % absolute tolerance
p1 = zeros(k,k);
flag=1;
iterk =0;
relres=1000; %%% To give a big value on relres
% Precondition
z1 = (M1\r1)/M2; %if M is not the identity matrix
rz1 = trace(r1'*z1);
rz2 = 1;
d1 = z1;
% CG iteration
m=1;
for m =1:maxit
    if m > 1
        beta = rz1/rz2;
        d1 = z1 + beta*d1;
    end
    w1=Jacobian_matrix0(d1,k,S2,S,c1,c2);%w = A(d);
    denom =  trace(d1'*w1);
    iterk =m;
    norm_z1=norm(z1,'fro');
    norm_z=norm_z1;
    
    norm_inf=max(max(abs(r1)));
    relres =norm_inf;              %relative residue =norm(z) / norm(b)
    if denom <= 0
        norm_d1=norm(d1,'fro');
        norm_d=norm_d1;
        p1 = d1/norm_d; % d is not a descent direction
        break % exit
    else
        alpha = rz1/denom;
        p1 = p1 + alpha*d1;
        r1 = r1 - alpha*w1;
    end
    z1= (M1\r1)/M2; %  z = M*r; if M is not the identity matrix ;
    norm_inf=max(max(abs(r1)));
    if norm_inf <= tolb % Exit if Hp=b solved within the relative tolerance
        iterk =m;
        relres =norm_inf;          %relative residue =norm(z) / norm(b)
        flag =0;
        break
    end
    rz2 = rz1;
    rz1 = trace(r1'*z1);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% end of pre_cg0.m %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the Jacobain product with x: F'(y,z)(x)
%%%%%% for the unconstrained case %%%%%%%
%%%%%%%

function Ax1= Jacobian_matrix0(x1,k,S2,S,c1,c2)
Ax1 =zeros(k,k);
Ax1 = x1*(S2'*S2/c1  + S'*S/c2) + S2*x1'*S2/c1 + S*x1'*S/c2 + x1 + x1';
% Ax1=c1^(-1)*Lambda2'*RTR*((Lambda2*x1)+(Lambda2*x1)')*RTR;
% Ax1=Ax1+ c2^(-1)*Lambda'*RTR*((Lambda*x1)+(Lambda*x1)')*RTR;
% Ax1=Ax1+ RTR*((x1)+(x1)')*RTR;
% Ax1=0.5*Ax1;

epsilon_add= 1.0e0;
Ax1 =Ax1 + epsilon_add*x1;


function A = AVAT1(Omega1,Omega3,Q1,Q3,S2,S,c1,c2,n,k,m1,m3,sigma)

A =zeros(k^2,k^2);


if m1 < n/2 && m1 >0
    Q1.llll = Q1.ll*Q1.ll';Q1.lrlr = Q1.lr*Q1.lr';
    S2P1ll = S2'*Q1.llll;S2P1lr = S2'*Q1.lr; S2Q1ll = Q1.ll'*S2;
elseif m1 <n
    Q1.llll = Q1.ll*Q1.ll';Q1.lrlr = Q1.lr*Q1.lr';
    S2P1lrlr = S2'*Q1.lrlr;S2P1lr = S2'*Q1.lr; S2Q1ll = Q1.ll'*S2;
end

if m3 < n/2 && m3 >0
    Q3.llll = Q3.ll*Q3.ll';Q3.lrlr = Q3.lr*Q3.lr';
elseif m3 <n
    Q3.llll = Q3.ll*Q3.ll';Q3.lrlr = Q3.lr*Q3.lr';
    
end

for j =1:k
    if m1 < n/2 && m1 >0
        S2P1llj = S2P1ll(j,:);
        S2P1lrj = S2P1lr(j,:);
        S2Q1llj = S2Q1ll(:,j);
    elseif m1 <n
        S2P1lrlrj = S2P1lrlr(j,:);
        S2P1lrj = S2P1lr(j,:);
        S2Q1llj = S2Q1ll(:,j);
    end
    
    A_21 = S(:,j);
    
    if m3 < n/2 && m3 >0
        Q3llj = Q3.llll(j,:);
        Q3lljt = Q3.ll(j,:);
        Q3lrj = Q3.lr(j,:);
    elseif m3 <n
        P3lrlrj = Q3.lrlr(j,:);
        
        Q3llj = Q3.ll(j,:);
    end
    
    
    for i=1:k
        if m1==0
            A1 =sparse(n,n);
        elseif m1 == n
            
            A1 = zeros(k,k);
            A1(:,i) =  S2(:,j);
            A1(i,:) = A1(i,:) + S2(:,j)';
            A1 = A1/(2*c1^0.5);
        elseif m1 <n/2
            Q1lli = Q1.ll(i,:);
            Q1lri = Q1.lr(i,:);
            A11 = Q1.llll(:,i)*S2P1llj;
            A11 = A11 + A11';
            
            A12 = Q1lli'*S2P1lrj;
            A12t = S2Q1llj* Q1lri;
            
            A12 = A12 + A12t;
            A12 = Omega1.*A12;
            A12 = Q1.ll * A12*Q1.lr';
            A1 = A11 + A12 + A12';
            A1  = A1/(2*c1^0.5);
            
        else
            A1org = zeros(k,k);
            
            A1org(:,i) =  S2(:,j);
            A1org(i,:) = A1org(i,:) + S2(:,j)';
            Q1lli = Q1.ll(i,:);
            
            A13 = Q1.lrlr(:,i) * S2P1lrlrj;
            A13 = A13 + A13';
            
            A12 = Q1lli'*S2P1lrj;
            A12t =  S2Q1llj*Q1.lr(i,:);
            A12 = A12 + A12t;
            A12 = Omega1.*A12;
            A12 = Q1.ll * A12*Q1.lr';
            
            A1 = A12 + A12' + A13;
            A1 = (A1org - A1)/(2*c1^0.5);
            
        end
        
        A2 = zeros(k,k);
        A2(:,i) =  A_21;
        A2(i,:) = A2(i,:) + A_21';
        A2 = A2/(2*c2^0.5);
        
        
        
        
        if m3==0
            A3 =sparse(n,n);
        elseif m3 == n
            A3 = zeros(k,k);
            A3(i,j) = 0.5;
            A3(j,i) = A3(j,i)+ 0.5;
        elseif m3 < n/2
            Q3lli = Q3.ll(i,:);
            Q3lri = Q3.lr(i,:);
            A31 = Q3.llll(:,i)*Q3llj';
            A31 = A31 + A31';
            
            A32 = Q3lli'*Q3lrj;
            A32t = Q3lljt'* Q3lri;
            
            A32 = A32 + A32t;
            A32 = Omega3.*A32;
            A32 = Q3.ll * A32*Q3.lr';
            A3 = A31 + A32 + A32';
            A3  = A3/2;
        else
            A3org = zeros(k,k);
            A3org(i,j) = 1;
            A3org(j,i) = A3org(j,i)+ 1;
            
            Q3lli = Q3.ll(i,:);
            
            A33 = Q3.lrlr(:,i) * P3lrlrj;
            A33 = A33 + A33';
            
            A32 = Q3lli'*Q3.lr(j,:);
            A32t =  Q3llj'*Q3.lr(i,:);
            A32 = A32 + A32t;
            A32 = Omega3.*A32;
            A32 = Q3.ll * A32*Q3.lr';
            
            A3 = A32 + A32' + A33;
            A3 = (A3org - A3)/2;
            
            
        end
        
        
        % A3 = Omega3.*A3;
        % A3 =(P3_1*A3)*P3';
        xx=c1^(-0.5)*A1*S2 + c2^(-0.5)*A2*S + A3;
        A(:,k*(j-1)+i) = reshape(xx,k^2,1);
    end
end
A = A + sigma*eye(k^2);


