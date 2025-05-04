
function [X, val_obj] = Correlation_Newton(G,C,m,n,tau,tau0,tol)

%%%%%%%%%%%%%%%%%%%                           %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This code is for the Local Correlation Stress Testing of the  Lagrangian Dual approach described in
%%%%%%  Houduo Qi and  Defeng Sun, "Correlation Stress Testing for Value-at-Risk: An Unconstrained
%%%%%%  Convex Optimization Approach", Department of Mathematics, National
%%%%%%  University of Singapore, March 2007 
%%%%%%%%%%%%%%%%%%%%%%%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 
% for solving 
%%  min  0.5 || X -  C ||^2 
%%  s.t.  X_ii = 1, i=1, ..., n
%%        X_ij = C_ij, i=1, ..., m; j =1, ..., n (Recall X = X^T)
%%        X- tau*I >=0 (symmetric and positive semi-definite)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Last updated on March 20, 2007 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ---Newton method starts--- ')
t0=cputime;


G =(G+G')/2; % make G symmetric
C =(C+C')/2;

if m<1
    tau =tau0;
else
    tau =min(tau, 0.90*tau0);
end

n_2=n-m;

if m>0
    rhs1 =ones(n,1);
    [P,D] = eig(G);
    d = diag(D);
    d =real(d);
    if min(d) < tau0
        fprintf('\n')
        fprintf('The smallest eignevalue of the initial G < %d \n', tau0)
        disp(' ---Start the Pre-processing--- ')
        [G,y]=CorNewton1(G,rhs1,tau0);
        C(1:m,1:n) = G(1:m,1:n); % update G to make sure  G is positive definite
        C((m+1):m,1:m) =  G((m+1):m,1:m);
        disp(' ---Pre-processing finished--- ')
    else
        fprintf('\n')
        fprintf('The smallest eignevalue of the initial matrix G(1:m,1:m) >= %d \n', tau0)
        disp(' ---No pre-processing needed--- ')
    end
end

C =C-tau*eye(n);  % reset C: C>= tau *I;

C_1=zeros(m, m);
C_2=zeros(m, n_2);

C_1=C(1:m, 1:m);
C_2=C(1:m, (m+1):n);

% Calculate Schur Complement
Y=C_2'*inv(C_1)*C_2;
dd=diag(Y);

C_3=zeros(n_2, n_2);
C_3=C((m+1):n, (m+1):n);


b=(1-tau)*ones(n_2,1)-dd;

C_3 = C_3-Y; % C_3 updated 
  
[Y0,y0] = CorNewton1(C_3,b,0); % updated C_3 is calibrated

C_3 = Y0+Y; % C_3 is updated back

X = [C_1,C_2; C_2',C_3]; 
G = X-C;
val_obj = sum(sum(G.*G))/2;
X = X + tau*eye(n); % optimal X


fprintf('Newton: Final function value %d \n', val_obj)

time_used = cputime -t0

