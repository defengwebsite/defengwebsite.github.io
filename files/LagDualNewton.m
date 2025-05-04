function [X,y,z,val_obj] = LagDualNewton(G,m,n,tau,tau0,tol)
%%%%%%%%%%%%%%%%%%%                           %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This code is for the Band Correlation Stress Testing of the  Lagrangian Dual approach described in
%%%%%%  Houduo Qi and  Defeng Sun, "Correlation Stress Testing for Value-at-Risk: An Unconstrained
%%%%%%  Convex Optimization Approach", Department of Mathematics, National
%%%%%%  University of Singapore, March 2007 
%%%%%%%%%%%%%%%%%%%%%%%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 
% for solving 
%%  min  0.5 || X -  G ||^2 
%%  s.t.  X_ii = 1, i=1, ..., n
%%        X_ij = C_ij, 1<= i < j <=m
%%        X- tau*I >=0 (symmetric and positive semi-definite)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %%%%%%%%%%%%%%%%%%%%%%%
% G: stressed correlation matrix (n by n)
% G(1:m,1:m) is fixed
% tau0: the threshold for positive definiteness of  G  
% tau: the lower bound for the smallest eigevalue of X   
% tol: stopping crierion for the KKT system: 1.0e-5 ~ 1.0e-6
% X: the calibrated correlation matrix
% y and z: dual solutions
% val_obj: final objective function value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Last updated on March 20, 2007 %%%%%%%%%%%%%%%%%%%%%%%%%%%%



disp(' ---Newton method starts--- ')
t0=cputime;
global b0
global c0


if m<2
    tau =tau0;
else
    tau =min(tau, 0.90*tau0);
end

if m>1
    I=(m-1)*ones(1,1);
    J=m*ones(1,1);
else
    I =[];
    J =[];
end
for i=(m-2):-1:1
    I=[i.*ones((m-i),1); I];
    J=[[(i+1):m]'; J];
end


[n_1, n_c]=size(I);

%[m,m_c] =size(I);


if m>1
    rhs1=ones(m,1);
    C1 = G(1:m,1:m);

    [P1,D1] = eig(C1);
    d1 = diag(D1);
    d1 =real(d1);
    if min(d1) <tau0
        fprintf('\n')
        fprintf('The smallest eignevalue of the initial G(1:m,1:m) < %d \n', tau0)
        disp(' ---Start the Pre-processing--- ')
        [C1,y]=CorNewton1(C1,rhs1,tau0);
        G(1:m,1:m) = C1(1:m,1:m); % update G to make sure  G(1:m,1:m) is positive definite
        disp(' ---Pre-processing finished--- ')
    else
        fprintf('\n')
        fprintf('The smallest eignevalue of the initial matrix G(1:m,1:m) >= %d \n', tau0)
        disp(' ---No pre-processing needed--- ')
    end
end

s =zeros(n_1,1);
for i=1:n_1
    s(i)=G(I(i),J(i));
end
 
b0 =ones(n,1);
c0 =s;


% disp('Start the ininial correlation test approximation')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% [X0,y0] = CorNewton1(G,b0,tau0);
% 
% 
% disp('The ininial correlation test approximation finished')
% 
% 
% 
% disp(' ---Pre-processing finished--- ')


if tau>0
   G = G-tau*eye(n); % reset G
   b0 = (1-tau)*b0; % reset b0   
end
 val_obj=0;

Res_b=zeros(300,1);

y=zeros(n,1);
z=zeros(n_1,1);
y=b0-diag(G);              %Initial point

Fy=zeros(n,1);
Fz=zeros(n_1,1);

k=0;
f_eval =0;

Iter_Whole=100;
Iter_inner =40; % Maximum number of Line Search in Newton method
%maxit =1000; %Maximum number of iterations in PCG
maxit =min(20000, max(10000, (n_1+n+1)));    % maximum step the CG method
Inner =0;


error_tol=1.0e-6; % termination tolerance
sigma_1=1.0e-4; %tolerance in the line search of the Newton method


x0_y=y;
x0_z=z;
M =eye(n,n); % Preconditioner if M is not I

dy =zeros(n,1);
dz =zeros(n_1,1);

Z =sparse(I,J,z,n,n);
Z =0.5*(Z+Z');
 C =G+diag(y)+Z;
 [P,D]=eig(C);
  P =real(P);
 lambda=diag(D);
 lambda= real(lambda);
 
 [f0,Fy,Fz] = gradientS(y,z,lambda,P,b0,c0,I,J);
 f_eval =f_eval + 1;
 by =b0-Fy;
 bz =c0-Fz;
 norm_b=sqrt(by'*by+bz'*bz);

 fprintf('Newton: Norm of Gradient %d \n',norm_b)
 Omega = omega_matrix(lambda);
 x0_y =y;
 x0_z =z;
 
 while (norm_b>error_tol & k< Iter_Whole)

 [dy,dz,flag,relres,iterk]=pre_cgS(by,bz,tol,maxit,M,Omega,P,I,J);

  fprintf('Newton: Number of CG Iterations %d \n', iterk)
  
  if (flag~=0); % if CG is unsuccessful, use the negative gradient direction
      dy =b0-Fy;
      dz =c0-Fz;
      fprintf('Newton:  CG method fails as flag== \n', flag)
  end
 slope =(Fy-b0)'*dy+(Fz-c0)'*dz;  %%% nabla f d
 

    y =x0_y+dy; 
    z =x0_z+dz; %temporary x0+d 
    
    Z =sparse(I,J,z,n,n);
    Z =0.5*(Z+Z');
     C = G +diag(y)+Z;
     [P,D]=eig(C); % Eig-decomposition: C =P*D*P^T
     P =real(P);
     lambda=diag(D);
     lambda=real(lambda);
     
     [f,Fy,Fz] = gradientS(y,z,lambda,P,b0,c0,I,J);

     k_inner=0;
     while(k_inner <=Iter_inner & f>f0+sigma_1*0.5^k_inner*slope)
         k_inner=k_inner+1;
         y = x0_y+0.5^k_inner*dy; 
         z = x0_z+0.5^k_inner*dz; % backtracking
         Z = sparse(I,J,z,n,n);
         Z = 0.5*(Z+Z');
         C = G +diag(y)+Z;
         [P,D]=eig(C); % Eig-decomposition: C =P*D*P^T
         P =real(P);
         lambda=diag(D);
         lambda =real(lambda);
         [f,Fy,Fz] = gradientS(y,z,lambda,P,b0,c0,I,J);
      end % loop for while
      f_eval =f_eval+k_inner+1;
         x0_y=y;
         x0_z=z;
     k=k+1;
     by=b0-Fy;
     bz=c0-Fz;
     norm_b=sqrt(by'*by+bz'*bz);
     fprintf('Newton: Steplength ==  %d \n',0.5^k_inner)
     fprintf('Newton: Norm of Gradient %d \n',norm_b)

     Res_b(k)=norm_b;
     Omega = omega_matrix(lambda);
 end %end loop for while
 
 
fprintf('Newton: function value %d \n',f0)
fprintf('Newton: Norm of Gradient %d \n',norm_b)
fprintf('Newton: Number of Iterations %d \n', k)
fprintf('Newton: Number of Function Evaluations %d \n', f_eval)
i=1;
C =P';
while (i<n+1)
    C(i,:) = max(0,lambda(i))*C(i,:);
    i=i+1;
end
X =P*C + tau*eye(n); % Optimal solution X* (Recall that tau*I is restored) 

for i=1:n
    for j=1:n
       val_obj=val_obj+(G(i,j)-X(i,j))^2;
    end
end
val_obj = val_obj/2;
fprintf('Newton: objective function value %d \n',val_obj)


time_used= cputime-t0

%%% end of the main program

function [f,Fy,Fz]= gradientS(y,z,lambda,P,b0,c0,I,J)

[n, n_c]=size(b0);
[m, m_c]=size(c0);

f=0.0;
Fy =zeros(n,1);
Fz =zeros(m,1);
%lambdap=max(0,lambda);
%H =diag(lambdap); %% H =P^T* diag(x) *P
%  H =H*P'; %%% Assign H*P' to H
 H=P';
 i=1;
 while (i<=n)
     H(i,:)=max(lambda(i),0)*H(i,:);
     i=i+1;
 end
 i=1;
 while (i<=n)
       Fy(i)=P(i,:)*H(:,i);
 i=i+1;     
 end
 i=1;
 while (i<=m)
       Fz(i)=P(I(i),:)*H(:,J(i));
       i=i+1;
 end
 i=1;
 while (i<=n)
     f =f+(max(lambda(i),0))^2;
     i=i+1;
 end
 
f =0.5*f -b0'*y -c0'*z;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% end of gradient.m %%%%%%

%%%%%%%%%%%%%%        %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% To generate the first -order difference of lambda
%%%%%%%
function omega = omega_matrix(lambda)

[n, n_c]=size(lambda);
omega =ones(n,n);
%Im =find(lambda<0);
%Ip =find(lamba>=0);
i=1;
while (i<=n)
    j=1;
    while (j<=n)
        if abs(lambda(i)-lambda(j))>1.0e-10
            omega(i,j) = (max(0,lambda(i))-max(0,lambda(j)))/(lambda(i)-lambda(j));
             elseif max(lambda(i),lambda(j))<=1.0e-15
                  omega(i,j)=0;
         end
      j=j+1;
    end
    i=i+1;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of omega_matrix.m %%%%%%%%%%
%%%%%% PCG method %%%%%%%
%%%%%%% This is exactly the algorithm by  Hestenes and Stiefel (1952)
%%%%%An iterative method to solve A(x) =b  
%%%%%The symmetric positive definite matrix M is a
%%%%%%%%% preconditioner for A. 
%%%%%%  See Pages 527 and 534 of Golub and va Loan (1996)

function [py,pz,flag,relres,iterk] = pre_cgS(by,bz,tol,maxit,M,Omega,P,I,J);

[n,n_c]=size(by);
[m,m_c]=size(bz);

% Initializations
ry = by;  
rz = bz; %We take the initial guess x0=0 to save time in calculating A(x0) 
n2b =sqrt(by'*by+bz'*bz);    % norm of b
tolb = tol * n2b;  % relative tolerance 
py = zeros(n,1);
pz = zeros(m,1);
flag=1;
iterk =0;
relres=1000; %%% To give a big value on relres
% Precondition 
zy = ry;
zz = rz;  %%%%% z = M\r; if M is not the identity matrix 
rz1 = ry'*zy + rz'*zz; 
rz2 = 1; 
dy = zy; 
dz = zz;
% CG iteration
for k = 1:maxit
   if k > 1
       beta = rz1/rz2;
       dy = zy + beta*dy;
       dz = zz + beta*dz;
   end
   [wy, wz]= Jacobian_matrixS(dy,dz,Omega,P,I,J); %w = A(d); 
   denom = dy'*wy + dz'*wz;
   iterk =k;
   relres =sqrt(zy'*zy+zz'*zz)/n2b;   %relative residue =norm(z) / norm(b)
   if denom <= 0 
       sssss=0;
       normd = sqrt(dy'*dy+dz'*dz);
       py = dy/normd; 
       pz = dz/normd; % d is not a descent direction
       break % exit
   else
       alpha = rz1/denom;
       py = py + alpha*dy;
       pz = pz + alpha*dz;
       ry = ry - alpha*wy;
       rz = rz - alpha*wz;
   end
   zy = ry; 
   zz = rz; %  z = M\r; if M is not the identity matrix ;
   normz = sqrt(zy'*zy + zz'*zz);
   if normz <= tolb % Exit if Hp=b solved within the relative tolerance
       iterk =k;
       relres =normz/n2b;          %relative residue =norm(z) / norm(b)
       flag =0;
       break
   end
   rz2 = rz1;
   rz1 = ry'*zy + rz'*zz;
end

return

%%%%%%%% %%%%%%%%%%%%%%%
%%% end of pre_cg.m%%%%%%%%%%%

%%%%%% To generate the Jacobain product with x: F'(y)(x) %%%%%%%
%%%%%%%

function [Ax_y, Ax_z]= Jacobian_matrixS(x_y,x_z,omega,P,I,J)

[n, n_c]=size(x_y);
[m, m_c]=size(x_z);

Ax_y =zeros(n,1);
Ax_z =zeros(m,1);

Z=sparse(I,J,x_z,n,n);
Z=0.5*(Z+Z');
Z =diag(x_y)+Z;

PHP=P'*Z*P;
H=omega.*PHP;
H=H*P'; %%% Assign H*P' to H= Omega o P^T*diag(x)*P)*P^T
 i=1;
 while (i<=n)
       Ax_y(i)=P(i,:)*H(:,i);
       Ax_y(i) = Ax_y(i) + 1.0e-8*x_y(i); % add a small perturbation
       i=i+1;
 end
 i=1;
 while (i<=m)
       Ax_z(i)= P(I(i),:)*H(:,J(i));
       Ax_z(i) = Ax_z(i) + 1.0e-8*x_z(i); % add a small perturbation
       i=i+1;
 end
 return
 

%%%%%%%%%%%%%%%
%end of Jacobian_matrix.m%%%
