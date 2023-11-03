function [x, logres]=pd_gmres(A,b, mPD, alphaP, alphaD,itermax,tol)
%
% Description:
%
% pd_GMRES is a modified implementation of the restarted Generalized
% Minimal Residual Error or GMRES(m) (Saad, 1986), performed by using a 
% proportional-derivative control-inspired law to update adaptively the 
% restarting parameter m before each restart.
%
%
% Input Parameters:
%
% A:        n-by-n matrix
%           left-hand side of the linear system Ax = b
%
% b:        n-by-1 vector
%           right-hand side of the linear system Ax = b
%
% m_PD:     int
%           restart parameter (similar to 'restart' in MATLAB)
%
% tol:      float
%           tolerance error threshold for relative residual norm
%           
% max_iter: int
%           maximum number of (inner?) iterations
%
% alphaP:   float
%           proportional coefficient from PD controller for 'm'
%
% alphaD:   float
%           derivative coefficient from PD controller for 'm'
%
% Output parameters:
%
% log_res:  (1 up to to max_iter)-by-1 vector
%           relative residual norms
%
% xx:       n-by-1 vector
%           approximate solution of the linear system 
%
% References:
%
% Nunez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A proportional-derivative 
% control strategy for restarting the GMRES (m) algorithm. Journal of 
% Computational and Applied Mathematics, 337, 209-224.
%

[s, n] = size(A);
x0 = zeros(size(b,1),1); % x0: First guess for solution vector 'x' 
mInitial = mPD;
mmin = 1;
mmax = n-1;
mstep=1;
maxit=itermax;
flag=0;
if (s~=n)
    error ('Matrix not square');
end
[i,j]=size (b);
if (s~=i)
    error ('Vector b does not match size of matrix A');
end
if (j~=1)
    error ('Vector is not a column vector')
end
if (size (b)~=size(x0))
    error('Incorrect size of initial guess vector x0');
end
restart=1;
r0=b-A*x0;
res(1,:)=norm(r0);
logres(1,:)=(norm(r0)/res(1,1));
iter(1,:)=restart;
mIteracion(1,1)=mInitial;
while flag==0
    if iter(size(iter,1),:) ~=1
        [miter]=pd_rule(m,mInitial,mmin,res,iter(size(iter,1),:),mstep, mmax,alphaP, alphaD); %cab
        m=miter(1,1);
        mInitial=miter(1,2);
    else
        m=mInitial;
    end
    mIteracion(iter(size(iter,1),:)+1,1)=m;
    v=zeros(n,m+1);
    w=zeros(n,m);
    r=b-A*x0;
    beta=norm(r);
    v(:,1)=r/beta; 
    h=zeros(m+1,m);
    for j=1:m                       % Modified gram schmidt--Arnoldi
        w(:,j)=A*v(:,j);
        for i=1:j
            h(i,j)=w(:,j)'*v(:,i);
            w(:,j)=w(:,j)-h(i,j)*v(:,i);
        end
        h(j+1,j)=norm(w(:,j));
        if h(j+1,j)==0
            m=j;
            h2=zeros(m+1,m);    % Comment by JCC: VERIFY!!! (why?)
            for k=1:m
                h2(:,k)=h(:,k);
            end
            h=h2;
        else        
        v(:,j+1)=w(:,j)/h(j+1,j);
        end
    end
    g=zeros(m+1,1);
    g(1,1)=beta;
    for j=1:m                       % Plane rotations (QR decompostion)
        P=eye(m+1);   
        sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
        cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));
        P(j,j)=cos;
        P(j+1,j+1)=cos;
        P(j,j+1)=sin;
        P(j+1,j)=-sin;
        h=P*h;
        g=P*g;
    end
    R=zeros(m,m);
    G=zeros(m,1);
    V=zeros(n,m);
    for k=1:m
        G(k)=g(k);
        V(:,k)=v(:,k);
        for i=1:m
            R(k,i)=h(k,i);
        end
    end
    minimizer=R\G;
    xm=x0+V*minimizer;
    res(restart+1,:)=abs(g(m+1,1));
    iter(restart+1,:)=restart+1;
    logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,1));
    if abs (g(m+1,1))/res(1,1) <tol   || size(logres,1)==maxit   % Using last component of g as residual
        flag=1;
        x = xm;
    else
        x0=xm ;                       %update and restart
        restart=restart+1;
    end
end