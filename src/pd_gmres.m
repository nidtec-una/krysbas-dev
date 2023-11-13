function [x, flag, relres, iter, resvec, time] = pd_gmres(A, b, ...
    mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD, varargin)
%PD-GMRES Proportional-Derivative GMRES(m)
% 
%   pd_gmres is a modified implementation of the restarted Generalized
%   Minimal Residual Error or GMRES(m) [1], performed by using
%   a proportional-derivative control-inspired law to update adaptively
%   the restarting parameter m before each restart. This implementation
%   follows closely the one presented in [2].
%
%   Signature:
%   ----------
% 
%   [x, flag, relres, iter, resvec, time] = pd_gmres(A, b, ...
%       mInitial, mMinMax, mStep, tol, maxit, xInitial, alphaPD, varargin)
%   
%
%   Input Parameters:
%   -----------------
%
%   A:        n-by-n matrix
%             Left-hand side of the linear system Ax = b.
%
%   b:        n-by-1 vector
%             Right-hand side of the linear system Ax = b.
%
%   mInitial: int, optional
%             Initial restart parameter (similar to 'restart' in MATLAB).
%             If empty, not given, or equal to n, then mInitial is not used
%             and the full unrestarted gmres algorithm with maxit = min(n,
%             10) is employed. Note that we require 1 <= mInitial <= n.
%
%   mMinMax:  2-by-1 vector, optional
%             Minimum and maximum values of the restart paramter m.
%             Default is [1; n-1]. Note that  we require 
%             1 <= mMinMax[1] < mMinMax[2] <= n.
%
%   mStep:    int, optional
%             Step size for incresing the restart parameter m between
%             cycles. Default is 1 if n <= 10 and 3 otherwise.
%
%   tol:      float, optional
%             Tolerance error threshold for the relative residual norm.
%             Default is 1e-6.
%           
%   maxit:    int, optional
%             Maximum number of iterations. 
%
%   xInitial: n-by-1 vector, optional
%             Vector of initial guess. Default is zeros(n, 1).
% 
%   alphaPD:  2-by-1 vector, optional
%             Proportional and derivative coefficients from the
%             proportional-derivative controller. Default is [-5, 3].
%
%   Output parameters:
%   ------------------
%
%   x:        n-by-1 vector
%             Approximate solution of the linear system.
%
%   flag:     boolean
%             1 if the algorithm converged, 0 otherwise.
%
%
%   log_res:  (1 up to to max_iter)-by-1 vector
%             relative residual norms
% 
%   References:
%   -----------
%
%   [1] Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal
%   residual algorithm for solving nonsymmetric linear systems. SIAM 
%   Journal on scientific and statistical computing, 7(3), 856-869.
%
%   [2] Nunez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A
%   proportional-derivative control strategy for restarting the GMRES(m)
%   algorithm. Journal of Computational and Applied Mathematics,
%   337, 209-224.
%

% ----> Sanity check on the number of input parameters
if nargin < 2
    error("Too few input parameters. Expected at least A and b.")
elseif nargin > 8
    error("Too many input parameters.")
end

% ----> Sanity checks on matrix A

% Check whether A is non-empty
if isempty(A)
    error("Matrix A cannot be empty.")
end

% Check whether A is square
[rowsA, colsA] = size(A);
if rowsA ~= colsA
    error("Matrix A must be square.")
end

n = rowsA;
clear rowsA colsA;

% ----> Sanity checks on vector b

% Check whether b is non-empty
if isempty(b)
    error("Vector b cannot be empty.")
end

% Check whether b is a column vector
[rowsb, colsb] = size(b);
if colsb ~= 1
    error("Vector b must be a column vector.")
end

% Check whether b has the same number of rows as b
if rowsb ~= n
    error("Dimension mismatch between matrix A and vector b.")
end

clear rowsb colsb;

% ----> Default values and sanity checks for parameter m0

% When the restart parameter is not specified explicitly as an input
% argument, pd_gmres() behaves identically to the unrestarted gmres()
% with a fixed number of iterations given by min(n, 10).
%
% There are three possibilities when this might happen: 
%   (1) If only two arguments are given i.e., A and b.
%   (2) If m0 equals the dimension of A i.e., m0 = n.
%   (3) If an empty matrix is given as restart parameter, i.e., m0 = [].
if (nargin < 3) || isempty(mInitial) || (mInitial == n)
    restarted = false;  % use unrestarted pd_gmres()
else
    restarted = true;  % use restarted pd_gmres()
end

if restarted
    if mInitial <= 0
        error("Restart parameter m0 must be a positive integer.")
    elseif mInitial > n
        error("Restart parameter m0 cannot be greater than n.")
    end
end

% ----> Default value and sanity checks for tol
if (nargin < 4) || isempty(tol)
    tol = 1e-6;
end

if tol < eps
    warning("Tolerance is too small and it will be changed to eps.")
    tol = eps;
elseif tol >= 1
    warning("Tolerance is too large and it will be changed to 1-eps.")
    tol = 1 - eps;
end

% ----> Default value for maxit
if (nargin < 5) || isempty(maxit)
    if restarted
        % If the restarted version of pd_gmres() version must be used
        % but maxit is not given, we take the min between n/m0 and 10.
        maxit = min(ceil(n/mInitial), 10);
    else
        % If the unrestarted version must ...
        maxit = min(n, 10);
    end
end
    
% ----> Default value and sanity checks for initial guess x0
if (nargin < 6) || isempty(xInitial) 
    xInitial = zeros(n, 1);
end

% Check whether x0 is a column vector
[rowsx0, colsx0] = size(xInitial);          
if colsx0 ~= 1
    error("Initial guess x0 is not a column vector.")
end

% Check whether x0 has the right dimension
if rowsx0 ~= n
    error("Initial guess x0 does not have the correct dimension.")
end

clear rowsx0 colsx0;

% ----> Default value for propotional parameter alphaP
% The default value for this parameter was taken from page 217 of
% https://doi.org/10.1016/j.cam.2018.01.009 
if (nargin < 7) || isempty(alphaP)
    alphaP = -3;
end

% ----> Default value for proportional parameter alphaD
% The default value for this parameter was taken from page 217 of
% https://doi.org/10.1016/j.cam.2018.01.009 
if (nargin < 8) || isempty(alphaD)
    alphaD = 5;
end

% ---> Algorithm starts here 

% Unrestarted version of the PD-GMRES
% This block calls the bult-in gmres() function but outputs slightly
% different variables. In particular, the resvec vector only shows the
% initial and the last residual (contrary to gmres() that stores all the
% residuals).
if ~restarted
    tic();
    [x, flag, relres, iter, resvec] = gmres(A, b, [], tol, maxit, [], [], xInitial);
    resvec = [resvec(1); resvec(end)];
    time = toc();
    return
end

% Restarted version of PD-GMRES

% Algorithm setup 
%mmin = 1;
mmax = n-1;
mstep =1 ;
flag=0;
restart=1;
r0=b-A*xInitial;
res(1,:)=norm(r0);
resvec(1,:)=(norm(r0)/res(1,1));
iter(1,:)=restart;
mIteration(1,1)=mInitial;

tic();  % start measuring CPU time

while flag==0
    if iter(size(iter,1),:) ~=1
        [miter]=pd_rule(m,mInitial,mMin,res,iter(size(iter,1),:),mstep, mmax,alphaP, alphaD);
        m=miter(1,1);
        mInitial=miter(1,2);
    else
        m=mInitial;
    end
    mIteration(iter(size(iter,1),:)+1,1)=m;
    v=zeros(n,m+1);
    w=zeros(n,m);
    r=b-A*xInitial;
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
    xm=xInitial+V*minimizer;
    res(restart+1,:)=abs(g(m+1,1));
    iter(restart+1,:)=restart+1;
    resvec(size(resvec,1)+1,:)=abs(g(m+1,1)/res(1,1));
    if abs(g(m+1, 1)) / res(1, 1) < tol || size(resvec, 1) == maxit % Using last component of g as residual
        flag=1;  % solution has converged
        x = xm;  % solution vector
    else
        xInitial=xm;  %update and restart
        restart=restart+1;
    end
    % Compute the relative residual
    relres = norm(b-A*xm)/norm(b);
end

time = toc();  % record CPU time

end
