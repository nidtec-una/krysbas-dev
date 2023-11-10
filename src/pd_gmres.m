function [x, flag, relres, iter, resvec, time] = pd_gmres(A, b, m0, tol, maxit, x0, alphaP, alphaD)
%PD-GMRES   Proportional-Derivate GMRES(m)
% 
%   pd_gmres is a modified implementation of the restarted Generalized
%   Minimal Residual Error or GMRES(m) (Saad, 1986), performed by using
%   a proportional-derivative control-inspired law to update adaptively
%   the restarting parameter m before each restart.
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
%   m0:       int
%             Restart parameter (similar to 'restart' in MATLAB).
%
%   tol:      float
%             Tolerance error threshold for relative residual norm.
%           
%   maxit:    int
%             Maximum number of iterations.
%
%   alphaP:   float
%             Proportional coefficient from the PD controller.
%
%   alphaD:   float
%             Derivative coefficient from the PD controller.
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
%   Nunez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A
%   proportional-derivative control strategy for restarting the GMRES (m)
%   algorithm. Journal of Computational and Applied Mathematics,
%   337, 209-224.
%
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

n = rowsA; delete rowsA colsA;

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

delete rowsb colsb;

% ----> Default values and sanity checks for parameter m0

% When the restart parameter is not specified explicitly as an input
% argument, pd_gmres() behaves identically to the unrestarted gmres()
% with a fixed number of iterations given by min(n, 10).
%
% There are three possibilities when this might happen: 
%   (1) If only two arguments are given i.e., A and b.
%   (2) If m0 equals the dimension of A i.e., m0 = n.
%   (3) If an empty matrix is given as restart parameter, i.e., m0 = [].
if isempty(m0) || m0 == n || nargin < 3
    restarted = false;  % use unrestarted pd_gmres()
else
    restarted = true;  % use restarted pd_gmres()
end

if restarted
    if m0 <= 0
        error("Restart parameter m0 must be a positive integer.")
    elseif m0 > n
        error("Restart parameter m0 cannot be greater than n.")
    end
end

% ----> Default value and sanity checks for tol
if isempty(tol) || nargin < 4
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
if isempty(maxit) || nargin < 5
    if restarted
        % If the restarted version of pd_gmres() version must be used
        % but maxit is not given, we take the min between n/m0 and 10.
        maxit = min(ceil(n/m0), 10);
    else
        % If the unrestarted version must ...
        maxit = min(n, 10);
    end
else
    
% ----> Default value and sanity checks for initial guess x0
if isempty(x0) || nargin < 6
    x0 = zeros(n, 1);
end

% Check whether x0 is a column vector
[~, colsx0] = size(x0);          
if colsx0 ~= 1
    error("Initial guess x0 is not a column vector.")
end

% Check whether x0 has the right dimension
if rowsb ~= n
    error("Initial guess x0 does not have the correct dimension.")
end

% ----> Default value for propotional parameter alphaP
if isempty(alphaP) || nargin < 7
    alphaP = -3;
end

% ----> Default value for proportional parameter alphaD
if isempty(alphaD) || nargin < 8
    alphaD = 5;
end

% ---> Algorithm starts here 

% Unrestarted version of the PD-GMRES
% This block calls the bult-in gmres() function but outputs slightly
% different variables. In particular, the resvec vector only shows the
% initial and the last residual (contrary to gmres() that stores all the
% residuals).
if ~restarted
    tic
    [x, flag, relres, iter, resvec] = gmres(A, b, [], tol, maxit, [], [], x0);
    resvec = [resvec(1); resvec(end)];
    time = toc();
    return
end

% Restarted version of PD-GMRES

% Algorithm setup 
mmin = 1;
mmax = n-1;
mstep =1 ;
flag=0;
restart=1;
r0=b-A*x0;
res(1,:)=norm(r0);
resvec(1,:)=(norm(r0)/res(1,1));
iter(1,:)=restart;
mIteration(1,1)=m0;

tic();  % start measuring CPU time

while flag==0
    if iter(size(iter,1),:) ~=1
        [miter]=pd_rule(m,m0,mmin,res,iter(size(iter,1),:),mstep, mmax,alphaP, alphaD); %cab
        m=miter(1,1);
        m0=miter(1,2);
    else
        m=m0;
    end
    mIteration(iter(size(iter,1),:)+1,1)=m;
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
    resvec(size(resvec,1)+1,:)=abs(g(m+1,1)/res(1,1));
    if abs(g(m+1, 1)) / res(1, 1) < tol || size(resvec, 1) == maxit % Using last component of g as residual
        flag=1;  % solution has converged
        x = xm;  % solution vector
    else
        x0=xm;  %update and restart
        restart=restart+1;
    end
    % Compute the relative residual
    relres = norm(b-A*xm)/norm(b);
end

time = toc();  % record CPU time

end
