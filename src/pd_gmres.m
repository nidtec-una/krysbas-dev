function [x, flag, relres, iter, resvec, time] = pd_gmres(A, b, m0, tol, maxit, x0, alphaP, alphaD)
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

%% Sanity checks and default values of parameters

% ----> Matrix A

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
delete rowsA colsA

% ----> Vector b

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
if rowsb ~= rowsA
    error("Dimension mismatch between matrix A and vector b.")
end

delete rowsb colsb

% ----> Restart parameter m0

% In the case where the restart parameter is not specified explicitly as an
% input argument, pd_gmres() behaves identically to the unrestarted gmres()
% with a fixed number of iterations given by min(n, 10).
%
% There are three possibilities where this might happen: 
%   (1) If only two arguments are given i.e., 'A' and 'b'.
%   (2) If 'm0' equals the dimension of 'A' i.e., 'm0' = 'n'.
%   (3) If an empty matrix is given as restart parameter, i.e., 'm0' = [].
%
% In cases (2) and (3), we still allow the user to give 'tol' and 'x0'.
% Since the number of iterations are fixed, a warning is raised if the
% input argument 'maxit' is given.
if isempty(m0) || m0 == n || nargs < 3
    % Tolerance specification
    if exist('tol', 'var') == 1
        if isempty(tol) == 1
            tol = 1e-6; % use default tolerance
        else
            ... % use user-specified tolerance
        end
    else
        tol = 1e-6; % use default tolerance
    end

    % Initial guess specification
    if exist('x0', 'var') == 1
        if isempty(x0) == 1
            x0 = zeros(n, 1);
        else
            ... % use user-specified initial guess
        end
    else
        x0 = zeros(n, 1);  % use default initial guess
    end

    % Warn user if 'maxit' is given
    if exist("maxit", "var") == 1
        warning("Number of maximum iterations will be overwritten to min(n, 10).")
    end
    
    % We are now ready to call the pd_gmres_unrestarted() routine   
    % Note: We could have called gmres() directly. However, pd_gmres()
    % there are subtle differences in the output. For instance, pd_gmres()
    % stores only the last residual in the vector resvec, whereas gmres()
    % stores all the residuals. In addition, gmres() does not provide the
    % cpu_time natively.
    [x, flag, relres, iter, resvec, time] = pd_gmres_unrestarted(A, b, tol, x0);
    return
end

% Check whether m0 is strictly less than the dimensionality of the system
if m0 > n
    error("Restart parameter m0 cannot be greater than n.")
end

% ----> Tolerance
if isempty(tol) || nargs < 4
    tol = 1e-6;
end

% ----> Maximum number of iterations maxit
if isempty(maxit) || nargs < 5
    maxit = min(round(n/m0), 10);
end

% ----> Initial guess x0
if isempty(x0) || nargs < 6
    x = zeros(n, 1);
end

% Check whether x0 is a column vector
rowsx0, colsx0 = size(x0);          
if colsx0 ~= 1
    error("Initial guess x0 is not a column vector.")
end

% Check whether x0 has the right dimension
if rowsb ~= n
    error("Initial guess x0 does not have the correct dimension.")
end

[s, n] = size(A);
x0 = zeros(size(b,1),1); % x0: First guess for solution vector 'x' 
mInitial = m0;
mmin = 1;
mmax = n-1;
mstep=1;
maxit=maxit;
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
resvec(1,:)=(norm(r0)/res(1,1));
iter(1,:)=restart;
mIteracion(1,1)=mInitial;

tic();  % start measuring CPU time

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


function [x, flag, relres, iter, resvec, time] = pd_gmres_unrestarted(A, b, tol, x0)
    % This is the so-called unrestarted version of the PDGMRES algorithm
    % In pratice, this routine calls the bult-in gmres() function but
    % outputs slightly different variables.
    
    % Unrestarted version of the PD-GMRES
    tic
    [x, flag, relres, iter, resvec] = gmres(A, b, [], tol, [], [], [], x0);
    resvec = [resvec(1); resvec(end)];
    time = toc();
end
