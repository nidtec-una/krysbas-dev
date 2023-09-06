% Solve a 1D-Poisson equation at [aStart, aEnd] = [0, 1]
% u''(x) = f(x), u(0) = g(0) = 0, u(1) = g(1) = 1;

% f(x) = 6*x; g(x) = x^3
f = @(x) 6*x;
g = @(x) x.^3;

NODES = 31;
aStart = 0;
aEnd = 1;
h = ( aEnd - aStart ) / (NODES - 1);
b = h*h*f( linspace(aStart, aEnd, NODES) )';
b(1,1) = g( aStart );
b(NODES,1) = g( aEnd );

A = 2*eye(NODES) - diag(ones(NODES-1,1), 1) - diag(ones(NODES-1, 1), -1);
A(1, 1:2) = [1 0];
A(NODES, NODES-1:NODES) = [0 1];

u1 = A\b;
u2 = g( linspace(aStart, aEnd, NODES) )';

rootFolder = fileparts(pwd); % go to the root folder
srcFolder = fullfile(rootFolder, 'src'); % enter the data folder
cd(srcFolder)
%     [time, logres_pd_gmres]=pd_gmres(A,b, mPD, alpha, delta,itermax); 
% Parameters
alpha=-3; %2
delta=5; %0.8
opts_tol=1e-9;
itermax=1000;
p = 1;

%%       %PD-GMRES(m)
% The original idea was to compute the average execution time,
% we may discuss if this is still necessary 
for i=1:p
    color_pd_gmres='b';
    mPD=30;
    %alpha=2;
    rootFolder = fileparts(pwd); % go to the root folder
    srcFolder = fullfile(rootFolder, 'src'); % enter the data folder
    cd(srcFolder)
    [time, logres_pd_gmres, u3]=pd_gmres(A,b, mPD, alpha, delta,itermax);
    %sol7(size(sol7,1)+1,:)= [time];
end