clear all
clc
sol7=[]; %PD-GMRES(m)

% This matrix is intended to be in /data
% Assistance on how to add path to /data would be helpful 
% This matrix is from SuiteSparse Matrix Collection
load sherman3.mat;

A=Problem.A;
b=Problem.b;

% Parametros
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
 [time, logres_pd_gmres]=pd_GMRES(A,b, mPD, alpha, delta,itermax);
 sol7(size(sol7,1)+1,:)= [time];
end