%%modificado por jccf 14 de febrero 2022
clear all
clc
%inicio=cputime;         %Control de tiempo de ejecucion
%count=1;
sol1=[]; %GMRES(m)
sol2=[]; %LGMRES(m,l)
sol3=[]; %GMRES-E(m,d)
sol4=[]; %LGMRES-E(m,d,l)
sol5=[]; %A-PD-LGMRES-E(m,d,l)
sol6=[]; %SLGMRES-E-s(m,d,l)
sol7=[]; %PD-GMRES(m)

%load add20.mat
%load circuit_2.mat;
%load raefsky1.mat;
%load raefsky2.mat;
%load fpga_trans_01.mat;
%load sherman4.mat;
% load orsreg_1.mat;
% load b_orsreg_1.mat;
%load cdde1.mat;
% %load b_cdde1.mat;
%  load orsirr_1.mat;
%  load b_orsirr_1.mat;
% load pde2961.mat;
% load b_pde2961.mat;
% load rdb2048.mat;
% load b_rdb2048.mat;
% load steam2.mat;
% load b_steam2.mat;
%  load wang2.mat;
%  load b_wang2.mat;
% load watt_1.mat;
% load b_watt_1.mat;
%load memplus.mat;
%load cavity05.mat
% load young3c.mat;
% load b_young3c.mat;
load sherman5.mat;
%load example2_morgan.mat;
%load matriz_16_05_3agrupado.mat;
%load memplus.mat;
%load cavity10.mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load mahindas.mat; %inluir b=b(:,1)
% %load AbD.mat  %propuesto por Diego
%  load ex40.mat;
%  load b_ex40.mat;
%load wang4.mat;
%load b_wang4.mat;
%load ex5.mat;
% load shl_400.mat;
% load b_shl_400.mat;
%load bips98_1142.mat; %muy dificil!..
%load raefsky3.mat;
A=Problem.A;
b=Problem.b;
%b=b(:,1);
[ss,n]=size (A);
%   bmax=max(max(A));
%   bmin=min(min(A));
%   b=bmin*ones(n,1)+(bmax-bmin)*rand(n,1);
% Parametros
        alpha=2; %
        delta=0.8; %
       opts_tol=1e-9;
       itermax=1000;
  %     [L,U] = ilu(A);

 %       [L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));
  %     M=L*U;
  %     A=M\A; b=M\b;

Name_Matrix= 'Sherman5';
p=1;
for c=1:1
%% GMRES(m)
%    for i=1:p
% 
%         color_gmres='b.-';
%         m=30;
%         [time,logres_gmres]=gmres_m(A,b,m,itermax);
%         sol1(size(sol1,1)+1,:)= [time];
%      end
%%
%% LGMRES(m,l)
%      for i=1:p
% 
%          color_lgmres='g';
%          mL=28;
%          lL=2;
%          [time, logres_lgmres]=lgmres(A,b,mL,lL,itermax);
%          sol2(size(sol2,1)+1,:)= [time];
%      end
%%
%%   GMRES-E(m,d)
%      for i=1:p
% 
%           color_gmrese='r';
%           mE=28;
%           dE=2;
%          [time,logres_gmrese]=GMRES_E(A,b,mE,dE,itermax);
%           sol3(size(sol3,1)+1,:)= [time];
%      end
%%
%% LGMRES-E(m,d,l)
%     for i=1:p
% 
%           color_lgmrese='c';
%           mLE=27;
%           dLE=2;
%           lLE=1;
%           [time, logres_lgmrese]=LGMRES_E_v2(A,b,mLE, lLE, dLE,itermax);
%           sol4(size(sol4,1)+1,:)= [time];
%      end
%%
%% A-PD-LGMRES-E(m,d,l)
%      for i=1:p
% 
%         color_pd_lgmrese='m';
%         mApd=28;
%         dApd=2;
%         lApd=2;
%        e=[0.01 0.03 0.05 0.1 0.15 0.2 0];
%          e=[0.2 0];
%         e=[0 1e-6 1e-3 0.01 0.1];
%         epsilon=e(1,c);
%         epsilon=0.01;
%         [time, logres_pd_lgmrese]=Adaptive_PD_lgmres_e(A,b,itermax,alpha,delta,mApd,dApd,lApd, epsilon);
%         sol5(size(sol5,1)+1,:)= [time];
%      end
%%
%% SLGMRES-E-s(m,d,l)
%     for i=1:p
% 
%           color_slgmrese='c+:';
%           mLE=28;
%           dLE=2;
%           lLE=2;
% %         e=[0.03 0.05];
%     %    e=[0 1e-6 1e-3 0.01 0.03 0.05 0.1];
% %         e=[0.01];
% %          epsilon=e(1,c);
%          epsilon=0.01;
%           [time, logres_slgmrese]=Adaptive_lgmres_e_switch(A,b,mLE, lLE, dLE,itermax,epsilon);
%           sol6(size(sol6,1)+1,:)= time;
%      end
%%
%%       %PD-GMRES(m)
      for i=1:p

         color_pd_gmres='b';
         mPD=30;
         %alpha=2;
         [time, logres_pd_gmres]=pd_GMRES(A,b, mPD, alpha, delta,itermax);
         sol7(size(sol7,1)+1,:)= [time];
     end
% % % % % %

 end
%tiempototal=cputime - inicio     %Imprime tiempo de ejecuci�n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure(1)
% subplot(1,1,1);
% xlabel('Number of restart cycles');ylabel('||r_m^{(j)}||/||r_{0}||');
%
% %['GMRES(30); cycle= ', num2str(sol(1,1)),'; t= ', num2str(sol(1,2)),' s.'], ...
%
%  legend(['GMRES(30); cycle= ', num2str(sol(1,1)),'; t= ', num2str(sol(1,2)),' s.'], ...
%  ['LGMRES(28,2); cycle= ', num2str(sol1(1,1)),'; t= ', num2str(sol1(1,2)),' s.'], ...
%  ['GMRES-E(28,2); cycle= ', num2str(sol2(1,1)),'; t= ', num2str(sol2(1,2)),' s.'], ...
%  ['LGMRES-E(27,2,1); cycle= ', num2str(sol3(1,1)),'; t= ', num2str(sol3(1,2)),' s.'], ...
%  ['A-SLGMRES-E(28,2,2); cycle= ', num2str(sol4(1,1)),'; t= ', num2str(sol4(1,2)),' s.'],'Location','SouthOutside');

% legend(['GMRES-E(26,8); cycle= ', num2str(sol2(1,1)),'; t= ', num2str(sol2(1,2)),' s.'], ...
%  ['SLGMRES-E(26,8,2); cycle= ', num2str(sol5(1,1)),'; t= ', num2str(sol5(1,2)),' s.'],'Location','SouthOutside');



 %hold on
% subplot(2,1,2);
% xlabel('Number of restart cycles'); ylabel('m, restart parameter')
%  legend(['A-SLGMRES-E(28,2,2); last cycle= ', num2str(sol4(1,1)),'; t= ', num2str(sol4(1,2)),'sec.'], ...
%['SLGMRES-E(28,2,2); last cycle= ', num2str(sol5(1,1)),'; t= ', num2str(sol5(1,2)),'sec.'],...
%['PD-GMRES(30); cycle= ', num2str(sol6(1,1)),'; t= ', num2str(sol6(1,2)),'s.'],'Location','SouthOutside');
%hold on
% legend(['LGMRES(28,2); last cycle= ', num2str(sol1(1,1)),'; t= ', num2str(sol1(1,2)),'sec.'], ...
% ['GMRES-E(28,2); last cycle= ', num2str(sol2(1,1)),'; t= ', num2str(sol2(1,2)),'sec.'], ...
% ['LGMRES-E(27,2,1); last cycle= ', num2str(sol4(1,1)),'; t= ', num2str(sol4(1,2)),'sec.'], ...
% ['Adaptive-SLGMRES-E(28,2,2); last cycle= ', num2str(sol3(1,1)),'; t= ', num2str(sol3(1,2)),'sec.'],'Location','SouthOutside');

%   legend(['A-SLGMRES-E(28,2,2); cycle= ', num2str(sol4(1,1)),'; t= ', num2str(sol4(1,2)),'s.'], ...
%   ['PD-GMRES(30); cycle= ', num2str(sol6(1,1)),'; t= ', num2str(sol6(1,2)),'s.'],'Location','SouthOutside');
% %

% legend(['LGMRES(28,2); last cycle= ', num2str(sol1(1,1)),'; t= ', num2str(sol1(1,2)),'sec.'], ...
% ['GMRES-E(28,2); last cycle= ', num2str(sol2(1,1)),'; t= ', num2str(sol2(1,2)),'sec.'], ...
% ['LGMRES-E(27,2,1); last cycle= ', num2str(sol4(1,1)),'; t= ', num2str(sol4(1,2)),'sec.'],'Location','SouthOutside');

% legend(['A-SLGMRES-E(28,2,2); last cycle= ', num2str(sol3(1,1)),'; t= ', num2str(sol3(1,2)),'sec.'], ...
% ['Switch-LGMRES-E(28,2,2); last cycle= ', num2str(sol5(1,1)),'; t= ', num2str(sol5(1,2)),'sec.'],'Location','SouthOutside');


%title('cavity10.mat')
%  hold on
%  figure(2)
% %subplot(2,1,2);
% title(Name_Matrix)
% xlabel('Number of restart cycles'); ylabel('m, restart parameter')
%  legend(['A-SLGMRES-E(28,2,2)'], ...
%  ['PD-GMRES(30)'],'Location','SouthOutside');


% semilogy(logres,color)

 %semilogy(logres_gmres,'r')
% ## semilogy(logres_gmres,color_gmres)
% hold on
% ## semilogy(logres_lgmres,color_lgmres)
% ## semilogy(logres_gmrese,color_gmrese)
% ## semilogy(logres_lgmrese,color_lgmrese)
%  semilogy(logres_pd_lgmrese,color_pd_lgmrese)
%  semilogy(logres_slgmrese,color_slgmrese)
%  hold on
%  semilogy(logres_pd_gmres,color_pd_gmres)
% title(Name_Matrix);
%  xlabel('Number of Restart Cycle');ylabel('|rj|/|r0|');
% legend(
% ##legend(['GMRES(30); cycle= ', num2str(sol1(1,1)),'; t= ', num2str(sol1(1,2)),' s.'], ...
% ##  ['LGMRES(28,2); cycle= ', num2str(sol2(1,1)),'; t= ', num2str(sol2(1,2)),' s.'], ...
% ##  ['GMRES-E(28,2); cycle= ', num2str(sol3(1,1)),'; t= ', num2str(sol3(1,2)),' s.'], ...
% ##  ['LGMRES-E(27,2,1); cycle= ', num2str(sol4(1,1)),'; t= ', num2str(sol4(1,2)),' s.'], ...
% legend( ['A-PD-LGMRES-E(28,2,2); cycle= ', num2str(sol5(1,1)),'; t= ', num2str(sol5(1,2)),' s.'], ...
%  legend(['SLGMRES-E-s(28,2,2); cycle= ', num2str(sol6(1,1)),'; t= ', num2str(sol6(1,2)),' s.'], ...
%  ['PD-GMRES(30;2;0.8); cycle= ', num2str(sol7(1,1)),'; t= ', num2str(sol7(1,2)),' s.'],'Location','SouthOutside');
%  hold on
