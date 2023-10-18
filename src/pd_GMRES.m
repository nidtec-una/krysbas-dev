function [tiempoC, logres, xx]=pd_GMRES(A,b, mPD, alpha, delta,itermax)
% Based on PD_GMRES_m_1.m by Juan Carlos Cabral
% Agregar tol? 
% 
% Description:
% 
% PD_GMRES_m_1 is a modified implementation of the restarted Generalized
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
%           irestarting parameter
%
% tol:      float
%           tolerance error threshold for relative residual norm
%           
% max_iter: int
%           maximum number of (inner?) iterations
%
%
% Output parameters:
%
% tiempo_C: float
%           Computational time of the algorithm
%
% log_res:  (1 up to to max_iter)-by-1 vector
%           relative residual norms
%
% References:
% 
% Nunez, R. C., Schaerer, C. E., & Bhaya, A. (2018). A proportional-derivative 
% control strategy for restarting the GMRES (m) algorithm. Journal of 
% Computational and Applied Mathematics, 337, 209-224.
% 

tic;         %Time Control 
%load wang4.mat; load b_wang4.mat;
tol=1e-9;
%A=Problem.A;
%b=Problem.b;
[s,n]=size (A);
% bmax=max(max(A));
% bmin=min(min(A));
% b=bmin*ones(n,1)+(bmax-bmin)*rand(n,1);

%b=b(:,1);
%b=ones(size(A,1),1);
x0=zeros(size(b,1),1);
minitial=mPD;
mmin=1;
mmax=n-1; %se puede considerar que no tiene cota superior, antes de las 1000 iteraciones  no logra alcanzar mmax con alpha=-3 y delta=5
mstep=1;
maxit=itermax;
%maxit=1000;

%Parameters 
alpha0=alpha;
delta0=delta;



flag=0;
%convergenica=0; 
if (s~=n)
    error ('Matrix not square');
end

[i,j]=size (b);

if (s~=i)
    error ('Vector does not match size of matrix A');
end

if (j~=1)
    error ('Vector is not a column vector')
end

if (size (b)~=size(x0))
    error('Initial guess not right size');
end

restart=1;
r0=b-A*x0;
res(1,:)=norm(r0);
%beta0=res(1,:);
logres(1,:)=(norm(r0)/res(1,1));
%logres=[];
%logres(1,:)=(norm(r0));

iter(1,:)=restart;
%COEF=[];
miteracion(1,1)=minitial;


while flag==0
    
    if iter(size(iter,1),:) ~=1
%         if Hs(1,:) < 0.1
%             alpha0= alpha0 +1;
%         end
        [miter]=pd_rule(m,minitial,mmin,res,iter(size(iter,1),:),mstep, mmax,alpha0, delta0); %cab
        m=miter(1,1);
        minitial=miter(1,2);
    else
        m=minitial;
    end
    miteracion(iter(size(iter,1),:)+1,1)=m;

    v=zeros(n,m+1);
    w=zeros(n,m);
    r=b-A*x0;
    beta=norm(r);
    v(:,1)=r/beta; 
    h=zeros(m+1,m);

    for j=1:m                       %modified gram schmidt--Arnoldi
        w(:,j)=A*v(:,j);
        for i=1:j
            h(i,j)=w(:,j)'*v(:,i);
            w(:,j)=w(:,j)-h(i,j)*v(:,i);
        end
        h(j+1,j)=norm(w(:,j));
        if h(j+1,j)==0
            m=j;
            h2=zeros(m+1,m);    %VERIFICAR!!!...
            for k=1:m
                h2(:,k)=h(:,k);
            end
            h=h2;
        else        
        v(:,j+1)=w(:,j)/h(j+1,j);
        end
    end
    %Hs=h;
    g=zeros(m+1,1);
    g(1,1)=beta;
 
    for j=1:m                       %plane rotations (QR decompostion)
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
%     %logres(size(logres,1)+1,:)=abs(g(m+1,1));
%Calculo de coeficientes de polinomiode GMRES(m)
%         n1=m;
%         C=[];
%         for i=1:n1
%             if i==1
%                 C(i,i)=1/beta;
%             end
%             
%             B1=[0; C(:,i)];
%             hs=[];
%             for j=1:i
%                 hs=Hs(1:j,j);
%             end
%             B2=[C*hs; 0];
%             C(i+1,i+1)=0;
%             C(:,i+1)=inv(Hs(i+1,i))*B1 - inv(Hs(i+1,i))*B2;
%         end
% %Calculos de cooeficientes
%         U=C;
%         C(:,m+1)=[];
%         C(m+1,:)=[];
%         coef=C*minimizer;
%         
% % if size(logres,1)==1
%         cf=size(coef,1);
%         for i=1:cf
%             CF(i,1)=-coef(cf+1-i,1);
%         end
%         
%         CF(cf+1,1)=1;
%         
%         
%         coef=[];
%         coef=CF';
%         COEF(size(COEF,1)+1,:)=coef;
% % Graficos de polinomios        
% %         x=linspace(0,5);
% % hold on 
% % grid on
% % for i=1:size(COEF,1)
% %      y=polyval(COEF(i,:),x);  
% %      plot(x,y)
% % end
% 
% 
%  COEF=[];
%     
    
    
    
    
    if abs (g(m+1,1))/res(1,1) <tol   || size(logres,1)==maxit   %empleando �ltima componente de g como residuo
    %if (abs (g(m+1,1))) <tol   || restart== 1000   %empleando �ltima componente de g como residuo
        flag=1;
        xx = xm;
 %       residuo= (abs (g(m+1,1)))/res(1,1);
    else
        x0=xm;                        %update and restart
        restart=restart+1;
    end
  
end
%if (abs (g(m+1,1)))/res(1,1) <tol      %Verification of convergence
%    convergencia=1;
%end

tiempo=toc;     %Imprime tiempo de ejecuci�n
color = 'r';
%semilogy(logres,color)
semilogy(iter,logres,'r')
%semilogy(logres,'g')
xlabel('Number of Restart Cycle');ylabel('|rj|/|r0|');
title('Residuals,  tol =10-8');

%tiempo=cputime - inicio;     %Imprime tiempo de ejecuci�n
figure(1)
subplot(1,1,1);
semilogy(logres,color)
lastcycle=size(logres,1);
tiempoC= [lastcycle tiempo];
legend(['PD-GMRES(30;2;0.8); cycle= ', num2str(lastcycle),'; t= ', num2str(tiempo),' s.'],'Location','SouthOutside');
 
title('sherman3');
 %title(color);
xlabel('Number of Restart Cycle');ylabel('|rj|/|r0|');
legend(['PD-GMRES(27,alpha_{P}=', num2str(alpha0),',alpha_{D}=', num2str(delta0),'), t= ', num2str(tiempo)],'Location','Best');
%title(['Example 2.2 - Complementary cycles of GMRES. Nl=', num2str(Nl),'; delta=', num2str(dl)])
hold on
figure(2)
subplot(1,1,1);
plot(miteracion,color)
xlabel('Number of restart cycles');ylabel('m, restart parameters');
hold on

%legend('Backer-GMRES(27,3)', 'Location','EastOutside');
%hold on



