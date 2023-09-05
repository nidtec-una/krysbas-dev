% Modificado por geem (mayo 2017)
% Gr�fico GMRES: norma residual normalizada |rj|/|r0| en funcion del numero
% de ciclos
% Notaci�n: Baker, Jessup & Manteuffel, 
% "Accelerating the Convergence of Restarted GMRES" ====> LGMRES(m,k)
% Comentarios: ...
%clear all
function [vec_sol]=mi_LGMRES(A, b, m1, k1, itermax, tol, color, print, Name_Matrix)
maxit=itermax;
x0=zeros(size(A,1),1);
m = m1;
k = k1;
s = m + k;
[~,n]=size (A);
flag=0;
restart=1; % Numero de iteraciones
r=b-A*x0;
norma_r0=norm(r);
logres(1,:)= 1; % Residuos normalizados |r(j)|/|r(0)| = 1 cuando j = 0s
iter(1,:)=restart;
miteracion(iter(size(iter,1),:),1)=m;
z=[];
%H1s=[];
tic;         %Time Control 
while flag==0
    beta=norm(r);
    v(:,1)=r/beta; 
    h=zeros(m+1,m); % Matriz de Hessenberg
    if size(logres,1)<=k; % Se ejecuta solamente GMRES(m), ej. GMRES(27)
        for j=1:m                       %modified gram schmidt--Arnoldi
            w(:,j)=A*v(:,j);
            for i=1:j
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));
            if h(j+1,j)==0
                m=j;
                h2=zeros(m+1,m);    
                for kp=1:m
                    h2(:,kp)=h(:,kp);
                end
                h=h2;
            else        
                v(:,j+1)=w(:,j)/h(j+1,j);
            end
        end
        %Hs=h; % Omitir?
        g=zeros(m+1,1);
        g(1,1)=beta;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%         O directamente: [Q1 R1] = qr(h) y copiar R <== R1
%         [Q1 R1] = qr(h)
%         R=zeros(m,m);
%         G=zeros(m,1);
%         V=zeros(n,m);
%         for kp=1:m
%             G(kp)=g(kp);
%             V(:,kp)=v(:,kp);
%             for i=1:m
%                 R(kp,i)=h(kp,i);
%             end
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R=zeros(m,m);
        G=zeros(m,1);
        V=zeros(n,m);
        for kp=1:m
            G(kp)=g(kp);
            V(:,kp)=v(:,kp);
            for i=1:m
                R(kp,i)=h(kp,i);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        minimizer=R\G;
        Z=V*minimizer; %error de aproximacion
        xm=x0 + Z;
        r=b-A*xm;
        %res(restart+1,:)=abs(g(m+1,1));
        iter(restart+1,:)=restart+1;
        miteracion(iter(size(iter,1),:),1)=m;
        logres(size(logres,1)+1,:)=abs(g(m+1,1)/norma_r0);

        if abs(g(m+1,1))/norma_r0 <tol %| iter(size(iter,1),:)==00
            flag=1;
        else
            x0=xm;                        %update and restart
            restart=restart+1;
        end
        
        %Calculo de z(k)
        z(:,size(z,2)+1)= Z;

    else % Se ejecuta el LGMRES(m,k); k != 0
        for j=1:s                       %modified gram schmidt--Arnoldi
                if j<=m
                    w(:,j)=A*v(:,j);
                else
                    w(:,j)=A*z(:,k-(j-m-1));
                end
            for i=1:j
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));
            if h(j+1,j)==0
                s=j;
                h2=zeros(s+1,s);    
                for kp=1:s
                    h2(:,kp)=h(:,kp);
                end
                h=h2;
            else        
            v(:,j+1)=w(:,j)/h(j+1,j);
            end
        end

        g=zeros(s+1,1);
        g(1,1)=beta;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j=1:s                       %plane rotations (QR decompostion)
            P=eye(s+1);   
            sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2)); 
            cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2)); 
            P(j,j)=cos; 
            P(j+1,j+1)=cos;
            P(j,j+1)=sin; 
            P(j+1,j)=-sin;
            h=P*h;
            g=P*g;
        end                              % Reemplazar?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        R=zeros(s,s);
        G=zeros(s,1);
        V=zeros(n,s);
        for kp=1:s
            G(kp)=g(kp);
            V(:,kp)=v(:,kp);
            for i=1:s
                R(kp,i)=h(kp,i);
            end
        end
        for kp=m+1:s
            V(:,kp)=z(:,k-(kp-m-1));
        end
        minimizer=R\G;
        xm=x0+V*minimizer;
        r=b-A*xm;
        %res(restart+1,:)=(norm(r));
        iter(restart+1,:)=restart+1;
        miteracion(iter(size(iter,1),:),1)=m;
        %logres(restart+1,:)=abs(norm(r)/norma_r0);
        logres(size(logres,1)+1,:)=abs(g(s+1,1)/norma_r0);
        aux=V*minimizer;
        Z=z;
        if k==1
            z=aux;
        else
            for j=2:k
                z(:,j-1)=Z(:,j);
            end
            z(:,k)=aux;
            
        end
   
        if abs(g(s+1,1))/norma_r0 <tol || size(logres,1)==maxit
            flag=1;
            residuo1= (abs (g(s+1,1)))/norma_r0;
        else
            x0=xm;                        %update and restart
            restart=restart+1;
        end
       
    end %if restarted
   
end  %while flag
t1 = toc;     %Imprime tiempo de ejecuci�n
if print == 1
    figure(3)
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Number of Restart Cycles');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
end

save('logres_LGMRES.mat','logres');
vec_sol = [t1 restart m*restart+k1*(restart-k1)];