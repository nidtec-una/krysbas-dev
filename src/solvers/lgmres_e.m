%modificado por jccf octurbre 2018
% Clear all variables from the workspace
clear all
% Clear the screen
clc        
%function [tiempoC]=LGMRES_E_v2(A,b,color,~,mLE, lLE, dLE,itermax)
%Time Control
tic; 

% Setting tolerance parameters and maximum number of iterations
tol=1e-9;   % Tolerance for convergence
opts.tol=  tol;   % Tolerance options for the algorithm
%opts.tol= 1e-9;
% opts.disp= 0;
MAX_ITERATION=1000;       % Tolerance options for the algorithm
%maxit=itermax;

% Load the problem matrix from the 'sherman4.mat' file
load sherman4.mat;
A=Problem.A;     % System matrix
b=Problem.b;     % Right-hand side vector of the system

% Matrix name and color for graphics (if needed)
Name_Matrix='Scherman'; color='r';
%b=ones(size(A,1),1);

% Initialization of the initial solution and parameters for the LGMRES-E method
x0=zeros(size(A,1),1);
% m=mLE;
% d=dLE;
% lL=lLE;
% l=lL;

m=27;  % Number of GMRES iterations
d=2;   % Dimension for harmonic Ritz vectors
lL=2;  % Length of storage for vectors
l=lL;  % Initialization of storage length

% Total subspace size
s = m + d +l;
ST=s;
n=size (A,1);   % Matrix size
flag=0;         % Flag for convergence control

restart=1;

% First GMRES iteration (Restart)
disp('Iniciando primera ejecución de GMRES (Arnoldi)...');  % Mensaje de inicio
r=b-A*x0;         % Calculation of the initial residual
res(1,:)=norm(r); % Norm of the initial residual
%beta2=res(1,:);  
beta= res(1,1);   % Initial beta value
logres(1,:)= 1; %because norm(r)/norm(ro)=1    % Initialization of log values of the residual norm
%iter(1,:)=restart;
miteracion(1,1)=m;     % Number of iterations
z=[];               % Initialization of vector z

%Preallocating for speed
% Memory preallocation to improve efficiency
w=zeros(n,m);  % Vector w
h=zeros(m+1,m); % Matrix h
%z=zeros(n,1);
ii=1; % Variable for matrix z


%%%%%%%%%%%First run 

% First run of GMRES (Arnoldi)
    v(:,1)=r/beta;   % Normalization of the first vector
        for j=1:m                       %modified gram schmidt--Arnoldi
            w(:,j)=A*v(:,j); % Matrix-vector multiplication
            for i=1:j        % Modified Gram-Schmidt orthogonalization
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));  % Calculation of the norm
            if h(j+1,j)==0          % Stopping condition if there is an Arnoldi breakdown
                m=j;
                h2=zeros(m+1,m);    % Adjust the size of the h matrix
                for k=1:m
                    h2(:,k)=h(:,k);
                end
                h=h2;
            else        
                v(:,j+1)=w(:,j)/h(j+1,j); % Update of vector v
            end
        end
        H0=h;
        %hs=h;
        % Solution of the reduced system using QR decomposition
        g=zeros(m+1,1);  % Vector of constants
        g(1,1)=beta;     % Initialization of beta value

        for j=1:m                       %plane rotations (QR decompostion)
            disp('Se realiza QR descomposition');
            P=eye(m+1);       % Identity matrix
            sin=h(j+1,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));  % Sin calculation
            cos=h(j,j)/(sqrt(h(j+1,j)^2 + h(j,j)^2));    % Cos calculation
            P(j,j)=cos;
            P(j+1,j+1)=cos;
            P(j,j+1)=sin;
            P(j+1,j)=-sin;
            h=P*h;          % Update of the h matrix
            g=P*g;          % Update of vector g
        end

      % Calculation of the minimum norm minimizer
      disp('Se calcula el minimizador del mínimo de norma');
        R=zeros(m,m);  % Initialization of R
        G=zeros(m,1);  % Initialization of G
        V=zeros(n,m);  % Initialization of V
        for k=1:m
            G(k)=g(k); % Update of G
            V(:,k)=v(:,k); % Update of V
            for i=1:m
                R(k,i)=h(k,i);  % Update of R
            end
        end
        minimizer=R\G;  % Solution of the system
        Z=V*minimizer;  % Update of the approximate solution
        xm=x0+Z;        % Update of the current solution

        disp('Primera ejecución de GMRES completa.');  % Mensaje de finalización
        % Update of the residual and restart
        x0=xm;
        logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,1));  % Update of the residual norm
%Calculo de z(k) para iteración 1  (errors approximations)
        z(:,1)= Z;
%%%%%% Calculo de vectores armonicos de Ritz
H0(m+1,:)=[];   % Removal of the last row
Fold=H0';       % Transposed Fold matrix
opts.v0=ones(size(Fold,1),1);  % Initialization of vector v0 for eigs

disp('Calculando vectores armónicos de Ritz...');  % Mensaje informativo

% Calculation of harmonic Ritz values and vectors
gf=R'*R;   % Product of R transpose and R
E=zeros(m,d);
D=zeros(d,1);
%Compute the harmonic Ritz vectors
        [E2,D2]=eigs(Fold,gf,d,'LM',opts);  % Calculation of largest eigenvalues
        for pq=1:d
            D(pq,1)=abs(D2(pq,pq));  % Take the absolute values
        end
        [~,I]=sort(D,1);  % Sort the eigenvalues
        for qq=1:d
            E(:,qq)=E2(:,I(qq,1));  % Extract the corresponding eigenvectors
        end
        dy0= V*E;  % Calculation of dy0 vectors

        % Handling of complex numbers in Ritz vectors
        dy=zeros(n,1);
    %%%% if dy0 is complex%%%%%%%%%%%%%%%%%%%%%%
            if isreal(dy0)==0  % If the vectors are not real
                disp('Manejando vectores de Ritz complejos...');
                ij=1; jj=0;
                while size(dy,2)<=d && ij <=d
                    if isreal(dy0(:,ij))==0 && norm(real(dy0(:,ij)))>0
                        dy(:,jj+1)= real(dy0(:,ij)); jj=size(dy,2);
                        if ij <=d
                            dy(:,jj+1)= abs(imag(dy0(:,ij)))*sqrt(1); jj=size(dy,2);
                            if ij <d
                                if dy0(:,ij)== conj(dy0(:,ij+1))
                                    ij=ij+2;
                                else
                                    ij=ij+1;
                                end
                            end
                        end
                    else
                        dy(:,jj+1)= dy0(:,ij);
                        ij=ij+1;
                        jj=size(dy,2);
                    end
                end
            else
                dy=dy0;  % If the vectors are real, do nothing
            end
            disp('Vectores armónicos de Ritz calculados correctamente.');  % Mensaje informativo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ev0=[]; %para comando opts.v0
while flag==0
%    v=[];
    if ii ~= lL
        disp('Se realiza lgmres-e sin la totalidad de los errores de aproximación');
    else
        disp('Se realiza LGMRESE(m,d,l)');
    end
    r=b-A*x0;
    beta=norm(r);
    v(:,1)=r/beta; 
    h=zeros(m+1,m);
    %s1=s;
    if size(z,2)<lL
        if ii<=lL
            l=ii;
            ii=ii+1;
        end
%         P=[];
%         h=[];
            s=m+d+l;
            for j=1:s                       %modified gram schmidt--Arnoldi
%                 if j<=m
%                     w(:,j)=A*v(:,j);
%                 else
%                     w(:,j)=A*dy(:,j-m);
%                 end
                if j<=m
                    w(:,j)=A*v(:,j);
                end
                if j>=m+1 &&  j<=m+d
                    w(:,j)=A*dy(:,j-m);
                end
                if j>=m+ d + 1 %&& j<=s -1
                     w(:,j)=A*z(:,l-(j-m-d-1));
                end
                for i=1:j
                    h(i,j)=w(:,j)'*v(:,i);
                    w(:,j)=w(:,j)-h(i,j)*v(:,i);
                end
                h(j+1,j)=norm(w(:,j));
                if h(j+1,j)==0
                    s=j;
                    h2=zeros(s+1,s);   
                    for k=1:s
                        h2(:,k)=h(:,k);
                    end
                    h=h2;
                else       
                    v(:,j+1)=w(:,j)/h(j+1,j);
                end
            end       
       %     H0=h;
       %     hs=h;
            g=zeros(s+1,1);
            g(1,1)=beta;
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
            end
            R=zeros(s,s);
            G=zeros(s,1);
            V=zeros(n,s);
            for k=1:s
                G(k)=g(k);
                V(:,k)=v(:,k);
                for i=1:s
                    R(k,i)=h(k,i);
                end
            end
            minimizer=R\G;
            
            
            for k=m+1:m+d       %se incluye vectores armónicos de Ritz
                V(:,k)=dy(:,k-m);
            end
            aux=V*minimizer;
            xm=x0+aux;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         Z=z;
%         if l==1
%             z=aux;
%         else
%             for j=2:size(logres,1)-1
%                 z(:,j-1)=Z(:,j);
%             end
%             z(:,l)=aux;
%         end
        if size(z,2)<lL
            z(:,size(z,2)+1)=aux;
        else
            for j=2:lL
                z(:,j-1)=Z(:,j);
            end
                z(:,lL)=aux;
        end     
            miteracion(size(miteracion,1)+1,1)=m;
            logres(size(logres,1)+1,:)=abs(g(s+1,1)/res(1,1));
            if (abs (g(s+1,1)))/res(1,1) <tol  || size(logres,1)==MAX_ITERATION   %empleando �ltima componente de g como residuo
            %if (res(size(res,1)))/res(1,1) <tol      %empleando �ltima componente de g como residuo
                flag=1;
            else
                x0=xm;                        %update and restart
                restart=restart+1;
            end
    %%%%%%%%%%%%%%%%%%%%Compute the harmonic Ritz vectors
    W=V(1:n,1:s);
    Fold=W'*A'*W;
    gf=R'*R;
    opts.v0=ones(size(Fold,1),1);
%     if size(Ev0)~=0
%         if s >= size(Ev0,1)
%             Ev1=zeros(s,1);
%             Ev1(1:size(Ev0,1),:)=Ev0;
%         end
%         if s < size(Ev0,1)
%             Ev0((s+1):size(Ev0,1),:)=[];
%             Ev1=Ev0;
%         end
%         Ev0=real(Ev1);
%         opts.v0= Ev0; %modif. octubre
%         Ev1=[];
%     end
    E=zeros(s,d);
    D=zeros(d,1);
%%%Compute the harmonic Ritz vectors
        %[E,D2]=eigs(df,d,'SM', opts);
        [E2,D2]=eigs(Fold,gf,d,'LM',opts);
        for pq=1:d
            D(pq,1)=abs(D2(pq,pq));
        end
        [~,I]=sort(D,1);
        for qq=1:d
            E(:,qq)=E2(:,I(qq,1));
        end
        dy0= V*E;
   dy=zeros(n,1);
        %%%% if dy0 is complex%%%%%%%%%%%%%%%%%%%%%%    
            if isreal(dy0)==0
                ij=1; jj=0;
                while size(dy,2)<=d && ij <=d
                    if isreal(dy0(:,ij))==0 && norm(real(dy0(:,ij)))>0
                        dy(:,jj+1)= real(dy0(:,ij)); jj=size(dy,2);
                        if ij <=d
                            dy(:,jj+1)= abs(imag(dy0(:,ij)))*sqrt(1); jj=size(dy,2);
                            if ij <d
                                if dy0(:,ij)== conj(dy0(:,ij+1))
                                    ij=ij+2;
                                else
                                    ij=ij+1;
                                end
                            end
                        end
                    else
                        dy(:,jj+1)= dy0(:,ij);
                        ij=ij+1;
                        jj=size(dy,2);
                    end
                end
            else
                dy=dy0;
            end
            
    else % if size(logres,1)<=l;
 %-----------------------------------
        s=ST;
        l=lL;
        for j=1:s                       %modified gram schmidt--Arnoldi
            if j<=m
                w(:,j)=A*v(:,j);
            end
            if j>=m+1 &&  j<=m+d
                w(:,j)=A*dy(:,j-m);
            end
            if j>=m+ d + 1 %&& j<=s -1
                 w(:,j)=A*z(:,l-(j-m-d-1));
            end
            for i=1:j
                h(i,j)=w(:,j)'*v(:,i);
                w(:,j)=w(:,j)-h(i,j)*v(:,i);
            end
            h(j+1,j)=norm(w(:,j));
            if h(j+1,j)==0
                s=j;
                h2=zeros(s+1,s);    
                for k=1:s
                    h2(:,k)=h(:,k);
                end
                h=h2;
            else        
            v(:,j+1)=w(:,j)/h(j+1,j);
            end
        end
%         H0=h;
%         hs=h;
        g=zeros(s+1,1);
        g(1,1)=beta;

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
        end
        R=zeros(s,s);
        G=zeros(s,1);
        V=zeros(n,s);
        for k=1:s
            G(k)=g(k);
            V(:,k)=v(:,k);
            for i=1:s
                R(k,i)=h(k,i);
            end
        end
        for k=m+1:m+d
            V(:,k)=dy(:,k-m);
        end
        for k=m+d+1:s
            V(:,k)=z(:,l-(k-m-d-1));
        end
        minimizer=R\G;
        aux=V*minimizer;
        xm=x0+ aux;
        
        Z=z;
%         if l==1
%             z=aux;
%         else
%             for j=2:l
%                 z(:,j-1)=Z(:,j);
%             end
%             z(:,l)=aux;
%             
%         end
        if size(z,2)<lL
            z(:,size(z,2)+1)=aux;
        else
            for j=2:lL
                z(:,j-1)=Z(:,j);
            end
                z(:,lL)=aux;
        end   

        %res(restart+1,:)=(g(m+1,1));
        %iter(restart+1,:)=restart+1;
             %iter(restart+1,:)=restart+1;
            miteracion(size(miteracion,1)+1,1)=m;
            logres(size(logres,1)+1,:)=abs(g(s+1,1)/res(1,1));

            if (abs (g(s+1,1)))/res(1,1) <tol  || size(logres,1)==MAX_ITERATION   %empleando �ltima componente de g como residuo
            %if (res(size(res,1)))/res(1,1) <tol      %empleando �ltima componente de g como residuo
                flag=1;
   %             residuo= (abs (g(s+1,1)))/res(1,1);
            else
                x0=xm;                        %update and restart
                restart=restart+1;
            end
        %%%%%%%%%%%%%%%%%%%%Compute the harmonic Ritz vectors
        W=V(1:n,1:s);
        Fold=W'*A'*W;
        gf=R'*R; 
        opts.v0=ones(size(Fold,1),1);
        %gf=hs'*hs;
        %gf=W'*A'*A*W;
        %df=Fold\gf;
        %if size(logres,1)>3
        %    opts.v0=abs(E2(:,3));
        %else
        %    opts.v0= rand(s,1);
        %end
        E=zeros(s,d);
        D=zeros(d,1);
        %[E,D2]=eigs(df,d,'SM', opts);
        [E2,D2]=eigs(Fold,gf,d,'LM',opts);
        for pq=1:d
            D(pq,1)=abs(D2(pq,pq));
        end
        [~,I]=sort(D,1);
        for qq=1:d
            E(:,qq)=E2(:,I(qq,1));
        end
        dy0= V*E;
        dy=zeros(n,1);
        %%%% if dy0 is complex%%%%%%%%%%%%%%%%%%%%%%    
            if isreal(dy0)==0
                ij=1; jj=0;
                while size(dy,2)<=d && ij <=d
                    if isreal(dy0(:,ij))==0 && norm(real(dy0(:,ij)))>0
                        dy(:,jj+1)= real(dy0(:,ij)); jj=size(dy,2);
                        if ij <=d
                            dy(:,jj+1)= abs(imag(dy0(:,ij)))*sqrt(1); jj=size(dy,2);
                            if ij <d
                                if dy0(:,ij)== conj(dy0(:,ij+1))
                                    ij=ij+2;
                                else
                                    ij=ij+1;
                                end
                            end
                        end
                    else
                        dy(:,jj+1)= dy0(:,ij);
                        ij=ij+1;
                        jj=size(dy,2);
                    end
                end
            else
                dy=dy0;
            end


    end %if restarted
end  %while flag
tiempo=toc;     %Imprime tiempo de ejecuci�n
%finalizamos con el método
%subplot(1,1,1);
semilogy(logres,color)
% title(Name_Matrix);
% title(color);
% xlabel('Number of Restart Cycle');ylabel('|rj|/|r0|');
% legend(['PD-GMRES(27,alpha_{P}=', num2str(alpha0),',alpha_{D}=', num2str(delta0),'), t= ', num2str(tiempo)],'Location','Best');
%title(['Example 2.2 - Complementary cycles of GMRES. Nl=', num2str(Nl),'; delta=', num2str(dl)])
% hold on
%  subplot(2,1,2);
%  plot(miteracion,color)
%  xlabel('Number of restart cycles');ylabel('m, restart parameters');
%   hold on
lastcycle=size(logres,1);
tiempoC= [lastcycle tiempo];