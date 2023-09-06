% Modificado por geem (julio 2017)
% Gráfico GMRES: norma residual normalizada |rj|/|r0| en funcion del numero
% de ciclos
% Notación: R. Morgan, "GMRES with Deflated Restarting"
% Inputs include A,x0,b,m,d,tol

% clear all
function [vec_sol] = gmres_e(A,b, m1, d1, itermax, tol, color,print, Name_Matrix)
tic;         %Time Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tol=1e-10;
%opts.tol=  opts_tol;
opts.tol=  eps; % ó 1e-6, 1e-9, 1e-12
%maxit = 1000
maxit=itermax;
x0=zeros(size(A,1),1);
m=m1;
%d=nritz;
d=d1;
s = m + d;
n = size(A,1);
flag=0;
restart=1;
iter(1,:)=restart;
miteracion(iter(size(iter,1),:),1)=m;
r0=b-A*x0;
beta= norm(r0);
res(1,:)= beta;
logres(1,:)=abs(norm(r0)/res(1,1));
%logres(1,:)= beta;
%logres=[];
%iter(1,:)=restart;
flag1=0;    %define restart 1
%opts.tol = 1e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%r=b-A*x0;
%beta=norm(r0);
v(:,1)=r0/beta;
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
            for k=1:m
                h2(:,k)=h(:,k);
            end
            h=h2;
        else       
            v(:,j+1)=w(:,j)/h(j+1,j);
        end
    end
    H0=h;
    hs=h;
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
    x0=xm;   
    %res(restart+1,:)=abs(g(m+1,1));
    iter(restart+1,:)=restart+1;
    miteracion(iter(size(iter,1),:),1)=m;
    logres(size(logres,1)+1,:)=abs(g(m+1,1)/res(1,1));
    %logres(size(logres,1)+1,:)=abs(g(m+1,1));
    restart= restart +1;
    %iter(restart+1,:)=restart+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H0(m+1,:)=[];
Fold=H0';
opts.v0=ones(size(Fold,1),1);
%gf=hs'*hs;
gf=R'*R;
%gf=V'*A'*A*V;
%df=Fold\gf;
E=zeros(m,d);
D=zeros(d,1);

%gf=hs'*hs;
gf=R'*R;
%gf=V'*A'*A*V;
%df=Fold\gf;
%Compute the harmonic Ritz vectors
        %[E,D2]=eigs(df,d,'SM', opts);
        [E2,D2]=eigs(Fold,gf,d,'LM',opts);
        for pq=1:d
            D(pq,1)=abs(D2(pq,pq));
        end
        [Y,I]=sort(D,1);
        for qq=1:d
            E(:,qq)=E2(:,I(qq,1));
        end
        dy0= V*E;

        dy=[];
    %%%% if dy0 is complex%%%%%%%%%%%%%%%%%%%%%%
    if isreal(dy0)==0
        ij=1;
        while size(dy,2)<=d && ij <=d
            if isreal(dy0(:,ij))==0
                dy(:,size(dy,2)+1)= real(dy0(:,ij));
                if ij <=d
                    dy(:,size(dy,2)+1)= abs(imag(dy0(:,ij))*sqrt(1));
                    if ij <d
                        if dy0(:,ij)== conj(dy0(:,ij+1))
                            ij=ij+2;
                        else
                            ij=ij+1;
                        end
                    end
                end
            else
                dy(:,size(dy,2)+1)= dy0(:,ij);
                ij=ij+1;
            end
        end
    else
        dy=dy0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E=[];

while flag==0
    r=b-A*x0;
    beta=norm(r);
    v(:,1)=r/beta;
    for j=1:s                       %modified gram schmidt--Arnoldi
        if j<=m
            w(:,j)=A*v(:,j);
        else
            w(:,j)=A*dy(:,j-m);
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
    H0=h(1:s,1:s);
    hs=h;
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
    for k=m+1:s
        V(:,k)=dy(:,k-m);
    end
    
    xm=x0+V*minimizer;
    %res(restart+1,:)=abs(g(s+1,1));
    iter(restart+1,:)=restart+1;
    miteracion(iter(size(iter,1),:),1)=m;
    logres(size(logres,1)+1,:)=abs(g(s+1,1)/res(1,1));
    %logres(size(logres,1)+1,:)=abs(g(s+1,1));
    %logres(restart+1,:)=abs(g(s+1,1));

    if (abs (g(s+1,1)))/res(1,1) <tol  || size(logres,1)==maxit   %empleando ï¿½ltima componente de g como residuo
    %if (abs (g(s+1,1))) <tol  || size(logres,1)==maxit    %empleando ï¿½ltima componente de g como residuo
        flag=1;
        %residuo= (abs (g(s+1,1)))/res(1,1);
        residuo= abs (g(s+1,1));        
    else
        x0=xm;                        %update and restart
        restart=restart+1;
    end
  
    %%%%%Compute the harmonic Ritz vectors
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
    
  %Compute the harmonic Ritz vectors
        %[E,D2]=eigs(df,d,'SM', opts);
        [E2,D2]=eigs(Fold,gf,d,'LM',opts);
        for pq=1:d
            D(pq,1)=abs(D2(pq,pq));
        end
        [Y,I]=sort(D,1);
        for qq=1:d
            E(:,qq)=E2(:,I(qq,1));
        end
        dy0= V*E;
        dy=[];
        
    %%%% if dy0 is complex%%%%%%%%%%%%%%%%%%%%%%
    if isreal(dy0)==0
        ij=1;
        while size(dy,2)<=d && ij <=d
            if isreal(dy0(:,ij))==0
                dy(:,size(dy,2)+1)= real(dy0(:,ij));
                if ij <=d
                    dy(:,size(dy,2)+1)= abs(imag(dy0(:,ij))*sqrt(1));
                    if ij <d
                        if dy0(:,ij)== conj(dy0(:,ij+1))
                            ij=ij+2;
                        else
                            ij=ij+1;
                        end
                    end
                end
            else
                dy(:,size(dy,2)+1)= dy0(:,ij);
                ij=ij+1;
            end
        end
    else
        dy=dy0;
    end
end  %while flag

t1 = toc;     %Imprime tiempo de ejecuciï¿½n
if print == 1
    figure(3)
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Number of Restart Cycles');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
end
save('logres_GMRESE.mat','logres');
vec_sol = [t1 restart m*restart+d*(restart-1)];
