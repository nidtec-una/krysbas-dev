% Modificado por geem (junio 2023)
function [vec_sol, vecnormy, xx]=mi_LGMRESE(A,b,mLE, dLE, lLE, ...
    max_ciclos, tol, alpha, delta, color, print, Name_Matrix) %Time Control 
n = size(A,1);
opts.tol= eps;
maxit=max_ciclos;
x0=zeros(n,1);
m=mLE;
d=dLE;
l=lLE;
s = m + d +l;
ST=s;
m_max = 96;
flag=0;
flag2=0;
restart=1; % Primer ciclo de reinicio
r=b-A*x0;
res(1,:)=norm(r);
beta= res(1,1);
logres(1,:)= 1; %norm(r)/norm(ro)=1 
iter(1,:)=restart;
log_de_m(iter(size(iter,1),:),1)=m;
log_de_s(iter(size(iter,1),:),1)=m; % inicialmente no hay d ni l
z=[];
% SOLVER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    tic; 
    v(:,1)=r/beta; 
    for j=1:m % gram schmidt--Arnoldi
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
    g=zeros(m+1,1);
    g(1,1)=beta;

    for j=1:m % min-cuadrados - plane rotations (QR decompostion)
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
%   H0=h;
%   hs=h;
     
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
    yj=R\G;
    norm_y = norm(yj);
    Z=V*yj;
    norm_z=norm(z);
    xm=x0+Z;
    x0=xm;
    
    if norm_y < delta  %para el caso de estancamiento en el primer ciclo.
        flag2=1; % Los vectores de error de aproximación no sirven
    else
        flag2=0;
    end  
    logres(size(logres,1)+1,:)=norm(b-A*xm)/res(1,1);
    lognormy(size(logres,1)+1,:)=norm_y;
    
%%%%% Calculo de z(k) para iteracion 1
    z(:,size(z,2)+1)= Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculo de vect armonicos de Ritz  %%%%%%%%%%%%%
H0(m+1,:)=[];
Fold=H0';
opts.v0=ones(size(Fold,1),1);
gf=R'*R;
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
%%%% if dy0 is complex%%%%
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
    E=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ev0=[]; %para comando opts.v0
while flag==0   %Mientras no se alcance convergencia.
    if flag2 == 1  % GMRES-E, no agrega vectores de error de aprox
          if m < m_max
            m = m + alpha;
          else
            m = m_max;
          end
          s = m + d;
          ST = s;
    end
    v=[];
    r = b-A*x0;
    beta = norm(r);
    v(:,1)=r/beta; 
    %h=zeros(m+1,m);
    %s1=s;
        P=[];
        h=[];
        W=[];
        V=[];
    if size(logres,1) <= l+1
       % P=[];
       % h=[];
        %W=[];
            s=s-l; % ese menos l
            for j=1:s  %modified gram schmidt--Arnoldi
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
      
            H0=h;
            hs=h;
            g=zeros(s+1,1);
            g(1,1)=beta;

            for j=1:s     %plane rotations (QR decompostion)
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
            yj=R\G;
            norm_y = norm(yj,inf);
            for k=m+1:s
                V(:,k)=dy(:,k-m);
            end
            
            aux=V*yj;
            xm=x0+aux;
       
       
        Z=z;
        if l==1
            z=aux;
        else
            for j=2:size(logres,1)-1
                z(:,j-1)=Z(:,j);
            end
            z(:,l)=aux;
            
        end

            iter(restart+1,:)=restart+1;
            log_de_m(iter(size(iter,1),:),1)=m;
            log_de_s(iter(size(iter,1),:),1)=s;
            logres(size(logres,1)+1,:)=norm(b-A*xm)/res(1,1);
            lognormy(size(logres,1)+1,:)=norm_y;
            
            if norm(b-A*xm)/res(1,1) <tol  || size(logres,1)== maxit 
            
                flag=1;
            else
                x0=xm;                        %update and restart
                restart=restart+1;
                if norm_y < delta
                  flag2=1;
                else
                  flag2=0;
                end 
            end



    %%%%%%%%%%%%%%%%%%%%Compute the harmonic Ritz vectors
    W=V(1:n,1:s);
    Fold=W'*A'*W;
    gf=R'*R;
    opts.v0=ones(size(Fold,1),1);

    if size(Ev0)~=0
        if s >= size(Ev0,1)
            Ev1=zeros(s,1);
            Ev1(1:size(Ev0,1),:)=Ev0;
        end
        if s < size(Ev0,1)
            Ev0((s+1):size(Ev0,1),:)=[];
            Ev1=Ev0;
        end
        Ev0=real(Ev1);
        opts.v0= Ev0; %modif. octubre
        Ev1=[];
    end
     %Compute the harmonic Ritz vectors
        %[E,D2]=eigs(df,d,'SM', opts);
        E=[];
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
        
    %%%% if dy0 is complex%%%%
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

  
    else % if size(logres,1)>=l;
        s=ST;
%         i=[];
%         j=[];
        for j=1:s   % Gram-Schmidt--Arnoldi
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
        H0=h;
        hs=h;

        g=zeros(s+1,1);
        g(1,1)=beta;

        for j=1:s   %plane rotations (QR decompostion)
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
        
        yj=R\G;
        norm_y = norm(yj);
        aux=V*yj;
        xm=x0+ aux;
        Z=z;
        if l == 1
            z = aux;
        else
            for j=2:l
                z(:,j-1)=Z(:,j);
            end
            z(:,l)=aux;
        end
        %res(restart+1,:)=(g(m+1,1));
        %iter(restart+1,:)=restart+1;
        iter(restart+1,:)=restart+1;
        log_de_m(iter(size(iter,1),:),1)=m;
        log_de_s(iter(size(iter,1),:),1)=s;
        %logres(size(logres,1)+1,:)=abs(g(s+1,1)/res(1,1));
        logres(size(logres,1)+1,:)=norm(b-A*xm)/res(1,1);
        lognormy(size(logres,1)+1,:)=norm_y;

            if norm(b-A*xm)/res(1,1) <tol  || size(logres,1)==maxit   
                %empleando ï¿½ltima componente de g como residuo
            %if (res(size(res,1)))/res(1,1) <tol      
                %empleando ï¿½ltima componente de g como residuo
                
                flag=1;
                residuo= norm(b-A*xm)/res(1,1);
            else
                x0=xm;                        %update and restart
                restart=restart+1;
                if norm_y < delta 
                  flag2=1;
                  %m = m  + alpha;
                  %s = m + d +l;
                else
                  flag2=0;
                end
            end
        %%%%%%%%%%%%%%%%%%%%Compute the harmonic Ritz vectors
        W=V(1:n,1:s);
        Fold=W'*A'*W; %F_old
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
        
        E=[];
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end %if restarted
    
end  %while flag
tiempo=toc; % Tiempo de ejecuciï¿½n

if print == 1
    figure(1)
    subplot(2,1,1);
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Cantidad de ciclos de reinicio');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(2,1,2);
    semilogy(lognormy(3:size(lognormy)),color)
    ylabel('||y_j||');
    xlabel('Cantidad de ciclos de reinicio')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    figure(2)
    subplot(2,1,1);
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Cantidad de ciclos de reinicio');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(2,1,2);
    plot(log_de_s, color)
    ylabel('s_j');
    xlabel('Cantidad de ciclos de reinicio')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    figure(3)
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    xlabel('Number of Restart Cycles');ylabel('||r_j||/||r_0||');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    figure(4)
    subplot(2,1,1);
    plot(log_de_s, color)
    ylabel('s_j');
    xlabel('Number of Restart Cycles')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
   
    subplot(2,1,2);
    semilogy(lognormy(3:size(lognormy)),color)
    ylabel('|| y_j ||');
    xlabel('Cantidad de ciclos de reinicio')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    figure(5)
    subplot(3,1,1);
    semilogy(logres,color,'LineWidth',2)
    title(Name_Matrix)
    ylabel('|| r_j ||_2 / || r_0 ||_2');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(3,1,2);
    plot(log_de_s, color)
    ylabel('m_j (s_j)');
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on
    
    subplot(3,1,3);
    semilogy(lognormy(3:size(lognormy)),color)
    ylabel('|| y_j ||');
    xlabel('Cantidad de ciclos de reinicio')
    %xlabel('Cantidad de ciclos de reinicio');ylabel('m (s)');
    hold on  
end
save('logres_LGMRESE.mat','logres');
sum_s = sum(log_de_s);
vec_sol = [tiempo restart sum_s];
vecnormy = lognormy;
xx = xm;