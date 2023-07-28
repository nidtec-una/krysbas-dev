function [tiempo_C, log_res] = PD_GMRES_m_1(A, b, m_PD, alpha, delta, ...
                                           iter_max)
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

tic; % Time Control 
tol = 1e-9;
[s, n] = size(A);
x_0 = zeros(size(b,1),1);
m_initial = m_PD;
m_min = 1;
m_max = n - 1; % se puede considerar que no tiene cota superior, 
% antes de las 1000 iteraciones no logra alcanzar m_max con 
% alpha = -3 y delta = 5
m_step = 1;
max_it = iter_max;
% OR
% max_it=1000;

%Parameters 
alpha_0 = alpha; % BORRAR?
delta_0 = delta; % BORRAR?

flag = 0; % FLAG DE CONVERGENCIA?
% convergencia = 0; % BORRAR?

if (s~=n)
    error ('Matrix not square');
end

[i, j] = size(b);

if (s~=i)
    error ('Vector does not match size of matrix A');
end

if (j~=1)
    error ('Vector is not a column vector')
end

if (size (b)~=size(x_0))
    error('Initial guess not right size');
end

restart = 1;
r0 = b - A *  x_0;
res(1, :) = norm(r0);
%beta0 = res(1,:);
% CREAR UN log_res DE TAMANHO MAX_ITER?
log_res(1, :) = (norm(r0)/res(1,1));
%log_res = [];
%log_res(1, :) = (norm(r0));

iter(1,:) = restart;
%COEF = [];
m_iteracion(1, 1) = m_initial;


while flag==0
    
    if iter(size(iter,1),:) ~=1
%         if Hs(1,:) < 0.1
%             alpha_0delta_0 alpha_0 + 1;
%         end
        [m_iter] = pd_rule(m, m_initial, m_min, res, iter(size(iter,1),:),...
                        m_step, m_max, alpha_0, delta_0); %cab
        m = m_iter(1,1);
        m_initial = m_iter(1,2);
    else
        m = m_initial;
    end
    m_iteracion(iter(size(iter, 1), :) + 1, 1) = m;

    v = zeros(n, m + 1);
    w = zeros(n, m);
    r = b - A * x_0;
    beta = norm(r);
    v(:, 1) = r / beta; 
    h = zeros(m + 1, m);

    for j = 1:m % Modified Gram Schmidt, or Arnoldi Method
        w(:, j) = A * v(:, j);
        for i = 1:j
            h(i, j) = w(:, j)' * v(:, i);
            w(:, j) = w(:, j) - h(i, j) *  v(:, i);
        end
        h(j + 1, j) = norm(w(:, j));
        if h(j + 1, j)==0
            m = j;
            h2 = zeros(m + 1, m);    %VERIFICAR!!!...
            for k = 1:m
                h2(:, k) = h(:, k);
            end
            h = h2;
        else        
        v(:, j + 1) = w(:, j)/h(j + 1, j);
        end
    end
    %Hs = h;
    g = zeros(m + 1, 1);
    g(1, 1) = beta;
 
    for j = 1:m   % Plane rotations (QR decompostion)
        P = eye(m + 1);   
        sin = h(j + 1, j)/(sqrt(h(j + 1, j)^2 + h(j, j)^2));
        cos = h(j, j)/(sqrt(h(j + 1, j)^2 + h(j, j)^2));
        P(j, j) = cos;
        P(j + 1, j + 1) = cos;
        P(j, j + 1) = sin;
        P(j + 1, j) = -sin;
        h = P *  h;
        g = P *  g;
    end
    R = zeros(m, m);
    G = zeros(m, 1);
    V = zeros(n, m);
    for k = 1:m
        G(k) = g(k);
        V(:, k) = v(:, k);
        for i = 1:m
            R(k, i) = h(k, i);
        end
    end
    minimizer = R\G;
    x_m = x_0 + V *  minimizer;
    res(restart + 1, :) = abs(g(m + 1, 1));
    iter(restart + 1, :) = restart + 1;
    log_res(size(log_res, 1) + 1, :) = abs(g(m + 1, 1)/res(1, 1));
%     %log_res(size(log_res, 1) + 1, :) = abs(g(m + 1, 1));
 

    if abs(g(m + 1, 1))/res(1, 1) < tol || size(log_res, 1)==max_it   %empleando ultima componente de g como residuo
    %if (abs (g(m + 1, 1))) <tol   || restart== 1000   %empleando ultima componente de g como residuo
        flag = 1;
 %       residuo = (abs (g(m + 1, 1)))/res(1, 1);
    else
        x_0 = x_m;                        %update and restart
        restart = restart + 1;
    end
  
end
%if (abs (g(m + 1, 1)))/res(1, 1) <tol      %Verification of convergence
%    convergencia = 1;                      % BORRAR?
%end

tiempo = toc;     %Imprime tiempo de ejecuci�n
last_cycle = size(log_res, 1);
tiempo_C = [last_cycle tiempo];