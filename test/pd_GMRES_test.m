load 'matrices\sherman5.mat';
A=Problem.A; 
b=Problem.b; 
% OR
% b=b(:,1);
% OR
% b=ones(size(A,1),1);



%Calculo de coeficientes de polinomio de GMRES(m)
%         n1 = m;
%         C = [];
%         for i = 1:n1
%             if i==1
%                 C(i, i) = 1/beta;
%             end
%             
%             B1 = [0; C(:, i)];
%             hs = [];
%             for j = 1:i
%                 hs = Hs(1:j, j);
%             end
%             B2 = [C *  hs; 0];
%             C(i + 1, i + 1) = 0;
%             C(:, i + 1) = inv(Hs(i + 1, i)) *  B1 - inv(Hs(i + 1, i)) *  B2;
%         end
% %Calculos de cooeficientes
%         U = C;
%         C(:, m + 1) = [];
%         C(m + 1, :) = [];
%         coef = C *  minimizer;
%         
% % if size(log_res, 1)==1
%         cf = size(coef, 1);
%         for i = 1:cf
%             CF(i, 1) = -coef(cf + 1 - i, 1);
%         end
%         
%         CF(cf + 1, 1) = 1;
%         
%         
%         coef = [];
%         coef = CF';
%         COEF(size(COEF, 1) + 1, :) = coef;
% % Graficos de polinomios        
% %         x = linspace(0, 5);
% % hold on 
% % grid on
% % for i = 1:size(COEF, 1)
% %      y = polyval(COEF(i, :), x);  
% %      plot(x, y)
% % end
% 
% 
%  COEF = [];
%    

%semilogy(log_res, color)
%semilogy(iter, log_res, 'r')
%semilogy(log_res, 'g')
%xlabel('Number of Restart Cycle');ylabel('|rj|/|r0|');
%title('Residuals,  tol = 10-8');

%tiempo = cputime - inicio;     %Imprime tiempo de ejecuci�n
%figure(1)
%subplot(1, 1, 1);
 %semilogy(log_res, color)

 %legend(['PD-GMRES(30;2;0.8); cycle =  ', num2str(last_cycle), ...
 %'; t =  ', num2str(tiempo), ' s.'],'Location','SouthOutside');
 
 %title(Name_Matrix);
 %title(color);
% % xlabel('Number of Restart Cycle');ylabel('|rj|/|r0|');
% % legend(['PD-GMRES(27,alpha_{P} = ', num2str(alpha_0),',alpha_{D} = ', num2str(delta_0),'), t =  ', num2str(tiempo)],'Location','Best');
% %title(['Example 2.2 - Complementary cycles of GMRES. Nl = ', num2str(Nl),'; delta = ', num2str(dl)])
 %hold on
% figure(2)
% %subplot(1,1,1);
 %plot(m_iteracion,color)
 %xlabel('Number of restart cycles');ylabel('m, restart parameters');
 %hold on

%legend('Baker-GMRES(27,3)', 'Location','EastOutside');
%hold on