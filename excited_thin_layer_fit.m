function solved_index = excited_thin_layer_fit(start_at_point,go_to_point,w,d,theta,r_exp,Eq_n0,n1_guess,k1_guess)

Tol = 1e-7;
options = optimset('LargeScale','off','TolX',Tol,'TolFun',Tol,'MaxFunEvals',Inf,'MaxIter',Inf);%,'PlotFcns',@optimplotfval);

count = 1;

colors = colormap(jet(go_to_point - start_at_point+1));
n1_guess = 20;
k1_guess = 100;


    for ii = start_at_point:go_to_point
        x = fminsearch(@(n)excited_thin_layer_residual(n,w(ii),d,theta, ...
                r_exp(ii),Eq_n0(ii)), [n1_guess k1_guess],options);
        if x(1) < 0  && x(2) < 0
            x(1) = -x(1);
            x(2) = -x(2);
        end
            
        results(count) = x(1) + 1i*x(2);
        
        n1_guess = x(1);
        k1_guess = x(2);
        
%          figure(333)
%          hold on;
%          
% %          plot(x(1),x(2),'LineStyle','none','Marker','*','Color',colors(count,:));
        
        count = count+1;
        
        
        
        
        
        %%% All this was to look at a colormap of to see where the fits
        %%% found the local minimum
        %%% 
        
%         residuals(count) = excited_residual([x(1), x(2)],w(ii),d,theta,r_exp(ii),Eq_n0(ii));
%         r_theory = test_r_coeff([x(1) x(2)],w(ii),d,theta,r_exp(ii),Eq_n0(ii));
%         
%         nr = linspace(-10,200,100);
%         ni = linspace(-20,200,100);
%         for k=1:length(ni)
%             
%             for j=1:length(nr)
%                 rtheory = test_r_coeff([nr(j) ni(k)],w(ii),d,theta,r_exp(ii),Eq_n0(ii));
%                 A = real(r_exp(ii)) - real(rtheory);
%                 B = imag(r_exp(ii)) - imag(rtheory);
%                 tot = A^2 + B^2;
%                 if j==1;
%                     Ztemp = tot;
%                 else
%                     Ztemp = [Ztemp,tot];
%                 end
%             end
%             
%             
%             if k == 1;
%                 Z = Ztemp;
%             else
%                 Z = [Z;Ztemp];
%             end
%         end
%         %[0,0.001,0.003,0.005,0.01,0.05,0.1,0.5,1]
%         
%         [X,Y] = meshgrid(nr,ni);
%         figure(303)
%         contour(X,Y,Z,[0,0.001,0.003,0.005,0.01,0.05,0.1,0.5,1]);
%         hold on
%         plot(x(1),x(2),'r','Marker','o','MarkerSize',15,'MarkerFaceColor','r');
%         hold off
%        
%         figure(88)
%         hold on;
%         plot(w(ii)/(2*pi),abs(x(1)),'b','Marker','*','MarkerSize',10);
%         title('Abs n_fit');
%         
%         figure(89)
%         hold on;
%         plot(w(ii)/(2*pi),angle(x(1)+1i*x(2)),'b','Marker','*','MarkerSize',10);
%         title('Phase of n_fit');
%         
%         
%         figure(90)
%         hold on
%         plot(x(1),x(2),'b','Marker','*','MarkerSize',10);
%         title('n_fit complex plane');
%         
%         figure(92)
%         hold on
%         plot(real(r_exp(ii)),imag(r_exp(ii)),'b','Marker','*','MarkerSize',10);
%         plot(real(r_theory),imag(r_theory),'r','Marker','*','MarkerSize',10);
%         title('Complex r_{exp} and r_{theory}');
        
        
        
        
    end

  
solved_index = results;


end

