function x = excited_thin_layer_residual(index, w, d, theta, r_exp, Eq_n0)

% For Thin film excited layer of thickness d 

PI = 3.141592653;
n1 = 1.0;                      %Vacuum
n2 = index(1) + 1i*index(2);   %Index of uniformly excited material


n3 = Eq_n0;                    %Index of equilibrium material
thetarad = PI*theta/180;       %Angle of incidence in radians
cost = cos(thetarad);
sint = sin(thetarad);


c = 2.99796e8;      % In m/s

%Reflection coefficient from vacumm to excited layer
r12 = ( n1*cost - sqrt(n2^2 - n1^2*sint^2 ))/( n1*cost + sqrt(n2^2 - n1^2*sint^2 ) );

%Reflection coefficient from excited layer to equilibrium layer
r23 = ( sqrt(n2^2 - n1^2*sint^2) - sqrt(n3^2 - n1^2*sint^2) )/( sqrt(n2^2 - n1^2*sint^2) + sqrt(n3^2 - n1^2*sint^2) );

%Phase shift after passing through excited layer
delta = (w*d/c)*sqrt(n2^2 - n1^2*sint^2);

r_theory = ( r12 + r23*exp(1i*2*delta) )/( 1 + r12*r23*exp(1i*2*delta) );

%{
if real(n2) < 0 || imag(n2) < 0
    r_theory = 10;
end
%}
M = real(r_theory) - real(r_exp);
A = imag(r_theory) - imag(r_exp);

%Returns value to be minimized
x = abs(M).^2 + abs(A).^2;


end

