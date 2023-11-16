%Specific Heat Analysis applied to La_2-xBa_xCuO_4(x=11.5%). Parameter values taken
% Solid State Communications, Vol. 74, No. 8, pp. 767-771, 199 (Specific
% Heat)
% Jpn. J. Applied Physics Vol. 26 L368-L370 1987 (Lattice dimensions)
% Energy units in mJ, distance units in cm.

A=9890000; %mJ/mol*K
D=129; %K
T_c = 17.5; %K
gamma_1 = 2.5; %mJ/mol*K^2
%gamma_2=10; %mJ/mol*K^2
beta    = 0.25; %mJ/mol*K^4

skin_depth = 400e-7; %[cm], 44 for approximated from optical conductivity data
r = 0.3;            %[cm] radius of pump beam
%fluence=1*10^-3; %mJ/cm^2
molar_volume = 59;   %118; %51 cm^3/mol
volume = skin_depth*(3.14159*r^2); % [cm^3] excited volume
R=0.15; % reflectivity, approximated from optical conductivity data


T_max=100;
T_initial=30;        % initial temperature
L=100;
num_T_Steps=1+(T_max-T_initial)*L;
F_Gradient= [0.050, 0.100, 0.200, 0.380, 0.760];  % Fluences to compute
T_final=zeros(length(F_Gradient),1);
absorbed_energy=( (1-R)*3.14159*r^2)*F_Gradient;
T=zeros(num_T_Steps,1);
C_e=zeros(num_T_Steps,1);
C_l=zeros(num_T_Steps,1);
E_e=zeros(num_T_Steps,1);
E_l=zeros(num_T_Steps,1);

for i=1:(num_T_Steps);
    T(i)=T_initial+((i-1)/L);
    C_l(i)=beta*(T(i)^3);
    C_e(i)=gamma_1*T(i); %(gamma_1*T(i)+(A)*exp(-D/T(i)))*heaviside(-(T(i)-T_c))+(gamma_2*T(i))*heaviside((T(i)-T_c));
    E_e(i)=(volume/molar_volume)*trapz(T(1:i),C_e(1:i),1);
    E_l(i)=(volume/molar_volume)*trapz(T(1:i),C_l(1:i),1);
    for j=1:length(F_Gradient);
        [val,T_f]=min(abs(E_e+E_l-absorbed_energy(j)));
        T_final(j)=T_initial+((T_f)-1)/L;
    end
    
    
end

% Plot Fluence vs final temperature
figure(99);
hold on;
plot(F_Gradient,T_final)