function nrandnk = Excited_thin_layer_reflection_THz(Eq_file_name,Ref_file_name,DeltaEfile,d,e_inf,min_freq,max_freq,color,deltaPhaseShift)

% Author: Kevin Cremin
% email: kacremin@physics.ucsd.edu
% 
% This program is set up to take THz TDS data and calculate the optical
% parameters in equilibrium as well as after photo-excitation.
% Due to the deep penetration depth of THz, the photo-excited material is
% modeled as a thin excited layer sitting on top of the bulk containing the
% equilibrium response. In the gradient model, this thin layer is divided
% into many thin layers and implemented with a transfer matrix method. 

% Function Input variables description
% DeltaEfile        = File name of deltaE TDS scan
% d                 = Penetration depth of pump [m]
% e00               = epsilon infinity
% min_freq          = minimum frequency to be evaluated [Hz]
% max_freq          = maximum frequency to be evaluated [Hz]
% color             = color used in plots
% deltaPhaseShift   = Amount to shift between two reference TDS scans at different Temps [ps]

% T_sc              = deltaE pump induced TDS scan
% T_sub             = Equilibrium TDS scan
% T_sub30K          = Equilibrium TDS scan at higher T where opt parameters are known 

%-- PHOTO-EXCITATION MODEL -----------------------------------------
ThinFilm    = 1;    % Model excited layer as uniform film on bulk
Gradient    = 0;    % Model excited layer as a gradient of many thin layers
%-------------------------------------------------------------------

%-- SAVING & PLOTTING PARAMETERS -----------------------------------
OptParamFig     = 888;     % Figure number to plot Opt parameters (i.e. R, Loss Function, sigma)
SaveParameters  = 0;        % Save optical parameters from fits to a file? (same directory as input file)
PlotxRange      = [0 7];    % Frequency range to plot
%---------------------------------------------------------------------


%-- CONSTANTS & EXPERIMENTAL PARAMETERS ------------------------------
e0              = 8.8542e-12;       % vacuum permitivity [m^-3 kg^-1 s^4 A^2] or [Ohm^-1 m^-1 s]
c               = 2.9979e8;         % speed of light in vacuum [m s^-1]
theta           = 30.;              % Angle of incidence for THz reflection setup
theta_rad       = theta*pi/180;     % Angle in radians
points_to_avg   = 10;               % simply the number of points to determine baseline avg
%ScanCutOff      = 200;
%---------------------------------------------------------------------


%-- LOAD EXPERIMENTAL DATA -------------------------------------------
%%% Data file names for Eq and reference scans
%Eq_file_name    = 'GaAs_TDS_Static_1.txt';
%Ref_file_name   = 'CopperRef_11.txt';

%%%% Load data files 
disp('Loading Delta E scan')    % Load delta E data file
Sig     = load(DeltaEfile);
T_sc    = Sig(:,2);           % Minus sign becauase sample and ref should have sign flip
t1      = Sig(:,1);

%%% Load Equilibrium THz scan
disp('Loading Equilibrium scan') 
Eq_file = load(Eq_file_name);
T_sub   = Eq_file(:,2);                   % Eq sample also gets sign flip
t       = Eq_file(:,1);                     
t       = t-min(t);                         

%%% Load Reference THz scan 
disp('Loading Reference Scan');
Ref_file = load(Ref_file_name);
T_ref    = Ref_file(:,2);

%%% Baseline subtraction for delta E scan %%%
tst1 = sum(T_sc(1:points_to_avg))/points_to_avg;    %Baseline Average of first few data points
T_sc = T_sc - tst1;                                 %Baseline Subtraction

%%% Baseline subtraction for the sample/ref scans %%%
tst2  = sum(T_sub(1:points_to_avg))/points_to_avg;
T_sub = T_sub - tst2;

tst3  = sum(T_ref(1:points_to_avg))/points_to_avg;
T_ref = T_ref - tst3;
%---------------------------------------------------------------------

 
%-- PHASE SHIFT ADJUSTMENT ------------------------------------------
% This is to correct for any temporal offset caused by referencing with 
% a gold mirror. DeltaPhaseshift is a variable that can be played with
% to give real physical values for optical parameters

t2 = t + deltaPhaseShift;  
t3 = linspace(min(t),max(t),length(t));
temp = interp1(t2,T_ref,t3,'spline');
T_ref = temp';
%--------------------------------------------------------------------


%-- TIME DOMAIN WINDOWING -------------------------------------------
% Window of "1" would be no windowing. Hanning windowing seems the most
% appropriate, although note that the matlab windowing functions are
% centered with respect to the data array, and may not be peaked where the 
% THz pulse is peaked in time. 

wind_sc         = tukeywin(length(T_sc),0.25);%hann(length(T_sc));
wind_sub        = tukeywin(length(T_sc),0.25);%hann(length(T_sc));
T_sc            = T_sc.*wind_sc;
T_sub           = T_sub.*wind_sub;
T_ref           = T_ref.*wind_sub;     % Reference pulse gets a sign flip, since we take the pulse reflected from gold 
                                        % as the same as the incinident pulse(reference) which would  
%--------------------------------------------------------------------


%-- PERFORM FOURIER TRANSFORM ---------------------------------------
disp('Initializing variables for FFT')

scan_length = length(T_sc);
zero_pad = 1024;%4096;                  % fixed value of 1024 seems to give desired spectral resolution
time_length = t(scan_length) - t(1);    %assumes in ps
TS = (zero_pad/scan_length)*time_length*1e-12; % length of zero padded scan in seconds
dt = TS/zero_pad;
time = 0:dt:(zero_pad-1)*dt;
dv = 1/TS;
v = 0:dv:(zero_pad-1)*dv;               % Frequency array in [THz]
w = v*2*pi;                             % Angular frequency
w = w';


%%% Perform FFT on TDS scans
disp('Doing FFT')
Tfft_ref    = fft(T_ref,zero_pad);
Tfft_sub    = fft(T_sub,zero_pad);%;*1.07;
Tfft_sc     = fft(T_sc,zero_pad);

%%% Complex Conjugate of fourier transform
%%% This is needed to get physical results for optical parameters
Tfft_ref    = real(Tfft_ref) - 1i*imag(Tfft_ref);
Tfft_sub    = real(Tfft_sub) - 1i*imag(Tfft_sub);
Tfft_sc     = real(Tfft_sc)  - 1i*imag(Tfft_sc);
%---------------------------------------------------------------------


%-- SOLVE FOR INDEX OF REFRACTION ------------------------------------
%%% Calc dE/E
dEE = Tfft_sc./Tfft_sub;       

%%% Calc Eq reflection coefficient r = E_reflected/E_incident
Eq_r0 = Tfft_sub./Tfft_ref;




%%% Solve for Equilibrium index (two solutions to sqrt)
Eq_n0_1 = sqrt( sin(theta_rad)^2 + cos(theta_rad)^2 * ((1 - Eq_r0)./(1+Eq_r0)).^2);
Eq_n0_2 = -sqrt( sin(theta_rad)^2 + cos(theta_rad)^2 * ((1 - Eq_r0)./(1+Eq_r0)).^2);

%%% Assume first solution (Checks sqrt solutions to obtain physical solution i.e. real(n)>0 )
Eq_n0 = Eq_n0_1;
for n=1:length(Eq_n0)
    if real(Eq_n0(n)) < 0
        Eq_n0(n) = Eq_n0_2(n);
    end
end

%%%% Trouble-shooting 
figure(101);hold on;
plot(v/1e12,real(Eq_r0),'r');
plot(v/1e12,imag(Eq_r0),'b');
xlim(PlotxRange);

figure(102);hold on;
plot(t,T_sub,'r');
plot(t,-T_ref,'k');

%%% Solve for the photoexcited reflecion coefficient
r_exp = Eq_r0 .* (dEE+1);

%%% Set initial guess for photoexcited n and k to the equilibrium values            
n1_guess_pumped = real(Eq_n0);
k1_guess_pumped = imag(Eq_n0);

%%% Time to solve for n and k 
start_at_freq = min_freq;   % Frequency to start at

start_at_point = round(start_at_freq/dv);
if start_at_point == 0
   start_at_point = 1;
end

go_to_freq = max_freq;              %solve for n and k out to this freq
go_to_point = round(go_to_freq/dv); %solve for n and k out to this point

%%% Extract index by fitting to model
if ThinFilm == 1
    n_fit = excited_thin_layer_fit(start_at_point,go_to_point,w,d,theta,r_exp,Eq_n0,n1_guess_pumped,k1_guess_pumped).';
end

if Gradient == 1
    n_fit = excited_ref_fit_JPR_multilayer(start_at_point,go_to_point,w,d,theta,r_exp,Eq_n0,n1_guess_pumped,k1_guess_pumped)';
end

nr = real(n_fit); % The index of refraction
nk = imag(n_fit); % The imaginary index - for some reason I get the right value ut it is negative

vs = v(start_at_point:go_to_point); %The frequency range of extraction
ws = vs'*2*pi;


tst         = nr+1i*nk;                         % Complex index 
tstt        = tst.*tst;                         % Epsilon
sigma       = -1i*e0*ws.*(tstt - e_inf)/100;   % Added /100 to get units of Ohm^-1 cm^-1
epsilonEq   = Eq_n0.*Eq_n0;                     % Equilibrium Epsilon

temp        = interp1(v,epsilonEq,vs,'spline'); % Make epsilonEq have the same dimensions as the epsilon vector
epsilonEq   = temp.';



sigmaEq     = -1i*e0*ws.*(epsilonEq - e_inf)/100;  % Equilibrium conductivity
epsilon     = tstt;
R           = abs((1 - tst)./(1+tst)).^2;
Req         = abs((1-sqrt(epsilonEq))./(1+sqrt(epsilonEq))).^2;

%--------------------------------------------------------------------------------


%-- SAVE PARAMETERS TO FILE -----------------------------------------------
[~,Root,ext] = fileparts(DeltaEfile); 
[~,RootEq,ext] = fileparts(Eq_file_name);

%%% Save extracted epsilon to file
if SaveParameters == 1;
    
    epsReal = real(epsilon);
    epsImag = imag(epsilon);
    eps = [vs;epsReal';epsImag'];
    eps_file_name = strcat(Root,'_eps',ext);
    fileID = fopen(eps_file_name,'w');
    fprintf(fileID,'%12s %12s %12s \r\n','% Freq [Hz]','Eps(real)','Eps(imag)');
    fprintf(fileID,'%12.4e %12.4e %12.4e\r\n',eps);
    fclose(fileID);
    
    epsEqReal = real(epsilonEq);
    epsEqImag = imag(epsilonEq);
    epsEq     = [vs;epsEqReal';epsEqImag'];
    epsEq_file_name = strcat(RootEq,'_eps',ext);
    fileID = fopen(epsEq_file_name,'w');
    fprintf(fileID,'%12s %12s %12s \r\n','% Freq [Hz]','Eps(real)','Eps(imag)');
    fprintf(fileID,'%12.4e %12.4e %12.4e\r\n',epsEq);
    fclose(fileID);
          
end
%------------------------------------------------------------------------------


%-- PLOT EVERYTHING -------------------------------------------------------
%%% Plot all time domain scans to compare temporal shifts
figure(OptParamFig+1)
hold on
plot(t,T_ref,'k','LineStyle','-');
plot(t,T_sub,'b','LineStyle','-');
plot(t,T_sc,'color',color,'LineStyle','-','Linewidth',2);
xlabel('Time (ps)')
ylabel('Amplitude')
title('TD THz traces')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot optical parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(OptParamFig)
hold on;
subplot(2,2,1);hold on;
plot(vs/1e12,R,'color',color);
plot(vs/1e12,Req,'k');hold on;
ylabel('R')
xlabel('Frequency (THz)');
xlim(PlotxRange);
ylim([0 1.2]);

%Loss function
figure(OptParamFig)
hold on;
subplot(2,2,2);
plot(vs/1e12,-imag(1./epsilon),'color',color);hold on;
plot(vs/1e12,-imag(1./epsilonEq),'k');
ylabel('LF -Im(1/\epsilon)');
xlabel('Frequency (THz)');
xlim(PlotxRange);


%Conductivity
figure(OptParamFig)
hold on;
subplot(2,2,3);
plot(vs/1e12,real(sigma),'color',color);hold on;
plot(vs/1e12,real(sigmaEq),'k');hold on;
ylabel('\sigma_1 (\Omega^{-1}cm^{-1})');
xlabel('Frequency (Thz)');
xlim(PlotxRange);
figure(OptParamFig)
hold on

figure(OptParamFig)
hold on;
subplot(2,2,4);
plot(vs/1e12,imag(sigma),'color',color);hold on;
plot(vs/1e12,imag(sigmaEq),'k');hold on;
ylabel('\sigma_2 (\Omega^{-1}cm^{-1})');
xlabel('Frequency (THz)');
xlim(PlotxRange);
figure(OptParamFig)
hold on

% Epsilon
% subplot(3,2,5);
% plot(vs/1e12,real(epsilon),'color',color);hold on;
% plot(vs/1e12,real(epsilonEq),'k');hold on;
% ylabel('\epsilon_1');
% xlabel('Frequency (THz)');
% xlim(PlotxRange);
% 
% subplot(3,2,6);
% plot(vs/1e12,imag(epsilon),'color',color);hold on;
% plot(vs/1e12,imag(epsilonEq),'k');hold on;
% ylabel('\epsilon_2');
% xlabel('Frequency (THz)');
% xlim(PlotxRange);

%----------------------------------------------------------------------


