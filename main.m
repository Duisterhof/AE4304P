%%Bart Duisterhof -- 4442695

%% Clear workspace 
clear

%% Retrieve user input preference
disp('Decide what to plot:');
disp('0): Root Locus plot');
disp('1): Time domain signals');
disp('2): Spectral analysis ');
disp('3): Variances' );
n= input('Your choice: ');
if n~=0
    k = input('Use pitch damper? Yes = 1, no = 0 ');
end

%% INPUT TURBULENCE- AND AIRCRAFT PARAMETERS
% AIRCRAFT FLIGHT CONDITION 'LANDING'.
V     = 51.4;
m     = 4556;
twmuc = 2*76;
KY2   = 0.980;
c     = 2.022;
S     = 24.2;
lh    = 5.5;   %horizontal tail length
g = 9.80665;
Kt = -0.117;

% TURBULENCE PARAMETERS
sigma = 1;
Lg    = 1500;

sigmaug_V = sigma/V;
sigmaag   = sigma/V;

%%%AIRCRAFT SYMMETRIC AERODYNAMIC DERIVATIVES : 
CX0 = 0.0000;     CZ0  =-1.1360;     Cm0  =  0.0000;
CXu =-0.2173;     CZu  =-2.2720;     Cmu  =  0.0000;
CXa = 0.4692;     CZa  =-5.13;       Cma  = -0.400;
CXq = 0.0000;     CZq  =-3.8400;     Cmq  = -7.3500;
CXd = 0.0000;     CZd  =-0.6238;     Cmd  = -1.5530;
CXfa= 0.0000;     CZfa =-1.4050;     Cmfa = -3.615;
                  CZfug= 0.0000;     Cmfug= -Cm0*lh/c;
                  CZfag= CZfa-CZq;   Cmfag=  Cmfa-Cmq;

% CALCULATION OF AIRCRAFT SYMMETRIC STABILITY DERIVATIVES
xu   = (V/c)*(CXu/twmuc);
xa   = (V/c)*(CXa/twmuc);
xt   = (V/c)*(CZ0/twmuc);
xq   = 0;
xd   = (V/c)*(CXd/twmuc);
xug  = xu;
xfug = 0;
xag  = xa;
xfag = 0;

zu   = (V/c)*( CZu/(twmuc-CZfa));
za   = (V/c)*( CZa/(twmuc-CZfa));
zt   = (V/c)*(-CX0/(twmuc-CZfa));
zq   = (V/c)*((CZq+twmuc)/(twmuc-CZfa));
zd   = (V/c)*( CZd/(twmuc-CZfa));
zug  = zu;
zfug = (V/c)*( CZfug/(twmuc-CZfa));
zag  = za;
zfag = (V/c)*( CZfag/(twmuc-CZfa));

mu   = (V/c)*(( Cmu+CZu*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
ma   = (V/c)*(( Cma+CZa*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mt   = (V/c)*((-CX0*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mq   = (V/c)*(Cmq+Cmfa*(twmuc+CZq)/(twmuc-CZfa))/(twmuc*KY2);
md   = (V/c)*((Cmd+CZd*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mug  = mu;
mfug = (V/c)*(Cmfug+CZfug*Cmfa/(twmuc-CZfa))/(twmuc*KY2);
mag  = ma;
mfag = (V/c)*(Cmfag+CZfag*Cmfa/(twmuc-CZfa))/(twmuc*KY2);

% STATE- AND INPUT MATRICES
A=[xu xa xt 0    xug                  xag       0 ;
   zu za zt zq   zug-zfug*V/Lg*(c/V)  zag       zfag*(c/V) ;
   0  0  0  V/c  0                    0         0 ;
   mu ma mt mq   mug-mfug*V/Lg*(c/V)  mag       mfag*(c/V) ;
   0  0  0  0   -V/Lg                 0         0 ;
   0  0  0  0    0                    0         1 ;
   0  0  0  0    0                   -(V/Lg)^2 -2*V/Lg ];

B= [xd 0                                 0;
  zd zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  0                                 0;
  md mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  sigmaug_V*sqrt(2*V/Lg)            0;
  0  0                                 sigmaag*sqrt(3*V/Lg);
  0  0                                 (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3)];

%define C and D matrix
%C-matrix with 4 standard aircraft states and load-factor
n_z_row = (V/g)*(A(3,:)-A(2,:));
C = [eye(4) zeros(4,3);
n_z_row]; 
D = [zeros(4,3);
    (V/g)*(B(3,:)-B(2,:))];

if n==0
    %TF's for gain determination
    sys_no_pd = ss(A,B,C,D);
    tfs = tf(sys_no_pd);
    tf_theta = tf(tfs.Numerator{3,1},tfs.Denominator{3,1});
    
    %compute RL plot
    [R,K_RL] = rlocus(tf_theta,-2:0.001:0);
    plot(real(R), imag(R),'.','color',[10/255,10/255,10/255],'markersize',1)
    title('Root Locus','fontsize',12);
    xlabel('Real Axis (seconds^{-1})')
    ylabel('Imaginary Axis (seconds^{-1})')
    return
end
%% Input gains and matrices for pitch damper
%Gains obtained through root-locus method
Ka = 0;
Kq = 0;
K = [0 Ka Kt Kq 0 0 0 ];

%if pitchdamer is enabled, modify A-matrix
if k==1
    A = A-B(:,1)*K;         % new A matrix = (A - BK) because of feedback
end

%computing system for pzmap with relevant eigenmodes
C_red = eye(7);
D_red = zeros(7,3);
sys = ss(A,B,C_red,D_red);

switch n
    case 1
        %% Time-domain analysis

        % TIME AXIS INPUT VECTOR DEFINITION
        rng(1)
        dt = 0.01;
        T  = 60; 
        t = 0:dt:T;
        N = length(t);

        % INPUT VECTOR DEFINITION
        nn = zeros(1,N);             % zero input elevator
        w1 = randn(1,N)/sqrt(dt);    % scaled input hor. turbulence,
                                     % note the sqrt(dt) because of lsim
        w3 = randn(1,N)/sqrt(dt);    % scaled input vert. turbulence,
                                     % note the sqrt(dt) because of lsim
        u  = [nn' nn' w3'];          % input vector definition (vertical
                                     % turbulence only, can be changed).

        % SIMULATION OF MOTION VARIABLES
        y = lsim(A,B,C,D,u,t);
        a_z_ref = diff(y(:,3)-y(:,2))*V/(g*dt);

        % PLOTTING RESULTS
        subplot(6,1,1);
        plot(t,y(:,1))
        xlabel('time [s]'); ylabel('$\frac{u}{V}$ [-]','interpreter','latex'); title('airspeed deviation');

        subplot(6,1,2);
        plot(t,y(:,2)*180/pi)
        xlabel('time [s]'); ylabel('$\alpha$ [deg]','interpreter','latex'); title('angle of attack');

        subplot(6,1,3);
        plot(t,y(:,3)*180/pi)
        xlabel('time [s]'); ylabel('$\theta$ [deg]','interpreter','latex'); title('pitch angle');

        subplot(6,1,4);
        plot(t,y(:,4)*180/pi)
        xlabel('time [s]'); ylabel('$\frac{qc}{V}$ [deg]','interpreter','latex'); title('pitch rate');

        subplot(6,1,5);
        plot(t,y(:,5))
        xlabel('time [s]'); ylabel('$n_z$ [-]','interpreter','latex'); title('Load factor from state space');

        subplot(6,1,6);
        N = length(t);
        plot( t(:,1:N-1),a_z_ref(1:N-1))
        xlabel('time [s]'); ylabel('$n_z$ [-]','interpreter','latex'); title('Load factor from time signals');
    case 2
        %% Spectral analysis
        %% Time-domain signals to be used

        % TIME AXIS INPUT VECTOR DEFINITION
        dt = 0.05;
        T  = 200; 
        t = [0:dt:T];
        N = length(t);
        Fs = 1/dt;

        % INPUT VECTOR DEFINITION
        nn = zeros(1,N);             % zero input elevator
        w1 = randn(1,N)/sqrt(dt);    % scaled input hor. turbulence,
                                     % note the sqrt(dt) because of lsim
        w3 = randn(1,N)/sqrt(dt);    % scaled input vert. turbulence,
                                     % note the sqrt(dt) because of lsim
        u  = [nn' nn' w3'];          % input vector definition (vertical
                                     % turbulence only, can be changed).

        % SIMULATION OF MOTION VARIABLES
        y = lsim(A,B,C,D,u,t);

        %% Analytical method
        % GET INPUT PARAMETERS
        Wc = sigma;
        Nf = N;

        % DEFINE NOISE INTENSITY
        W  = Wc/dt;    % discrete time covariance, remember?
      
        % DEFINE FREQUENCY AXIS
        f_1 = logspace(-2,2,Nf);

        % COMPUTE FREQUENCY RESPONSE
        mag = bode(A,B,C,D,3,f_1);

        % COMPUTE POWER SPECTRA
        Suu_1 = mag(:,1).^2;
        Saa_1 = mag(:,2).^2;
        Stt_1 = mag(:,3).^2;
        Sqq_1 = mag(:,4).^2;
        Snn_1 = mag(:,5).^2;

        % PLOT POWER SPECTRA
        figure(1)
        subplot(2,3,1); 
        loglog(f_1,Suu_1); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Suu [$\frac{1}{rad/s}$]','interpreter','latex');
        subplot(2,3,2); 
        loglog(f_1,Saa_1); xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('Saa [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,3); 
        loglog(f_1,Stt_1); xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('Stt [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,4); 
        loglog(f_1,Sqq_1);xlabel('$\omega$ [rad/sec]','interpreter','latex');  ylabel('Sqq [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,5); 
        loglog(f_1,Snn_1); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Snn [$\frac{1}{rad/s}$]','interpreter','latex');

        %% experimental method, fft.m

        X = [y(:,1),y(:,2),y(:,3),y(:,4),y(:,5)];
        mag = dt*fft(X);
        
        Suu = mag(1:N/2+1,1);
        Saa = mag(1:N/2+1,2);
        Stt = mag(1:N/2+1,3);
        Sqq = mag(1:N/2+1,4);
        Snn = mag(1:N/2+1,5);

        %convert to PSD
        %this is different to the equation in the report, due to matlab
        %defenitions
        Suu_2 = (1/T)*abs(Suu).^2;
        Saa_2 = (1/T)*abs(Saa).^2;
        Stt_2 = (1/T)*abs(Stt).^2;
        Sqq_2 = (1/T)*abs(Sqq).^2;
        Snn_2 = (1/T)*abs(Snn).^2;
        
        f_2 = 2*pi*Fs*(0:(N/2))/N;          %frequency array in rad/s
        N = length(f_2);

        % PLOT POWER SPECTRA
        figure(2)
        subplot(2,3,1); 
        loglog(f_2,Suu_2); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Suu [$\frac{1}{rad/s}$]','interpreter','latex');
        subplot(2,3,2); 
        loglog(f_2,Saa_2);xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Saa [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,3); 
        loglog(f_2,Stt_2);xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Stt [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,4); 
        loglog(f_2,Sqq_2); xlabel('$\omega$ [rad/sec]','interpreter','latex');ylabel('Sqq [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,5); 
        loglog(f_2,Snn_2); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Snn [$\frac{1}{rad/s}$]','interpreter','latex');
        
 
        %% Experimental PSD, pwelch.m
        X = [y(:,1),y(:,2),y(:,3),y(:,4),y(:,5)];
        f = (Fs/N)*(0:(N/2)); %frequency input array in Hz
        [pxx,f] = pwelch(X,[],[],f,Fs);
        
        %converting everything to rad
        f_3 = f*2*pi;
        Suu_3 = pxx(:,1);
        Saa_3 = pxx(:,2);
        Stt_3 = pxx(:,3);
        Sqq_3 = pxx(:,4);
        Snn_3 = pxx(:,5);
                
        % PLOT POWER SPECTRA
        figure(3)
        subplot(2,3,1); 
        loglog(f_3,Suu_3); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Suu [$\frac{1}{rad/s}$]','interpreter','latex');
        subplot(2,3,2); 
        loglog(f_3,Saa_3);xlabel('$\omega$ [rad/sec]','interpreter','latex');ylabel('Saa [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,3); 
        loglog(f_3,Stt_3); xlabel('$\omega$ [rad/sec]','interpreter','latex');ylabel('Stt [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,4); 
        loglog(f_3,Sqq_3); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Sqq [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,5); 
        loglog(f_3,Snn_3); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Snn [$\frac{1}{rad/s}$]','interpreter','latex');
        
        
        figure(4)

        subplot(2,3,1); 
        loglog(f_1,Suu_1,f_2,Suu_2,f_3,Suu_3); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Suu [$\frac{1}{rad/s}$]','interpreter','latex');legend('analytical','fft.m','pwelch.m');
        subplot(2,3,2); 
        loglog(f_1,Saa_1,f_2,Saa_2,f_3,Saa_3); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Saa [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,3); 
        loglog(f_1,Stt_1,f_2,Stt_2,f_3,Stt_3); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Stt [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,4); 
        loglog(f_1,Sqq_1,f_2,Sqq_2,f_3,Sqq_3); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Sqq [$\frac{rad^2}{rad/s}$]','interpreter','latex');
        subplot(2,3,5); 
        loglog(f_1,Snn_1,f_2,Snn_2,f_3,Snn_3); xlabel('$\omega$ [rad/sec]','interpreter','latex'); ylabel('Snn [$\frac{1}{rad/s}$]','interpreter','latex');

    case 3
        %% Variances
        %compute time simulation

        % TIME AXIS INPUT VECTOR DEFINITION
        dt = 0.01;
        T  = 30000; 
        t = [0:dt:T];
        N = length(t);
        Fs = N/T;

        % INPUT VECTOR DEFINITION
        nn = zeros(1,N);             % zero input elevator
        w1 = randn(1,N)/sqrt(dt);    % scaled input hor. turbulence,
                                     % note the sqrt(dt) because of lsim
        w3 = randn(1,N)/sqrt(dt);    % scaled input vert. turbulence,
                                     % note the sqrt(dt) because of lsim
        u  = [nn' nn' w3'];          % input vector definition (vertical
                                     % turbulence only, can be changed).

        % SIMULATION OF MOTION VARIABLES
        y = lsim(A,B,C,D,u,t);
            
        %% analytical variance
        
        Bin = B(:,3);
        L   = lyap(A,Bin*sigma*Bin');
        % take only the part that belongs to the 4 states
        L   = L(1:4,1:4);   
        var_arr = [L(1,1),L(2,2),L(3,3),L(4,4)];
        disp('Using lyapunov');
        disp('1.0e-06 *')
        disp(10^6*var_arr)
        
        
        %% NUMERICAL INTEGRATION OF PSD's
        
        %compute analytical PSD
     
        % DEFINE FREQUENCY AXIS
        omega = logspace(-2,2,N);
        
        %% bod() routine 
        % COMPUTE FREQUENCY RESPONSE
        mag = bode(A,B,C,D,3,omega);

        % COMPUTE POWER SPECTRA
        Suu = mag(:,1).^2;
        Saa = mag(:,2).^2;
        Stt = mag(:,3).^2;
        Sqq = mag(:,4).^2;
        Snn = mag(:,5).^2;
   
        % NUMERICAL INTEGRATION OF PSD's
        do = diff(omega)';   % compute "difference vector" in omega
                             % i.e., omega(k+1)-omega(k);
        % then perform (very crude) integration
        var_arr(1) = sum(do.*Suu(1:N-1));
        var_arr(2) = sum(do.*Saa(1:N-1));
        var_arr(3) = sum(do.*Stt(1:N-1));
        var_arr(4) = sum(do.*Sqq(1:N-1));    
        var_arr(5) = sum(do.*Snn(1:N-1));
        disp('Integrating analytical PSDs');
        disp('1.0e-06 *')
        disp(10^6*var_arr/pi)
               

        
        %% using fft.m
              
        % generate plots with fft
        
        X = [y(:,1),y(:,2),y(:,3),y(:,4),y(:,5)];
        mag = dt*fft(X);
        
        %cut axis in half, as it's a fourier transform
        Suu = mag(1:N/2+1,1);
        Saa = mag(1:N/2+1,2);
        Stt = mag(1:N/2+1,3);
        Sqq = mag(1:N/2+1,4);
        Snn = mag(1:N/2+1,5);

        %square and scale
        Suu = (1/T)*abs(Suu).^2;
        Saa = (1/T)*abs(Saa).^2;
        Stt = (1/T)*abs(Stt).^2;
        Sqq = (1/T)*abs(Sqq).^2;
        Snn = (1/T)*abs(Snn).^2;
        
        f = 2*pi*Fs*(0:(N/2))/N;          %frequency array in rad/s
        N = length(f);
        
        %integrate
        do = diff(f)';
        var_arr = zeros(1,5);
        var_arr(1) = sum(do.*Suu(1:N-1));
        var_arr(2) = sum(do.*Saa(1:N-1));
        var_arr(3) = sum(do.*Stt(1:N-1));
        var_arr(4) = sum(do.*Sqq(1:N-1));    
        var_arr(5) = sum(do.*Snn(1:N-1));
        disp('Integrating fft.m PSDs:');
        disp('1.0e-06 *')
        disp(10^6*var_arr/pi)     
        
        %% using pwelch.m
        %call function
        f = (Fs/N)*(0:(N/2)); %frequency input array in Hz
        [pxx,f] = pwelch(X,[],[],f,Fs);
        Suu = pxx(:,1);
        Saa = pxx(:,2);
        Stt = pxx(:,3);
        Sqq = pxx(:,4);
        Snn = pxx(:,5);
        
        %integration
        f = (Fs/N)*(0:(N/2))*2*pi; %frequency array in rad/s
        do = diff(f)';
        N = length(f);
        
        var_arr(1) = sum(do.*Suu(1:N-1));
        var_arr(2) = sum(do.*Saa(1:N-1));
        var_arr(3) = sum(do.*Stt(1:N-1));
        var_arr(4) = sum(do.*Sqq(1:N-1));    
        var_arr(5) = sum(do.*Snn(1:N-1));
        disp('Integrating pwelch.m PSDs:');
        disp('1.0e-06 *')
        disp(10^6*var_arr/pi)
        
        %% using var.m
        var_arr = var(y);
        disp('Using var.m routine:')
        disp('1.0e-06 *')
        disp(10^6*var_arr)
        
    otherwise 
        disp('Invalid input');
end



