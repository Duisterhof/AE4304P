% Filename:  cit2s.m
%https://www.linkedin.com/in/janneke-blok-183416140?miniProfileUrn=urn%3Ali%3Afs_miniProfile%3AACoAACI7qnMBrb-HCX-umMUyA2F5WK10_OY8gx4
% Calculation of state matrix and input matrix for calculation
% of symmetric aircraft response to atmospheric turbulence.
% The system model is in the form
%
%       .
%       x = Ax + Bu
%       -    -    -
%
% with 
%     x = [u/V alpha theta qc/V u_g/V alpha_g alpha_g*]' 
% and 
%     u = [delta_e w_1 w_3]'.
%
% The turbulence filters are according to Dryden.

%
%  Cessna Citation Ce-500, landing (1)
%
clear

%% Get user input preference
disp('Decide what to plot:');
disp('1): Time domain signals');
disp('2): Spectral analysis ');
disp('3): Variances' );
n= input('Your choice:');

%% INPUT TURBULENCE- AND AIRCRAFT PARAMETERS
% AIRCRAFT FLIGHT CONDITION 'LANDING'.
V     = 51.4;
m     = 4556;
twmuc = 2*76;
KY2   = 0.980;
c     = 2.022;
S     = 24.2;
lh    = 5.5; %horizontal tail length
g = 9.80665;
Kt = -0.117;

% TURBULENCE PARAMETERS
sigma = 1;
Lg    = 1500;

sigmaug_V = sigma/V;
sigmaag   = sigma/V;

% AIRCRAFT SYMMETRIC AERODYNAMIC DERIVATIVES : 
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

B=...
 [xd 0                                 0;
  zd zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  0                                 0;
  md mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  sigmaug_V*sqrt(2*V/Lg)            0;
  0  0                                 sigmaag*sqrt(3*V/Lg);
  0  0                                 (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3)];

%% Determine gain for desired Phugoid performance
Ka = 0;
Kq = 0;
Kt = -0.117;
K = [0 Ka Kt Kq 0 0 0 ];
At = A-B(:,1)*K;         % new A matrix = (A - BK) because of feedback
C = [eye(7) ;
   -(V/g)*zu -(V/g)*za -(V/g)*zt (V/c)-(V/g)*zq   -(V/g)*(zug-zfug*V/Lg*(c/V))  -(V/g)*zag       -(V/g)*zfag*(c/V) 
    ];


D = zeros (8 , 3 ) ;
sys_no_pd = ss(A,B,C,D);
sys_pd = ss(At,B,C,D);

tfs = tf(sys_no_pd);
tf_theta = tf(tfs.Numerator{8,1},tfs.Denominator{8,1});
%sys = linmod('pitch_damper');
%sys = ss(sys.a,sys.b,sys.c,sys.d);

switch n
    case 1
        %% Time-domain analysis

        % TIME AXIS INPUT VECTOR DEFINITION
        dt = 0.01;
        T  = 3000; 
        t = [0:dt:T];
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
        y = lsim(At,B,C,D,u,t);

        a_z = a_z(V,c,g,zu,za,zt,zq,zug,zfug,Lg,zag,zfag,y(:,1),y(:,2),y(:,3),y(:,4),y(:,5),y(:,6),y(:,7),zd,sigmaug_V,sigmaag,u(:,1),u(:,2),u(:,3));
        a_z_ref = diff(y(:,3)-y(:,2))*V/(g*dt);

        % PLOTTING RESULTS
        subplot(2,3,1);
        plot(t,y(:,1)*180/pi)
        xlabel('time [s]'); ylabel('u/V [-]'); title('airspeed deviation');

        subplot(2,3,2);
        plot(t,y(:,2)*180/pi)
        xlabel('time [s]'); ylabel('alpha [deg]'); title('angle of attack');

        subplot(2,3,3);
        plot(t,y(:,3)*180/pi)
        xlabel('time [s]'); ylabel('theta [deg]'); title('pitch angle');

        subplot(2,3,4);
        plot(t,y(:,4)*180/pi)
        xlabel('time [s]'); ylabel('qc/V [deg]'); title('pitch rate');

        subplot(2,3,5);
        plot(t,y(:,8))
        xlabel('time [s]'); ylabel('load_factor'); title('Load factor from state space');

        subplot(2,3,6);
        N = length(t);
        plot( t(:,1:N-1),a_z(1:N-1))
        xlabel('time [s]'); ylabel('load_factor'); title('Load factor from time signals');
    case 2
        %% Spectral analysis

        %% Analytical method
        % GET INPUT PARAMETERS
        Wc = sigma;
        Nf = N;

        % DEFINE NOISE INTENSITY
        W  = Wc/dt;    % discrete time covariance, remember?

        % DEFINE TIME AXIS
        t  = [0:dt:T]; N = length(t);

        % Input vector u=[0 0 w3]'.
        D = zeros(5,3);
        C = [eye(4) zeros(4,3);
           -(V/g)*zu -(V/g)*za -(V/g)*zt (V/c)-(V/g)*zq   -(V/g)*(zug-zfug*V/Lg*(c/V))  -(V/g)*zag       -(V/g)*zfag*(c/V) ]; 

        % DEFINE FREQUENCY AXIS
        omega = logspace(-2,2,Nf);

        % COMPUTE FREQUENCY RESPONSE
        mag = bode(A,B,C,D,3,omega);

        % COMPUTE POWER SPECTRA
        Suu = mag(:,1).^2;
        Saa = mag(:,2).^2;
        Stt = mag(:,3).^2;
        Sqq = mag(:,4).^2;
        Snn = mag(:,5).^2;

        % PLOT POWER SPECTRA
        subplot(2,3,1); 
        loglog(omega,Suu); xlabel('omega [rad/sec]'); ylabel('Suu [rad^2]');
        subplot(2,3,2); 
        loglog(omega,Saa); xlabel('omega [rad/sec]'); ylabel('Saa [rad^2]');
        subplot(2,3,3); 
        loglog(omega,Stt); xlabel('omega [rad/sec]'); ylabel('Stt [rad^2]');
        subplot(2,3,4); 
        loglog(omega,Sqq); xlabel('omega [rad/sec]'); ylabel('Sqq [rad^2]');
        subplot(2,3,5); 
        loglog(omega,Snn); xlabel('omega [rad/sec]'); ylabel('Snn [rad^2]');


        %% experimental method, fft.m

        X = [y(:,1),y(:,2),y(:,3),y(:,4),y(:,8)];
        mag = fft(X);

        %double-sided
        Suu = abs(mag(1:N/2+1,1)/N);
        Saa = abs(mag(1:N/2+1,2)/N);
        Stt = abs(mag(1:N/2+1,3)/N);
        Sqq = abs(mag(1:N/2+1,4)/N);
        Snn = abs(mag(1:N/2+1,5)/N);

        %convert to single sided
        Suu(2:end-1) = 2*Suu(2:end-1);
        Saa(2:end-1) = 2*Saa(2:end-1);
        Stt(2:end-1) = 2*Stt(2:end-1);
        Sqq(2:end-1) = 2*Sqq(2:end-1);
        Snn(2:end-1) = 2*Snn(2:end-1);

        Fs = N/T;
        f = (Fs/N)*(0:(N/2));

        % PLOT POWER SPECTRA
        subplot(2,3,1); 
        loglog(f,Suu); xlabel('omega [rad/sec]'); ylabel('Suu [rad^2]');
        subplot(2,3,2); 
        loglog(f,Saa); xlabel('omega [rad/sec]'); ylabel('Saa [rad^2]');
        subplot(2,3,3); 
        loglog(f,Stt); xlabel('omega [rad/sec]'); ylabel('Stt [rad^2]');
        subplot(2,3,4); 
        loglog(f,Sqq); xlabel('omega [rad/sec]'); ylabel('Sqq [rad^2]');
        subplot(2,3,5); 
        loglog(f,Snn); xlabel('omega [rad/sec]'); ylabel('Snn [rad^2]');

        %% Experimental PSD, pwelch.m
        X = [y(:,1),y(:,2),y(:,3),y(:,4),y(:,8)];
        [pxx,f] = pwelch(X,[],[],f,Fs);
        Suu = pxx(:,1);
        Saa = pxx(:,2);
        Stt = pxx(:,3);
        Sqq = pxx(:,4);
        Snn = pxx(:,5);

        % PLOT POWER SPECTRA
        subplot(2,3,1); 
        loglog(f,Suu); xlabel('omega [rad/sec]'); ylabel('Suu [rad^2]');
        subplot(2,3,2); 
        loglog(f,Saa); xlabel('omega [rad/sec]'); ylabel('Saa [rad^2]');
        subplot(2,3,3); 
        loglog(f,Stt); xlabel('omega [rad/sec]'); ylabel('Stt [rad^2]');
        subplot(2,3,4); 
        loglog(f,Sqq); xlabel('omega [rad/sec]'); ylabel('Sqq [rad^2]');
        subplot(2,3,5); 
        loglog(f,Snn); xlabel('omega [rad/sec]'); ylabel('Snn [rad^2]');
    case 3
        %% Variances

        %analytical variance
        % NUMERICAL INTEGRATION OF PSD's
        do = diff(f)';   % compute "difference vector" in omega
                             % i.e., omega(k+1)-omega(k);
        Ndo = length(do)
        % then perform (very crude) integration
        var(1) = sum(do.*Suu(1:Ndo));
        var(2) = sum(do.*Saa(1:Ndo));
        var(3) = sum(do.*Stt(1:Ndo));
        var(4) = sum(do.*Sqq(1:Ndo));
    otherwise 
        disp('Invalid input');
end



