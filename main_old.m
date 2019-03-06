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
% INPUT TURBULENCE- AND AIRCRAFT PARAMETERS
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
dt = 1;
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

A11 = xu;
A12 = xa ;
A13 = xt ;
A15 = xug ;
A16 = xag ;
A21 = zu;
A22 = za;
A23 = zt;
A24 = zq;
A25 = zug-zfug*V/Lg*(c/V);
A26 = mag;
A27 = mfag*(c/V);
A34 = V/c;
A41 = mu ;
A42 = ma ;
A43 = mt ;
A44 = mq ;
A45 = mug-mfug*V/Lg*(c/V);
A46 = mag;
A47 = mfag*(c/V);
A55 = -V/Lg;
A67 = 1;
A76 = -(V/Lg)^2;
A77 = -2*V/Lg;
B11 = xd;
B21 = zd;
B22 = zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg);
B23 = zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
B41 = md;
B42 = mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg);
B43 = mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
B52 = sigmaug_V*sqrt(2*V/Lg);
B63 = sigmaag*sqrt(3*V/Lg);
B73 = (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3);
A81 = A34*A41-A21*A11-A22*A21-A24*A41;
A82 = A34*A42-A21*A12-A22*A22-A24*A42;
A83 = A34*A43-A21*A13-A22*A23-A24*A43;
A84 = A34*A44-A22*A24-A23*A34-A24*A44-B21*Kt*A34;
A85 = A34*A45-A21*A15-A22*A25-A24*A45-A25*A55;
A86 = A34*A46-A21*A16-A22*A26-A24*A46-A27*A76;
A87 = A34*A47-A24*A47-A26*A67-A27*A77;
A88 = A34*A47-A24*A47-A26*A67-A27*A77;

B81 = B41*A34-A21*B11-A22*B21-A24*B41;
B82 = B42*A34-A22*B22-A24*B42-A25*B52-B22/dt;
B83 = B43*A34-A22*B23-A24*B43-A26*B63-A27*B73-B23/dt;

% STATE- AND INPUT MATRICES
A=[xu xa xt 0    xug                  xag       0 0;
   zu za zt zq   zug-zfug*V/Lg*(c/V)  zag       zfag*(c/V) 0;
   0  0  0  V/c  0                    0         0 0;
   mu ma mt mq   mug-mfug*V/Lg*(c/V)  mag       mfag*(c/V) 0;
   0  0  0  0   -V/Lg                 0         0 0;
   0  0  0  0    0                    0         1 0;
   0  0  0  0    0                   -(V/Lg)^2 -2*V/Lg 0;
   A81 A82 A83 A84 A85 A86 A87 0];

B=...
 [xd 0                                 0;
  zd zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  0                                 0;
  md mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  sigmaug_V*sqrt(2*V/Lg)            0;
  0  0                                 sigmaag*sqrt(3*V/Lg);
  0  0                                 (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3);
  B81 B82 B83];

%% Determine gain for desired Phugoid performance
C = eye (8) ; D = zeros (8 , 3 ) ;
Ka = 0;
Kq = 0;
Kt = -0.117;
K = [0 Ka Kt Kq 0 0 0 0];
At = A-B(:,1)*K;         % new A matrix = (A - BK) because of feedback


sys_no_pd = ss(A,B,C,D);
sys_pd = ss(At,B,C,D);