%{
Define bode plot settings
%}
opts = bodeoptions;
opts.Grid = 'on';
opts.Xlim = [10e-3, 10e3];
%opts.Ylim = [-80, 80];

%{
Design a control system using frequency shaping for the following system:
%}
s = tf('s');
G = 10 / (s+1)^2;

%{
Performance requirements:
    0 steady-state error to a unit step
    -40dB attenuation in [0.01:0.1] rad/s
    -40dB attenuation in [100:1000] rad/s
    10 rad/s bandwidth
    30 deg phase margin
%}

%{
Notes on W(s):
    Performance guarantee is given by
    |W*S| < 1 for all freq
    
%}

% Designing W:
%{
zeta_z = .707;
wn_z = 0.1;
numerator_w = (s^2 + 2*zeta_z * wn_z*s + wn_z);

zeta_p = .707;
wn_p = 1000;
denominator_w = (s^2 + 2*zeta_p*wn_p*s + wn_p);
gain_w = 100;
LAMBDA = gain_w * (numerator_w) / denominator_w;
W = LAMBDA^-1;
%}
load('W_current.mat')
LAMBDA
W = LAMBDA^-1

load('K_current.mat');
K = K_des

L = G*K;
S = 1/(1 + L);
T = L / (1 + L);

T
stepinfo(T)

%{
margin(T)

figure(1);
f = bodeplot(LAMBDA, opts);
title('Bode Plot of \Lambda');

figure(2);
a = bodeplot(S, LAMBDA, opts);
legend('S', 'Lambda');
title('Bode Plot of Sensitivity');

figure(3);
h = bodeplot(L, opts);
title('Bode Plot of G(s) * K(s)');

figure(4)
g = bodeplot(S*W, opts);
title('Bode Plot of S(s) * W(s)');

figure(5)
b = bodeplot(G, opts);
title('Bode Plot of G(s)');

figure(6)
step(T);
title("Closed Loop Response");
%}