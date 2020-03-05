




%This script is to illustrate how one can design H-infinity controllers in
%Matlab. The example shows some different ways to synthesize a controller
%for the mixed S/KS weighted sensitivity problem. The plant and the weights
%were found in (Skogestad and Postlethwaite, 1996 ed.1 p.60) and the
%weights are not necessarily the "best".
%---------------------------------
%Defining the subsystems
%---------------------------------
%Plant G=200/((10s+1)(0.05s+1)^2)
%Alternative 1, mu-tools directly:
G = nd2sys(1,conv([10,1],conv([0.05 1],[0.05 1])),200);
%Alternative 2, indirectly via cst:
%s = tf('s');
%Gcst = 200/((10*s+1)*(0.05*s+1)^2);
%[a,b,c,d] = ssdata( balreal(Gcst) );
%G = pck(a,b,c,d);
G,pause
%Weights Ws = (s/M+w0)/(s+w0*A), Wks=1
M = 1.5; w0 = 10; A=1.e-4;
Ws = nd2sys([1/M w0],[1 w0*A]);
Wks = 1;
%---------------------------------
%Creating the generalized plant P
%---------------------------------
%Alternative 0, direct approach:
% /z1\ /Ws -Ws*G\ /r\
% |z2| = |0 Wks | | |
% \ v/ \I -G / \u/
%Transfer matrix representation
Z1 = sbs(Ws,mmult(-1,Ws,G));
Z2 = sbs(0,Wks);
V = sbs(1,mmult(-1,G));
P0 = abv(Z1,Z2,V);
%P0 is generally not a minimal realization, so we use reduction methods
[a,b,c,d] = unpck(P0);
[ab,bb,cb,db] = ssdata( balreal( minreal( ss(a,b,c,d) ) ) );
P0 = pck(ab,bb,cb,db); %now we have a System description
%---------------------------------
%Creating the generalized plant P
%---------------------------------
%Alternative 1, direct approach:
% /z1\ /W1 -W1*G\ /r\
% |z2| = |0 W2 | | |
% \ v/ \I -G / \u/
%ss realization of the subsystems:
[A,B,C,D] = unpck(G);
[A1,B1,C1,D1] = unpck(Ws);
[A2,B2,C2,D2] = unpck(Wks);
%number of inputs and outputs to the different subsystems:
n1 = size(A1,1); [q1, p1] = size(D1);
n2 = size(A2,1); [q2, p2] = size(D2);
n = size(A,1) ; [p , q ] = size(D);
%ss realization of the whole thing:
Ap = [ A1 , zeros(n1,n2) , -B1*C ;
zeros(n2,n1) , A2 , zeros(n2,n) ;
zeros(n ,n1) , zeros(n ,n2) , A ];
Bp = [ B1 ,-B1*D;
zeros(n2,p) , B2 ;
zeros(n ,p) , B ];
Cp = [ C1 , zeros(q1,n2) , -D1*C ;
zeros(q2,n1), C2 , zeros(q2,n) ;
zeros(q ,n1), zeros(q ,n2) , -C ];
Dp = [ D1 , -D1*D;
zeros(q2,p ), D2 ;
eye(p) , -D ];
%making a balanced realization reduces the likelihood of numerical problems
[Apb,Bpb,Cpb,Dpb] = ssdata( balreal( ss(Ap,Bp,Cp,Dp) ) );
P1 = pck(Apb,Bpb,Cpb,Dpb);
%---------------------------------
%Creating the generalized plant P
%---------------------------------
%Alternative 2, using sysic:
systemnames = 'G Ws Wks';
inputvar = '[r(1); u(1)]'; %all inputs are scalar, r(2) would be a 2dim signal
outputvar = '[Ws; Wks; r-G]';
input_to_G = '[u]';
input_to_Ws = '[r-G]';
input_to_Wks = '[u]';
sysoutname = 'P2';
cleanupsysic = 'yes';
sysic
%---------------------------------
%Synthesizing the controller
%---------------------------------
%All the methods presented here use the System
%matrix representation of the generalized plant P
%Choose your favourite method and your favourite P
%Choose plant
P = P1; %(0-4)
%then some parameters (number of measurements and inputs,bounds on gamma )
nmeas = 1; nu = 1; gmn=0.5; gmx=20; tol = 0.001;
%uncomment your favourite controller
[K1,CL,gopt] = hinfsyn(P,nmeas,nu,gmn,gmx,tol);
%[gopt,K1] = hinflmi(P,[nmeas, nu],0,tol); CL = starp(P,K,nmeas,nu);
%[gopt,K1] = hinfric(P,[nmeas, nu],gmn,gmx); CL = starp(P,K,nmeas,nu);
K=K1;
K, pause
%Alternative for RCT v3.0.1.
%Normally you would of course not do all the transfers between the system
%representations, but rather do everything using standard ss objects
%[a,b,c,d] = unpck(G); Gcst = ss(a,b,c,d);
%[a,b,c,d] = unpck(Ws); Wscst = ss(a,b,c,d);
%[a,b,c,d] = unpck(Wks); Wkscst = ss(a,b,c,d);
%[K4,CL,gopt] = mixsyn(Gcst,Wscst,Wkscst,[]);
%[a,b,c,d] = ssdata( balreal(K) ); K = pck(a,b,c,d);
%[a,b,c,d] = ssdata( balreal(CL) ); CL = pck(a,b,c,d);
%K=K4;
%K, pause
%---------------------------------
%Analysis of the result
%---------------------------------
%plot singular values of (weighted) closed loop system
w = logspace(-4,6,50);
CLw = vsvd(frsp(CL,w));
figure(1); vplot('liv,m',CLw);
title('singular values of weighted closed loop system');
pause
%generate typical transfer matrices
[type,out,in,n] = minfo(G);
I = eye(out);
S = minv(madd(I,mmult(G,K))); %sensitivity
T = msub(I,S); %complementary sensitivity
KS = mmult(K,S); %input to G
GK = mmult(G,K); %loop transfer function
%singular values as a function of frequency
Sw = vsvd(frsp(S,w));
Tw = vsvd(frsp(T,w));
Kw = vsvd(frsp(K,w));
KSw = vsvd(frsp(KS,w));
GKw = vsvd(frsp(GK,w));
%Plot singular value plots
%Note: if desired, you can change vplot to plot the amplitude in dB. Type
%edit vplot and uncomment the appropriate lines in the code
figure(2); vplot('liv,lm',Sw,'-',Tw,'--',GKw,'-.');
title('\sigma(S(jw)) (solid), \sigma(T(jw)) (dashed) and \sigma(GK(jw)) (dashdot)');
xlabel('Frequency [rad/sec]'); ylabel('Amplitude')
figure(3); vplot('liv,lm',Kw);
title('\sigma(K(jw))');
xlabel('Frequency [rad/sec]’); ylabel(’Amplitude')
pause
%Did we get what we asked for?
Sd = minv(Ws); Sdw = vsvd(frsp(Sd,w)); %"desired" sensitivity
KSd = minv(Wks); KSdw = vsvd(frsp(KSd,w)); %"desired" output
figure(4); vplot('liv,lm',Sw,'-',Sdw,'--');
title('\sigma(S(jw)) (solid) and \sigma(Ws^{-1}(jw)) (dashed)');
xlabel('Frequency [rad/sec]'); ylabel('Amplitude')
figure(5); vplot('liv,lm',KSw,'-',KSdw,'--')
title('\sigma(KS(jw)) (solid) and \sigma(Wks^{-1}(jw)) (dashed)');
xlabel('Frequency [rad/sec]'); ylabel('Amplitude')
pause
%Finally the step response
reference = 1; tfinal = 1; step = 0.01;
y = trsp(T,reference,tfinal,step);
u = trsp(KS,reference,tfinal,step);
figure(6); subplot(2,1,1); vplot('iv,d',y);
title('Step response'); ylabel('y');
subplot(2,1,2); vplot('iv,d',u);
ylabel('u'); xlabel('time');

