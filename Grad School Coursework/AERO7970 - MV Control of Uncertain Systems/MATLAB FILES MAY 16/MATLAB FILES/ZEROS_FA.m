%  ZEROS AND TIME RESPONSE
%   SCRIPT FILE TO ALANYZE THE INFLUENCE
%   OF ZEROS IN THE TIME RESPONSE
s=tf('s')
%   POINT MASS MODEL NUMBER 1
g1=(s-1)/(s^2)
pause
%   STATE SPACE REPRESENTATION
%   x1: velocity
%   x2: position
%   input: mass acceleration
%   output: x1 - x2
A=[0 0;1 0];
B=[1 0]';
C=[1 -1];
%   FIND A CONTROLLER K(S) SUCH THAT POSITION ERROR 
%   TO A UNIT STEP IS EQUAL TO ZERO
%
%   STRUCTURAL PROPERTIES
%
CONT=rank(ctrb(A,B))
OBSV=rank(obsv(A,C))
pause
%THE SYSTEM IS UNSTABLE
impulse(g1),grid, pause
%   POSITION AND VELOCITY FROM OUTPUT
%   x1 = u - ydot
%   x2 = u - ydot - y
%
%   NOTE: TRANSFER FUNCTION BETWEEN U AND X2 IS 1/S^2
%
%   UNITY FEEDBACK TO STABILIZE THE SYSTEM
%   K1(s)=-0.045(1+10s)/(1+0.2s)
%
k1=-0.045*(1+10*s)/(1+0.2*s);
gcl1=k1*g1/(1+k1*g1)
%rhStabilityCriterion
%
pause
ltiview(k1*g1/(1+k1*g1)),grid
pause
%
%   STEP RESPONSE Y, X1, X2
%
y=gcl1;
x1=(k1-s*k1*g1)/(1+k1*g1);
x2=x1/s;
ltiview(y,x1,x2),grid
pause
%
%   VERIFY THAT A PI CONTROLLER SATISFIES THE REQUIREMENTS
%
k2=-0.95*(1+2*s)/s
pause
x2fin=k2*x2/(1+k2*x2)
x1fin=s*x2fin
yfin=x1fin-x2fin
ltiview(x1fin,x2fin,yfin)
pause
%   CHANGE SENSOR, MEASURE POSITION
%
pos=1/(s^2)
k=0.1*(1+10*s)/(1+0.2*s)
poscl=k*pos/(1+k*pos);
velcl=s*poscl;
ycl=velcl-poscl;
ltiview(poscl,velcl,ycl)
%



