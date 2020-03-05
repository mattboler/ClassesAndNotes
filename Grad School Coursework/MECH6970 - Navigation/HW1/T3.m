clc; clear all; close all;

f = 100*10^6;
df1 = -33.1679; 
df2 = -33.1711;
df3 = -33.1743;

f1 = f + df1;
f2 = f + df2;
f3 = f + df3;

vs = 3*10^8;

dr1 = vs * (1 - f1/f);
dr2 = vs * (1 - f2/f);
dr3 = vs * (1 - f3/f);
