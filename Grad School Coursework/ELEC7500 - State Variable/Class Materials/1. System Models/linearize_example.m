% Example of extracting state equation from Simulink model
%  John Hung, Nov 20, 2013
%
% Simulink model is named 'test_model'
%  The model includes two nonlinear blocks (cubic)
%
% Example shows how to get a different model by re-ordering
%  the variables.

clear all  % clear all variables
clc    % clear command window also

% Get default model from 'linearize' command
g1 = linearize('test_model')

% Get names of the state variables in model g1
g1_state=g1.statename
pause

% Change order of the state variables
% Use parentheses here, not brackets!
order = {'test_model/variable 3';
         'test_model/variable 2';
         'test_model/variable 1'}
pause

% New model with revised ordering of state variables
% Use 'StateOrder' option to specify the order
g2 = linearize('test_model','StateOrder',order)

% Get names of the state variables in model g2
g2_state=g1.statename
pause

% Extract the A,B,C,D matrices for each model
clc

A1 = g1.a
A2 = g2.a

B1 = g1.b
B2 = g2.b

C1 = g1.c
C2 = g2.c

D1 = g1.d
D2 = g2.d