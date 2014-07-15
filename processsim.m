%close all; clear all;
% Polystyrene: Deterministic model
% this file was transformed from script file to function file for the
% purpose of integration of unscented kalman filter.
% created by santosh on 10 may 2014
% for the purpose of maintaing uniform variables same as in literature here
% we have changed the variable 'y' to 'x'.

function xresplant = processsim (xinit, tinit, sampletime);

global f R kp0 kd0 kt0 kf0 kfs0 Ep Ed Et Ef Efs deltaHr UA rhoCp rhoCpc Qi Qs Qm Qc V Vc Cif Cmf Tf Tcf Cpdot Csf dinit 

% Parameters
f = 0.6; %initiator efficiency
R = 8.314; % Universal gas constant
kp0 = 1.06*10e7; kd0 = 5.95*10e13; kt0 = 1.25*10e9; kf0 = 53020.29; kfs0 = 91457.3;
Ep = 29572.898; Ed = 123853.658; Et = 7017.27; Ef = 53029.29; Efs = 91457.3;
deltaHr = 69919.56; % Heat of reaction
UA = 293.076; 
rhoCp = 1507.248; rhoCpc = 4045.7048; 
%Qi = 0.03;%Initiator flowrate
%Qs = 0.1275; % Solvent flowrate
%Qm = 0.105; % Monomer flowrate
%Qc = 0.131;  % Cooling jacket flowrate   original value 0.131 l/s
V = 3000; %Reactor volume
Vc = 3312.4; %Cooling jacket volume
Cif = 0.5888; % Initial feed concentration
Cmf = 8.6981; % Monomer feed concentration
Csf = 13; % ????? check this value later. Also Cmf etc given above. Look suspect.
Tf = 330; % Reactor feed temperature
Tcf = 295; % Cooling jacket feed temperature

%xauginit = [xinit;wnoise1; wnoise2; wnoise3; wnoise4; wnoise5; wnoise6;vnoise];

% Initial conditions:
%xinit = [0.5888 8.6981 13  330 295 0 0 0 ]';
%tfinal = tinit + tsample;
% Call to ode45
%T = 0.1;
%N = 5;
[t, xresplant] = ode45(@odes_plant,[tinit sampletime],xinit);   % ????? plz chech the function everytime if sumthing is changed
    %Xres = xres(N/T +1,:);
%xpredict = xres;
% Result analysis
