% 17 May 2014 Saturday

% Kalman Filter-based state estimatioin for a polystyrene CSTR
global f R kp0 kd0 kt0 kf0 kfs0 Ep Ed Et Ef Efs deltaHr UA rhoCp rhoCpc Qi Qs Qm Qc V Vc Cif Cmf Tf Tcf Cpdot Csf dinit 

% -------------------------------------
% Initialization
% -------------------------------------

f = 0.6; %initiator efficiency
R = 8.314; % Universal gas constant
kp0 = 1.06*10e7; kd0 = 5.95*10e13; kt0 = 1.25*10e9; kf0 = 53020.29; kfs0 = 91457.3;
Ep = 29572.898; Ed = 123853.658; Et = 7017.27; Ef = 53029.29; Efs = 91457.3;
deltaHr = 69919.56; % Heat of reaction
UA = 293.076; 
rhoCp = 1507.248; rhoCpc = 4045.7048; 
Qi = 0.03;%Initiator flowrate
Qs = 0.1275; % Solvent flowrate
Qm = 0.105; % Monomer flowrate
Qc = 0.131;  % Cooling jacket flowrate   original value 0.131 l/s
V = 3000; %Reactor volume
Vc = 3312.4; %Cooling jacket volume
Cif = 0.5888; % Initial feed concentration
Cmf = 8.6981; % Monomer feed concentration
Csf = 13; % ????? check this value later. Also Cmf etc given above. Look suspect.
Tf = 330; % Reactor feed temperature
Tcf = 295; % Cooling jacket feed temperature



dinit=zeros(6,53);

% Parameters of process model

C=[0 0 1 0 0 0 0 0;
   0 0 0 1 0 0 0 0;
   0 0 0 0 1 0 0 0;
   0 0 0 0 0 1 0 0;
   0 0 0 0 0 0 1 0;
   0 0 0 0 0 0 0 1];
processnoise =[ 0.000162 0.0001 0.0001 9.564e-8 0.0009784 0.001170]';
xinit = [0.06677 3.167 6.3816 308.89 299.95 0.00030926 0.3128 474.65 ]';
% ----------------------------------------------------
% Parameters of the discretized estimation model
tinit = 0;
sampletime = 300;
stepfinal = 400;
xmeanprev = [0.06677 3.167 6.3816 308.89 299.95 0.00030926 0.3128 474.65 ]';
xwmeanprev = [0 0 0 0 0 0]'; 
%Pold =diag ([0.005 0.4 0.4 0.001 0.001 0.0001 0.1 1.0 0.001 0.1 0.01 0.001 0.001 0.006 ]);
res = [0,0,0,0,0,0,0, xmeanprev']

xcovprev = diag([0.00005 0.000004 0.00004 0.00001 0.000001 0.0000001 0.000001 0.00001]);
xwcovprev = diag([0.000001 0.000001 0.00001 0.000001 0.0001 0.0006 ]);



% -------------------------------------
% Main loop
% ------------------------------------
for stepno = 0:stepfinal;
   if stepno <= 200;
        Qc =0.131;
        Qs=0.1275;
        Qm=0.105;
        Qi=0.03;
   else Qc = 0.1441;
        Qs=0.14025;
        Qm=0.1155
        Qi=0.033
   end

    % -------------------------------------
    % Simulate process
    % ------------------------------------
    xresplant = processsim (xinit, tinit, sampletime);
    [row,col] = size(xresplant);
    ymeas = C*xresplant(row,:)';
    xinit = xresplant(row,:)';
    ymeasured = ymeas + diag(processnoise*randn(6,1)');
    
    % -------------------------------------
    % Kalman Filter Estimation
    % ------------------------------------
    [xmeannew,xwmeannew,xcovnew,xwcovnew] = ukfsigma1(xmeanprev,xwmeanprev,xcovprev,xwcovprev,ymeasured);
    presenttime = sampletime*stepno;
    resrow = [presenttime, ymeasured', xmeannew'];
    res = [res ; resrow];
    xmeanprev = xmeannew;
    xcovprev = xcovnew;
    xwmeanprev = xwmeannew;
   xwcovprev = xwcovnew;
   %Pold=Pupdated;
end

% -------------------------------------
% Post processing


% ------------------------------------
figure
plot(res(:,1), res(:,8),'-r');
title('Inititor Conc  vs time');
xlabel('Time'),ylabel('Inititor Conc');


figure
plot(res(:,1), res(:,3),'-r', res(:,1), res(:,11));
title('Temperature  vs time');
xlabel('Time'),ylabel('Temperature');
legend('ymeasured','Estimated Temp');

figure
plot(res(:,1), res(:,9));
title('Monomer conc Vs Time');
xlabel ('Time'), ylabel('Monomer Conc');

figure
plot(res(:,1), res(:,2),'-r', res(:,1), res(:,10));
title('Solvent Conc  vs time');
xlabel('Time'),ylabel('Solvent Concentraion');
legend('ymeasured','Estimated Solvent Conc');

figure
plot(res(:,1), res(:,4),'-r', res(:,1), res(:,12));
title('Cooling Jacket Temp  vs time');
xlabel('Time'),ylabel('Cooling Jacket Temperature');
legend('ymeasured','Estimated Cooling Jacket Temp');

figure
plot(res(:,1), res(:,5),'-r', res(:,1), res(:,13));
title('Zeroth Moment  vs time');
xlabel('Time'),ylabel('Zeroth Moment');
legend('ymeasured','Estimated Zeroth Moment');

figure
plot(res(:,1), res(:,6),'-r', res(:,1), res(:,14));
title('First Moment  vs time');
xlabel('Time'),ylabel('First Moment');
legend('ymeasured','Estimated First Moment');

figure
plot(res(:,1), res(:,7),'-r', res(:,1), res(:,15));
title('Second Moment  vs time');
xlabel('Time'),ylabel('Second Moment');
legend('ymeasured','Estimated Second Moment');


