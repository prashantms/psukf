%ODEs of the polystyrene model
%in this file variabel 'y' is being changed to 'x'.
%created by santosh on 10 may 2014.
function dxdt=odes_plant(t,x)
global f R kp0 kd0 kt0 kf0 kfs0 Ep Ed Et Ef Efs deltaHr UA rhoCp rhoCpc Qi Qs Qm Qc V Vc Cif Cmf Tf Tcf Cpdot Csf dinit 
% state variables - readable names:
Ci = x(1); Cm = x(2); Cs = x(3); T = x(4); Tc = x(5); lambda0 = x(6); lambda1 = x(7); lambda2 = x(8);


kd = kd0*exp(-Ed/(R*T));
kp = kp0*exp(-Ep/(R*T));
ktc = kt0*exp(-Et/(R*T));
kf = kf0*exp(-Ef/(R*T));
kfs = kfs0*exp(-Efs/(R*T));

Cpdot = (2*f*kd*Ci/ktc)^(1/2);   %Cpdot=Cgp ----> growing polymer concentration

alpha = (kp*Cm)/(kp*Cm + ktc*Cpdot);
Qout = Qi + Qm + Qs;



dxdt(1) = ((Qi*Cif-Qout*Ci)/V)-kd*Ci;
dxdt(2) = ((Qm*Cmf-Qout*Cm)/V)-kp*Cm*Cpdot;
dxdt(3) = ((Qs*Csf + Qi*Cif-Qout*Cs)/V);
dxdt(4) = (Qout*(Tf-T)/V)+((-deltaHr)*kp*Cm*Cpdot/(rhoCp))-(UA*(T-Tc)/(rhoCp*V));
dxdt(5) = (Qc*(Tcf-Tc)/Vc)+ (UA*(T-Tc)/(rhoCpc*Vc));
dxdt(6) = (kfs*Cs*kf*Cm)*alpha*Cpdot+ ((ktc*(Cpdot)^2)/2)- Qout*lambda0/V ;
dxdt(7) = (((Cpdot)^2)/(1-alpha)*ktc) - Qout*lambda1/V ;
dxdt(8) = ((3*(Cpdot)^2)*ktc/(1-alpha)^2)-Qout*lambda2/V ;

dxdt = dxdt';
