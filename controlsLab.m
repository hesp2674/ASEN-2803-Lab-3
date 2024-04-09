%% Housekeeping
clc; close all; clear

%% loading in constants
Kg = 33.3;
Km = .0401;
Rm = 19.2;
J = 5e-4 + (0.2*0.2794^2) + 0.0015;
fc = 1.8;
% important as fuck
Kpt = 1;
Kdt = 1;

% data = readmatrix('Position Control Constants 2022 - Aerospace Modules.xlsx')

%% modeling closed loop behaviour of the system
bigW_squared = omega_n_square(Kpt, Kg, Km, J, Rm);
bigZeta = zeta(Kpt, Kg, Km, J, Rm, Kdt);

num = [bigW_squared];
denom = [1 2*bigZeta*sqrt(bigW_squared) bigW_squared];
cltf = tf(num, denom);
[x, t] = step(cltf);





% zeta: dampening ratio omega_n = natural frequency response








%% Equations 17/18
% 18
omega_n_square = @(Kpt, Kg, Km, J, Rm) (Kpt*Kg*Km)/(J*Rm);
zeta = @(Kpt, Kg, Km, J, Rm, Kdt) (Kg^2*Km^2 + Kdt*Kg*Km)/(2*sqrt(...
        Kpt*Kg*Km*J*Rm));

