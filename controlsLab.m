%% Housekeeping
clc; close all; clear

%% loading in constants
Kg = 33.3;
Km = .0401;
Rm = 19.2;
J = 5e-4 + (0.2*0.2794^2) + 0.0015;
fc = 1.8;
% important as fuck
KptVec = [5 10 20 10 10 10];
KdtVec = [0 0 0 1 -1 -0.5];
Kpt = KptVec(5);
Kdt = KdtVec(6);

% data = readmatrix('Position Control Constants 2022 - Aerospace Modules.xlsx')

%% modeling closed loop behaviour of the system

% Equations 18
omega_n_square = @(Kpt, Kg, Km, J, Rm) (Kpt*Kg*Km)/(J*Rm);
zeta = @(Kpt, Kg, Km, J, Rm, Kdt) (Kg^2*Km^2 + Kdt*Kg*Km)/(2*sqrt(...
        Kpt*Kg*Km*J*Rm));

% calling anonymous func handles
bigW_squared = omega_n_square(Kpt, Kg, Km, J, Rm);
bigZeta = zeta(Kpt, Kg, Km, J, Rm, Kdt);

% [numVec denomVec] = dtsg(KptVec, KdtVec, omega_n_square, zeta);

num = [bigW_squared];
denom = [1 2*bigZeta*sqrt(bigW_squared) bigW_squared];
cltf = tf(num, denom);

% creating input function
t = 0:0.01:20; % Time from 0 to 20 seconds with a step size of 0.01

% Create a vector of zeros
frequency = 0.1; % Frequency of the square wave (in Hz)
amplitude = 0.5; % Amplitude of the square wave

% Define time span
t = 0:0.01:20; % Time from 0 to 20 seconds with a step size of 0.01

% Generate square wave
input = amplitude * square(2 * pi * frequency * t);

[x, t] = lsim(cltf, input, t);


%% plot you fuck
figure; hold on;
plot(t, x, 'LineWidth',1.2)
plot(t, input, LineWidth=1.5)
% yline(1.1, 'k--')
% yline(0.9, 'k--')
xlabel('time [s]')
ylabel('theta [rad]')
title('Kp = 10, Kd = 0')
grid minor;





% zeta: dampening ratio omega_n = natural frequency response

% Function Library
function [numVec denomVec] = dtsg(KptVec, KdtVec, omega_n_square, zeta)
    numVec = zeros(length(KptVec));
    denomVec = zeros(length(KptVec));
    
    for i = 1:length(KptVec)
        numVec(i) = omega_n_square(KptVec(i), Kg, Km, J, Rm);
        denomVec(i) = zeta(KptVec(i), Kg, Km, J, Rm, KdtVec(i));
    end    
end
