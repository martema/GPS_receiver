%Acquisition stage - Parallel acquisition in time domain
%
% February 2017

%% Step 2 - Parallel acquisition in time domain
clc
clear all
close all

% General parameters
fs = 12e6; % Sampling frequency (12 MHz)
f_IF = 4.152e6; % Carrier frequency (4.152 MHz)
PRN = 21; % SV number
t_Coh = 1e-3; % Coherent integration time (1 ms)
n_Blocks = 1; % Number of blocks generated
block_Length = 1e-3; % Block length (1ms)
n_Samples = fs * block_Length; % Number of samples per block

% Received signal
file_Id=fopen('signalrx_1.bin','rb');
signal_Rx = fread(file_Id,'double');

% Gold code
load('CA_Code.mat');
code_In = CA_Code(:,PRN);

% Code generation
chip_Fraction_In = 1e-7;
chip_Index_In = 1;
chip_Rate=1.023e6; % 1.023 MHz

code_Out = []; % Generated code with duration 1ms times n_Blocks

for k = 1:n_Blocks
    [code_Out_Block, chip_Fraction_Out] = SampleCode( fs, code_In, chip_Index_In, chip_Rate, chip_Fraction_In);
    code_Out = [code_Out code_Out_Block];
    chip_Fraction_In = chip_Fraction_Out;
end

% Carrier generation and CAF evaluation
delta_Fd = 2/(3*t_Coh); % Resolution in the doppler domain
fd_Min = -3000; % Minimum doppler (-3 MHz)
fd_Max = 3000; % Maximum doppler (+3 MHz)

carrier_Fraction_In = 0;
phase_In = 0;
carrier = [];
carrier_Cos = [];
carrier_Sin = [];

fd = fd_Min:delta_Fd:fd_Max;

n_Bins_Fd = length(fd);
L = length(code_Out);

caf_Cos = zeros(n_Bins_Fd, L);
caf_Sin = zeros(n_Bins_Fd, L);
caf = zeros(n_Bins_Fd, L);

tic
for k = 1:n_Bins_Fd
    [carrier, ~, carrier_Cos, carrier_Sin, ~ ] = GenerateCarriers( fs, f_IF + fd(k), phase_In, n_Samples, carrier_Fraction_In); % the carrier_Fraction_Out and the phase_Out are not needed since we are generating just one block of the carrier
    
    signal_BB = signal_Rx(1:L)' .* carrier;
    caf(k,:) = abs(CirCorrFFT(signal_BB, code_Out)).^2;

end
time = toc

[max_Caf_fd fd_Est] = max(max(caf,[],2)) % Estimated doppler frequency
[max_Caf_tau tau_Est] = max(max(caf,[],1)) % Estimated code delay

% Plots
% 3D CAF
figure
mesh((1:L),fd,caf);
set(gca, 'FontSize',12);
xlabel('$\bar{\tau}$','Interpreter','Latex','FontSize',20);
ylabel('${\bar{f_d}}$','Interpreter','Latex','FontSize',20)
zlabel('$\left|R(\bar{\tau},\bar{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
axis([1 L fd(1) fd(end) 0 1]);
title('CAF','FontSize',18)
saveas(gcf, 'Figures/PCAF3D', 'fig');
saveas(gcf, 'Figures/PCAF3D', 'png');

% 2D CAF - Code delay
figure
plot(caf(fd_Est,:));
set(gca, 'FontSize',12);
xlabel('$\bar{\tau}$','Interpreter','Latex','FontSize',20);
ylabel('$\left|R(\bar{\tau},\hat{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
title('Code delay domain 2D CAF','FontSize',18);
saveas(gcf, 'Figures/PCAFCode', 'fig');
saveas(gcf, 'Figures/PCAFCode', 'png');

% 2D CAF - Doppler frequency
figure
plot(fd,caf(:,tau_Est),'-*');
set(gca, 'FontSize',12);
xlabel('$\bar{f_d}$','Interpreter','Latex','FontSize',20);
ylabel('$\left|R(\hat{\tau},\bar{f_d})\right|$','Interpreter','Latex','FontSize',20)
title('Doppler domain 2D CAF','FontSize',18);
saveas(gcf, 'Figures/PCAFDoppler', 'fig');
saveas(gcf, 'Figures/PCAFDoppler', 'png');


