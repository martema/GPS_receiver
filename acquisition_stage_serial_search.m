% Acquisition stage - Serial search


clc
clear all
close all

%% Received signal and Gold sequence
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

%% Code generation
chip_Fraction_In = 1e-7;
chip_Index_In = 1;
chip_Rate=1.023e6; % 1.023 MHz

code_Out = []; % Generated code with duration 1ms times n_Blocks

for k = 1:n_Blocks
    [code_Out_Block, chip_Fraction_Out] = SampleCode( fs, code_In, chip_Index_In, chip_Rate, chip_Fraction_In);
    code_Out = [code_Out code_Out_Block];
    %chip_Fraction_In = chip_Fraction_Out;
end

%% Carrier generation and CAF evaluation
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
    for l = 1:L
        local_Replica_Cos = code_Out.*carrier_Cos;
        local_Replica_Sin = code_Out.*carrier_Sin;
        local_Replica = code_Out.*carrier;
        
        caf_Cos(k,l) = local_Replica_Cos * signal_Rx(l:l+L-1);
        caf_Sin(k,l) = local_Replica_Sin * signal_Rx(l:l+L-1);
        caf(k,l) = local_Replica * signal_Rx(l:l+L-1);
        
    end
end
caf = abs(caf/L).^2;
time = toc

[max_Caf fd_Est] = max(max(caf,[],2)) % Estimated doppler frequency
[max_Caf tau_Est] = max(max(caf,[],1)) % Estimated code delay

%% Plots
% 3D CAF
figure
mesh((1:L),fd,caf);
set(gca, 'FontSize',12);
xlabel('$\bar{\tau}$','Interpreter','Latex','FontSize',20);
ylabel('${\bar{f_d}}$','Interpreter','Latex','FontSize',20)
zlabel('$\left|R(\bar{\tau},\bar{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
axis([1 L fd(1) fd(end) 0 1]);
title('CAF','FontSize',18)
saveas(gcf, 'Figures/CAF3D', 'fig');
saveas(gcf, 'Figures/CAF3D', 'png');

% 2D CAF - Code delay
figure
plot(caf(fd_Est,:));
set(gca, 'FontSize',12);
xlabel('$\bar{\tau}$','Interpreter','Latex','FontSize',20);
ylabel('$\left|R(\bar{\tau},\hat{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
title('Code delay domain 2D CAF','FontSize',18);
saveas(gcf, 'Figures/CAFCode', 'fig');
saveas(gcf, 'Figures/CAFCode', 'png');

% 2D CAF - Doppler frequency
figure
plot(fd,caf(:,tau_Est),'-o','LineWidth',1.5);
set(gca, 'FontSize',12);
xlabel('$\bar{f_d}$','Interpreter','Latex','FontSize',20);
ylabel('$\left|R(\hat{\tau},\bar{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
title('Doppler domain 2D CAF','FontSize',18);
saveas(gcf, 'Figures/CAFDoppler', 'fig');
saveas(gcf, 'Figures/CAFDoppler', 'png');
