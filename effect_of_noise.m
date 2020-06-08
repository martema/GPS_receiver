%Acquisition stage - Effect of noise


%% Step 3 - Effect of noise
clc
clear all
close all

% General parameters
fs = 12e6; % Sampling frequency (12 MHz)
f_IF = 4.152e6; % Carrier frequency (4.152 MHz)
t_Coh = 1e-3; % Coherent integration time (1 ms)
block_Length = 1e-3; % Block length (1ms)

% Received signal 
file_Id=fopen('signalrx_2.bin','rb');
signal_Rx = fread(file_Id,'double');

% Gold code
load('CA_Code.mat');

p_Fa_Cell = 1e-9;
PRN_1 = 3; 
PRN_2 = 16; 
PRN_3 = 22;

code_In_1 = CA_Code(:,PRN_1);
code_In_2 = CA_Code(:,PRN_2);
code_In_3 = CA_Code(:,PRN_3);

% Code generation
chip_Fraction_In = 1e-7;
chip_Index_In = 1;
chip_Rate=1.023e6; % 1.023 MHz

code_Out_1 = []; % Generated code with duration 1ms times n_Blocks
code_Out_2 = []; % Generated code with duration 1ms times n_Blocks
code_Out_3 = []; % Generated code with duration 1ms times n_Blocks

n_Blocks = t_Coh / block_Length; % Number of blocks generated
n_Samples = fs * block_Length; % Number of samples per block

for k = 1:n_Blocks
    [code_Out_Block_1, chip_Fraction_Out] = SampleCode( fs, code_In_1, chip_Index_In, chip_Rate, chip_Fraction_In);
    code_Out_1 = [code_Out_1 code_Out_Block_1];
    chip_Fraction_In = chip_Fraction_Out;
end

chip_Fraction_In = 1e-7;

for k = 1:n_Blocks
    [code_Out_Block_2, chip_Fraction_Out] = SampleCode( fs, code_In_2, chip_Index_In, chip_Rate, chip_Fraction_In);
    code_Out_2 = [code_Out_2 code_Out_Block_2];
    chip_Fraction_In = chip_Fraction_Out;
end

chip_Fraction_In = 1e-7;

for k = 1:n_Blocks
    [code_Out_Block_3, chip_Fraction_Out] = SampleCode( fs, code_In_3, chip_Index_In, chip_Rate, chip_Fraction_In);
    code_Out_3 = [code_Out_3 code_Out_Block_3];
    chip_Fraction_In = chip_Fraction_Out;
end


% Carrier generation and CAF evaluation signal 1
delta_Fd = 2/(3*t_Coh); % Resolution in the doppler domain
fd_Min = -3000; % Minimum doppler (-3 MHz)
fd_Max = 3000; % Maximum doppler (+3 MHz)

carrier_Fraction_In = 0;
phase_In = 0;
carrier = [];

fd = fd_Min:delta_Fd:fd_Max;

n_Bins_Fd = length(fd);
L = length(code_Out_1);

caf_Cos = zeros(n_Bins_Fd, L);
caf_Sin = zeros(n_Bins_Fd, L);
caf = zeros(n_Bins_Fd, L);

for k = 1:n_Bins_Fd
    [carrier, ~, ~, ~, ~ ] = GenerateCarriers( fs, f_IF + fd(k), phase_In, n_Samples, carrier_Fraction_In); 
    signal_BB = signal_Rx(1:L)' .* carrier;
    caf_1(k,:) = abs(CirCorrFFT(signal_BB, code_Out_1)).^2;

end

[max_Caf_1 fd_Est_1] = max(max(caf_1,[],2)); % Estimated doppler frequency
[max_Caf_1 tau_Est_1] = max(max(caf_1,[],1)); % Estimated code delay


% 3D CAF
figure
mesh(1:L, fd, caf_1);
set(gca, 'FontSize',12);
set(gca, 'YLim',[fd(1) fd(end)]);
set(gca, 'XLim',[1 L]);
xlabel('$\bar{\tau}$','Interpreter','Latex','FontSize',20);
ylabel('${\bar{f_d}}$','Interpreter','Latex','FontSize',20)
zlabel('$\left|R(\bar{\tau},\bar{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
title('CAF 45 dbHz','FontSize',18)
saveas(gcf, 'Figures/CAF_Noise_45', 'fig');
saveas(gcf, 'Figures/CAF_Noise_45', 'png');

% Carrier generation and CAF evaluation signal 2
delta_Fd=2/(3*t_Coh); % Resolution in the doppler domain
fd_Min = -3000; % Minimum doppler (-3 MHz)
fd_Max = 3000; % Maximum doppler (+3 MHz)

carrier_Fraction_In = 0;
phase_In = 0;
carrier = [];

fd = fd_Min:delta_Fd:fd_Max;

n_Bins_Fd = length(fd);
L = length(code_Out_2);

caf_Cos = zeros(n_Bins_Fd, L);
caf_Sin = zeros(n_Bins_Fd, L);
caf = zeros(n_Bins_Fd, L);

for k = 1:n_Bins_Fd
    [carrier, ~, ~, ~, ~ ] = GenerateCarriers( fs, f_IF + fd(k), phase_In, n_Samples, carrier_Fraction_In); 
    signal_BB = signal_Rx(1:L)' .* carrier;
    caf_2(k,:) = abs(CirCorrFFT(signal_BB, code_Out_2)).^2;

end

[max_Caf_2 fd_Est_2] = max(max(caf_2,[],2)); % Estimated doppler frequency
[max_Caf_2 tau_Est_2] = max(max(caf_2,[],1)); % Estimated code delay

% Plots
% 3D CAF
figure
mesh(1:L, fd,caf_2);
set(gca, 'FontSize',12);
set(gca, 'YLim',[fd(1) fd(end)]);
set(gca, 'XLim',[1 L]);
xlabel('$\bar{\tau}$','Interpreter','Latex','FontSize',20);
ylabel('${\bar{f_d}}$','Interpreter','Latex','FontSize',20)
zlabel('$\left|R(\bar{\tau},\bar{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
title('CAF 50 dbHz','FontSize',18)
saveas(gcf, 'Figures/CAF_Noise_50', 'fig');
saveas(gcf, 'Figures/CAF_Noise_50', 'png');

% Carrier generation and CAF evaluation signal 3
delta_Fd=2/(3*t_Coh); % Resolution in the doppler domain
fd_Min = -3000; % Minimum doppler (-3 MHz)
fd_Max = 3000; % Maximum doppler (+3 MHz)

carrier_Fraction_In = 0;
phase_In = 0;
carrier = [];
carrier_Cos = [];
carrier_Sin = [];

fd = fd_Min:delta_Fd:fd_Max;

n_Bins_Fd = length(fd);
L = length(code_Out_3);

caf_Cos = zeros(n_Bins_Fd, L);
caf_Sin = zeros(n_Bins_Fd, L);
caf = zeros(n_Bins_Fd, L);

for k = 1:n_Bins_Fd
    [carrier, ~, ~, ~, ~ ] = GenerateCarriers( fs, f_IF + fd(k), phase_In, n_Samples, carrier_Fraction_In); 
    signal_BB = signal_Rx(1:L)' .* carrier;
    caf_3(k,:) = abs(CirCorrFFT(signal_BB, code_Out_3)).^2;

end

[max_Caf_3 fd_Est_3] = max(max(caf_3,[],2)); % Estimated doppler frequency
[max_Caf_3 tau_Est_3] = max(max(caf_3,[],1)); % Estimated code delay

% Plots
% 3D CAF
figure
mesh(1:L, fd,caf_3);
set(gca, 'FontSize',12);
set(gca, 'YLim',[fd(1) fd(end)]);
set(gca, 'XLim',[1 L]);
xlabel('$\bar{\tau}$','Interpreter','Latex','FontSize',20);
ylabel('${\bar{f_d}}$','Interpreter','Latex','FontSize',20)
zlabel('$\left|R(\bar{\tau},\bar{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
title('CAF 60 dbHz','FontSize',18)
saveas(gcf, 'Figures/CAF_Noise_60', 'fig');
saveas(gcf, 'Figures/CAF_Noise_60', 'png');
