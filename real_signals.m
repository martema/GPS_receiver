%Acquisition stage - Paralell acquisition real signals


%% Step 4 - Real GNSS signal
clc
close all
clear all

% General parameters
fs = 16.3676e6; % Sampling frequency (12 MHz)
f_IF = 4.1304e6; % Carrier frequency (4.152 MHz)
t_Coh = 1e-3; % Coherent integration time (1 ms)
n_Blocks = 1; % Number of blocks generated
block_Length = 1e-3; % Block length (1ms)
n_Samples = fs * block_Length; % Number of samples per block

% Received signal 
file_Id=fopen('signalrx_real.bin','rb');
signal_Rx = fread(file_Id,'int8');

% Gold code
load('CA_Code.mat');
num_Code = size(CA_Code,2);

% Code generation
chip_Fraction_In = 1e-7;
chip_Index_In = 1;
chip_Rate=1.023e6; % 1.023 MHz

code_Out = [];
for k = 1:num_Code
    [code_Out_K, ~] = SampleCode( fs, CA_Code(:,k), chip_Index_In, chip_Rate, chip_Fraction_In);
    code_Out = [code_Out; code_Out_K];
end

% Carrier generation and CAF evaluation 
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
L = size(code_Out,2);

caf = zeros(n_Bins_Fd, L);
max_Caf = zeros(num_Code,1);
fd_Est = zeros(num_Code,1);
tau_Est = zeros(num_Code,1);

for n = 1:num_Code
    for k = 1:n_Bins_Fd
        [carrier, ~, carrier_Cos, carrier_Sin, ~ ] = GenerateCarriers( fs, f_IF + fd(k), phase_In, n_Samples, carrier_Fraction_In); % the carrier_Fraction_Out and the phase_Out are not needed since we are generating just one block of the carrier

        signal_BB = signal_Rx(1:L)' .* carrier;
        caf(k,:) = abs(CirCorrFFT(signal_BB, code_Out(n,:))).^2;

    end

    [max_Caf(n) fd_Est(n)] = max(max(caf,[],2)); % Estimated doppler frequency
    [max_Caf(n) tau_Est(n)] = max(max(caf,[],1)); % Estimated code delay

    figure
    mesh((1:L),fd,caf);
    set(gca, 'FontSize',12,'XLim',[1 L],'Ylim',[fd(1) fd(end)])
    xlabel('$\bar{\tau}$','Interpreter','Latex','FontSize',20);
    ylabel('${\bar{f_d}}$','Interpreter','Latex','FontSize',20)
    zlabel('$\left|R(\bar{\tau},\bar{f_d})\right|^2$','Interpreter','Latex','FontSize',20)
    title('CAF','FontSize',18)

end
