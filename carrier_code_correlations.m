% Generation of carrier and code, correlations
%
% February 2017

clc
clear all
close all

%% Step 1 - Generation of the Gold codes used by GPS 

[CA_Code] = GoldCodeGenerator();
save('CA_Code.mat','CA_Code');

% Autocorrelation of C/A code 1
PRN_1 = 1;
code_1 = CA_Code(:,PRN_1)';
p = length(code_1);
auto_Corr =  CirCorr(code_1, [code_1(floor(p/2):end) code_1(1:(floor(p/2-1)))]);

figure
plot(auto_Corr);
set(gca, 'XLim', [0 p]);
set(gca, 'Fontsize',12);
xlabel('Chip', 'Fontsize', 16, 'Fontweight', 'bold');
ylabel('R[l]', 'Fontsize', 16, 'Fontweight', 'bold');
title('Autocorrelation function','FontSize', 18);
saveas(gcf,'figures/Autocorrelation','fig');
saveas(gcf,'figures/Autocorrelation','png');

% Correlation of C/A code 
PRN_1 = 1;
PRN_2 = 2;
code_1 = CA_Code(:,PRN_1)';
code_2 = CA_Code(:,PRN_2)';
p = length(code_1);
cross_Corr =  CirCorr(code_1, code_2);

figure
plot(cross_Corr);
set(gca, 'Fontsize', 12);
set(gca, 'XLim', [0 p]);
xlabel('Chip', 'Fontsize', 16, 'Fontweight', 'bold');
ylabel('R[l]', 'Fontsize', 16, 'Fontweight', 'bold');
title('Crosscorrelation function','FontSize', 18);
saveas(gcf,'figures/Crosscorrelation','fig');
saveas(gcf,'figures/Crosscorrelation','png');

%% Step 2 - Generation of the local replica of the IF carrier

% Common parameters
signal_Length = 20e-3; % Length of the signal (20 ms)
block_Length = 1e-3; % Length of the signal block (1 ms)
initial_Phase = 0.1; % Initial phase of the signal (rad)
n_Blocks = signal_Length/block_Length; % Number of blocks of the signal

% Parameter for the case INTEGER #SAMPLES
fs_1 = 16.368e6; % Sampling frequency (16.368 MHz)
f_IF_1 = 4.092e6; % Carrier frequency (4.092 MHz)
n_Samples_1 = fs_1*block_Length; % Number of samples per block


% Generation of the carrier with INTEGER #SAMPLES
phase_In = initial_Phase;
carrier_Fraction_In = 0;
% carrier_1 = zeros(1,n_Samples_1*n_Blocks);
% carrier_Cos_1 = zeros(1,n_Samples_1*n_Blocks);
% carrier_Sin_1 = zeros(1,n_Samples_1*n_Blocks);
carrier_1 = [];
carrier_Cos_1 = [];
carrier_Sin_1 = [];

for k = 1:n_Blocks
%     index_1 = (k-1) * n_Samples_1 + 1;
%     index_2 = k*n_Samples_1;
%     [carrier_1(index_1:index_2), phase_Out, carrier_Cos_1(index_1:index_2), carrier_Sin_1(index_1:index_2)]=...
%         GenerateCarriers(sampling_Rate_1, carrier_Frequency_1, phase_In, n_Samples_1);
    [carrier_1_Block, phase_Out, carrier_Cos_1_Block, carrier_Sin_1_Block, carrier_Fraction_Out]=...
        GenerateCarriers(fs_1, f_IF_1, phase_In, n_Samples_1, carrier_Fraction_In);
    carrier_1 = [carrier_1 carrier_1_Block];
    carrier_Cos_1 = [carrier_Cos_1 carrier_Cos_1_Block];
    carrier_Sin_1 = [carrier_Sin_1 carrier_Sin_1_Block];
    phase_In = phase_Out;
    carrier_Fraction_In = carrier_Fraction_Out;
end 

% Parameter for the case FRACTIONAL #SAMPLES

fs_2 = 16.3676e6; % Sampling frequency (16.3676 MHz)
f_IF_2 = 4.1304e6; % Carrier frequency (4.1304 MHz)
n_Samples_2 = fs_2*block_Length; % Number of samples per block

% Generation of the carrier with FRACTIONAL #SAMPLES
phase_In = initial_Phase;
carrier_2 = [];
carrier_Cos_2 = [];
carrier_Sin_2 = [];

for k = 1:n_Blocks
    [carrier_2_Block, phase_Out, carrier_Cos_2_Block, carrier_Sin_2_Block, carrier_Fraction_Out]=...
        GenerateCarriers(fs_2, f_IF_2, phase_In, n_Samples_2, carrier_Fraction_In);
    carrier_2 = [carrier_2 carrier_2_Block];
    carrier_Cos_2 = [carrier_Cos_2 carrier_Cos_2_Block];
    carrier_Sin_2 = [carrier_Sin_2 carrier_Sin_2_Block];
    phase_In = phase_Out;
    carrier_Fraction_In = carrier_Fraction_Out;
end 

% Generation of a carrier with fractional #SAMPLES without phase update
% between one block and the other
phase_In = initial_Phase;
carrier_Cos_3 = [];
for k = 1:n_Blocks
    [carrier_3_Block, phase_Out, carrier_Cos_3_Block, ~, carrier_Fraction_Out]=...
        GenerateCarriers(fs_2, f_IF_2, phase_In, n_Samples_2, carrier_Fraction_In);
    carrier_Cos_3 = [carrier_Cos_3 carrier_Cos_3_Block];
    carrier_Fraction_In = carrier_Fraction_Out;
end 


%% Step 3 - Generation of the local replica of the despreading code

% Common parameters
signal_Length = 20e-3; % Length of the signal (20 ms)
block_Length = 1e-3; % Length of the signal block (1 ms)
n_Blocks = signal_Length/block_Length; % Number of blocks of the signal
chip_Rate = 1.023e6; % Chip rate (1.023 MHz)
PRN_Number = 1;

chip_Index_In = 1; % Index #the starting chip 
chip_Fraction_In = 1e-7 % ~0;
code_In = CA_Code(:,PRN_Number);

% Parameter for the case INTEGER #SAMPLES
fs_1 = 16.368e6; % Sampling frequency (16.3686 MHz)

% Generation of the code with INTEGER #SAMPLES
code_Out_1 = [];
code_In = CA_Code(:,PRN_Number);

for k = 1:n_Blocks
    [code_Out_1_Block, chip_Fraction_Out] = SampleCode(fs_1, code_In, chip_Index_In, chip_Rate, chip_Fraction_In);
    code_Out_1 = [code_Out_1 code_Out_1_Block];
    chip_Fraction_In = chip_Fraction_Out;
end

% Parameter for the case FRACTIONAL #SAMPLES
fs_2 = 16.3676e6; % Sampling frequency (16.3676 MHz)
% Generation of the code with FRACTIONAL #SAMPLES
code_Out_2 = [];

for k = 1:n_Blocks
    [code_Out_2_Block, chip_Fraction_Out] = SampleCode(fs_2, code_In, chip_Index_In, chip_Rate, chip_Fraction_In);
    code_Out_2 = [code_Out_2 code_Out_2_Block];
    chip_Fraction_In = chip_Fraction_Out;
end

%% Step 4 - Plots
% % Analog code integer number of samples
% mul = 10;
% fs_Anal = fs_1 * mul; % Sampling frequency of the analog code
% n_Samples_Anal = fs_Anal * block_Length;
% 
% code_Anal = [];
% for k = 1:n_Blocks
%     [code_Anal_Block, chip_Fraction_Out] = SampleCode(fs_Anal, code_In, chip_Index_In, chip_Rate, chip_Fraction_In);
%     code_Anal = [code_Anal code_Anal_Block];
%     chip_Fraction_In = chip_Fraction_Out;
% end
% 
% t_Anal = (0:length(code_Anal)-1)*1/(fs_Anal)*1e3; % Time scale of the analog code

% Code plots
n_Chips_Plot = 8; % Number of plotted chips
n_Chips_Code = 1023; % Number of chips per code

t_1 = (0:length(code_Out_1)-1)/fs_1*1e3; % Time scale generated code 1(ms)
t_1 = (0:length(code_Out_1)-1)/fs_2*1e3; % Time scale generated code 2(ms)

figure,
subplot(2,1,1), plot(t_1(n_Samples_1-n_Samples_1/n_Chips_Code*n_Chips_Plot:n_Samples_1+n_Samples_1/n_Chips_Code*n_Chips_Plot), code_Out_1(n_Samples_1-n_Samples_1/n_Chips_Code*n_Chips_Plot:n_Samples_1+n_Samples_1/n_Chips_Code*n_Chips_Plot));
xlabel('time (ms)','FontWeight','bold','FontSize', 12)
title('Code with integer number of samples','FontSize',18)
subplot(2,1,2), plot(t_1(floor(n_Samples_2-n_Samples_2/n_Chips_Code*n_Chips_Plot):floor(n_Samples_2+n_Samples_2/n_Chips_Code*n_Chips_Plot)), code_Out_2(floor(n_Samples_2-n_Samples_2/n_Chips_Code*n_Chips_Plot):floor(n_Samples_2+n_Samples_2/n_Chips_Code*n_Chips_Plot)));
xlabel('time (ms)','FontWeight','bold','FontSize', 12)
title('Code with fractional number of samples','FontSize',18)
saveas(gcf,'figures/CodeInteger','fig');
saveas(gcf,'figures/CodeInteger','png');

figure,
subplot(2,1,1), plot(code_Out_1(1:n_Samples_1/n_Chips_Code*n_Chips_Plot),'-o');
xlabel('Sample','FontWeight','bold','FontSize', 12)
title('Code with integer number of samples','FontSize',18)
set(gca,'XTick',[1:5:150])
subplot(2,1,2), plot(code_Out_2(1:floor(n_Samples_2/n_Chips_Code*n_Chips_Plot)),'-o');
xlabel('Sample','FontWeight','bold','FontSize', 12)
title('Code with fractional number of samples','FontSize',18)
set(gca,'XTick',[1:5:150])
saveas(gcf,'figures/CodeFractional','fig');
saveas(gcf,'figures/CodeFractional','png');

% Carrier plots
% Carrier plots integer number of samples t = (0:length(code_Out_1)-1)/fs_1*1e-3; % Time scale generated code 1(ms)
mul = 10;
n_Period_Plot = 4;
fs_Anal = fs_1*mul; % Sampling frequency of the analog sinusoid
t = (0:length(code_Out_1)-1)/fs_1*1e3; % Time scale generated code 1(ms)
t_Anal = (0:length(code_Out_1)-1)*1/fs_Anal*1e3; % Times vector for the analog sinusoid

figure, plot(t(n_Samples_1-(fs_1/f_IF_1*n_Period_Plot):n_Samples_1+(fs_1/f_IF_1*n_Period_Plot)),carrier_Cos_1(n_Samples_1-(fs_1/f_IF_1*n_Period_Plot):n_Samples_1+(fs_1/f_IF_1*n_Period_Plot)),'o', 'MarkerFaceColor','b','MarkerSize',4);
hold on
plot(t_Anal(((n_Samples_1-(fs_1/f_IF_1*n_Period_Plot))-1)*mul:(n_Samples_1+(fs_1/f_IF_1*n_Period_Plot))*mul),sqrt(2)*cos(2*pi*f_IF_1*t_Anal(((n_Samples_1-(fs_1/f_IF_1*n_Period_Plot))-1)*mul:(n_Samples_1+(fs_1/f_IF_1*n_Period_Plot))*mul)*1e-3+initial_Phase),'r');
axis('tight')
xlabel('time (ms)','FontWeight','bold','FontSize', 12)
title('Carrier with integer number of samples','FontSize',18)
saveas(gcf,'figures/CarrierInteger','fig');
saveas(gcf,'figures/CarrierInteger','png');

% Carrier plots fractional number of samples
t = (0:length(code_Out_1)-1)/fs_2*1e3; % Time scale generated code 1(ms)
fs_Anal = fs_2*mul; % Sampling frequency of the analog sinusoid
t_Anal = (0:length(code_Out_1)-1)*1/fs_Anal*1e3; % Times vector for the analog sinusoid

figure, plot(t(n_Samples_1-(fs_1/f_IF_1*n_Period_Plot):n_Samples_1+(fs_1/f_IF_1*n_Period_Plot)),carrier_Cos_2(n_Samples_1-(fs_1/f_IF_1*n_Period_Plot):n_Samples_1+(fs_1/f_IF_1*n_Period_Plot)),'o', 'MarkerFaceColor','b','MarkerSize',4);
hold on
plot(t_Anal(((n_Samples_1-(fs_1/f_IF_1*n_Period_Plot))-1)*mul:(n_Samples_1+(fs_1/f_IF_1*n_Period_Plot))*mul),sqrt(2)*cos(2*pi*f_IF_2*t_Anal(((n_Samples_1-(fs_1/f_IF_1*n_Period_Plot))-1)*mul:(n_Samples_1+(fs_1/f_IF_1*n_Period_Plot))*mul)*1e-3+initial_Phase),'r');
axis('tight')
xlabel('time (ms)','FontWeight','bold','FontSize', 12)
title('Carrier with fractinalnumber of samples','FontSize',18)
saveas(gcf,'figures/CarrierFractional','fig');
saveas(gcf,'figures/CarrierFractional','png');

% Carrier plots fractional number of samples without phase update
n_Period_Plot = 2;

t = (0:length(code_Out_1)-1)/fs_2*1e3; % Time scale generated code 1(ms)
fs_Anal = fs_2*mul; % Sampling frequency of the analog sinusoid
t_Anal = (0:length(code_Out_1)-1)*1/fs_Anal*1e3; % Times vector for the analog sinusoid

figure, plot(t(3*n_Samples_1-(fs_1/f_IF_1*n_Period_Plot):3*n_Samples_1+(fs_1/f_IF_1*n_Period_Plot)),carrier_Cos_3(3*n_Samples_1-(fs_1/f_IF_1*n_Period_Plot):3*n_Samples_1+(fs_1/f_IF_1*n_Period_Plot)),'-o', 'MarkerFaceColor','b','MarkerSize',4);
axis('tight')
xlabel('time (ms)','FontWeight','bold','FontSize', 12)
title('Carrier without phase update','FontSize',18)
saveas(gcf,'figures/CarrierFractionalNoPhase','fig');
saveas(gcf,'figures/CarrierFractionalNoPhase','png');

% Carrier times code integer number of samples
t = (0:length(code_Out_1)-1)/fs_1*1e3; % Time scale generated code 1(ms)
figure
plot(t(1:64),carrier_Cos_1(1,1:64).*code_Out_1(1:64),'-o');
axis('tight')
xlabel('time (ms)','FontWeight','bold','FontSize', 12)
title('Modulated code with integer number of samples','FontSize',18)
saveas(gcf,'figures/ModulatedInteger','fig');
saveas(gcf,'figures/ModulatedInteger','png');

% Carrier times code fractional number of samples
t = (0:length(code_Out_1)-1)/fs_2*1e3; % Time scale generated code 1(ms)
figure
plot(t(16368-8:16368+8),carrier_Cos_2(16368-8:16368+8).*code_Out_2(16368-8:16368+8),'-o');
axis('tight')
xlabel('time (ms)','FontWeight','bold','FontSize', 12)
title('Modulated code with fractional number of samples','FontSize',18)
saveas(gcf,'figures/ModulatedFractional','fig');
saveas(gcf,'figures/ModulatedFractional','png');

% Spectra integer number of samples
Nw = 1000;
NoOverlap = 900;
win = rectwin(Nw);

[pdf_code f1] = pwelch(code_Out_1,win,NoOverlap,Nw,1);
[pdf_carrier f2] = pwelch(carrier_Cos_1(1,:),win,NoOverlap,Nw,1);
[pdf_modulated f3] = pwelch(code_Out_1.*carrier_Cos_1,win,NoOverlap,Nw,1);

figure,
subplot(3,1,1), plot(f1,pdf_code)
xlabel('f_n','FontWeight','bold','FontSize', 12)
title('Code spectrum','FontSize',18)
subplot(3,1,2), plot(f2,pdf_carrier)
xlabel('f_n','FontWeight','bold','FontSize', 12)
title('Carrier spectrum','FontSize',18)
subplot(3,1,3), plot(f3,pdf_modulated)
xlabel('f_n','FontWeight','bold','FontSize', 12)
title('Modulated code spectrum','FontSize',18)
saveas(gcf,'figures/SpectraInteger','fig');
saveas(gcf,'figures/SpectraInteger','png');

% Spectra fractional number of samples
Nw = 1000;
NoOverlap = 900;
win = rectwin(Nw);

[pdf_code f1] = pwelch(code_Out_2,win,NoOverlap,Nw,1);
[pdf_carrier f2] = pwelch(carrier_Cos_2(1,:),win,NoOverlap,Nw,1);
[pdf_modulated f3] = pwelch(code_Out_2.*carrier_Cos_2,win,NoOverlap,Nw,1);

figure,
subplot(3,1,1), plot(f1,pdf_code)
xlabel('f_n','FontWeight','bold','FontSize', 12)
title('Code spectrum','FontSize',18)
subplot(3,1,2), plot(f2,pdf_carrier)
xlabel('f_n','FontWeight','bold','FontSize', 12)
title('Carrier spectrum','FontSize',18)
subplot(3,1,3), plot(f3,pdf_modulated)
xlabel('f_n','FontWeight','bold','FontSize', 12)
title('Modulated code spectrum','FontSize',18)
saveas(gcf,'figures/SpectraFractional','fig');
saveas(gcf,'figures/SpectraFractional','png');

%% Step 5 - Correlation functions
% a) Auto-correlation function

% Common parameters
fs = 16.368e6; % Sampling frequency (16.3686 MHz)
chip_Rate = 1.023e6; % Chip rate (1.023 MHz)
block_Length = 1e-3; % Length of the signal block (1 ms)

% Parameter signal 1
PRN_1 = 1; % SV number for signal 1
code_In_1 = CA_Code(:,PRN_1);
signal_Length_1 = 1e-3; % Length of the first signal (1ms)
chip_Index_In_1 = 1; % Signal 1 starts from chip 1
chip_Fraction_In_1 = 1e-7;
n_Blocks_1 = signal_Length_1/block_Length; % Number of blocks of the signal 1


% Parameter signal 2
PRN_2 = 1; % SV number for signal 2
code_In_2 = CA_Code(:,PRN_2);
signal_Length_2 = 2e-3; % Length of the second signa (2ms)
chip_Index_In_2 = 800;% Signal 2 starts from chip 800
chip_Fraction_In_2 = 1e-7;
n_Blocks_2 = signal_Length_2/block_Length; % Number of blocks of the signal 2

% Signal 1 generation
code_Out_1 = [];
for k = 1:n_Blocks_1
    [code_Out_1_Block, chip_Fraction_Out_1] = SampleCode(fs, code_In_1, chip_Index_In_1, chip_Rate, chip_Fraction_In_1);
    code_Out_1 = [code_Out_1 code_Out_1_Block];
    chip_Fraction_In_1 = chip_Fraction_Out_1;
end

% Signal 2 generation
code_Out_2 = [];
for k = 1:n_Blocks_2
    [code_Out_2_Block, chip_Fraction_Out_2] = SampleCode(fs, code_In_2, chip_Index_In_2, chip_Rate, chip_Fraction_In_2);
    code_Out_2 = [code_Out_2 code_Out_2_Block];
    chip_Fraction_In_2 = chip_Fraction_Out_2;
end

l_Corr_1 = LinCorr(code_Out_1, code_Out_2);
c_Corr_1 = CirCorr(code_Out_1, code_Out_2);
c_Corr_FFT_1 = CirCorrFFT(code_Out_2, code_Out_1);

figure
subplot(3,1,1), plot(l_Corr_1);
xlabel('Sample','FontWeight','bold','FontSize', 12)
title('Linear correlation','FontSize',18)
subplot(3,1,2), plot(c_Corr_1);
xlabel('Sample','FontWeight','bold','FontSize', 12)
title('Circular correlation','FontSize',18)
subplot(3,1,3), plot(c_Corr_FFT_1);
xlabel('Sample','FontWeight','bold','FontSize', 12)
title('Circular correlation using FFT','FontSize',18)
saveas(gcf,'figures/Autocorr3tipi','fig');
saveas(gcf,'figures/Autocorr3tipi','png');

% a) Cross-correlation function

% Common parameters
fs = 16.3676e6; % Sampling frequency (16.3686 MHz)
chip_Rate = 1.023e6; % Chip rate (1.023 MHz)
block_Length = 1e-3; % Length of the signal block (1 ms)

% Parameter signal 1
PRN_1 = 1; % SV number for signal 1
code_In_1 = CA_Code(:,PRN_1);
signal_Length_1 = 1e-3; % Length of the first signal (1ms)
chip_Index_In_1 = 1; % Signal 1 starts from chip 1
chip_Fraction_In_1 = 1e-7;
n_Blocks_1 = signal_Length_1/block_Length; % Number of blocks of the signal 1


% Parameter signal 2
PRN_2 = 2; % SV number for signal 2
code_In_2 = CA_Code(:,PRN_2);
signal_Length_2 = 2e-3; % Length of the second signa (2ms)
chip_Index_In_2 = 800;% Signal 2 starts from chip 800
chip_Fraction_In_2 = 1e-7;
n_Blocks_2 = signal_Length_2/block_Length; % Number of blocks of the signal 2

% Signal 1 generation
code_Out_1 = [];
for k = 1:n_Blocks_1
    [code_Out_1_Block, chip_Fraction_Out_1] = SampleCode(fs, code_In_1, chip_Index_In_1, chip_Rate, chip_Fraction_In_1);
    code_Out_1 = [code_Out_1 code_Out_1_Block];
    chip_Fraction_In_1 = chip_Fraction_Out_1;
end

% Signal 2 generation
code_Out_2 = [];
for k = 1:n_Blocks_2
    [code_Out_2_Block, chip_Fraction_Out_2] = SampleCode(fs, code_In_2, chip_Index_In_2, chip_Rate, chip_Fraction_In_2);
    code_Out_2 = [code_Out_2 code_Out_2_Block];
    chip_Fraction_In_2 = chip_Fraction_Out_2;
end

l_Corr_1 = LinCorr(code_Out_1, code_Out_2);
tic
c_Corr_1 = CirCorr(code_Out_1, code_Out_2);
t1 = toc;
tic
c_Corr_FFT_1 = CirCorrFFT(code_Out_2, code_Out_1);
t2 = toc;

figure
subplot(3,1,1), plot(l_Corr_1);
xlabel('Sample','FontWeight','bold','FontSize', 12)
title('Linear correlation','FontSize',18)
subplot(3,1,2), plot(c_Corr_1);
xlabel('Sample','FontWeight','bold','FontSize', 12)
title('Circular correlation','FontSize',18)
subplot(3,1,3), plot(c_Corr_FFT_1);
xlabel('Sample','FontWeight','bold','FontSize', 12)
title('Circular correlation using FFT','FontSize',16)
saveas(gcf,'figures/Crosscorr3tipi','fig');
saveas(gcf,'figures/Crosscorr3tipi','png');

