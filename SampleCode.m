function [ code_Out, chip_Fraction_Out ] = SampleCode( sampling_Rate, code_In, chip_Index_In, chip_Rate, chip_Fraction_In)
%SampleCode - Generate a local replica of the despreading code sampled at
%the sampling frequency of the RF Front-End.
%
% This function generates a block of duration 1ms of a PRN code
% Syntax:  [ code_Out, chip_Fraction_Out ] = SampleCode( sampling_Rate, code_In, chip_Index_In, chip_Rate, chip_Fraction_In)
%
% Inputs:
%    sampling_Rate - Sampling frequency of the RF Front-End (Hz)
%    code_In - Vector containing a period of a PRN code sampled at one sample per chip
%    chip_Index_In - Index of the chip to which the first generated sample
%                    belongs. Assumes values in [1,L] where L is the number
%                    of code chip
%    chip_Rate -  CHop rate of the despreading code
%    chip_Fraction_In - Fraction part of the number of samples from the
%                       previous block
%
% Outputs:
%    code_Out - Generated sequence of samples of the PRN code code_In
%    chip_Fraction_Out - Fractional part of the number of samples for the
%                        next block
%
% Author: Cilia Martina, Arena Remo
% February 2017


n_Samples_Chip = sampling_Rate / chip_Rate; % #samples per chip (Fractional)
code_In = circshift(code_In, chip_Index_In-1); % Vector of the PRN code starting from chip_Index_In

n_Chip_Code = 1023; % #chip of a PRN code
code_Out = [];

for k = 1:n_Chip_Code
    n_Samples_Eff = floor(n_Samples_Chip + chip_Fraction_In); % Actual #samples for the k-th chip
    chip_Fraction_In = abs(n_Samples_Chip + chip_Fraction_In - n_Samples_Eff); % Fraction of samples for the (k+1)-th chip
    chip_Up = code_In(k) + zeros(1,n_Samples_Eff); % Upsampled version of k-th chip
    code_Out = [code_Out chip_Up];
end

chip_Fraction_Out = chip_Fraction_In;

