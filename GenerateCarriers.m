function [ carrier, phase_Out, carrier_Cos, carrier_Sin, carrier_Fraction_Out ] = GenerateCarriers( sampling_Rate, carrier_Frequency, phase_In, n_Samples, carrier_Fraction_In)
%GenerateCarriers - Generates the carrier used to perform the
%down-conversion of a GNSS signal from the Intermediate Frequency to the
%Base Band
%
% This function generates a block of duration 1ms of a carrier

% Syntax:  [carrier, phase_Out, carrier_Cos, carrier_Sin, carrier_Fraction_Out] = GenerateCarriers(sampling_Rate, carrier_Frequency, phase_In, n_Samples, carrier_FractionIn)
%
% Inputs:
%    sampling_Rate - Sampling frequency of the RF Front-End (Hz)
%                    carrier_Frequency - Intermediate frequency (Hz)
%    phase_In - Phase of the carrier for the first sample (rad)
%    n_Samples - Number of samples to be generated in the block
%                (Fractional)
%    carrier_Fraction_In - Fractional part of the number of samples from
%                          the previous block
%
% Outputs:
%    carrier - Generated sequence of complex samples of the carrier
%    phase_Out - Phase of the carrier for the (nSamples+1)-th sample (rad)
%    carrier_Cos - Generated sequence of real samples of the cosine part of
%                  the carrier
%    carrier_Sin - Generated sequence of real samples of the sine part of
%                  the carrier
%    carrier_Fraction_Out - Fractional part of the number of samples for
%                           the next block
%
% Author: Cilia Martina, Arena Remo
% February 2017

n_Samples_Eff=floor(n_Samples+carrier_Fraction_In); % Integer part of the number of samples

n = (0:1:(n_Samples_Eff))./sampling_Rate; % Digital time instant
arg = 2*pi*carrier_Frequency*n + phase_In; % Argument of the sin and cos functions

carrier_Sin = sqrt(2) * sin(arg(1:end-1));
carrier_Cos = sqrt(2) * cos(arg(1:end-1));
carrier = carrier_Cos + 1i*carrier_Sin;

phase_Out = mod(arg(end), 2*pi); % Phase of the carrier for the (nSamples+1)-th sample
carrier_Fraction_Out=abs(n_Samples - n_Samples_Eff + carrier_Fraction_In); % Fractional part of the number of samples for the next block
end
