function [ c_Corr_FFT ] = CirCorrFFT( input_1, input_2 )
%CIRCORR - Returns the circular correlation function between two input
%sequences
%
% This function uses the FFT approach
% R[tau] = IFFT(FFT(input_1)*conj(FFT(input_2)))
%
% Syntax:  [c_Corr_FFT] = CirCorrFFT(input_1, input_2);
%
% Inputs:
%    input_1 - Input sequence 1, dimensions 1*N1
%    input_2 - Input sequence 2, dimensions 1*N2
%
% Outputs:
%    l_Corr - Linear correlation function between sequence input_1 and
%    sequence input_2, dimensions 1*(length of the longest sequence);
%
% Author: Cilia Martina


N1 = length(input_1); % Length of sequence 1
N2 = length(input_2); % Length of sequence 2

N = max(N1,N2); % Length of the longest sequence

X1 = fft(input_1, N); % FFT(input_1)
X2 = fft(input_2, N); % FFT(input_2)

% Normalized Correlation
c_Corr_FFT = ifft(X1.*conj(X2))/N;
end

