function [ l_Corr] = LinCorr( input_1, input_2 )
%LinCorr - Returns the linear correlation function between two input
%sequences
%
% Syntax:  [l_Corr] = LinCorr(input_1, input_2);
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

% Zero padding of the longest sequence
if (N1 < N2)
    long = [input_2 zeros(1,N1-1)];
    short = input_1;
    N = N1;
    N_Corr = N2;
else
    long = [input_1 zeros(1,N2-1)];
    short = input_2;
    N = N2;
    N_Corr = N1;
end

% Normalized Correlation
l_Corr = zeros(1,N_Corr);

for m = 1:N_Corr
    l_Corr(m) = short*(long(m:(m+N-1)))';
end

l_Corr = l_Corr/N_Corr;
