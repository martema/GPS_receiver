function [ c_Corr ] = CirCorr( input_1, input_2)
%CIRCORR - Returns the circular correlation function between two input
%sequences
%
% Syntax:  [l_Corr] = CirCorr(input_1, input_2);
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

% Circular insertion on the longest sequence
if (N1 < N2)
    long = [input_2 input_2(1:N1-1)];
    short = input_1;
    N = N1; % Short sequence length
    N_Corr = N2; % Long sequence length
    
else
    long = [input_1 input_1(1:N2-1)];
    short = input_2;
    N = N2; % Short sequence length
    N_Corr = N1; % Long sequence length
end

% Normalized Correlation
c_Corr = zeros(1,N_Corr);

for m = 1:N_Corr
    c_Corr(m) = short*(long(m:(m+N-1)))';
end

c_Corr = c_Corr/N_Corr;
