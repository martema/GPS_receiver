function [CA_Code] = GoldCodeGenerator()
%GoldCodeGenerator - Generates the Gold sequences used by GPS
%
% Syntax:  [CA_Code] = GoldCodeGenerator();
%
% Parameters:
%    N - Number of shift register states 
%    p - Codes length
%    sv_Num - Number of GPS SVs
%
% Outputs:
%    CA_Code - p*sv_Num matrix that contains the code of length p for each
%    of the GPS satellites
%
% Author: Cilia Martina, Arena Remo
% February 2017

%% Parameters
N = 10;
p = 2^N-1;
sv_Num = 33;

% GPS SVs spreading codes 
phaseSelector = [2 6;3 7;4 8;5 9;1 9;2 10; 1 8;2 9;3 10;2 3;3 4;5 6;6 7;7 8;8 9;9 10;...
    1 4;2 5;3 6;4 7;5 8;6 9;1 3;4 6;5 7;6 8;7 9;8 10;1 6;2 7;3 8;4 9;5 10];

%% Initialization of the shift registers
G1_Shift_Reg = -ones(1,N); % Shift register for G1
G2_Shift_Reg = -ones(1,N); % Shift register for G2

%% C/A code generation
G1 = zeros(p,sv_Num);
G2 = zeros(p,sv_Num);

for ind = 1:p
    %Generation of the (ind)-th bit for all the satellites
    G1(ind,:) = G1_Shift_Reg(end)*ones(1,sv_Num); % Value of the last cell
    G2_Cells = G2_Shift_Reg(phaseSelector); % Value of the cells phaseSelector
    G2(ind,:) = prod(G2_Cells,2);
    
    % Shift register G1 update
    G1_Value = G1_Shift_Reg(3)*G1_Shift_Reg(10);
    G1_Shift_Reg = circshift(G1_Shift_Reg,1);
    G1_Shift_Reg(1) = G1_Value;
    
    % Shift register G2 update
    G2_Value = G2_Shift_Reg(2)*G2_Shift_Reg(3)*G2_Shift_Reg(6)*G2_Shift_Reg(8)*G2_Shift_Reg(9)*G2_Shift_Reg(10);
    G2_Shift_Reg = circshift(G2_Shift_Reg,1);
    G2_Shift_Reg(1) = G2_Value;
end

    CA_Code = G1.*G2;
end

