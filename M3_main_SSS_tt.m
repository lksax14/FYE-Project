function M3_main_SSS_tt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132
% Program Description
% The purpose of this function is to determine the initial velocities of
% each concentration of the different substrates.
%
% Function Call
% M2_Algorithm_014_24(data)
%
% Input Arguments
% N / A
%
% Output Arguments
% N / A
%
% Assignment Information
%   Assignment:     M02, Problem 1
%   Team member:    William Albright, albrighw@purdue.edu [repeat for each person]
%                   Leo Sakuma, lsakuma@purdue.edu
%                   Luke Sachse, lsachse@purdue.edu
%                   Alec Liu, ajliu@purdue.edu
%   Team ID:        014-24
%   Academic Integrity:
%     [] We worked with one or more peers but our collaboration
%        maintained academic integrity.
%     Peers we worked with: N/A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = readmatrix('Data_PGOX50_enzyme.csv'); %Imports the file to be used
[xDim, yDim] = size(data);
M2_Algorithm_014_24(data, xDim, yDim);