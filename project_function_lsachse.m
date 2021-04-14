function [v_naught] = M2_Algorithm_014_24(data)
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
% data - The data that contains the information needed on the different
% substrates
%
% Output Arguments
% v_naught - The array that will be returned with each initial velocity for
% each concentration of substrate (m/s)
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

%% ____________________
%% INITIALIZATION
data = readmatrix('Data_nextGen_KEtesting_allresults.csv'); %Imports the file to be used
data_x = data(3:end, 1); %Time data (seconds)
vnaught = zeros([1, 50]); %Pre-allocates the array of vnaughts
z = 1; %Vnaught index counter
col = 2; %Data column index counter


%%____________________
%% CALCULATIONS
%These statements will fill an array of vnaughts respective to their
%concentrations. A simple moving average (SMA) is used to try and smooth
%the noisy data as much as possible. An average of the test data and the
%replicate test data is then taken and the initial slope of that line
%(calculated within the first 21 seconds) and that value is placed into the
%vnaught array. This is repeated for each concentration to add up to a
%total of 50 vnaughts.
for ct = 1:50
    j = 1; %SMA index counter
    data_y = data(3:end, col); %Test data ([P] (uM))
    data_yrep = data(3:end, col + 10); %Replicate test data ([P] (uM))
    SMA = zeros([1, 7482]); %Pre-allocates an array of zeros for use. This array represents the Simple Moving Average values
    SMAyrep = zeros([1, 7482]); %Pre-allocates an array of zeros for use. This array represents the Simple Moving Average values for the replicate test data
    col = col + 1;
    for i = 3:numel(data_y)
        for k = 1:-1:0
            SMA(j) = SMA(j) + (data_y(i - k));
            SMAyrep(j) = SMAyrep(j) + (data_yrep(i - k));
        end
        SMA(j) = SMA(j) / 2;
        SMAyrep(j) = SMAyrep(j) / 2;
        j = j + 1;
    end
    data_SMAave = (SMA(1:end) + SMAyrep(1:end)) / 2; %This is the average of the SMA and SMArep values to find the best fitting and smoothed line that represents the data accurately
    coeffs = polyfit(data_x(11:21), data_SMAave(1:11), 1);
    vnaught(z) = coeffs(1); %This is Vnaught, the initial slope of the line, and will be placed into the array that is returned
    z = z + 1;
    if col == 11
        col = 21;
    end
    if col == 31
        col = 41;
    end
    if col == 51
        col = 61;
    end
    if col == 71
        col = 81;
    end
end

v_naught = vnaught;
%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS


%% ____________________
%% RESULTS


%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.
