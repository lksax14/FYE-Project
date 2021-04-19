function [v_naught, vMax, Km] = M2_Algorithm_014_24(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132
% Program Description
% The purpose of this function is to determine the initial velocities of
% each concentration of the different substrates. It is also to return the
% paramaters of vMax and Km for each enzyme.
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
% vMax - The array that will be returned with each vMax for each enzyme
% (m/s)
%
% Km - The array that will be returned with each Km for each enzyme ([S] (uM))
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
data = readmatrix('Data_nextGen_KEtesting_allresults.csv');
data_x = data(3:end, 1); %Time data (seconds)
v_naught = zeros([5, 10]); %Pre-allocates the array of vnaughts (5x10)
Km = zeros([1, 5]); %Pre-allocates the array of Km values (1x5)
Vrow = 1; %V_naught array row counter
Vcol = 1; %V_naught array column counter
col = 2; %Data column index counter
Sub_Concentrations = data(1, 2:11); %Creates array of the substrate concentrations for use ([S] (uM))

%%____________________
%% CALCULATIONS
%All of these statements will fill an array of vnaughts respective to their
%concentrations. A simple moving average (SMA) is used to try and smooth
%the noisy data as much as possible. An average of the test data and the
%replicate test data is then taken and the initial slope of that data
%(calculated within the first 20 seconds) is found. That value is placed 
%into the vnaught array. This is repeated for each concentration to add up
%to a total of 50 vnaughts.
for ct = 1:50
    j = 1; %SMA index counter
    data_y = data(3:end, col); %Test data ([P] (uM))
    data_yrep = data(3:end, col + 10); %Replicate test data ([P] (uM))
    SMA = zeros([1, 7482]); %Pre-allocates an array of zeros for use
    SMAyrep = zeros([1, 7482]); %Pre-allocates an array of zeros for use
    col = col + 1;
    %These statements smooth both the original test data and the replicate
    %test data
    for i = 3:numel(data_y)
        for k = 1:-1:0
            SMA(j) = SMA(j) + (data_y(i - k));
            SMAyrep(j) = SMAyrep(j) + (data_yrep(i - k));
        end
        SMA(j) = SMA(j) / 2;
        SMAyrep(j) = SMAyrep(j) / 2;
        j = j + 1;
    end
    data_SMAave = (SMA(1:end) + SMAyrep(1:end)) / 2; %Takes the average of the original and replicate test data in their smoothed form
    coeffs = polyfit(data_x(2:20), data_SMAave(1:19), 1); %Finds the slope of that averaged data within the first 20 seconds
    v_naught(Vrow, Vcol) = coeffs(1); %Places that value into the array as needed
    Vcol = Vcol + 1;
    if Vcol == 11
        Vcol = 1;
        Vrow = Vrow + 1;
    end
    %If the column index reaches the highest substrate concentration, it
    %will jump over the replicate test data to get to the next enzyme that
    %needs to be analyzed
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

vMax = max(v_naught, [], 2); %Finds vMax for each enzyme

%These statements use the Lineweaver-Burk method (double reciprocals) to
%try and find Km for each enzyme
for ct = 1:5
    V_recip = 1 ./ v_naught(ct, :);
    S_recip = 1 ./ Sub_Concentrations;
    m = polyfit(S_recip, V_recip, 1);
    Km(ct) = m(1) * vMax(ct);
end

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
