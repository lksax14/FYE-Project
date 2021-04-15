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
%Vi=Vmax−Km×Vi/[S] %this is the equation we need to use.
% "Vmax is the highest vnaught value of the ones you calculated using different concentrations"
% Vmax might be the last value in each row?
%% ____________________
%% INITIALIZATION
data = readmatrix('Data_nextGen_KEtesting_allresults.csv'); %Imports the file to be used
data_x = data(3:end, 1); %Time data (seconds)
vnaught = zeros([5, 10]); %Pre-allocates the array of vnaughts
Vrow = 1; %Vnaught row counter
Vcol = 1; %Vnaught column counter
col = 2; %Data column index counter
vMax =  zeros([5, 1]); % Pre-allocates an array of vmaxes


%%____________________
%% CALCULATIONS
%These statements will fill an array of vnaughts respective to their
%concentrations. A simple moving average (SMA) is used to try and smooth
%the noisy data as much as possible. An average of the test data and the
%replicate test data is then taken and the initial slope of that line
%(calculated within the first 20 seconds) and that value is placed into the
%vnaught array. This is repeated for each concentration to add up to a
%total of 50 vnaughts.
for ct = 1:50
    j = 1; %SMA index counter
    data_y = data(3:end, col); %Test data ([P] (uM))
    data_yrep = data(3:end, col + 10); %Replicate test data ([P] (uM))
    SMA = zeros([1, 7482]); %Pre-allocates an array of zeros for use
    SMAyrep = zeros([1, 7482]); %Pre-allocates an array of zeros for use
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
    data_SMAave = (SMA(1:end) + SMAyrep(1:end)) / 2;
    coeffs = polyfit(data_x(2:20), data_SMAave(1:19), 1);
    vnaught(Vrow, Vcol) = coeffs(1);
    Vcol = Vcol + 1;
    if Vcol == 11
        Vcol = 1;
        Vrow = Vrow + 1;
    end
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

vMax = max(v_naught, [], 2);


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
