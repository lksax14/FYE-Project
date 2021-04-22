function [v_naught, vMax, Km] = M3_MM_PGOX50_014_24(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENGR 132
% Program Description
% The purpose of this function is to determine the initial velocities of
% each substrate of the different enzymes. It is also to return the
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
% each substrate of the enzymes (m/s)
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
data = readmatrix('Data_PGOX50_enzyme.csv');
data_x = data(5:end, 1); %Time data (seconds)
v_naught = zeros([1, 10]); %Pre-allocates the array of vnaughts (1x10)
Km = zeros([1, 1]); %Pre-allocates the array of Km values (1x1)
Vrow = 1; %V_naught array row counter
Vcol = 1; %V_naught array column counter
col = 2; %Data column index counter
Sub_Concentrations = data(2, 2:11); %Creates array of the substrate concentrations for use ([S] (uM))
data_y1 = data(5:end, 2);
data_y2 = data(5:end, 3);
data_y3 = data(5:end, 4);
data_y4 = data(5:end, 5);
data_y5 = data(5:end, 6);
data_y6 = data(5:end, 7);
data_y7 = data(5:end, 8);
data_y8 = data(5:end, 9);
data_y9 = data(5:end, 10);
data_y10 = data(5:end, 11);
v0line1 = v_naught(1) * data_x;
v0line2 = v_naught(2) * data_x;
v0line3 = v_naught(3) * data_x;
v0line4 = v_naught(4) * data_x;
v0line5 = v_naught(5) * data_x;
v0line6 = v_naught(6) * data_x;
v0line7 = v_naught(7) * data_x;
v0line8 = v_naught(8) * data_x;
v0line9 = v_naught(9) * data_x;
v0line10 = v_naught(10) * data_x;

%%____________________
%% CALCULATIONS
%All of these statements will fill an array of vnaughts respective to their
%concentrations. A simple moving average (SMA) is used to try and smooth
%the noisy data as much as possible. An average of the test data and the
%replicate test data is then taken and the initial slope of that data
%(calculated within the first 20 seconds) is found. That value is placed 
%into the vnaught array. This is repeated for each concentration for each 
%enzyme to add up to a total of 50 vnaughts.
for ct = 1:10
    j = 1; %SMA index counter
    data_y = data(3:end, col); %Test data ([P] (uM))
    SMA = zeros([1, 1222]); %Pre-allocates an array of zeros for use (original data)
    col = col + 1;
    %These statements smooth both the original test data and the replicate
    %test data
    for i = 3:numel(data_y)
        for k = 1:-1:0
            SMA(j) = SMA(j) + (data_y(i - k));
        end
        SMA(j) = SMA(j) / 2;
        j = j + 1;
    end
    coeffs = polyfit(data_x(2:20), SMA(1:19), 1); %Finds the slope of that averaged data within the first 20 seconds
    v_naught(Vrow, Vcol) = coeffs(1); %Places that value into the array as needed
    Vcol = Vcol + 1;
    if Vcol == 11
        Vcol = 1;
        Vrow = Vrow + 1;
    end
end

vMax = max(v_naught, [], 2); %Finds vMax for each enzyme

%These statements use the Lineweaver-Burk method (double reciprocals) to
%find the Km value for each enzyme
V_recip = 1 ./ v_naught(1, :);
S_recip = 1 ./ Sub_Concentrations;
m = polyfit(S_recip, V_recip, 1);
Km = m(1) * vMax;

michaelis_menten = vMax * Sub_Concentrations ./ (Km + Sub_Concentrations);

%SSE Calculation
SSE = sum((v_naught - michaelis_menten) .^ 2);

%% ____________________
%% FORMATTED TEXT/FIGURE DISPLAYS
figure(1);

subplot(2, 5, 1); %first row, first column subplot of 2x5 grid
plot(data_x, data_y1, 'r-');
grid on;
hold on; 
plot(data_x, v0line1, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM))");
title("PGO-X50 Values (3.75 uM)");
legend('Enzyme Data','Initial Velocity','Location','best');
ylim([0 10]);
grid on;
hold off 


subplot(2, 5, 2); %first row, second column subplot of 2x5 grid %%10
plot(data_x, data_y2, 'r-');
hold on;
plot(data_x, v0line2, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (7.5 uM)");
ylim([0 20]);
grid on;
hold off 


subplot(2, 5, 3); %first row, third column subplot of 2x5 grid %%20
plot(data_x, data_y3, 'r-');
hold on;
plot(data_x, v0line3, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (15 uM)");
ylim([1 40]);
grid on;
hold off 


subplot(2, 5, 4); %first row, fourth column subplot of 2x5 grid 
plot(data_x, data_y4, 'r-');
hold on;
plot(data_x, v0line4, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (30 uM)");
ylim([0 100]);
grid on;
hold off 


subplot(2, 5, 5); %first row, fifth column subplot of 2x5 grid
plot(data_x, data_y5, 'r-');
hold on;
plot(data_x, v0line5, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (65 uM)");
ylim([0 250]);
grid on;
hold off 


subplot(2, 5, 6); %second row, first column subplot of 2x5 grid
plot(data_x, data_y6, 'r-');
hold on;
plot(data_x, v0line6, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (125 uM)");
ylim([0 450]);
grid on;
hold off 


subplot(2, 5, 7); %second row, second column subplot of 2x5 grid
plot(data_x, data_y7, 'r-');
hold on;
plot(data_x, v0line7, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (250 uM)");
ylim([0 700]);
grid on;
hold off 


subplot(2, 5, 8); %second row, third column subplot of 2x5 grid
plot(data_x, data_y8, 'r-');
hold on;
plot(data_x, v0line8, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (500 uM)");
ylim([0 850]);
grid on;
hold off 


subplot(2, 5, 9); %second row, fourth column subplot of 2x5 grid
plot(data_x, data_y9, 'r-');
hold on;
plot(data_x, v0line9, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (1000 uM)");
ylim([0 1000]);
grid on;
hold off 


subplot(2, 5, 10); %second row, fifth column subplot of 2x5
plot(data_x, data_y10, 'r-');
hold on;
plot(data_x, v0line10, 'k-');
xlabel("Time (sec)");
ylabel("Concentration of Product (uM)");
title("PGO-X50 Values (2000 uM)");
ylim([0 1400]);
grid on;
hold off 

sgtitle("PGOX50 Enzyme Data and Initial Velocities")

figure(2)
plot(Sub_Concentrations, v_naught, 'k.');
hold on
plot(Sub_Concentrations, michaelis_menten, 'r-');
title("PGO-X50 Michaelis-Menten Plot");
xlabel("Substrate Concentration (uM)");
ylabel("Velocities (m/s)");
legend("Initial Velocities", "Michaelis Menten Model", "Location", "best");
grid on;
hold off


%% ____________________
%% RESULTS


%% ____________________
%% ACADEMIC INTEGRITY STATEMENT
% We have not used source code obtained from any other unauthorized
% source, either modified or unmodified. Neither have we provided
% access to my code to another. The program we are submitting
% is our own original work.
