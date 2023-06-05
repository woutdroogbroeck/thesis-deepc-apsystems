% clear; close all;
% 
% 
% % Physiologic evaluation of factors controlling glucose tolerance in man: 
% % measurement of insulin sensitivity and beta-cell glucose sensitivity 
% % from the response to intravenous glucose. Bergman et al.
% 

% Define the data for each patient
data = [
    80 0 7 1.36e-2 3.41e-2 17.3e-6 0.22 .12*80 1.6*80;
    82 0 15 1.52e-2 3.13e-2 9.7e-6 0.22 .12*80 1.6*80;
    87 0 14 2.17e-2 2.92e-2 19.1e-6 0.11 .12*80 1.6*80;
];

% Define the column names
colNames = {'Gb', 'Xb', 'Ib','p1', 'p2', 'p3', 'n', 'Vi', 'Vg'};

% Create a table to store the data
minimal_patients = array2table(data, 'VariableNames', colNames);

% Save the table to a .mat file
save('minimal_patients.mat', 'minimal_patients');
