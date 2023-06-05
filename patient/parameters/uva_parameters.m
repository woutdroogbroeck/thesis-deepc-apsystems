clear; close all; clc;

maximal_patients = readtable('parameters/patient.csv');

adult_idx = find(contains(maximal_patients.Name, 'adult'));
adults = maximal_patients(adult_idx,:);
save('adults.mat','adults')

child_idx = find(contains(maximal_patients.Name, 'child'));
children = maximal_patients(child_idx,:);
save('chlidren.mat','children')

adolescent_idx = find(contains(maximal_patients.Name, 'adolescent'));
adolescents = maximal_patients(adolescent_idx,:);
save('adolescents.mat','adolescents')

save('maximal_patients.mat', 'maximal_patients');



% 
% load('patient/parameters/adolescents.mat')
% load('patient/parameters/adults.mat')
% load('patient/parameters/children.mat')

% param = adolescents;
% 
% for i = 1:10
%     all_param(i,:) = param(i,3:end).Variables;
% end
% 
% for i = 1:length(all_param)
%     mean_param(i) = mean(all_param(:,i));
% end
% 
% save('mean_adolescents','mean_param');

% load('patient/parameters/maximal_patients.mat')
% 
% headers = maximal_patients.Properties.VariableNames(3:end);  % Get the headers
% 
% load('mean_adolescents.mat')
% baba = array2table(mean_param, 'VariableNames', headers);
% maximal_patients(31,3:end) = baba;
% maximal_patients(31,1) = {'mean_adolescents'};
% 
% load('mean_adults.mat')
% baba = array2table(mean_param, 'VariableNames', headers);
% maximal_patients(32,3:end) = baba;
% maximal_patients(32,1) = {'mean_adults'};
% 
% load('mean_children.mat')
% baba = array2table(mean_param, 'VariableNames', headers);
% maximal_patients(33,3:end) = baba;
% maximal_patients(33,1) = {'mean_children'};
% 
% 
% save('maximal_patienst','maximal_patients');