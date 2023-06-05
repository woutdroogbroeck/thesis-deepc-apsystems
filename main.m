%% main simulator script 
clear; close all; clc;

% load functions/classes from folders and subfolders
setup_path;

%--- create patient ---%
dt = 1;              % integration step [min]
Ts_controller = 5;   % controller sample time
model = "minimal";   % "minimal" or "uvapadova" simulation model
linear = false;      % simulations on linearized model or not
patient_id = 1;      % select patient -> "minimal" 1-3 | "uvapadova" 1-30
patient = load_patient(model, dt, patient_id, linear); 

%--- CGM and insulin pump ---%
sensor_id = 2;
seed = 0;

% patient.add_sensor(sensor_id,seed,Ts_controller); % add noise to measurements
% 
% pump_id = 1;
% patient.add_pump(pump_id); % add insulin pump for discrete insulin delivery

%--- simulation time ---%
begin_time = 6;
end_time = 18;
t = begin_time*60:dt:end_time*60-dt;
t_h = t/60;
simulation_steps = length(t);

%--- meals ---%
meals = zeros(1,simulation_steps);
% idx_meal1 = find(t_h==8);
% idx_meal2 = find(t_h==13);
% idx_meal3 = find(t_h==18);
% meals(idx_meal1) = 60;
% meals(idx_meal2) = 70;
% meals(idx_meal3) = 80;

% meals = fisher_disturbance(70,dt,meals,10); disturbance for minimal model

meal_announcement = repmat(struct('size',[],'minutes',[]),1,simulation_steps); 
meal1.size = 40;
meal1.minutes= 30;
% meal_announcement(idx_meal1-6) = meal1;

%--- controllers ---%
ref = 0; % set reference for tracking (controller tuning)
main_controllers;

%--- simulation ---%
simulate = @(patient, control_policy) simulation(patient, simulation_steps, control_policy, ...
    Ts_controller, meals, meal_announcement, false, false);

result_basal = simulate(patient, basal_policy);
result_deepc = simulate(patient, deepc_policy);

%%

LQR_COLOR = '#0072BD';
MPC_COLOR = '#D95319';
D_COLOR = 	'#77AC30';
DEEPC_COLOR = "#EDB120";

lw = 2;
fontsize_axis = 16;
fontsize_legend = 12;
fontsize_title = 16;

width = 800;
height = 300;

f = figure;
subplot(121)
plot(t_h,result_basal.G, 'Color',D_COLOR, LineWidth=lw)

hold on
plot(t_h,result_deepc.G, 'Color',DEEPC_COLOR,LineWidth=lw);
xlabel('Time [h]','Interpreter','latex','FontSize',fontsize_axis)
ylabel('$y$ [mg/dL]','Interpreter','latex','FontSize',fontsize_axis)
legend('basal','MPC$_5$','DeePC$_5$','DeePC$_1$', 'Interpreter', 'latex','FontSize',fontsize_legend);
title(sprintf('Glucose Concentration'), 'Interpreter', 'latex','FontSize',fontsize_title);
axis square
xlim([t_h(1) t_h(end)])


subplot(122)
hold on
stairs(t_h,1e3*result_deepc.u_p,'Color',DEEPC_COLOR,LineWidth=lw)
xlabel('Time [h]','Interpreter','latex','FontSize',fontsize_axis)
ylabel('$u$ [mU/min]','Interpreter','latex','FontSize',fontsize_axis)
legend('$u_{MPC}$','$u_{DeePC_5}$','$u_{DeePC_1}$','$u_{MPCI}$', 'Interpreter', 'latex','FontSize',fontsize_legend);
title(sprintf('Insulin Infusion Rate'), 'Interpreter', 'latex','FontSize',fontsize_title);
axis square
xlim([t_h(1) t_h(end)])
ylim([0 110])
f.Color = "w";
f.Position = [300 300 width height];
