%% script to define control policies

% create patient object with controller sample time, for the correct
% discretization of the system matrices
ctrl_patient = load_patient(model, Ts_controller, patient_id, linear);
sys = ctrl_patient.sys_ctrl_meal; % linear system matrices
basal = patient.basal; % basal insulin flow for equilibrium

%% -- basal -- %%

basal_policy = @(x, d) basal; 

%% -- LQR -- %% 

% kalman filter parameters
Qk = 1e-9*eye(size(sys.A));
% Qk(1,1) = 10;
Rk = 1e-9;
P0 = 1e-9*eye(size(sys.A));
y0 = patient.G0;

Q = 100;
R = 1;

u_min = 0;
u_max = 0.1;

lqr_controller = LQR(sys, Q, R, ref, u_min, u_max, basal, Qk, Rk, P0);
lqr_policy = @(y, d) lqr_controller.solve(y-y0, d);

%% -- MPC -- %%

% kalman filter parameters
Qk = 1e-9*eye(size(sys.A));
Rk = 1e-9;
P0 = 1e-9*eye(size(sys.A));
y0 = patient.G0;

N = 20;
Q = 10;
R = 10^5;
Qt = 0*Q;

u_min = 0;
u_max = 0.1;
y_min = -100;

mpc_controller = MPC(N, Q, Qt, R, sys, ref, u_min, u_max, y_min, basal, Qk, Rk, P0);
mpc_policy = @(y, d) mpc_controller.solve(y-y0, d);

%% -- DeePC -- %%

deepc_data_patient = load_patient(model, dt, patient_id, linear); % same patient but set different seed for CGM
% deepc_data_patient.add_sensor(sensor_id,5,Ts_controller);
% deepc_data_patient.add_pump(pump_id); 

Tini = 6;
Tf = 24;
n = size(sys.A,1);
Td = max(2*(Tini+Tf+n)-1, 3*(Tini+Tf)-1);

u_offline = 0.04*prbs(6,Td);
y_offline = offline_data(deepc_data_patient, u_offline, Ts_controller);

lambda_y = 10^5;
lambda_g = 00;
lambda_proj = 10^4;
Q = 20;
R = 10^4;
dR = 0;
deepc_ref = patient.G0 + ref;
u_min = 0;
u_max = 0.1;
y_min = 0;
basal_meal_response = zeros(1,200); % update with estimate of meal response 

deepc_controller = DeePC(Tini, Tf, lambda_y, lambda_g, lambda_proj, ...
    Q, R, dR, deepc_ref, u_min, u_max, y_min, u_offline, y_offline, basal, patient.G0, Ts_controller,basal_meal_response,false);
deepc_policy = @(y, meal) deepc_controller.solve(y, meal);

lambda_y = 10^6;
R = 10^5;
lambda_proj = 00;
lambda_g = 10^3;

deepci_controller = DeePCI(Tini, Tf, lambda_y, lambda_g, lambda_proj, ...
    Q, R, dR, deepc_ref, u_min, u_max, y_min, u_offline, y_offline, basal, patient.G0, Ts_controller,basal_meal_response,false);
deepci_policy = @(y, meal) deepci_controller.solve(y, meal);


