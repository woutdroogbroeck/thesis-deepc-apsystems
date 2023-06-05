function results = simulation(patient, simulation_steps, control_policy, Ts_controller, meal_disturbance, meal_announcement, full_state, mean_G)
%SIMULATION function that performs meal simulation
%   patient            : patient object
%   simulation_steps   : total number of simulation steps
%   control_policy     : function that takes as argument the current state/output and the current meal announcement 
%   Ts_controller      : controller sample time, rate at which input is updated
%   meal_disturbance   : array with meals in grams carbohydrates
%   meal_announcements : array with meals (+- estimation error) a certain time beforehand, all zeros if no meal announcement 
%   full_state         : boolean to indicate of full state is used for controller
patient.set_initial_state;
t = zeros(1,simulation_steps);
x = zeros(size(patient.x,1),simulation_steps);
u = zeros(1,simulation_steps);
G = zeros(1,simulation_steps);
Gm = zeros(1,simulation_steps);
Gmean = zeros(1,simulation_steps);
% if meal_announcement == 0
%     meal_announcement = zeros(1,simulation_steps);
% end
if isa(control_policy,'function_handle')
    if full_state
        for k = 1:simulation_steps
            x(:,k) = patient.x;
            t(k) = patient.t;
            G(k) = patient.G;
            Gm(k) = patient.Gm;
            if mod(t(k),Ts_controller) == 0
                u(k) = control_policy(x(:,k), meal_announcement(k));
            else
                u(k) = u(k-1);
            end
            patient.update_state(u(k),meal_disturbance(k));
        end
    else
        if ~ isempty(patient.sensor) 
            for k = 1:simulation_steps
                x(:,k) = patient.x;
                t(k) = patient.t;
                G(k) = patient.G;
                Gm(k) = patient.Gm;
                Gmean(k) = patient.Gmean;
                if mod(t(k),Ts_controller) == 0
                    if mean_G 
                        u(k) = control_policy(Gmean(k), meal_announcement(k));
                    else
                        u(k) = control_policy(Gm(k), meal_announcement(k));
                    end                    
                else
                    u(k) = u(k-1);
                end
                patient.update_state(u(k),meal_disturbance(k));            
            end 
        else
            for k = 1:simulation_steps
                x(:,k) = patient.x;
                t(k) = patient.t;
                G(k) = patient.G;
                if mod(t(k),Ts_controller) == 0
                    u(k) = control_policy(G(k), meal_announcement(k));
                else
                    u(k) = u(k-1);
                end
                patient.update_state(u(k),meal_disturbance(k));            
            end 
        end
    
    
    end
else
    for k = 1:simulation_steps
        x(:,k) = patient.x;
        t(k) = patient.t;
        G(k) = patient.G;
        u(k) = control_policy(k);
        patient.update_state(u(k),meal_disturbance(k));            
    end 
end
results.t = t;
results.x = x;
results.u = u;
results.G = G;
results.Gm = Gm;
results.Gmean = Gmean;
results.u_p = patient.u(1,:);
results.d = patient.u(2,:);
[results.lbgi, results.hbgi, results.li, results.hi] = evaluate_bg(G);
end

