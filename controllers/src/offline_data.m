function y_offline = offline_data(patient, u_offline, Ts_controller)
patient.set_initial_state;
Td = length(u_offline);
u = zeros(1,Td*Ts_controller/patient.Ts);
y_offline = zeros(1,Td);
i = 1;
if isempty(patient.sensor)
    for k = 1:Td*Ts_controller/patient.Ts
        if mod(patient.t,Ts_controller) == 0
            u(k) = u_offline(i);
            y_offline(i) = patient.G; 
            i = i+1;
        else
            u(k) = u(k-1);
        end
        patient.update_state(u(k),0);
    end
else
    for k = 1:Td*Ts_controller/patient.Ts
        if mod(patient.t,Ts_controller) == 0
            u(k) = u_offline(i);
            y_offline(i) = patient.Gm; 
            i = i+1;
        else
            u(k) = u(k-1);
        end
        patient.update_state(u(k),0);
    end
end

% u_offline = patient.u(1,:);
patient.set_initial_state;
end