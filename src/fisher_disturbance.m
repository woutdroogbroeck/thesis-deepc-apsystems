function meal_disturbance = fisher_disturbance(meal_size, Ts, meal_disturbance, k_start)
    % compute fisher disturbance of form d = A*exp(-b t)
    eating_time = 120; % 120 minutes for meal to 
    b = 0.05; 

    A = meal_size*b; 
    t = 0:Ts:eating_time;
    d = A*exp(-b.*t);
    meal_disturbance(k_start:k_start+length(d)-1) = d;

end