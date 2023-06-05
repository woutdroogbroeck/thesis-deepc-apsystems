classdef Patient < handle
    
    properties
        param
        Ts
        linear
        model
    end
    properties
        t
        x
        xb
        G
        Gm = -1
        Gmean = 0
        sys_ctrl
        sys_ctrl_meal
        Ad
        Bd
        previous_portion
        last_Qsto
        last_foodtaken
        meal_not_consumed        
        eating
        linear_system
        all_portions
        EAT_RATE
        minimal_state
        basal
        x0
        G0
        sensor = []
        pump = []
        u = []
    end
    methods
        function obj = Patient(params, Ts, model, linear)
            obj.param = params;
            obj.Ts = Ts;
            obj.model = model;
            obj.linear = linear;
            obj.EAT_RATE = 5;           
            if model == "minimal"
                [control_system ,linear_system] = linearize_minimal_model(params);
            elseif model == "uvapadova"
                [control_system, linear_system] = linearize_uvapadova(params);                
            end
            obj.linear_system = linear_system;
            sys = c2d(linear_system,Ts);
            sys_ctrl = ss(control_system.A, control_system.B(:,1), control_system.C, 0);
            obj.sys_ctrl = c2d(sys_ctrl,Ts);
            obj.sys_ctrl_meal = c2d(control_system,Ts);
            obj.Ad = sys.A;
            obj.Bd = sys.B;
            obj = obj.set_initial_state;
            obj.x0 = obj.x;
            obj.G0 = obj.G;
        end

        function obj = set_initial_state(obj)
            obj.t = 0;
            if obj.model == "minimal"
                if obj.linear
                    obj.x = zeros(3,1);
                    obj.G = obj.x(1);
                    obj.xb = obj.x;
                    obj.previous_portion = 0;
                    obj.meal_not_consumed = 0;
                    obj.all_portions = [0];
                    obj.minimal_state = obj.x;
                    obj.basal = 0;
                else
                    obj.x = obj.param(1,1:3).Variables';
                    obj.G = obj.x(1);
                    obj.xb = obj.x;
                    obj.previous_portion = 0;
                    obj.meal_not_consumed = 0;
                    obj.all_portions = [0];
                    obj.minimal_state = obj.x;
                    obj.basal = 1e-3*obj.param.n*obj.param.Ib*obj.param.Vi;
                end                
            elseif obj.model == "uvapadova"
                if obj.linear
                    obj.x = zeros(13,1);
                    obj.G = obj.x(13)/obj.param.Vg;
                    obj.previous_portion = 0;
                    obj.meal_not_consumed = 0;
                    obj.xb = obj.x;
                    obj.all_portions = [0];
                    obj.minimal_state = [obj.x(13)/obj.param.Vg;
                                         obj.x(7);
                                         obj.x(6)];
                    obj.basal = 0;
                else
                    obj.x = obj.param(1,3:15).Variables';
                    obj.G = obj.x(13)/obj.param.Vg;
                    obj.last_Qsto = obj.x(1) + obj.x(2);
                    obj.previous_portion = 0;
                    obj.eating = false;
                    obj.meal_not_consumed = 0;
                    obj.xb = obj.x;
                    obj.all_portions = [0];
                    obj.minimal_state = [obj.x(13)/obj.param.Vg;
                                         obj.x(7);
                                         obj.x(6)];
                    obj.basal = obj.param.u2ss * obj.param.BW/6000;
                end   
            end
            if ~isempty(obj.sensor)
                obj.sensor.reset;
                obj.Gm = obj.sensor.measure(obj);   
                obj.Gmean = mean(obj.sensor.lastN_measurements);
            end
            obj.u = [];
        end

        function obj = add_sensor(obj, sensor_id, seed, N)
            sensors = load('sensor/parameters/sensors.mat').sensors;
            sensor_param = sensors(sensor_id,:);
            sensor_param.PACF = 0.7;
            obj.sensor = CGM(sensor_param, 1, seed, N);
            obj.sensor.reset;
            obj.Gm = obj.sensor.measure(obj); 
            obj.Gmean = mean(obj.sensor.lastN_measurements);
        end

        function obj = add_pump(obj, pump_id)
            pumps = load('pump/parameters/pumps.mat').pumps;
            pump_param = pumps(pump_id,:);
            pump_param.inc_basal = 0.005;
            obj.pump = InsulinPump(pump_param, obj.basal);
        end
        
        function obj = update_state(obj, insulin, meal)
            current_portion = obj.set_portion(meal);
            obj.all_portions = [obj.all_portions, current_portion];
            if ~ isempty(obj.pump)
                insulin = obj.pump.discrete_insulin(insulin);
            end
            current_u = [insulin; current_portion/obj.Ts];
            obj.u = [obj.u, current_u];
            if obj.model == "uvapadova"
                
                if obj.linear
                    obj.x = obj.Ad*obj.x + obj.Bd*[insulin; current_portion/obj.Ts];   
                    obj.G = obj.x(13)/obj.param.Vg; 
                    obj.minimal_state = [obj.x(13)/obj.param.Vg;
                                         obj.x(7);
                                         obj.x(6)];
                               
                else
                    if (current_portion>0 && obj.previous_portion<=0)
                        obj.last_Qsto = obj.x(1) + obj.x(2);
                        obj.last_foodtaken = 0;
                        obj.eating = true;
                    end
        
                    if obj.eating
                        obj.last_foodtaken = obj.last_foodtaken + current_portion;
                    end
        
                    if current_portion <=0 && obj.previous_portion>0
                        obj.eating = false;
                    end
        
                    f = @(x) uvapadova_fast(x,obj.param,insulin,current_portion/obj.Ts,obj.last_Qsto,obj.last_foodtaken);
                    obj.x = runge_kutta4(obj.x, f, obj.Ts);
                    obj.G = obj.x(13)/obj.param.Vg;
                    obj.minimal_state = [obj.x(13)/obj.param.Vg;
                                         obj.x(7);
                                         obj.x(6)];
                    
                end
                

            elseif obj.model == "minimal"
                if obj.linear
                     obj.x = obj.Ad*obj.x + obj.Bd*[insulin; meal/obj.Ts];    
                     obj.G = obj.x(1);
                     obj.minimal_state = obj.x;
                else
                    f = @(x) minimal_model(x, obj.param, insulin, meal/obj.Ts);
                    obj.x = runge_kutta4(obj.x, f, obj.Ts);           
                    obj.G = obj.x(1);   
                    obj.minimal_state = obj.x;
                end                
                
            end
            obj.t = obj.t + obj.Ts;
            if ~ isempty(obj.sensor)
                obj.Gm = obj.sensor.measure(obj);
                obj.Gmean = mean(obj.sensor.lastN_measurements);
            end
            obj.previous_portion = current_portion;
        end

        function [current_portion, obj] = set_portion(obj, meal)
            obj.meal_not_consumed = obj.meal_not_consumed + meal;
            if obj.meal_not_consumed>0
                current_portion = min(obj.EAT_RATE*obj.Ts, obj.meal_not_consumed);
                obj.meal_not_consumed = obj.meal_not_consumed-current_portion;
                obj.meal_not_consumed = max(0,obj.meal_not_consumed);
            else
                current_portion = 0;
            end
        end

    end
end