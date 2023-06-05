classdef DeePC < handle
    % a class that implements a Data-Enabled Predictive Control algorithm for a
    % single input/single output system (glucose-insulin)
    properties
        Tini
        Tf
        lambda_y
        lambda_g
        lambda_proj
        Q
        R
        dR
        ref
        u_ref
        EAT_RATE = 5
        U_MAX = 0.1
        Y_MIN = 70
        U_MIN
        Ts
    end
    properties
        Up
        Yp
        Uf
        Yf
        H
        N
        uini
        yini
        problem
        Hpinv  
        u1 = 0
        d
        basal_yd
        steps_to_meal = 0;
    end
    
    methods

        function obj = DeePC(Tini, Tf, lambda_y, lambda_g, lambda_proj, Q, R, dR, ref, u_min, u_max, y_min, u_offline, y_offline, u0, y0, Ts, basal_yd, page, sigma)
            obj.Tini = Tini;
            obj.Tf = Tf;
            obj.lambda_y = lambda_y;
            obj.lambda_g = lambda_g;
            obj.lambda_proj = lambda_proj;
            obj.Ts = Ts;
            obj.Q = Q;
            obj.R = R;  
            obj.dR = dR;
            obj.ref = ref;
            obj.u_ref = u0;
            obj.d = zeros(Tf,1);
            obj.basal_yd = basal_yd;

            % create hankel/page matrices
            if page
                % sigma for SVD cutoff
                [obj.Up, obj.Uf] = page_matrix(u_offline,Tini,Tf,sigma);
                [obj.Yp, obj.Yf] = page_matrix(y_offline,Tini,Tf,sigma);
            else
                [obj.Up, obj.Uf] = hankel_matrix(u_offline,Tini,Tf);
                [obj.Yp, obj.Yf] = hankel_matrix(y_offline,Tini,Tf);
            end            
            obj.H = [obj.Up; obj.Yp; obj.Uf; obj.Yf];
            obj.Hpinv = pinv(obj.H);
            
            obj.N = size(obj.H,2); % size of vector g

            obj.uini = ones(Tini,1)*u0;
            obj.yini = ones(Tini,1)*y0;

            obj.U_MIN = u_min;
            obj.U_MAX = u_max;
            obj.Y_MIN = y_min;
            obj.problem = obj.build_problem;

        end
               
        function problem = build_problem(obj)
            % setup the optimization problem
            yalmip('clear')
            
            % initialize optimization variables
            u = sdpvar(obj.Tf,1);
            y = sdpvar(obj.Tf,1);
            g = sdpvar(obj.N,1);

            % past input/output measurements
            u_ini = sdpvar(obj.Tini,1);  
            y_ini = sdpvar(obj.Tini,1); 
            
            % slack variable for initial condition
            y_slack = sdpvar(obj.Tini,1);

            % slack variable for output constraint
            y_constr_slack = sdpvar(1,1);
            
            % predicted output disturbance
            y_d = sdpvar(obj.Tf,1);

            % system "dynamics" constraints
            b = [u_ini; y_ini+y_slack; u; y-y_d];
            constraints = [obj.H*g == b];
            
            % initialize objective function
            objective = 0; 
            for k = 1:obj.Tf
                % penalty on output and input
                objective = objective + obj.Q*(y(k)-obj.ref)^2 + obj.R*(u(k)-obj.u_ref)^2;
                if k>=2
                    % penalty on input rate of change
                    objective = objective + obj.dR*(u(k)-u(k-1))^2;
                end
                % input and output constraints
                constraints = [constraints, obj.U_MIN <= u(k)<= obj.U_MAX];
                constraints = [constraints, obj.Y_MIN  <= y(k) + y_constr_slack ];
            end

            % high penalty on output constraint slack variable, its only
            % function is to avoid feasibility issues
            objective = objective + 10^9*y_constr_slack^2;
            
            % reference vector to motivate [uini; yini; u; y] == [u_ref ; y_ref; u_ref; u_ref]            
            stacked_ref = [ones(obj.Tini,1)*obj.u_ref; ones(obj.Tini,1)*obj.ref; ones(obj.Tf,1)*obj.u_ref; ones(obj.Tf,1)*obj.ref];
            % value of g if system perfectly at reference
            gr = obj.Hpinv*stacked_ref; 
            
            % penalty on deviation from "ideal" g
            if obj.lambda_proj > eps
                objective = objective + obj.lambda_proj*norm(gr-g,2);
            end
            
            % penalty on output slack variable
            if obj.lambda_y > eps
                objective = objective + obj.lambda_y*norm(y_slack,2);
            else 
                % y_slack = 0 if lambda_y is not specified
                constraints = [constraints, norm(y_slack,2) <= eps];
            end
            
            % g regularization
            if obj.lambda_g > eps
                objective = objective + obj.lambda_g*norm(g,2);
            end
        
            ops = sdpsettings('solver','mosek');
            ops.verbose = 0;
            
            % build the problem: input last Tini u and y measurements, 
            % predicted output disturbance and current y. result = g*
            problem = optimizer(constraints,objective,ops,{u_ini,y_ini, y_d, y(1)},g);

        end

        function [u, obj] = solve(obj, currenty, meal)
            % calculates the optimal control action
            % args
            %    currenty: Current output of the system
            %    meal: struct with meal information: size (in gram) and 
            %          steps until meal is consumed
            % returns
            %    u: The calculated control action

            % if meal announcement, predict output disturbance
                
            % minimal model meal disturbance estimate
%             y_d = zeros(obj.Tf,1);
%             if meal.size > 0
%                 fisher = fisher_disturbance(meal.size, obj.Ts, y_d, meal.steps);
%                 obj.d = fisher*1000/120;
%                 for i = 2:length(obj.d)
%                     moving_avg(i) = obj.d(i)+obj.d(i-1)/2;
%                 end
%                 obj.d = moving_avg';
%             end
%             y_d(1:obj.Tf) = obj.d(1:obj.Tf);

            y_d = zeros(obj.Tf,1);
            if meal.size > 0
                steps = round(meal.minutes/obj.Ts);
                obj.steps_to_meal = steps+round(20/obj.Ts); % 20 minutes longer
                obj.d(steps:steps+obj.Tf-1) = .5*meal.size*obj.basal_yd(1:obj.Tf)/70; 
            end
            y_d(1:obj.Tf) = obj.d(1:obj.Tf);
            
            % compute g*
            g_optimal = obj.problem(obj.uini, obj.yini, y_d, currenty);
            % compute u* = Uf g*
            u_optimal = obj.Uf*g_optimal;
            % use the first value of the sequence
            u = u_optimal(1);
            
            % update past Tini measurements
            obj.uini = [obj.uini; u];            
            obj.yini = [obj.yini; currenty];
            obj.uini(1) = [];
            obj.yini(1) = [];
            
            % update upcoming output disturbance
            obj.d = [obj.d; 0];
            obj.d(1) = [];

        end
        
    end
end