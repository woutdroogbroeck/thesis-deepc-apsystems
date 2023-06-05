classdef DeePCI < handle
    % a class that implements a Data-Enabled Predictive Control algorithm with integral action for a
    % single input/single output system (glucose-insulin) similar to MPCI
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
        U_MAX = .1
        Y_MIN
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
        Td
        uini
        yini
        problem
        Hpinv    
        u1
        y1
        d
        basal_yd
    end
    
    methods

        function obj = DeePCI(Tini, Tf, lambda_y, lambda_g, lambda_proj, Q, R, dR, ref, u_min, u_max, y_min, u_offline, y_offline, u0, y0, Ts, basal_yd, page, sigma)
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
            obj.basal_yd = basal_yd;

            for k = 2:length(u_offline)
                du_offline(k) = u_offline(k)-u_offline(k-1); 
                dy_offline(k) = y_offline(k)-y_offline(k-1);
            end
            if page
                [obj.Up, obj.Uf] = page_matrix(du_offline,Tini,Tf,sigma);
                [obj.Yp, obj.Yf] = page_matrix(dy_offline,Tini,Tf,sigma);
            else
                [obj.Up, obj.Uf] = hankel_matrix(du_offline,Tini,Tf);
                [obj.Yp, obj.Yf] = hankel_matrix(dy_offline,Tini,Tf);
            end
            obj.H = [obj.Up; obj.Yp; obj.Uf; obj.Yf];
            obj.Hpinv = pinv(obj.H);

            obj.N = size(obj.H,2);

            obj.uini = zeros(Tini,1);
            obj.yini = zeros(Tini,1);
            obj.u1 = u0;
            obj.y1 = y0;
            
            obj.d = zeros(Tf,1);
            obj.U_MIN = u_min;
            obj.U_MAX = u_max;
            obj.Y_MIN = y_min;
            obj.problem = obj.build_problem;

        end
               
        function problem = build_problem(obj)
            yalmip('clear')

            du = sdpvar(obj.Tf,1);
            dy = sdpvar(obj.Tf,1);
            y = sdpvar(obj.Tf,1);
            u = sdpvar(obj.Tf,1);
            g = sdpvar(obj.N,1);

            y_slack = sdpvar(obj.Tini,1);

            y_constr_slack = sdpvar(1,1);
            
            du_ini = sdpvar(obj.Tini,1);  
            dy_ini = sdpvar(obj.Tini,1); 

            dy_m = sdpvar(obj.Tf,1);
            
            b = [du_ini; dy_ini+y_slack; du; dy-dy_m];

            constraints = [obj.H*g == b];
            objective = 0; 
            
            A = tril(ones(obj.Tf));
            sdpvar y0 u0;

            y0_vec = y0*ones(obj.Tf,1);
            u0_vec = u0*ones(obj.Tf,1);
            
            constraints = [constraints, y == A*dy + y0_vec];
            constraints = [constraints, u == A*du + u0_vec];
            for k = 1:obj.Tf
                objective = objective + obj.Q*(y(k)-obj.ref)^2 + obj.R*du(k)^2;
                constraints = [constraints, obj.U_MIN <= u(k)<= obj.U_MAX];
                constraints = [constraints, obj.Y_MIN  <= y(k) + y_constr_slack ];
            end

            % high penalty on output constraint slack variable, its only
            % function is to avoid feasibility issues
            objective = objective + 10^9*y_constr_slack^2;         

            if obj.lambda_y > eps
                objective = objective + obj.lambda_y*norm(y_slack,2);
            else 
                constraints = [constraints, norm(y_slack,2) <= eps];
            end

            if obj.lambda_g > eps
                objective = objective + obj.lambda_g*norm(g,2);
            end

        
            ops = sdpsettings('solver','mosek');
            ops.verbose = 0;

            problem = optimizer(constraints,objective,ops,{du_ini, dy_ini, y(1), y0, u0, dy_m}, g);

        end

        function [u, obj] = solve(obj, currenty, meal)

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

            dd = zeros(obj.Tf,1);
            for i = 2:obj.Tf
                dd(i) = y_d(i) - y_d(i-1);
            end

            % Solve optimization problem
            g_optimal = obj.problem(obj.uini, obj.yini, currenty, obj.y1, obj.u1, dd);
            du_optimal = obj.Uf*g_optimal;
            du = du_optimal(1);
            u = du + obj.u1;

            obj.uini = [obj.uini; du];            
            obj.yini = [obj.yini; currenty-obj.y1];
            obj.uini(1) = [];
            obj.yini(1) = [];

            obj.u1 = u;
            obj.y1 = currenty;

            obj.d = [obj.d; 0];
            obj.d(1) = [];
        end
        
    end
end