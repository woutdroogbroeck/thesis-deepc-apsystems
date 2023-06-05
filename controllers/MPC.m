classdef MPC < handle
     % a class that implements a Model Predictive Control algorithm for a
     % single input/single output system (glucose-insulin)
    properties
        N
        Q
        Qt
        R
        system
        ref
        EAT_RATE = 5
        U_MAX = 0.1
        U_MIN
        Y_MIN
    end
    properties
        problem
        Ts
        kf = []
        xhat 
        u1
        d
        basal
    end
    
    methods
        function obj = MPC(N, Q, Qt, R, system, ref, u_min, u_max, y_min, basal, Qk, Rk, P0)
            obj.N = N;
            obj.Q = Q;
            obj.Qt = Qt;
            obj.R = R;
            obj.ref = ref;
            obj.u1 = 0;
            obj.system = system;
            obj.basal = basal;
            obj.U_MIN = u_min;
            obj.U_MAX = u_max - basal;
            obj.Y_MIN = y_min;

            obj.problem = obj.build_problem;
            obj.Ts = system.Ts;
            x0 = zeros(size(system.A,1),1);
            obj.xhat = x0;
            
            obj.d = zeros(N,1);
            % if kalman parameters are given, assume state estimation
            if nargin > 10
                obj.kf = KalmanFilter(system.A, system.B, system.C, Qk, Rk, x0, P0);
            end
        end

        function problem = build_problem(obj)
            % setup the optimization problem, clear yalmip 
            yalmip('clear')
            
            % system matrices for prediction model
            A = obj.system.A;
            B = obj.system.B;
            C = obj.system.C;

            % initialize optimization variables
            u = sdpvar(repmat(1,1,obj.N),repmat(1,1,obj.N));
            x = sdpvar(repmat(size(A,1),1,obj.N+1),repmat(1,1,obj.N+1));

            % copy of x for output constraint, ensure that even if the meal
            % prediction was wrong, output stays above bounds
            x2 = sdpvar(repmat(size(A,1),1,obj.N+1),repmat(1,1,obj.N+1));

            % slack variable for output constraint
            y_constr_slack = sdpvar(1,1);
            
            % predicted meal disturbance vector
            d_m = sdpvar(obj.N,1);

            constraints = [];
            constraints = [constraints, x2{1} == x{1}];
            objective = 0;
            % penalty on output and input
            for k = 1:obj.N
                objective = objective + (C*x{k}-obj.ref)'*obj.Q*(C*x{k}-obj.ref)+u{k}'*obj.R*u{k};
                % system dynamics
                constraints = [constraints, x{k+1} == A*x{k}+B*[u{k}; d_m(k) ]];
                constraints = [constraints, x2{k+1} == A*x2{k}+B*[u{k}; 0 ]];
                % input and output constraints
                constraints = [constraints, -obj.basal <= u{k}<= obj.U_MAX];
                constraints = [constraints, obj.Y_MIN <= C*x2{k}  + y_constr_slack];
            end
            % terminal constraint if Qt is nonzero
            objective = objective + (C*x{obj.N}-obj.ref)'*obj.Qt*(C*x{obj.N}-obj.ref);

             % high penalty on output constraint slack variable, its only
            % function is to avoid feasibility issues
            objective = objective + 10^9*y_constr_slack^2;

            ops = sdpsettings('verbose',0,'solver','mosek');
            
            % build the problem: input current state and estimated
            % disturbance vector, returns the optimal input u_1
            problem = optimizer(constraints,objective,ops,{x{1},d_m},u{1});

        end
        
        function [u, obj] = solve(obj, x, meal)
         % calculates the optimal control action
            % args
            %    x: current state or current output of the system
            %    meal: struct with meal information: size (in gram) and 
            %          steps until meal is consumed
            % returns
            %    u: The calculated control action

            d_m = zeros(obj.N,1);
            if meal.size > 0
                steps = meal.minutes/obj.Ts;
                meal_steps = round(meal.size/obj.EAT_RATE);
                meal_d = obj.EAT_RATE*ones(meal_steps,1);                
                obj.d(steps:steps+meal_steps-1) = meal_d;
            end
            d_m(1:obj.N) = obj.d(1:obj.N);

            if ~ isempty(obj.kf)
                x_hat = obj.kf.estimate([obj.u1; 0],x);
                u = obj.problem(x_hat, d_m);
                obj.xhat = [obj.xhat, x_hat];
            else
                u = obj.problem(x, d_m);
            end

            u = max(min(obj.basal+u, obj.U_MAX+obj.basal),obj.U_MIN);
            obj.u1 = u - obj.basal;
            obj.d = [obj.d; 0];
            obj.d(1) = [];

        end
    end
end

