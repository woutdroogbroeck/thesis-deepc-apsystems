classdef MPCI < handle
    % a class that implements a Model Predictive Control algorithm with integral action
    % for a single input/single output system (glucose-insulin)
    % adapted from CACSD class % Oscar Mauricio Agudelo (c) 2019 %
    properties
        N
        Q
        Qt
        R
        system
        ref
        EAT_RATE = 5
        U_MAX = .1
        U_MIN
        Y_MIN
    end
    properties
        problem
        Ts
        kf = []
        xhat 
        d
        d1 = 0
        u1 = 0
        y1 = 0
        x1
        sys
        basal
    end
    
    methods
        function obj = MPCI(N, Q, Qt, R, system, ref, u_min, u_max, y_min, basal, Qk, Rk, P0)
            obj.N = N;
            obj.Q = Q;
            obj.Qt = Qt;
            obj.R = R;
            obj.ref = ref;
            obj.system = system;
            obj.basal = basal;
            obj.U_MIN = u_min;
            obj.U_MAX = u_max-obj.basal;
            obj.Y_MIN = y_min;

            obj.problem = obj.build_problem;
            obj.Ts = system.Ts;
            x0 = zeros(size(system.A,1),1);
            obj.xhat = x0;
            
            obj.x1 = x0;
            % if kalman parameters are given, assume state estimation
            if nargin > 10
                obj.kf = KalmanFilter(system.A, system.B, system.C, Qk, Rk, x0, P0);
            end
            obj.d = zeros(N,1);
        end

        function problem = build_problem(obj)
        % setup the optimization problem, clear yalmip 
            yalmip('clear')
            
            sys_lqr = ss(obj.system.A,obj.system.B(:,1),obj.system.C,0,obj.system.Ts);
            obj.sys = sys_lqr;
            A = obj.system.A;
            B = obj.system.B;
            C = obj.system.C;

            n = size(A,1);  %Number of states
            ny = size(C,1); %Number of outputs
            nu = size(B,2); %Number of Inputs

            % creating the augmented system
            
            A_bar =[A zeros(n,ny)
                C eye(ny,ny)];
            
            B_bar = [B
                zeros(ny,nu)];
            
            C_bar = [C eye(ny,ny)];

            u = sdpvar(repmat(1,1,obj.N+1),repmat(1,1,obj.N+1));
            x = sdpvar(repmat(size(A,1)+ny,1,obj.N+1),repmat(1,1,obj.N+1));
            x2 = sdpvar(repmat(size(A,1)+ny,1,obj.N+1),repmat(1,1,obj.N+1));
            
            y_constr_slack = sdpvar(1,1);
            
            d_m = sdpvar(obj.N+1,1);

            constraints = [];
            constraints = [constraints, x2{1} == x{1}];
            objective = 0;
            for k = 2:obj.N+1
                objective = objective + (C_bar*x{k-1}-obj.ref)'*obj.Q*(C_bar*x{k-1}-obj.ref)+(u{k}-u{k-1})'*obj.R*(u{k}-u{k-1});
                constraints = [constraints, x{k} == A_bar*x{k-1}+B_bar*[u{k}-u{k-1}; d_m(k)-d_m(k-1)]];
                constraints = [constraints, x2{k} == A_bar*x2{k-1}+B_bar*[u{k}-u{k-1}; 0]];
                constraints = [constraints, -obj.basal <= u{k}<= obj.U_MAX];
                constraints = [constraints, obj.Y_MIN <= C_bar*x2{k} + y_constr_slack ];
            end

%             objective = objective + (C_bar*x{obj.N}-obj.ref)'*obj.Qt*(C_bar*x{obj.N}-obj.ref);
            objective = objective + 10^8*y_constr_slack^2;

            ops = sdpsettings('verbose',0,'solver','mosek');
            problem = optimizer(constraints,objective,ops,{x{1},u{1}, d_m},u{2});

        end
        
        function [u, obj] = solve(obj, x, meal)
             d_m = zeros(obj.N,1);
            if meal.size > 0
                steps = meal.minutes/obj.Ts;
                meal_steps = round(meal.size/obj.EAT_RATE);
                meal_d = obj.EAT_RATE*ones(meal_steps,1);                
                obj.d(steps:steps+meal_steps-1) = meal_d;
            end
            d_m(1:obj.N) = obj.d(1:obj.N);
            d_m = [obj.d1; d_m];

            if ~ isempty(obj.kf)
                x_hat = obj.kf.estimate([obj.u1; 0], x);
                x_bar = [x_hat-obj.x1; obj.y1];
                u = obj.problem(x_bar, obj.u1, d_m);
                obj.xhat = [obj.xhat, x_hat];
                obj.x1 = x_hat;
                obj.y1 = x;
                u = max(min(obj.basal+u, obj.U_MAX+obj.basal),obj.U_MIN);
                obj.u1 = u-obj.basal;

            else
                x_bar = [x-obj.x1; obj.y1];
                u = obj.problem(x_bar, obj.u1, d_m);
                obj.x1 = x;
                obj.y1 = obj.system.C*x;
                
                u = max(min(obj.basal+u, obj.U_MAX+obj.basal),obj.U_MIN);
                obj.u1 = u-obj.basal;
            end
            
            obj.d = [obj.d; 0];
            obj.d1 = obj.d(1);
            obj.d(1) = [];
        end
    end
end