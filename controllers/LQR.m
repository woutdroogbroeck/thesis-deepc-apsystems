classdef LQR < handle
    % simple LQR class
    
    properties
        K
        ref
        U_MAX = .1
        U_MIN
        basal
    end
    properties
        kf = []
        xhat
        Nx
        Nu
        d = 0
        u = 0
        steps_until_meal = 0
        Ts
    end
    
    methods
        function obj = LQR(system, Q, R, ref, u_min, u_max, basal, Qk, Rk, P0)

        sys_lqr = ss(system.A,system.B(:,1),system.C,0,system.Ts);
        obj.K = lqr(sys_lqr,Q,R);
        obj.ref = ref;
        obj.basal = basal;
        
        obj.U_MIN = u_min;
        obj.U_MAX = u_max;

        % full state feedback matrices
        [obj.Nx, obj.Nu] = fsf_matrices(sys_lqr);       
        obj.u = 0;
        obj.Ts = system.Ts;
        x0 = zeros(size(system.A,1),1);
        obj.xhat = x0;
        if nargin > 7
            obj.kf = KalmanFilter(system.A, system.B, system.C, Qk, Rk, x0, P0);
        end

        end
        
        function [u, obj] = solve(obj, x, meal)
           if meal.size > 0
                obj.steps_until_meal = 2*meal.minutes/obj.Ts;
            end
            if obj.steps_until_meal > 0
                current_ref = -10;
                obj.steps_until_meal = obj.steps_until_meal - 1;
            else 
                current_ref = obj.ref;
            end
            if ~ isempty(obj.kf)
                x_hat = obj.kf.estimate([obj.u; 0], x);
                u = max(min(obj.basal + obj.Nu*current_ref  -obj.K*(x_hat-obj.Nx*current_ref), obj.U_MAX),obj.U_MIN);

                obj.xhat = [obj.xhat, x_hat];
            else

                u = max(min(obj.basal + obj.Nu*current_ref  -obj.K*(x-obj.Nx*current_ref), obj.U_MAX),obj.U_MIN);
            end
            obj.u = u - obj.basal;        
        end
    end
end

