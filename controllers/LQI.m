classdef LQI < handle
    % simple lqr controller with integral action
    
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
        d = 0;
        u
        sys
        integral = 0;
        prev_error = 0;
        steps_until_meal = 0
        Ki
        Ts;
    end
    
    methods
        function obj = LQI(system, Q, R, ref, u_min, u_max, basal, Qk, Rk, P0)
    
        sys_lqr = ss(system.A,system.B(:,1),system.C,0,system.Ts);
        obj.sys = sys_lqr;

        n = size(sys_lqr.A,1);  %Number of states 
        ny = size(sys_lqr.C,1); %Number of outputs
        nu = size(sys_lqr.B,2); %Number of Inputs

        obj.U_MIN = u_min;
        obj.U_MAX = u_max;
        obj.Ts = system.Ts;

        Adi = [eye(ny) sys_lqr.C;
            zeros(n,1) sys_lqr.A];

        Bdi = [0; sys_lqr.B];
        Cdi = [0 sys_lqr.C];
        sysi = ss(Adi,Bdi,Cdi,sys_lqr.D,system.Ts);
        
        Ki_total = lqr(sysi,Q,R);
        
        obj.K = Ki_total(2:end);
        obj.Ki = Ki_total(1);

        obj.ref = ref;
        obj.basal = basal;

        % full state feedback matrices
        [obj.Nx, obj.Nu] = fsf_matrices(sys_lqr);       
        obj.u = 0;

        x0 = zeros(size(system.A,1),1);
        obj.xhat = x0;
        if nargin > 5
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
                x_hat = obj.kf.estimate([obj.u;0], x);
                error = x - obj.ref;
                obj.integral = obj.integral + obj.prev_error;
                obj.prev_error = error;
                ui = -obj.Ki*obj.integral;


                u = max(min(obj.basal -obj.K*(x_hat-obj.Nx*current_ref)+ ui, obj.U_MAX),obj.U_MIN);
                obj.xhat = [obj.xhat, x_hat];

            else
                    error = obj.sys.C*x - obj.ref;
                    obj.integral = obj.integral + obj.prev_error;
                    obj.prev_error = error;
                    ui = -obj.Ki*obj.integral;
                    u = max(-obj.K*(x-obj.Nx*current_ref) + ui,0);
                    u = min(obj.basal + u, obj.U_MAX);

            end
            obj.u = u - obj.basal;        
        end
    end
end

