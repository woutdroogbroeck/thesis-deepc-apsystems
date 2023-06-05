classdef LQRm < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        K
        ref
        U_MAX = .1
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
        function obj = LQRm(system, Q, R, ref, basal, Qk, Rk, P0)
        %UNTITLED Construct an instance of this class
        %   Detailed explanation goes here

        sys_lqr = ss(system.A,system.B(:,1),system.C,0,system.Ts);
        obj.K = lqr(sys_lqr,Q,R);
        obj.ref = ref;
        obj.basal = basal;
        % full state feedback matrices
        [obj.Nx, obj.Nu] = fsf_matrices(sys_lqr);       
        obj.u = 0;
        obj.Ts = sys_lqr.Ts;
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
                x_hat = obj.kf.estimate([obj.u; 0], x);
%                   x_hat = x;
%                 u =  min(max(-obj.K*(x_hat-obj.Nx*obj.ref),0), obj.U_MAX);
%                 if obj.fsf
%                     u = max(-obj.K*(x_hat-obj.Nx*obj.ref),0);
%                     
%                     u = min(obj.basal + obj.Nu*obj.ref + u, obj.U_MAX);

                    u = max(min(obj.basal + obj.Nu*current_ref  -obj.K*(x_hat-obj.Nx*current_ref), obj.U_MAX),0);
%                 else
%                     u = max(-obj.K*(x_hat),0);
%                     u = min(obj.basal + u, obj.U_MAX);
%                 end
% %                 u = u + obj.Nu*obj.ref;
%                 obj.kf.predict([u; 0], x_hat);
                obj.xhat = [obj.xhat, x_hat];
%                 obj.u = u;
            else
%                 u = min(max(-obj.K*(x-obj.Nx*obj.ref),0), obj.U_MAX);
%                 if obj.fsf
                    u = max(-obj.K*(x-obj.Nx*obj.ref),0);
                    u = min(obj.basal + obj.Nu*obj.ref + u, obj.U_MAX);
%                 else
%                     u = max(-obj.K*(x),0);
%                     u = min(obj.basal + u, obj.U_MAX);
%                 end
%                  u = u + obj.Nu*obj.ref;
            end
            obj.u = u - obj.basal;        
        end
    end
end

