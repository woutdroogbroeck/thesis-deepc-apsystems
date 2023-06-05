classdef KalmanFilter < handle
    properties
        A
        B
        C
        Q
        R
        P
        x_hat
        x0
    end

    methods
        function obj = KalmanFilter(A, B, C, Q, R, x0, P0)
            obj.A = A;
            obj.B = B;
            obj.C = C;
            obj.Q = Q;
            obj.R = R;
            obj.x_hat = zeros(size(A,1),1);
            obj.P = P0;
            obj.x0 = x0;
        end

        function [x_hat, obj] = estimate(obj, u1, y)

            P_minus = obj.P;

            x_hat_minus = obj.A * obj.x_hat + obj.B * u1;
            L = P_minus * obj.C' * pinv(obj.C * P_minus * obj.C' + obj.R);
            x_hat = x_hat_minus + L * (y - obj.C * x_hat_minus);

            obj.x_hat = x_hat;
            Pk = P_minus - L*obj.C*P_minus;
            obj.P = obj.A * Pk * obj.A' + obj.Q;
            

        end

        function obj = predict(obj, u, x_hat)
            obj.x_hat = obj.A * x_hat + obj.B * u;
        end
    end
end
