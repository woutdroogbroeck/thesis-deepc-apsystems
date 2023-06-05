function [sys_control, sys_sim] = linearize_minimal_model(p)

    % linearization of minimal model
    A = [-p.p1 -p.Gb  0;
          0  -p.p2  p.p3;
          0   0  -p.n];
%     B = [0;
%          0;
%          1e3/p.Vi];
    B = [0 1e3/p.Vg;
         0 0;
         1e3/p.Vi 0];
    C = [1 0 0];
    D = [0, 0];

    sys_control = ss(A,B,C,D);

    % state space for simulink 
    B = [0 1e3/p.Vg;
         0 0;
         1e3/p.Vi 0];
    C = eye(3);
    D = [0 0;
         0 0;
         0 0];

    sys_sim = ss(A,B,C,D);

end