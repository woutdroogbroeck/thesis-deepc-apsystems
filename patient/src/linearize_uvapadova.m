function [sys_control, sys_sim] = linearize_uvapadova(p)
    kempt = (p.kmax-p.kmin)/2; 
%     kempt = p.kmax;
    A = zeros(13);
    % glucose intestinal absorption
    % Qsto1
    A(1,1) = -p.kmax;
    % Qsto2
    A(2,1) = p.kmax;
    A(2,2) = -kempt; 
    % Qgut
    A(3,2) = kempt; 
    A(3,3) = -p.kabs;
    
    % glucose subsystem
    % Gp
    Et = p.ke1; % *x(4)
    A(4,3) = p.f*p.kabs/p.BW;
    A(4,4) = -p.kp2-p.k1 -Et;
    A(4,5) = p.k2;
    A(4,9) = -p.kp3;
    % Gt 
    A(5,4) = p.k1;
    A(5,5) = -p.k2 - p.Vm0*p.Km0/(p.Km0+p.Gtb)^2;
    A(5,7) = -p.Vmx*p.Gtb/(p.Km0+p.Gtb);
    
    % insulin subsystem
    % Ip 
    A(6,6)  = -(p.m2+p.m4);
    A(6,10) = p.m1;
    A(6,11) = p.ka1;
    A(6,12) = p.ka2;
    % X 
    A(7,6) = p.p2u/p.Vi;
    A(7,7) = -p.p2u;
    % I1
    A(8,6) = p.ki/p.Vi;
    A(8,8) = -p.ki;
    % Id 
    A(9,8) = p.ki;
    A(9,9) = -p.ki;
    % Il 
    A(10,6) = p.m2;
    A(10,10) = -(p.m1+p.m30);
    % S1
    A(11,11) = -(p.kd+p.ka1);
    % S2
    A(12,11) = p.kd;
    A(12,12) = -p.ka2;
    % Gm
    A(13,4) = p.ksc;
    A(13,13) = -p.ksc;
    
%     B = zeros(13,1);
%     B(11) = 6000/p.BW;
    B = zeros(13,2);
    B(11,1) = 6000/p.BW;
%     Bsim(11,1) = 1/p.BW;
    B(1,2) = 1000;
    
    C = zeros(1,13);
    C(13) = 1/p.Vg;
    
    D = 0;

    sys_control = ss(A,B,C,D);

    Bsim = zeros(13,2);
    Bsim(11,1) = 6000/p.BW;
%     Bsim(11,1) = 1/p.BW;
    Bsim(1,2) = 1000;

    Dsim = [0, 0];

    sys_sim = ss(A,Bsim,C,Dsim);

end