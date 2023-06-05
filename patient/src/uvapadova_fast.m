function dx = uvapadova_fast(x, param, insulin, CHO, last_Qsto, last_foodtaken)
    % model equations from https://doi.org/10.1177%2F1932296813514502
    % not efficiently implemented, but easier to understand
    dx = zeros(13,1);
    d = CHO*1000;  % g -> mg
    u = insulin*6000/param.BW;  % U/min -> pmol/kg/min
    
%     Qsto1  = x(1);
%     Qsto2  = x(2);
%     Qgut   = x(3);
%     Gp     = x(4);
%     Gt     = x(5);
%     Ip     = x(6);
%     X      = x(7);
%     I1     = x(8);
%     Id     = x(9);
%     Il     = x(10);
%     S1     = x(11);
%     S2     = x(12);
%     Gs     = x(13);

    %% glucose rate of appearance 
    % stomach
    Qsto = x(1) + x(2);

    Dbar = last_Qsto + last_foodtaken * 1000;  % mg
    % rate constant of gastric emptying kempt
    if Dbar > 0
        alpha = 5/(2*(1-param.b)*Dbar);
        beta = 5/(2*param.d*Dbar);
        kempt = param.kmin + (param.kmax - param.kmin)/2 * ( ...
            tanh(alpha*(Qsto - param.b*Dbar)) - ...
            tanh(beta*(Qsto - param.d*Dbar))+2);
    else
        kempt = param.kmax;
    end

    % solid
    dx(1) = -param.kmax*x(1) + d;
%     dx(1) = dx(1);

    % liquid
    dx(2) = -kempt*x(2) + param.kmax*x(1);
%     dx(2) = dx(2);

    % intestine
    dx(3) = -param.kabs*x(3) + kempt*x(2);
%     dx(3) = dx(3);

    % glucose appearance rate in plasma
    Ra = param.f*param.kabs*x(3)/param.BW;

    %% endogenous glucose production

    EGP = param.kp1 - param.kp2*x(4) - param.kp3*x(9);  
   
    %% glucose subsystem
    % glucose utilization
    Uii = param.Fsnc;

    % excretion by kidney
    if x(4)>param.ke2
        E = param.ke1*(x(4)-param.ke2);
    else
        E = 0;
    end
    % plasma glucose
    dx(4) = max(EGP,0) + Ra - Uii - E - param.k1*x(4) + param.k2*x(5);
    dx(4) = (x(4)>0)*dx(4);

    Vm = param.Vm0 + param.Vmx*x(7);
    Km = param.Km0;

    % insulin-dependent utilization
    Uid = Vm*x(5)/(Km+x(5));
    
    % tissue glucose
    dx(5) = -Uid + param.k1*x(4) - param.k2*x(5);
    dx(5) = (x(5)>0)*dx(5);

    %% insulin subsystem
    
    % insulin plasma
    dx(6) = -(param.m2+param.m4)*x(6) + param.m1*x(10) + param.ka1*x(11) + param.ka2*x(12);
    I = x(6)/param.Vi;    
    dx(6) = (x(6)>0)*dx(6);

    % insulin in interstitial fluid
    dx(7) = -param.p2u*x(7) + param.p2u*(I-param.Ib);
%     dx(7) = dx(7);

    dx(8) = -param.ki*(x(8)-I);
%     dx(8) = dx(8);

    dx(9) = -param.ki*(x(9)-x(8));
%     dx(9) = dx(9);

    % insulin liver
    dx(10) = -(param.m1+param.m30)*x(10) + param.m2*x(6);
    dx(10) = (x(10)>0)*dx(10);

    % subcutaneous insulin kinetics
    dx(11) = u - (param.kd+param.ka1)*x(11);
    dx(11) = (x(11)>0)*dx(11);
    
    dx(12) = param.kd*x(11) - param.ka2*x(12);
    dx(12) = (x(12)>0)*dx(12);

    % subcutaneous glucose kinetics
    dx(13) = -param.ksc*x(13) + param.ksc*x(4);
    dx(13) = (x(13)>0)*dx(13);    


end