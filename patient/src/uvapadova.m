function dx = uvapadova(x, param, insulin, CHO, last_Qsto, last_foodtaken)
    % model equations from https://doi.org/10.1177%2F1932296813514502
    % not efficiently implemented, but easier to understand
    dx = zeros(13,1);
    d = CHO*1000;  % g -> mg
    u = insulin*6000/param.BW;  % U/min -> pmol/kg/min
    
    Qsto1  = x(1);
    Qsto2  = x(2);
    Qgut   = x(3);
    Gp     = x(4);
    Gt     = x(5);
    Ip     = x(6);
    X      = x(7);
    I1     = x(8);
    Id     = x(9);
    Il     = x(10);
    S1     = x(11);
    S2     = x(12);
    Gs     = x(13);

    %% glucose rate of appearance 
    % stomach
    Qsto = Qsto1 + Qsto2;

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
    dQsto1 = -param.kmax*Qsto1 + d;
    dx(1) = dQsto1;

    % liquid
    dQsto2 = -kempt*Qsto2 + param.kmax*Qsto1;
    dx(2) = dQsto2;

    % intestine
    dQgut = -param.kabs*Qgut + kempt*Qsto2;
    dx(3) = dQgut;

    % glucose appearance rate in plasma
    Ra = param.f*param.kabs*Qgut/param.BW;

    %% endogenous glucose production

    EGP = param.kp1 - param.kp2*Gp - param.kp3*Id;  
   
    %% glucose subsystem
    % glucose utilization
    Uii = param.Fsnc;

    % excretion by kidney
    if Gp>param.ke2
        E = param.ke1*(Gp-param.ke2);
    else
        E = 0;
    end
    % plasma glucose
    dGp = max(EGP,0) + Ra - Uii - E - param.k1*Gp + param.k2*Gt;
    dx(4) = (Gp>0)*dGp;

    Vm = param.Vm0 + param.Vmx*X;
    Km = param.Km0;

    % insulin-dependent utilization
    Uid = Vm*Gt/(Km+Gt);
    
    % tissue glucose
    dGt = -Uid + param.k1*Gp - param.k2*Gt;
    dx(5) = (Gt>0)*dGt;

    %% insulin subsystem
    
    % insulin plasma
    dIp = -(param.m2+param.m4)*Ip + param.m1*Il + param.ka1*S1 + param.ka2*S2;
    I = Ip/param.Vi;    
    dx(6) = (Ip>0)*dIp;

    % insulin in interstitial fluid
    dX = -param.p2u*X + param.p2u*(I-param.Ib);
    dx(7) = dX;

    dI1 = -param.ki*(I1-I);
    dx(8) = dI1;

    dId = -param.ki*(Id-I1);
    dx(9) = dId;

    % insulin liver
    dIl = -(param.m1+param.m30)*Il + param.m2*Ip;
    dx(10) = (Il>0)*dIl;

    % subcutaneous insulin kinetics
    dS1 = u - (param.kd+param.ka1)*S1;
    dx(11) = (S1>0)*dS1;
    
    dS2 = param.kd*S1 - param.ka2*S2;
    dx(12) = (S2>0)*dS2;

    % subcutaneous glucose kinetics
    dGs = -param.ksc*Gs + param.ksc*Gp;
    dx(13) = (Gs>0)*dGs;    


end