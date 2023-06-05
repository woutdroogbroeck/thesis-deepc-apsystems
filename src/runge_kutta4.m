function xn = runge_kutta4(x,f,Ts)
    s1 = f(x);
    s2 = f(x + 0.5*Ts*s1);
    s3 = f(x + 0.5*Ts*s2);
    s4 = f(x + Ts * s3);
    xn = x + Ts/6 * (s1 + 2 * s2 + 2* s3 + s4);
end