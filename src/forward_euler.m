function xn = forward_euler(x,f,Ts)
    xn = x + Ts * f(x); 
end