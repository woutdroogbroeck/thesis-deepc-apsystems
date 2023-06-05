function [lbgi, hbgi, li, hi] = evaluate_bg(G)
     g = 1.509;
     a = 1.084;
     b = 5.381;
    num_samples = length(G);

    rh = 10*max(g*(log(G).^a-b),0).^2;
    rl = 10*min(g*(log(G).^a-b),0).^2;
        
    hbgi = cumsum(rh)./(1:num_samples);
    lbgi = cumsum(rl)./(1:num_samples);

    hi = hbgi(end);
    li = lbgi(end);

end


 