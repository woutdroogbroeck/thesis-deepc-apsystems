function y = johnson_transform(xi, lam, gamma, delta, x)
    y = xi + lam * sinh((x - gamma)/delta);
end

