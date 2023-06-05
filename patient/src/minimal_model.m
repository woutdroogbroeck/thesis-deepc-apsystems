function xdot = minimal_model(x,p,u,d)


Gdot = -p.p1*(x(1)-p.Gb) - x(1)*x(2) + 1e3*d/p.Vg;
Xdot = -p.p2*(x(2)-p.Xb) + p.p3*(x(3)-p.Ib);
Idot = -p.n*x(3) + 1e3*u/p.Vi;

xdot = [Gdot; Xdot; Idot];