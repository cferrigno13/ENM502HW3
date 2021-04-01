function [U,c,NORM] = newt(u0, lam, N, tol)
nm = tol + 10;
u = u0;
c = 1;

while nm > tol
    [J,R] = Jmaker(u,lam, N);
    du = J\(-R);
    nm = norm(du);
    u = u + du;
    c = c+1;
    if c > 16
        error('too many steps')
    end
end
[U] = u;
NORM = nm;
end