function [u] = disc(A,n,m)
N = 29;
h = 1/N;
i = 1;
u = zeros((N+1)^2,1);

for y = 0:h:1
    for x = 0:h:1
        u(i) = A * sin(n*pi*x) * sin(m*pi*y);
        i = i+1;
    end
end