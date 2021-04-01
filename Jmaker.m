function [J,R] = Jmaker(u,lam, N)
n = (N+1)^2;
R = zeros(n,1);
J = zeros(n,n);
b = ones(n,1);

for i = 1:n
    if i <= N
        J(i,i) = 1;
        R(i) = u(i);
        b(i) = 0;
    end
    for j = 1:N+1
        if i == j*(N+1)
            J(i,i) = 1;
            R(i) = u(i);
            b(i) = 0;
        end
        if i == j*(N+1)+1
            J(i,i) = 1;
            R(i) = u(i);
            b(i) = 0;
        end
    end
    if i <= n && i > (n-(N+1))
        J(i,i) = 1;
        R(i) = u(i);
        b(i) = 0;
    end
    if b(i) ~= 0
        J(i,i-(N+1)) = N^2;
        J(i,i-1) = N^2;
        J(i,i) = -4*(N^2) + lam*(1+2*u(i));
        J(i,i+1) = N^2;
        J(i,i+(N+1)) = N^2;
        R(i) = N^2 * (u(i+1)-4*u(i)+u(i-1)+u(i+(N+1))+u(i-(N+1))) + lam*u(i)*(1+u(i));
    end
end
[J] = sparse(J);
%[J] = J;
[R] = R;
end