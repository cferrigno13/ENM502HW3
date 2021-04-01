function [Jh, Rh] = AJmaker(ua, ub, lama, lamb)

ds = ((lamb - lama)^2 + (norm(ub - ua))^2)^(1/2);

[J, R] = Jmaker(ub, lamb, 29);

DrDlamb = ub.*(1+ub);

DadaDub = -2 * (ub - ua);

DadaDlamb = -2 * (lamb - lama);

Jh = [J DrDlamb;DadaDub' DadaDlamb];

ada = abs(ds)^2 - norm(ub - ua)^2 - abs(lamb - lama)^2;

Rh = [R;ada];