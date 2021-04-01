function [u,lam] = FullNewt(ua, ub0, lama, lamb0, Ltol, Utol)
Unorm = Utol + 10;
Lnorm = Ltol + 10;
c = 1;
ub = ub0;
lamb = lamb0;

while ((Unorm > Utol) || (Lnorm > Ltol))
    [Jhb,Rhb] = AJmaker(ua, ub, lama, lamb);
    Q = Jhb\(-Rhb);
    du = Q(1:900);
    dl = Q(901);
    Unorm = norm(du);
    Lnorm = abs(dl);
    ub = ub + du;
    lamb = lamb + dl;
    c = c + 1;
    if c > 16
        error('too many steps')
    end
end
u = ub;
lam = lamb;
