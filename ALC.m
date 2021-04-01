function [uc0, lamc0] = ALC(ua, ub, lama, lamb)

ds = ((lamb - lama)^2 + (norm(ub - ua))^2)^(1/2);

% Now build augmented Jacobian

Jh = AJmaker(ua, ub, lama, lamb);

% Now find dRhat/ds

DRDs = zeros(900,1);

DRhDs = [DRDs;2*ds];

% Now solve for du/ds and dlambda/ds

Z = Jh\(-DRhDs);

DubDs = Z(1:900);

DlambDs = Z(901);


% Find initial guesses for the third point

uc0 = ub + ds * DubDs;

lamc0 = lamb + ds * DlambDs;
end
