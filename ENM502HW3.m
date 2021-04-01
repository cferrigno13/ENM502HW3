%% ENM 502 HW 3



%% Finding initial guess for first eigenvalue using Newtons Method

% First doing n = 1, m = 1, so we are near lambda = 2 pi^2, lets start with
% slightly less than 2 pi^2

u00 = disc(1,1,1);
lam = 19;
N = 29;
tol = 10^-6;

[u0,c,NORM] = newt(u00,lam,N,tol);

%% Now one step of analytic continuation to get u1 at lambda = 18

% We have u0 at lambda = 19

% Find J
J = Jmaker(u0,lam,N);

% Find dr/dlambda

drdlam = u0.*(1+u0);

% Solve for du/dlamda

dudlam = J\(-drdlam);

% Initial guess for u1

u10 = u0 - dudlam;

% Newtons to converge to find u1 at lambda = 18

u1 = newt(u10, 18, N, tol);

%% Now Arc Length Continuation

% My function ALC takes two pairs of u solutions and lambdas and returns
% intelligent guesses for the next point using arc length continuation

[u20, lam20] = ALC(u0,u1,19,18);

% Now my full newton function takes a known initial U lambda combo and an
% initial guess at a second U lambda combo and runs the full newtons method
% with an augmented residual and jacobian to return the next U lambda
% solution.

Utol = 10^-6;
Ltol = 10^-5;

[u2, lam2] = FullNewt(u1, u20, 18, lam20, Ltol, Utol);

%%  Getting the rest of the curve for lambda less than 2 pi^2

lambda = [19;18;lam2];

U = [u0 u1 u2];

%UN = [norm(u0);norm(u1);norm(u2)];

for i = 3:29
    [ug,lg] = ALC(U(:,i-1),U(:,i),lambda(i-1),lambda(i));
    [U(:,i+1),lambda(i+1)] = FullNewt(U(:,i),ug,lambda(i),lg,Ltol,Utol);
end

UN = zeros(size(U,2),1);

for i = 1:size(U,2)
    UN(i) = norm(U(:,i));
end

%% Now lets jump on the curve to the right of lambda = 2 pi^2

u00R = disc(-1,1,1);
LR0 = 20;
uR0 = newt(u00R,LR0,N,tol);

% This puts us on the curve at lambda = 20

% One step of analytic continuation:

J = Jmaker(uR0,LR0,N);
drdlam = uR0.*(1+uR0);
dudlam = J\(-drdlam);
uR10 = uR0 + dudlam;
uR1 = newt(uR10, LR0+1, N, tol);

% Now take it away with arc length continuation

lambdaR = [LR0; LR0+1];
UR = [uR0,uR1];

for i = 2:19
    [ug,lg] = ALC(UR(:,i-1),UR(:,i),lambdaR(i-1),lambdaR(i));
    [UR(:,i+1),lambdaR(i+1)] = FullNewt(UR(:,i),ug,lambdaR(i),lg,Ltol,Utol);
end

UNR = zeros(size(UR,2),1);

for i = 1:size(UR,2)
    UNR(i) = norm(UR(:,i));
end

% We now have lambda's and U's to the right of lambda = 2 pi^2

%% The next solution curve: lambda = 5 pi^2, first with n = 2 and m = 1 and positive A

u_00 = disc(1,2,1);
lam_0 = 50;

u_0 = newt(u_00,lam_0,N,tol);

% We are on the curve for lambda = 50

% One step of analytic continuation:
J = Jmaker(u_0,lam_0,N);
drdlam = u_0.*(1+u_0);
dudlam = J\(-drdlam);
u_10 = u_0 + dudlam;
u_1 = newt(u_10, lam_0+1, N, tol);

% Finish it with arc length
lambda_ = [lam_0; lam_0+1];
U_ = [u_0,u_1];

for i = 2:9
    [ug,lg] = ALC(U_(:,i-1),U_(:,i),lambda_(i-1),lambda_(i));
    [U_(:,i+1),lambda_(i+1)] = FullNewt(U_(:,i),ug,lambda_(i),lg,Ltol,Utol);
end

UN_ = zeros(size(U_,2),1);

for i = 1:size(U_,2)
    UN_(i) = norm(U_(:,i));
end

%% Same n = 2, m = 1, but now with negative A to trace the 'valley' instead of the 'hill'

u_00V = disc(-1,2,1);
lam_0V = 50;

u_0V = newt(u_00V,lam_0V,N,tol);

% We are on the curve for lambda = 50

% One step of analytic continuation:
J = Jmaker(u_0V,lam_0V,N);
drdlam = u_0V.*(1+u_0V);
dudlam = J\(-drdlam);
u_10V = u_0V + dudlam;
u_1V = newt(u_10V, lam_0V+1, N, tol);

% Finish it with arc length
lambda_V = [lam_0V; lam_0V+1];
U_V = [u_0V,u_1V];

for i = 2:9
    [ug,lg] = ALC(U_V(:,i-1),U_V(:,i),lambda_V(i-1),lambda_V(i));
    [U_V(:,i+1),lambda_V(i+1)] = FullNewt(U_V(:,i),ug,lambda_V(i),lg,Ltol,Utol);
end

UN_V = zeros(size(U_V,2),1);

for i = 1:size(U_V,2)
    UN_V(i) = norm(U_V(:,i));
end

%% For lambda near 5 pi^2, n = 1, m = 2, and positive A

u_00_ = disc(1,1,2);
lam_0_ = 50;

u_0_ = newt(u_00_,lam_0_,N,tol);

% We are on the curve for lambda = 50

% One step of analytic continuation:
J = Jmaker(u_0_,lam_0_,N);
drdlam = u_0_.*(1+u_0_);
dudlam = J\(-drdlam);
u_10_ = u_0_ + dudlam;
u_1_ = newt(u_10_, lam_0_+1, N, tol);

% Finish it with arc length
lambda__ = [lam_0_; lam_0_+1];
U__ = [u_0_,u_1_];

for i = 2:9
    [ug,lg] = ALC(U__(:,i-1),U__(:,i),lambda__(i-1),lambda__(i));
    [U__(:,i+1),lambda__(i+1)] = FullNewt(U__(:,i),ug,lambda__(i),lg,Ltol,Utol);
end

UN__ = zeros(size(U__,2),1);

for i = 1:size(U__,2)
    UN__(i) = norm(U__(:,i));
end

%% Now for negative A, so a valley

u_00_V = disc(-1,1,2);
lam_0_V = 50;

u_0_V = newt(u_00_V,lam_0_V,N,tol);

% We are on the curve for lambda = 50

% One step of analytic continuation:
J = Jmaker(u_0_V,lam_0_V,N);
drdlam = u_0_V.*(1+u_0_V);
dudlam = J\(-drdlam);
u_10_V = u_0_V + dudlam;
u_1_V = newt(u_10_V, lam_0_V+1, N, tol);

% Finish it with arc length
lambda__V = [lam_0_V; lam_0_V+1];
U__V = [u_0_V,u_1_V];

for i = 2:9
    [ug,lg] = ALC(U__V(:,i-1),U__V(:,i),lambda__V(i-1),lambda__V(i));
    [U__V(:,i+1),lambda__V(i+1)] = FullNewt(U__V(:,i),ug,lambda__V(i),lg,Ltol,Utol);
end

UN__V = zeros(size(U__V,2),1);

for i = 1:size(U__V,2)
    UN__V(i) = norm(U__V(:,i));
end

%% Combining all of these into a plot for norm of u versus lambda

l = 1:62;
t = zeros(1,62);


uvlam = figure(1);

plot(lambda,UN, 'b--')
hold on
plot(lambdaR,-UNR,'r--')
plot(lambda_,UN_,'bo')
plot(lambda_V,-UN_V,'ro')
plot(lambda__,UN__,'b')
plot(lambda__V,-UN__V,'r')
plot(l,t,'k-')
hold off

legend('A = 1, m=n=1, \lambda_0=19','A = -1, m=n=1, \lambda_0=20','A = 1, n=2,m=1, \lambda_0=50','A = -1, n=2,m=1, \lambda_0=50','A = 1, n=1,m=2, \lambda_0=50','A = -1, n=1,m=2, \lambda_0=50','Trivial Solution')
xlabel('\lambda')
ylabel('||U||')
title('Solution Trace Plot')
saveas(uvlam,'Solution_Trace.png')


%% What the real solution curve looks like without counting 'valleys' as negative
uvlamTrue = figure(2);

plot(lambda,UN, 'b')
hold on
plot(lambdaR,UNR,'b')
plot(lambda_,UN_,'b')
plot(lambda_V,UN_V,'b')
plot(lambda__,UN__,'b')
plot(lambda__V,UN__V,'b')
plot(l,t,'b')
hold off

%legend('A = 1, m=n=1, \lambda_0=19','A = -1, m=n=1, \lambda_0=20','A = 1, n=2,m=1, \lambda_0=50','A = -1, n=2,m=1, \lambda_0=50','A = 1, n=1,m=2, \lambda_0=50','A = -1, n=1,m=2, \lambda_0=50','Trivial Solution')
xlabel('\lambda')
ylabel('||U||')
title('True Solution Trace Plot')
saveas(uvlamTrue,'True_Solution_Trace.png')


%% Plots for different U's

% Lambda about 7, from curve started with A=n=m=1 at lamda0 = 19
% Mountain
cont7 = figure(3);
surf(reshape(U(:,25),30,30))
title('Solution for \lambda \approx 7.07')
saveas(cont7,'sol_7.png')

% Lambda about 15 from curve started with A=n=m=1 at lambda0 = 19
% Mountain
cont15 = figure(4);
surf(reshape(U(:,5),30,30))
title('Solution for \lambda \approx 15.36')
saveas(cont15,'sol_15.png')

% Lambda about 23 from curve started with A=-1,n=m=1 and lambda0 = 20
% Valley
cont23 = figure(5);
surf(reshape(UR(:,4),30,30))
title('Solution for \lambda \approx 23.13')
saveas(cont23,'sol_23.png')

% Lambda about 40 from curve started with A=-1,n=m=1 and lambda0 = 20
% Valley
cont41 = figure(6);
surf(reshape(UR(:,18),30,30))
title('Solution for \lambda \approx 40.71')
saveas(cont41,'sol_41.png')

% Lambda about 55 from curve started with A=1,n=2,m=1 and lambda0 = 50
% Mountain
cont55m = figure(7);
surf(reshape(U_(:,5),30,30))
title('Mountain Solution for \lambda \approx 54.73')
saveas(cont55m,'sol_55m.png')

% Lambda about 55 from curve started with A=-1,n=2,m=1 and lambda0 = 50
% Valley
cont55v = figure(8);
surf(reshape(U_V(:,5),30,30))
title('Valley Solution for \lambda \approx 54.73')
saveas(cont55v,'sol_55v.png')

% Lambda about 59 from curve started with A=1,n=1,m=2 and lambda0 = 50
% Mountain
cont59m = figure(9);
surf(reshape(U__(:,8),30,30))
title('Mountain Solution for \lambda \approx 58.85')
saveas(cont59m,'sol_59m.png')

% Lambda about 59 from curve started with A=-1,n=1,m=2 and lambda0 = 50
% Valley
cont59v = figure(10);
surf(reshape(U__V(:,8),30,30))
title('Valley Solution for \lambda \approx 58.85')
saveas(cont59v,'sol_59v.png')