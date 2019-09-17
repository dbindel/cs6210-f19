
tau = zeros(10,1);
kappa_ref = 1;
kappa     = hw3p1(3, tau);
fprintf('Rel err when tau = 0: %e\n', abs(kappa_ref-kappa)/abs(kappa));

tau(4:10) = randn(7,1);
kappa_ref = hw3p1_ref(3, tau);
kappa     = hw3p1(3, tau);
fprintf('Rel err when tau = rand: %e\n', abs(kappa_ref-kappa)/abs(kappa));

% TODO: It may make sense to add a test here for tau large (~10^16)
%       This may require some care, as the reference calculation is not
%       guaranteed to obtain high relative accuracy in this case!

