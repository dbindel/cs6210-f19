% === Test things out on some sample functions ===

fprintf("\n=== Squared exp ===\n");
bump = @(x, y) exp(-(x.^2+y.^2)/2);
[X, Y, fpiv] = gecp_fapprox(bump, 20);
print_diagnostics(X, Y, fpiv);

fprintf("\n=== Ackley ===\n");
ackley = @(x, y) (20 - 20 * exp(-0.2*sqrt((x.^2+y.^2)/2)) ...
                  - exp((cos(2*pi*x)+cos(2*pi*y))/2));
[X, Y, fpiv] = gecp_fapprox(ackley, 20);
print_diagnostics(X, Y, fpiv);
