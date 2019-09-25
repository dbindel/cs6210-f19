% Print diagnostics
%
function print_diagnostics(X, Y, fpiv)
    for k = 1:length(X)
        fprintf("%d: err(%g, %g) = %e\n", k, X(k), Y(k), fpiv(k));
    end
end
