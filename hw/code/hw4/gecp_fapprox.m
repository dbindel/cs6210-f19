% [X, Y, fpiv] = gecp_fapprox(f, m=20, rtol=1e-8, atol=1e-8)
%
% Approximate a function f on [-1,1]^2 via (at most) m steps of
% Gaussian elimination with complete pivoting.
%
% Arguments:
%     f: Function to be approximated (should be vectorized)
%     m: Maximum number of terms (default 20)
%     rtol (keyword, 1e-8): Relative tolerance for termination
%     atol (keyword, 1e-8): Absolute tolerance for termination
%
% Returns:
%     X: array of X coordinates of pivot points
%     Y: array of Y coordinates of pivot points
%     fpiv: values of the error function (at the points of max modulus)
%
function [X, Y, fpiv] = gecp_fapprox(f, m, rtol, atol)

    % DSB: In my implementation, I maintained an LU factorization
    %      of the matrix with entries f(X[i], Y[j]), stored in packed
    %      form.  I also kept the G and H matrices, computed one column/row
    %      at a time.  I only searched for X and Y over a grid of points;
    %      something more robust might start with a grid search and then
    %      refine with Newton iteration (for example).
 
    % Fill in default arguments
    if nargin < 2, m = 20;      end
    if nargin < 3, rtol = 1e-8; end
    if nargin < 4, atol = 1e-8; end
    
    fpiv = zeros(m, 1);
    X = zeros(m, 1);
    Y = zeros(m, 1);
    LU_fXY = zeros(m, m);

    % Set up sample mesh and evaluate f on grid
    nmesh = 101;
    xmesh = linspace(-1, 1, nmesh);

    [XX, YY] = meshgrid(xmesh, xmesh);
    S = f(XX, YY);
    G = zeros(nmesh, m);
    H = zeros(m, nmesh);
    
    % GECP main loop
    for k = 1:m

        % TODO: Update X, Y, LU_fXY, G, H, and S; quit early if
        %       tolerances are satisfied.
        
    end

end
