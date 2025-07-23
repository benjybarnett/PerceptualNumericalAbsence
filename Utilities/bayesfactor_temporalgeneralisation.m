function [bf, bf_complement] = bayesfactor_temporalgeneralisation(X, varargin)
    % Bayes factor analysis for temporal generalisation matrices (T x T x N)
    %
    % Inputs:
    %   X: T x T x N matrix (time x time x subjects)
    % Optional Name-Value pairs:
    %   'nullvalue' - value under H0 (default: 0)
    %   'interval'  - prior interval (e.g., [0.5 Inf] for directional test)
    %   'rscale'    - Cauchy prior scale (default: 0.7071)
    %
    % Outputs:
    %   bf, bf_complement: T x T matrices of Bayes factors

    % --- Default options ---
    opt = struct();
    opt.nullvalue = 0;
    opt.interval = [-Inf Inf];
    opt.rscale = 0.7071;

    % --- Parse input ---
    fnames = varargin(1:2:end);
    fvalues = varargin(2:2:end);
    for f = 1:numel(fnames)
        assert(isfield(opt, fnames{f}), 'Invalid input: %s', fnames{f});
        opt.(fnames{f}) = fvalues{f};
    end
    assert(isinf(opt.interval(2)), 'Upper bound of interval must be infinite');

    % --- Prep ---
    [T1, T2, N] = size(X);
    X = X - opt.nullvalue;
    r = opt.rscale;
    r_scaled = r * sqrt(N);
    prior_full = @(x) tpdf(x / r_scaled, 1) / r_scaled;
    likelihood_fun = @(theta, t) nctpdf(t, N - 1, theta);
    LIM = @(T) 10 * sqrt((2 / N) + ((T^2) / (4 * N)));

    % --- Reshape data ---
    Xr = reshape(X, T1 * T2, N);
    mu = mean(Xr, 2);
    sd = std(Xr, 0, 2);
    se = sd / sqrt(N);
    tvals = mu ./ se;

    % --- Precompute H0 likelihoods ---
    M0 = likelihood_fun(0, tvals);

    % --- Initialize outputs ---
    bf = zeros(T1*T2, 1);
    bf_complement = zeros(T1*T2, 1);

    % --- Prior definitions ---
    if all(isinf(abs(opt.interval)))
        use_interval = false;
    else
        use_interval = true;
        lowerbound = opt.interval(1) * sqrt(N);
        upperbound = Inf;

        K_a = integral(prior_full, lowerbound, upperbound, 'RelTol',1e-4,'AbsTol',1e-6);
        K_b = 1 - K_a;
        prior_a = @(x) prior_full(x) / K_a;
        prior_b = @(x) prior_full(x) / K_b;
    end

    % --- Parallel computation ---
    parfor i = 1:T1*T2
        t = tvals(i);
        lim = LIM(t);

        if ~use_interval
            M1 = integral(@(x) likelihood_fun(x, t) .* prior_full(x), -lim, lim, 'RelTol',1e-4,'AbsTol',1e-6);
            bf(i) = M1 / M0(i);
            bf_complement(i) = bf(i);
        else
            M1a = integral(@(x) likelihood_fun(x, t) .* prior_a(x), lowerbound, lim, 'RelTol',1e-4,'AbsTol',1e-6);
            M1b = integral(@(x) likelihood_fun(x, t) .* prior_b(x), -lim, lowerbound, 'RelTol',1e-4,'AbsTol',1e-6);
            bf(i) = M1a / M0(i);
            bf_complement(i) = M1b / M0(i);
        end
    end

    % --- Reshape back to time x time ---
    bf = reshape(bf, T1, T2);
    bf_complement = reshape(bf_complement, T1, T2);
end
