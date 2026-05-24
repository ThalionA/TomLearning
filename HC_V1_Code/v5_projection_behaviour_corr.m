function out = v5_projection_behaviour_corr(u, behaviour, n_surrogates, seed)
%V5_PROJECTION_BEHAVIOUR_CORR  Correlate projected CCA activity with behaviour.
%
%   OUT = V5_PROJECTION_BEHAVIOUR_CORR(U, BEHAVIOUR) computes the linear
%   correlation between a projected canonical variate u(trial, bin) and
%   each behavioural variable provided. Follows the Gonzalez et al.
%   procedure: project subspace weights through local activity, correlate
%   with position / speed / trial number, contrast against a circular-
%   shift surrogate distribution.
%
%   Inputs
%       u            (n_trials x n_pos_bins) projected canonical variate.
%                    Typically v5_project(P_x_3d, x_mean, A(:, d)) for one
%                    canonical dimension d.
%       behaviour    struct with one or more of:
%                       .position    (1 x n_pos_bins) or
%                                    (n_trials x n_pos_bins) — position
%                                    index or value per bin.
%                       .speed       (n_trials x n_pos_bins) — speed per
%                                    (trial, bin). Optional.
%                       .trial_number  (n_trials x 1) — trial index. Will
%                                    be broadcast across bins.
%                    Any subset of these may be present; only present
%                    fields are correlated.
%       n_surrogates integer, default 100.
%       seed         RNG seed for the circular-shift surrogate. Default 0.
%
%   Output
%       out          struct with one sub-struct per behaviour variable
%                    present in input. Each sub-struct has fields:
%                       .r          Pearson correlation
%                       .r2         r^2
%                       .p          p-value from corrcoef
%                       .null_r2    (n_surrogates x 1) surrogate r^2
%                       .p_perm     fraction of null_r2 >= observed r2
%
%   Surrogate
%       Circularly shift u along the flattened-sample axis by a random
%       amount. This breaks the temporal alignment with the behavioural
%       variables while preserving the marginal distribution of u and
%       the autocorrelation structure within u. Matches the shuffle
%       described in the paper for "spatial information of pCCA
%       projected activity" (their procedure shuffles weights; we
%       shuffle the projection along the time axis — both produce a
%       null where the behaviour-projection link is broken while
%       respecting the marginal distribution of u).

    if nargin < 3 || isempty(n_surrogates), n_surrogates = 100; end
    if nargin < 4 || isempty(seed),         seed         = 0;   end

    if ~ismatrix(u)
        error('v5_projection_behaviour_corr:badShape', 'u must be 2D');
    end
    [n_trials, n_pos_bins] = size(u);
    if ~isstruct(behaviour)
        error('v5_projection_behaviour_corr:badShape', ...
            'behaviour must be a struct');
    end

    u_flat = u(:);
    n_samples = n_trials * n_pos_bins;

    rng_state = rng(seed, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    shift_amounts = randi(n_samples, n_surrogates, 1) - 1;

    out = struct();
    field_names = fieldnames(behaviour);
    for f = 1:numel(field_names)
        name = field_names{f};
        x = expand_behaviour(behaviour.(name), n_trials, n_pos_bins, name);
        x_flat = x(:);

        % Skip degenerate columns (constant behaviour).
        if std(x_flat) < 1e-15 || std(u_flat) < 1e-15
            out.(name) = struct( ...
                'r', NaN, 'r2', NaN, 'p', NaN, ...
                'null_r2', NaN(n_surrogates, 1), 'p_perm', NaN);
            continue;
        end

        [r_obs, p_obs] = corrcoef_safe(u_flat, x_flat);
        r2_obs = r_obs ^ 2;

        null_r2 = zeros(n_surrogates, 1);
        for s = 1:n_surrogates
            u_shift = circshift(u_flat, shift_amounts(s));
            [r_s, ~] = corrcoef_safe(u_shift, x_flat);
            null_r2(s) = r_s ^ 2;
        end
        p_perm = mean(null_r2 >= r2_obs);

        out.(name) = struct( ...
            'r', r_obs, 'r2', r2_obs, 'p', p_obs, ...
            'null_r2', null_r2, 'p_perm', p_perm);
    end
end


% -------------------------------------------------------------------- %
function x = expand_behaviour(raw, n_trials, n_pos_bins, name)
% EXPAND_BEHAVIOUR  Broadcast a behavioural variable to (n_trials x n_pos_bins).
    sz = size(raw);
    if isequal(sz, [n_trials, n_pos_bins])
        x = double(raw);
    elseif isvector(raw) && numel(raw) == n_pos_bins
        % Per-bin (e.g. position index) → broadcast across trials.
        x = repmat(double(raw(:)).', n_trials, 1);
    elseif isvector(raw) && numel(raw) == n_trials
        % Per-trial (e.g. trial number) → broadcast across bins.
        x = repmat(double(raw(:)), 1, n_pos_bins);
    else
        error('v5_projection_behaviour_corr:badShape', ...
            'behaviour.%s has shape %s; expected %s or vector of %d or %d', ...
            name, mat2str(sz), mat2str([n_trials n_pos_bins]), n_pos_bins, n_trials);
    end
end


function [r, p] = corrcoef_safe(a, b)
% CORRCOEF_SAFE  Pearson r and two-sided p, robust to constant inputs.
    a = a(:); b = b(:);
    if std(a) < 1e-15 || std(b) < 1e-15
        r = 0; p = 1;
        return;
    end
    [R, P] = corrcoef(a, b);
    r = R(1, 2);
    p = P(1, 2);
end
