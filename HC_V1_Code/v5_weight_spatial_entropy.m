function [entropy_per_dim, entropy_null] = v5_weight_spatial_entropy( ...
        W_cells, place_tuning, n_spatial_bins, n_surrogates, seed)
%V5_WEIGHT_SPATIAL_ENTROPY  Spatial entropy of CCA weight magnitudes.
%
%   [ENTROPY, ENTROPY_NULL] = V5_WEIGHT_SPATIAL_ENTROPY(W_CELLS, ...
%       PLACE_TUNING, N_SPATIAL_BINS, N_SURROGATES) quantifies how
%   uniformly each canonical dimension's per-cell weight magnitude is
%   distributed along the track, using each cell's preferred position
%   as its location.
%
%   Following Gonzalez et al.:
%       1. For each cell, find its preferred position (argmax of place
%          tuning curve).
%       2. Bin cells by preferred position into N_SPATIAL_BINS equal-
%          width bins across the position axis.
%       3. For each canonical dimension, sum |W_cells| within each
%          spatial bin. Normalise to a probability distribution.
%       4. Shannon entropy = -sum_b p_b log(p_b).
%
%   Surrogate: circularly shift each cell's place tuning curve by a
%   random amount, recompute preferred position, recompute entropy.
%   N_SURROGATES repetitions per call.
%
%   Inputs
%       W_cells          (n_units x r) cell-space weights from
%                        v5_cell_weights.
%       place_tuning     (n_units x n_position_bins). Per-cell tuning
%                        curve (e.g. trial-averaged firing rate per
%                        position bin, or the residualised PSTH
%                        averaged over trials).
%       n_spatial_bins   integer, default 5 (paper convention).
%       n_surrogates     integer, default 100.
%       seed             RNG seed for surrogate generation. Default 0.
%
%   Outputs
%       entropy_per_dim  (r x 1) Shannon entropy (nats) per canonical
%                        dimension.
%       entropy_null     (n_surrogates x r) surrogate entropies for
%                        per-dim significance testing.
%
%   Notes
%       - Maximum entropy = log(n_spatial_bins) (uniform distribution).
%       - Minimum entropy = 0 (all mass in one bin).
%       - The Gini coefficient is a related but distinct quantity not
%         computed here; see v5_weight_gini if added.

    if nargin < 3 || isempty(n_spatial_bins), n_spatial_bins = 5;   end
    if nargin < 4 || isempty(n_surrogates),   n_surrogates   = 100; end
    if nargin < 5 || isempty(seed),           seed           = 0;   end

    if ~ismatrix(W_cells)
        error('v5_weight_spatial_entropy:badShape', 'W_cells must be 2D');
    end
    if ~ismatrix(place_tuning)
        error('v5_weight_spatial_entropy:badShape', 'place_tuning must be 2D');
    end
    [n_units, r] = size(W_cells);
    if size(place_tuning, 1) ~= n_units
        error('v5_weight_spatial_entropy:badShape', ...
            'place_tuning has %d rows; W_cells has %d', ...
            size(place_tuning, 1), n_units);
    end
    n_pos = size(place_tuning, 2);
    if n_spatial_bins < 2
        error('v5_weight_spatial_entropy:badShape', ...
            'n_spatial_bins must be >= 2; got %d', n_spatial_bins);
    end
    if n_spatial_bins > n_pos
        error('v5_weight_spatial_entropy:badShape', ...
            'n_spatial_bins=%d > n_position_bins=%d', n_spatial_bins, n_pos);
    end

    abs_W = abs(W_cells);

    % Map each cell's preferred position to one of n_spatial_bins.
    pref_pos = preferred_position(place_tuning);              % (n_units x 1)
    spatial_bin = position_to_spatial_bin(pref_pos, n_pos, n_spatial_bins);

    entropy_per_dim = compute_entropy(abs_W, spatial_bin, n_spatial_bins);

    % Surrogates: circular shift each cell's tuning curve, recompute pref.
    rng_state = rng(seed, 'twister');
    cleanup = onCleanup(@() rng(rng_state));  %#ok<NASGU>
    entropy_null = zeros(n_surrogates, r);
    shifts = randi(n_pos, n_units, n_surrogates) - 1;          % uniform in 0..n_pos-1
    for s = 1:n_surrogates
        tuning_shifted = zeros(n_units, n_pos);
        for u = 1:n_units
            tuning_shifted(u, :) = circshift(place_tuning(u, :), shifts(u, s), 2);
        end
        pref_shifted = preferred_position(tuning_shifted);
        bin_shifted = position_to_spatial_bin(pref_shifted, n_pos, n_spatial_bins);
        entropy_null(s, :) = compute_entropy(abs_W, bin_shifted, n_spatial_bins);
    end
end


% -------------------------------------------------------------------- %
function pref = preferred_position(tuning)
% PREFERRED_POSITION  Argmax along the position axis. Ties → first.
%   tuning: (n_units x n_pos). Returns column vector of indices.
    [~, pref] = max(tuning, [], 2);
    pref = pref(:);
end


function bin_idx = position_to_spatial_bin(pos_idx, n_pos, n_spatial_bins)
% POSITION_TO_SPATIAL_BIN  Map (1..n_pos) -> (1..n_spatial_bins).
%   Equal-width bins. Indices in 1..n_pos map deterministically; the last
%   spatial bin captures the upper edge.
    edges = linspace(1, n_pos + 1, n_spatial_bins + 1);
    % discretize: idx i goes into bin b iff edges(b) <= i < edges(b+1).
    % MATLAB's `discretize` returns NaN for out-of-range; pos_idx is
    % guaranteed in 1..n_pos.
    bin_idx = discretize(pos_idx, edges);
    bin_idx(isnan(bin_idx)) = n_spatial_bins;  % defensive
end


function ent = compute_entropy(abs_W, spatial_bin, n_spatial_bins)
% COMPUTE_ENTROPY  Per-dim Shannon entropy of binned |W|.
%   abs_W: (n_units x r). spatial_bin: (n_units x 1).
    r = size(abs_W, 2);
    ent = zeros(r, 1);
    for d = 1:r
        bin_mass = zeros(n_spatial_bins, 1);
        for b = 1:n_spatial_bins
            m = spatial_bin == b;
            bin_mass(b) = sum(abs_W(m, d));
        end
        total = sum(bin_mass);
        if total <= 0
            ent(d) = 0;
            continue;
        end
        p = bin_mass / total;
        p_pos = p(p > 0);
        ent(d) = -sum(p_pos .* log(p_pos));
    end
end
