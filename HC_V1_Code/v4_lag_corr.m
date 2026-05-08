function r_lags = v4_lag_corr(u, v, lags)
%V4_LAG_CORR  In-window lagged correlation of two 1D signals.
%   r_lags = v4_lag_corr(u, v, lags)
%
%   For each L in lags:
%     L >= 0:  corr(u(1:end-L),     v(1+L:end))    -> u leads v at +L
%     L <  0:  corr(u(1-L:end),     v(1:end+L))    -> v leads u at -L
%
%   Returns NaN at any lag where the surviving overlap has <= 5 samples,
%   or where either segment has zero variance.
%
%   This is the v3 local helper promoted to a standalone file so v4 and
%   the v4 tests can share it.

    u = u(:); v = v(:);
    n = numel(u);
    r_lags = nan(1, numel(lags));
    for il = 1:numel(lags)
        L = lags(il);
        if L >= 0
            a = u(1:(n - L));
            b = v((1 + L):n);
        else
            a = u((1 - L):n);
            b = v(1:(n + L));
        end
        if numel(a) > 5 && std(a) > 0 && std(b) > 0
            r_lags(il) = corr(a, b, 'rows', 'complete');
        end
    end
end
