function ifi = v4_ifi_from_lags(r_lags, lags)
%V4_IFI_FROM_LAGS  Information Flow Index from lagged correlations.
%   ifi = v4_ifi_from_lags(r_lags, lags)
%
%   IFI = (|r_neg| - |r_pos|) / (|r_neg| + |r_pos|), where
%     r_neg = mean(r_lags(lags < 0))   (v leads u)
%     r_pos = mean(r_lags(lags > 0))   (u leads v)
%
%   Conventions (matching v3, but bounded):
%     u leads v -> |r_pos| > |r_neg| -> IFI < 0
%     v leads u -> |r_neg| > |r_pos| -> IFI > 0
%
%   This is the canonical lag-asymmetry index used in
%   cca_animal_sanity_check.m. Taking absolute values guarantees IFI in
%   [-1, 1] regardless of the sign of r_neg or r_pos. The earlier draft
%   (without abs) could produce |IFI| > 1 when r_neg and r_pos had
%   opposite signs — see commit history / PLAN_v4_fixes.md.
%
%   Returns NaN if both sides are essentially zero correlation
%   (|r_neg| + |r_pos| < eps_tol).

    eps_tol = 1e-3;
    r_neg = abs(mean(r_lags(lags < 0), 'omitnan'));
    r_pos = abs(mean(r_lags(lags > 0), 'omitnan'));
    if isnan(r_neg) || isnan(r_pos)
        ifi = NaN;
        return;
    end
    denom = r_neg + r_pos;
    if denom > eps_tol
        ifi = (r_neg - r_pos) / denom;
    else
        ifi = NaN;
    end
end
