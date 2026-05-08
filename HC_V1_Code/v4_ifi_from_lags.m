function ifi = v4_ifi_from_lags(r_lags, lags)
%V4_IFI_FROM_LAGS  Information Flow Index from lagged correlations.
%   ifi = v4_ifi_from_lags(r_lags, lags)
%
%   IFI = (r_neg - r_pos) / (r_neg + r_pos), where
%     r_neg = mean(r_lags(lags < 0))   (v leads u)
%     r_pos = mean(r_lags(lags > 0))   (u leads v)
%
%   Conventions (matching v3):
%     u leads v -> r_pos > r_neg -> IFI < 0
%     v leads u -> r_neg > r_pos -> IFI > 0
%
%   Returns NaN if (r_neg + r_pos) is too small to give a stable ratio.

    r_neg = mean(r_lags(lags < 0), 'omitnan');
    r_pos = mean(r_lags(lags > 0), 'omitnan');
    if (r_neg + r_pos) > 0.001
        ifi = (r_neg - r_pos) / (r_neg + r_pos);
    else
        ifi = NaN;
    end
end
