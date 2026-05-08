function [X_seg, Y_seg, ok, info] = v4_rz_align_pre(X_trial, Y_trial, pos_trial, is_valid_trial, rz_entry_cm, n_rz_bins)
%V4_RZ_ALIGN_PRE  Extract the n_rz_bins immediately *before* RZ entry.
%   [X_seg, Y_seg, ok, info] = v4_rz_align_pre(X_trial, Y_trial, pos_trial, is_valid_trial, rz_entry_cm, n_rz_bins)
%
%   Locates the first bin where pos_trial >= rz_entry_cm (call it i_rz),
%   and returns the segment (i_rz - n_rz_bins) : (i_rz - 1).
%   Segment ends one bin BEFORE RZ entry (we don't want the entry bin
%   itself, since the animal stops to consume reward).
%
%   The segment is rejected (ok=false, X_seg=Y_seg=[]) if:
%     - position never crosses rz_entry_cm
%     - the segment would extend before bin 1 (i_rz <= n_rz_bins)
%     - any bin in the segment is invalid (vel-drop or non-cued)
%
%   info.i_rz : index of RZ entry bin (NaN if no crossing)
%   info.reason : 'ok' / 'no_crossing' / 'too_early' / 'invalid_in_segment'

    info = struct('i_rz', NaN, 'reason', '');
    X_seg = []; Y_seg = []; ok = false;

    above = find(pos_trial(:) >= rz_entry_cm, 1, 'first');
    if isempty(above)
        info.reason = 'no_crossing'; return;
    end
    info.i_rz = above;
    s = above - n_rz_bins;
    if s < 1
        info.reason = 'too_early'; return;
    end
    seg = s:(above - 1);
    if ~all(is_valid_trial(seg))
        info.reason = 'invalid_in_segment'; return;
    end
    X_seg = X_trial(seg, :);
    Y_seg = Y_trial(seg, :);
    ok = true; info.reason = 'ok';
end
