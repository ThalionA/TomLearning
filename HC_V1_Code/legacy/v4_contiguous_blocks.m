function blocks = v4_contiguous_blocks(is_valid, min_block)
%V4_CONTIGUOUS_BLOCKS  Find contiguous runs of true with length >= min_block.
%   blocks = v4_contiguous_blocks(is_valid, min_block)
%
%   is_valid : logical column or row vector
%   min_block: scalar integer, minimum block length to keep (in bins)
%
%   Returns blocks as an Nx2 matrix of [start, stop] indices (1-based,
%   inclusive). Empty (0x2) if no run satisfies the length filter.

    is_valid = is_valid(:);
    n = numel(is_valid);
    blocks = zeros(0, 2);
    if n == 0, return; end

    % Pad with false, find rising/falling edges via diff.
    padded = [false; is_valid; false];
    d = diff(double(padded));
    starts = find(d ==  1);   % indices into is_valid, 1-based
    stops  = find(d == -1) - 1;
    lens   = stops - starts + 1;
    keep   = lens >= min_block;
    blocks = [starts(keep), stops(keep)];
end
