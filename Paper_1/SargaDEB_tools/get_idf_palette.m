function cmap = get_idf_palette(n)
% Base 5 RGB colors
    base = [
         0,  40,  60;     % dark navy blue
        0, 113, 145 ; % 007191
        98, 200, 211;    % 62C8D3
        244, 122,  0;     % F47A00
        211,  31,  17   % D31F11      
         56,   2,  59; % 38023B
    ] / 255;

    n_base = size(base, 1);

    if n <= n_base
        cmap = base(1:n, :);  % just subset if ≤5
        return;
    end

    % We want n total colors, placed evenly along the segments between the 5 base ones
    % That means: interpolate over the whole base gradient (which has 4 segments)
    cmap = zeros(n, 3);
    x_base = linspace(0, 1, n_base);
    x_interp = linspace(0, 1, n);  % full length

    for ch = 1:3  % R, G, B channels
        cmap(:, ch) = interp1(x_base, base(:, ch), x_interp, 'pchip');
    end
end