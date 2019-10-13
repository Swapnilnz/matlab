function [spill] = f_spill(K, L, h, h_w)
%F_SPILL Calculates Spill for a given NI or SI dam
%   Inputs K, a constant, L, length of spillway weir, and h_w, height of the
%   weir (different for NI and SI). Calculates average flow of spill for a
%   given h (height of water in dam)
if h > h_w
    spill = K.*L.*((h - h_w).^1.5);
else
    spill = 0;
end
end

