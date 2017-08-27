function [idx1, idx0, nr1, riv_msk_fil] = riv_mask_fil(riv_msk)
% function [idx1, idx0, nr1, riv_msk_fil] = riv_mask_fil(riv_msk)

%create a filled channel
riv_msk_fil = imfill(riv_msk, 'holes');
%index outside channel (0) and inside channel (1) and number of channel
%pixels
idx1 = find(riv_msk_fil == 1); idx0 = find(riv_msk_fil == 0); nr1 = length(idx1);