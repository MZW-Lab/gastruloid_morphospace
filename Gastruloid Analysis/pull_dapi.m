%% Function Description
%Returns the image in the DAPI channel

%% Pull DAPI Image
function [im_DAPI] = pull_dapi(loaded_im, C1)
    im_DAPI = loaded_im(:,:,C1); %Invert so nuclei are darker
    im_DAPI = rescale(im_DAPI);