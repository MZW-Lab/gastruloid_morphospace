%% Function Description
%Perform morphological dilation and erosion with strel "structured element"
%denoted as se. Element size is crucial here and should be approx. equal to
%nuclear radius

%% Region Segmentation
function [seg_im] = reg_seg(nuc_diam, im_DAPI, bin_mask)
    se = strel('disk', nuc_diam/2);                      %Generate secondary structural element
    se_dapi = imopen(im_DAPI, se);                                %opening on binary image using structural element
    dp_binmask = se_dapi.*bin_mask;                        %Combine with binary colonymask via dot product to further enhance image
    seg_im = imregionalmax(dp_binmask );                      %Return the binary image that identifies the regional maxiam in I
