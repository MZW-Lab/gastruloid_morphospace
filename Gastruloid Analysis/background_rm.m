%% Function Description
%using structural element to first identify elements followed by 
%connected components finder to general a mask and remove background

%% Remove Background
function [cm1, cm2, cm3, cm4, cm5, cm4_CC] = background_rm(nuc_diam, im_DAPI, erode_diam)
    se = strel('disk', nuc_diam);                   %Generates a structured disk element in a disk shape
    se_er = strel('disk', erode_diam);
    cm1 = imbinarize(im_DAPI, 'global');            %computes mean intensity around a pixel and takes threshold
    cm2 = imfill(cm1, 'holes');                     %Fills the holes in image
    cm2 = imerode(cm2, se_er);
    cm3 = imdilate(cm2, se);                        %Dilates image using structural element
    cm4 = imfill(cm3, 'holes');                     %further improves image graininess via filling
    cm4_CC = bwconncomp(cm4);                       %Finds and counts connected components (i.e. nuclei)
    %Count all the connected elements
    colIDX = 1;
    for ii = 1:(numel(cm4_CC.PixelIdxList)-1)
        if numel(cm4_CC.PixelIdxList{colIDX}) < numel(cm4_CC.PixelIdxList{ii+1})
            colIDX = ii+1;
        end
    end

    cm5 = zeros(cm4_CC.ImageSize);                  %Create empty mask
    cm5(cm4_CC.PixelIdxList{colIDX})=1;             %Fill empty mask in indices = 1