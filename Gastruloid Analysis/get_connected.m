%% Function Description
% uses regionally segemented image to generated connected components. After defining a lower and upper
% bound nuclei size, iterate through the indexes of found connected components, pull number of pixels
% and compare it to the upper and lower bounds. Finally store the nuclei that meet that criteria into list

%% Get Connected Components
function[im_props, ri_CC, nuc_info, num_objects] = get_connected(reg_im, nuc_diam)
    ri_CC = bwconncomp(reg_im);                                 %Get connected components
    im_props = regionprops(ri_CC, 'Centroid');                  %Returns measurements for the set of properties for a connected object in BW
    num_objects = 0;
    k = 0; to_remove = [];
    lowerbound = pi*(nuc_diam/2).^2/2;                      %Establish the upper and lower limits of nuclei size
    upperbound = pi*(nuc_diam/2).^2*2;
%Filter out large and small objects
    for ii = 1:numel(ri_CC.PixelIdxList)
        nucsize = numel(ri_CC.PixelIdxList{ii});
        if nucsize > lowerbound && nucsize < upperbound %Check Size (weird typo in original)
            num_objects = num_objects + 1;                            
            nuc_info(num_objects).PixelIdxList = ri_CC.PixelIdxList{ii};    %Assign Nuclei object to array
            nuc_info(num_objects).Centroid = im_props(ii).Centroid;         %Assign centroid object to array
            idxin(num_objects) = ii; 
        else
            k = k+1;
            to_remove(k) = ii;
        end
    end

%Update CC with correct number of objects so that nuclei we removed are
%not counted. 
    ri_CC.PixelIdxList(to_remove) = [];
    ri_CC.NumObjects = num_objects;