%% Function Description
% using the previously generated location data, distance to edge data, we propogate the function into each
% of the image channels and store them into each channels information space
%% Multichannel Image Generation
function [col_name, sub_struct, chan_idx] = image_gen_dapi_plus(distances, num_objects,tif_file, samplenames, Cnames, drug_name, ri_CC, i, border)       
    C1 = 1;
    C2 = 2;
    C3 = 3;
    C4 = 4;
    channels = {C1,C2,C3,C4};
    sub_struct = struct;
    chan_idx = cellfun(@(a) a ~= 0, channels);                                 %Generate Variable name for each Channel
    chan_var = channels(chan_idx);

    fprintf('\nNuclei Identified');
    %Iterate through each channel and run distance formula for nuclei
    for ii = 1:sum(chan_idx)
        fprintf(['\nProcessing sample: ', samplenames{i}, ', Channel: ', Cnames{ii}]);
        im_CurChan = tif_file(:,:,chan_var{ii});                                %Pulls image of specific channel chan_var
        im_props = regionprops(ri_CC,im_CurChan, 'Centroid', 'MeanIntensity');
        for iii = 1:num_objects
            %disp(sprintf("\n%d percent completed",round(100*iii/NumObjects)))
            nuc_centroids = im_props(iii).Centroid;                         %Setting Centroid
            for iiii = 1:size(border,1)
                distances(iiii) = sqrt(sum((border(iiii,:) - nuc_centroids) .^ 2));   %Distance from Centroid
            end
            [dist_min,nuc_min_dist_idx] = min(distances);
            im_props(iii).Dist2Edge = dist_min;
            im_props(iii).ClosestEdge = border(nuc_min_dist_idx,:);
        end
        col_name = [genvarname(drug_name{i})];      %Generates variable name from string
        sub_struct.(col_name).(Cnames{ii}) = im_props;                   %Setting S varname to regionprops

    end