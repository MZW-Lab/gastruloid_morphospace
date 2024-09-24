%% Function Description
% using bwboundaries, we find the edges of each respective nuclei and find the distance between the edges of 
% the nuclei and the edge of the colony. We then take the distances from the edge as well as the
% closest edge and store it into nuc_info

%% Find Nuclei Distance from the Edge
function [mask_edges, border, distances, nuc_info] = find_nuc_2_edge(colonymask5, nuc_info, num_objects)    
    mask_edges = bwboundaries(colonymask5);                                %Finds the region for a boundary, returns row and colony coordinates of boundary pixes
    mask_edges = mask_edges{:}; border = []; distances = [];
    border(:,1) = mask_edges(:,2);
    border(:,2) = mask_edges(:,1);
    %Iterate through the list of nuclei and find distance from the edge
    for ii = 1:num_objects
        nuc_Centroid = nuc_info(ii).Centroid;
        for iii = 1:size(border,1)
            distances(iii) = sqrt(sum((border(iii,:) - nuc_Centroid) .^ 2));
        end
        [min_dist, min_dist_idx] = min(distances);
        nuc_info(ii).Dist2Edge = min_dist;
        nuc_info(ii).ClosestEdge = border(min_dist_idx,:);
    end