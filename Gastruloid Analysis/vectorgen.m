function info_vec = vectorgen(num_objects, col_name, sub_struct, cm5)
%% Part 1 Description
%Find the colony farthest from edge, set as radius, and find the density of
%colony, robust to account for different colony sizes
%% Finding colony density    
    pixel_to_real = 0.65;
    col_area = sum(cm5, 'all')*pixel_to_real^2;                                %Radius of col defined as nuclei with farthest nuclei from edge
    col_density = num_objects/col_area;
%% Part 2 Description
%First using histcounts to find edges of each intensity, then use that and
%evenly spaced distance edges to find 2D histcount of each colony. Still
%need post processing
%% Data Establishment
tex_C2 = [sub_struct.(col_name).C2(:).MeanIntensity];
tex_C3 = [sub_struct.(col_name).C3(:).MeanIntensity];
tex_C4 = [sub_struct.(col_name).C4(:).MeanIntensity];
%% Binning for Underlying Distribution Via Hist, 25 Bins
    edge25 = linspace(0, 390, 26);               %Array for edges for a bin size of 25
    cell_DE = [sub_struct.(col_name).C2(:).Dist2Edge];
    G_DE_S25 = []; B_DE_S25 = []; S_DE_S25 = [];
    [cell_DE_bin_c25, cell_edges25, cell_DE_bins25] = histcounts(cell_DE, edge25);
%     for i = 1:50
%         indx = find(cell_DE_bins25 == i);
%         G_DE_S25(i) = sum([sub_struct.(col_name).C2(indx).MeanIntensity]);
%         B_DE_S25(i) = sum([sub_struct.(col_name).C3(indx).MeanIntensity]);
%         S_DE_S25(i) = sum([sub_struct.(col_name).C4(indx).MeanIntensity]);
%     end
%     G_MI50_um = G_DE_S25./cell_DE_bin_c25;
%     B_MI50_um = B_DE_S25./cell_DE_bin_c25;
%     S_MI50_um = S_DE_S25./cell_DE_bin_c25;
%% Binning for Underlying Distribution Via Hist, 50 Bins, With Intensity Normalization
    edge50 = linspace(0, 390, 51);
    [cell_DE_bin_c50, cell_edges50, cell_DE_bins50] = histcounts(cell_DE, edge50);
    for i = 1:50
        indx = find(cell_DE_bins50 == i);
        ex_C2 = [sub_struct.(col_name).C2(indx).MeanIntensity];
    	ex_C3 = [sub_struct.(col_name).C3(indx).MeanIntensity];
        ex_C4 = [sub_struct.(col_name).C4(indx).MeanIntensity];
        mean_C2 = (ex_C2-min(tex_C2))/(max(tex_C2)-min(tex_C2));
        mean_C3 = (ex_C3-min(tex_C3))/(max(tex_C3)-min(tex_C3));
    	mean_C4 = (ex_C4-min(tex_C4))/(max(tex_C4)-min(tex_C4));
        G_DE_S50(i) = sum(mean_C2); %#ok<AGROW>
        B_DE_S50(i) = sum(mean_C3); %#ok<AGROW>
        S_DE_S50(i) = sum(mean_C4); %#ok<AGROW>
    end
    G_MI50 = G_DE_S50./cell_DE_bin_c50;
    B_MI50 = B_DE_S50./cell_DE_bin_c50;
    S_MI50 = S_DE_S50./cell_DE_bin_c50;
    [g_corr,s_corr,b_corr] = chi_corr(sub_struct, col_name);
%% Binning for Underyling Distribution Via Hist, 25 Bins    
%     edge25 = linspace(0, 390, 26);
%     [cell_DE_bin_c25, cell_edges25, cell_DE_bins2] = histcounts(cell_DE, edge25);
%     for i = 1:25
%         indx = find(cell_DE_bins2 == i);
%         G_DE_S25(i) = sum([sub_struct.(col_name).C2(indx).MeanIntensity]);
%         B_DE_S25(i) = sum([sub_struct.(col_name).C3(indx).MeanIntensity]);
%         S_DE_S25(i) = sum([sub_struct.(col_name).C4(indx).MeanIntensity]);
%     end
%     G_MI25 = G_DE_S25./cell_DE_bin_c25;
%     B_MI25 = B_DE_S25./cell_DE_bin_c25;
%     S_MI25 = S_DE_S25./cell_DE_bin_c25;
%     [g_corr25,s_corr25,b_corr25] = chi_corr(sub_struct, col_name);
%% Final vector
%     info_vec = [num_objects,col_density, G_MI50, B_MI50, S_MI50, G_MI25, B_MI25, S_MI25, g_corr,s_corr,b_corr,g_corr25,s_corr25,b_corr25];
    info_vec = [num_objects,col_density, G_MI50, B_MI50, S_MI50, g_corr, s_corr, b_corr];