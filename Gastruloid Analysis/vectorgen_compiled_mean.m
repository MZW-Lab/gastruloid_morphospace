function vectorgen_compiled_mean(bra_sheet, gata_sheet, sox_sheet, dapi_sheet, dist_sheet, num_cells, cm5_master,C2_max, C3_max, C4_max, C2_min, C3_min, C4_min, plate_str)
    master_vec = []; g_corr =[]; s_corr=[]; b_corr=[];
    pixel_to_real = 0.65;    
    dapi_mean = mean((cell2mat([dapi_sheet{:}])));
%% Part 2 Description
%First using histcounts to find edges of each intensity, then use that and
%evenly spaced distance edges to find 2D histcount of each colony. Still
%need post processing
%% Binning for Underlying Distribution Via Hist, 50 Bins, With Intensity Normalization
    edge50 = linspace(0, 390, 51);
    c_area = 0.65*pi*edge50.^2;
    ring_areas = (c_area(2:end)-c_area(1:end-1));
%     c_area = pi*edge50.^2;
%     ring_areas = (c_area(2:end)-c_area(1:end-1));
%     gs = cell2mat([gata_sheet{:}]);
%     bs = cell2mat([bra_sheet{:}]);
%     ss = cell2mat([sox_sheet{:}]);
%     c2_med = median(gs, 'all');
%     c3_med = median(bs, 'all');
%     c4_med = median(ss, 'all');
    for i = 1:numel(bra_sheet)
        [cell_DE_bin_c50, cell_edges50, cell_DE_bins50] = histcounts(cell2mat(dist_sheet{i}), edge50); %#ok<ASGLU>
        col_area = sum((cm5_master{i}), 'all')*pixel_to_real^2;                                %Radius of col defined as nuclei with farthest nuclei from edge
        col_density =  num_cells(i)/col_area;
        for ii = 1:50
            indx = find(cell_DE_bins50 == ii);
            ex_C2 = cell2mat(gata_sheet{i}(indx));
            ex_C3 = cell2mat(bra_sheet{i}(indx));
            ex_C4 = cell2mat(sox_sheet{i}(indx));
            c_count = length(ex_C2);
            c_density = c_count/ring_areas(ii);
%             mean_C2 = (ex_C2-dapi_mean)/(C2_max - C2_min);
%             mean_C3 = (ex_C3-dapi_mean)/(C3_max - C3_min);
%             mean_C4 = (ex_C4-dapi_mean)/(C4_max - C4_min);
            mean_C2 = ex_C2 - dapi_mean; 
            mean_C3 = ex_C3 - dapi_mean; 
            mean_C4 = ex_C4 - dapi_mean; 
            mean_C2(mean_C2 < 0) = 0;
            mean_C3(mean_C3 < 0) = 0;
            mean_C4(mean_C4 < 0) = 0;
            c_densities(ii) = c_density;  %#ok<AGROW>
            G_DE_S50(ii) = sum(mean_C2); %#ok<AGROW>
            B_DE_S50(ii) = sum(mean_C3); %#ok<AGROW>
            S_DE_S50(ii) = sum(mean_C4); %#ok<AGROW>
            %[gata_corr,sox_corr,bra_corr] = chi_corr_compiled(ex_C2, ex_C3, ex_C4);
            [gata_corr,sox_corr,bra_corr] = standard_corr_compiled(ex_C2, ex_C3, ex_C4);
            g_corr(ii) = gata_corr; %#ok<AGROW>
            s_corr(ii) = sox_corr; %#ok<AGROW>
            b_corr(ii) = bra_corr; %#ok<AGROW>
        end
        G_MI50 = G_DE_S50./cell_DE_bin_c50;
        B_MI50 = B_DE_S50./cell_DE_bin_c50;
        S_MI50 = S_DE_S50./cell_DE_bin_c50;
        current_vec = [num_cells(i),col_density,G_MI50,B_MI50,S_MI50,g_corr,s_corr,b_corr,c_densities];
    master_vec = [master_vec; current_vec]; %#ok<AGROW>
    end
writematrix(master_vec, append('vector_info_subdapi_denP', plate_str, '.csv'));    