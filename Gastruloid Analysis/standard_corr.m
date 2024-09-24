%% A New Correlation Coefficient
%Using the method developed by Sourav Chatterjee to find the correlation
%between each of the channels.
%The method follows the following steps
%1. Given Vector X and Y, sort X from greatest to least.
%2.
%% Chi Correlation
function [G_DE_CORR50,B_DE_CORR50,S_DE_CORR50] = standard_corr(sub_struct, col_name)
    edge50 = linspace(0, 390, 51);                          %Array for edges for a bin size of 25
    cell_DE = [sub_struct.(col_name).C2(:).Dist2Edge];      %Pulling the Distance from Edge for each individual Cell
    G_DE_CORR50 = []; B_DE_CORR50 = []; S_DE_CORR50 = [];
    [cell_DE_bin_c50, cell_edges50, cell_DE_bins50] = histcounts(cell_DE, edge50);
    for i = 1:50
        indx = find(cell_DE_bins50 == i);                   %Find indices of each cell in each respective bind
        G_DE_S50 = [sub_struct.(col_name).C2(indx).MeanIntensity];  %Using index to find intensity, returns array
        B_DE_S50 = [sub_struct.(col_name).C3(indx).MeanIntensity];
        S_DE_S50 = [sub_struct.(col_name).C4(indx).MeanIntensity];
        bin_count = length(G_DE_S50);                       %Get Length of each index
        G_B_corr = corrcoef(G_DE_S50, B_DE_S50);
        G_S_corr = corrcoef(G_DE_S50, S_DE_S50);
        S_B_corr = corrcoef(S_DE_S50, B_DE_S50);
        G_DE_CORR50(i) = G_B_corr(0, 1); 
        B_DE_CORR50(i) = G_S_corr(0, 1); 
        S_DE_CORR50(i) = S_B_corr(0, 1);
    end
