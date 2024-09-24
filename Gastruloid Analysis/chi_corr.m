%% A New Correlation Coefficient
%Using the method developed by Sourav Chatterjee to find the correlation
%between each of the channels.
%The method follows the following steps
%1. Given Vector X and Y, sort X from greatest to least.
%2.
%% Chi Correlation
function [G_DE_CORR50,B_DE_CORR50,S_DE_CORR50] = chi_corr(sub_struct, col_name)
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

        [G_sorted,g_ord] = sort(G_DE_S50, 'ascend');        %Sorted X for each Channel
        [B_sorted,b_ord] = sort(B_DE_S50, 'ascend');
        [S_sorted,s_ord] = sort(S_DE_S50, 'ascend');
        
        g_s_sorted = S_DE_S50(g_ord);                       %Match Indices of each other channel to sorted channel
        g_b_sorted = B_DE_S50(g_ord);
        s_b_sorted = B_DE_S50(s_ord);

        [GS_sorted,gS_ord] = sort(g_s_sorted, 'ascend');   %Find the Rank of each item in channel
        [GB_sorted,gB_ord] = sort(g_b_sorted, 'ascend');
        [SB_sorted,sB_ord] = sort(s_b_sorted, 'ascend');
        g_rank = 1:length(G_DE_S50);                        %Generate list of values
        b_rank = 1:length(B_DE_S50);                        
        s_rank = 1:length(S_DE_S50);   
        g_rank(gS_ord) = g_rank;                            %return rank by indexing the indexing values. 
        b_rank(gB_ord) = b_rank;
        s_rank(sB_ord) = s_rank;
        g_corr = sum(abs(g_rank(2:bin_count) - g_rank(1:bin_count-1)));
        s_corr = sum(abs(b_rank(2:bin_count) - b_rank(1:bin_count-1)));
        b_corr = sum(abs(s_rank(2:bin_count) - s_rank(1:bin_count-1)));
        g_corr = (1 - ((3*g_corr)/(bin_count^2-1)));
        s_corr = (1 - ((3*s_corr)/(bin_count^2-1)));
        b_corr = (1 - ((3*b_corr)/(bin_count^2-1)));
        G_DE_CORR50(i) = g_corr; 
        B_DE_CORR50(i) = s_corr; 
        S_DE_CORR50(i) = b_corr;
    end
