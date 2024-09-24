%MAIN
%2021 MZW CQ

%% Initialize
clear all; close all;
plate_id = 'E';
% Define file locations
root = append(pwd, "/",plate_id,"/");

%Define Nuclear Marker index (e.g. Dapi channel number or H2B Channel number). 
C1 = 1; %DAPI 
C2 = 2; %GATA3 561
C3 = 3; %BRA 640
C4 = 4; %SOX2 488

% Define Channel Names. No spaces or special    characters. 
Cnames = {'C2','C3','C4'};
Cnames_dapi = {'C1','C2','C3','C4'};
%Approximate the diameter of a nucleus in pixels on the short axis
nuc_diam = 10;

%Set display Flag. 1 to display. 0 to hide. 
flag = 0;
%% Set up files for processing
filenames = dir(root); %Get filenames

% filenames = filenames(~[filenames.isdir]);
% [~,idx] = sort([filenames.datenum]);
% filenames = filenames(idx);
filenames = natsortfiles(filenames);
[samplenames, drug_name, drug_str, drug_str_pdf, filenames_mod] = return_drug_str(filenames); cm5_master = {};
information_matrix = []; dapi_sheet = {}; gata_sheet = {}; bra_sheet = {}; sox_sheet = {}; dist_sheet = {}; num_cells = [];
t_C2max_cl = []; t_C3max_cl = []; t_C4max_cl = []; 
t_C2min_cl = []; t_C3min_cl = []; t_C4min_cl = [];
for i = 1:numel(samplenames)
%% Load nuclear image image
    tif_file = loadTiffStack([root+filenames_mod(i).name]);
    im_DAPI = pull_dapi(tif_file,C1);
%Plotting Function
    if flag == 1
        figure(1)
        inv_dapi = imshow(-im_DAPI,[]);
        %imwrite(inv_dapi, 'Dapi_sample.jpeg');
    end

%% Background Removal
    [colonymask, colonymask2, colonymask3, colonymask4, colonymask5, cm4_CC] = background_rm(nuc_diam, im_DAPI);
%Plotting Function
    if flag == 1
        figure(1)
        %dapi_img = imshow(colonymask5);
        dapi_img = imshowpair(colonymask, colonymask5, 'falsecolor');
        imsave(dapi_img);
    end
    
%% Do region-based segmentation
    BW = reg_seg(nuc_diam, im_DAPI, colonymask5);
%Plotting Function
    if flag == 1
        figure(1)
        dapi = imshow(-BW);
        %dapi = imshowpair(-im_DAPI, BW, 'blen');
        %imsave(dapi);  
    end
    
%% Get connected components
    [im_props, ri_CC, nuc_info, num_objects] = get_connected(BW, nuc_diam);

%% Calculate distance from edge
    [mask_edges, border, distances, nuc_info] = find_nuc_2_edge(colonymask5, nuc_info, num_objects);

%% Get imaging data
    %[col_name, sub_struct, chan_idx] = image_gen(distances, num_objects, tif_file, samplenames, Cnames, drug_name, ri_CC, i, border);
    [col_name, sub_struct, chan_idx] = image_gen_dapi_plus(distances, num_objects, tif_file, samplenames, Cnames_dapi, drug_name, ri_CC, i, border);
    
%% Generate Information Vector
%Need to put multi-value generation into a for loop
    [info_vec] = vectorgen(num_objects, col_name, sub_struct, colonymask5);
    info_vec = info_vec';
    information_matrix = [information_matrix, info_vec]; 
    
 %% Generate information of each channel
    dapi_vec = {sub_struct.(col_name).C1(:).MeanIntensity};
    gata_vec = {sub_struct.(col_name).C2(:).MeanIntensity};
    bra_vec = {sub_struct.(col_name).C3(:).MeanIntensity};
    sox_vec = {sub_struct.(col_name).C4(:).MeanIntensity};
    dapi_sheet = [dapi_sheet; {dapi_vec}];
    gata_sheet = [gata_sheet; {gata_vec}];
    bra_sheet =  [bra_sheet; {bra_vec}];
    sox_sheet = [sox_sheet; {sox_vec}];
    dist_vec = {sub_struct.(col_name).C2(:).Dist2Edge};
    dist_sheet = [dist_sheet; {dist_vec}]; 
    dist_1 = [sub_struct.(col_name).C2(:).Dist2Edge];
    num_cells = [num_cells, num_objects];
    cm5_master = [cm5_master, {colonymask5}];
    texmax_C2 = max([sub_struct.(col_name).C2(:).MeanIntensity]);
    texmax_C3 = max([sub_struct.(col_name).C3(:).MeanIntensity]);
    texmax_C4 = max([sub_struct.(col_name).C4(:).MeanIntensity]);
    texmin_C2 = min([sub_struct.(col_name).C2(:).MeanIntensity]);
    texmin_C3 = min([sub_struct.(col_name).C3(:).MeanIntensity]);
    texmin_C4 = min([sub_struct.(col_name).C4(:).MeanIntensity]);
    texmin_C2 = min([sub_struct.(col_name).C2(:).MeanIntensity]);
    texmin_C3 = min([sub_struct.(col_name).C3(:).MeanIntensity]);
    texmin_C4 = min([sub_struct.(col_name).C4(:).MeanIntensity]);
    t_C2max_cl = [t_C2max_cl, texmax_C2];
    t_C3max_cl = [t_C3max_cl, texmax_C3];
    t_C4max_cl = [t_C4max_cl, texmax_C4];
    C2_max = max(t_C2max_cl);
    C3_max = max(t_C3max_cl);
    C4_max = max(t_C4max_cl);
    t_C2min_cl = [t_C2min_cl, texmin_C2];
    t_C3min_cl = [t_C3min_cl, texmin_C3];
    t_C4min_cl = [t_C4min_cl, texmin_C4];
    C2_min = min(t_C2min_cl);
    C3_min = min(t_C3min_cl);
    C4_min = min(t_C4min_cl);
    disp([C2_max, C3_max, C4_max, C2_min, C3_min, C4_min]);
 %% Plot Results

% L = labelmatrix(ri_CC);
% RGB2 = label2rgb(L,'spring','k','noshuffle'); 
% 
% figure(2)
% subplot(1,3,1)
% 
% imshow(RGB2)
% hold on
% for ii = 1:num_objects
%     edge = nuc_info(ii).ClosestEdge;
%     cent = nuc_info(ii).Centroid;
%     hline = line([edge(1), cent(1)], [edge(2), cent(2)]); 
%     hline.LineStyle = '-';
%     hline.Color = 'c';
%     shadow = plot(border(:,1),border(:,2),'w','LineWidth',3);
%     %saveas(shadow)
% end
% hold off
% 
% ylim([0 1024])
% xlim([0 1024])
% set(gca,'XTick',[], 'YTick', [])
% 
% 
% subplot(1,3,2)
% hold on
% plot([sub_struct.(col_name).C2(:).Dist2Edge],[sub_struct.(col_name).C2(:).MeanIntensity],'o','color','b')
% plot([sub_struct.(col_name).C3(:).Dist2Edge],[sub_struct.(col_name).C3(:).MeanIntensity],'o','color','r')
% plot([sub_struct.(col_name).C4(:).Dist2Edge],[sub_struct.(col_name).C4(:).MeanIntensity],'o','color','g')
% plot([sub_struct.(col_name).C2(:).Dist2Edge],[sub_struct.(col_name).C2(:).MeanIntensity],'o','color','b')
% plot([sub_struct.(col_name).C3(:).Dist2Edge],[sub_struct.(col_name).C3(:).MeanIntensity],'o','color','r')
% plot([sub_struct.(col_name).C4(:).Dist2Edge],[sub_struct.(col_name).C4(:).MeanIntensity],'o','color','g')
% 
% xlabel('Distance to Edge')
% ylabel('Mean Intensity')
% %legend(Cnames(chan_idx))
% hold off
% 
% 
% 
% subplot(1,3,3)
% view(3) 
% c = dist_1;
% S = repmat(15,numel(dist_1),1);
% scatter3([sub_struct.(col_name).C2(:).MeanIntensity],[sub_struct.(col_name).C3(:).MeanIntensity],[sub_struct.(col_name).C4(:).MeanIntensity],S,c);
% xlabel('Mean Intensity C2')
% ylabel('Mean Intensity C3')
% zlabel("Mean Intensity C4")
% colormap(jet(size(dist_1,2)))
% hold off
% % % 
%  sgtitle(sprintf(drug_name{i}))
%  set(gcf, 'Position',  [100, 100, 1000, 300])
%  %savefig(gcf, drug_str{i});        
%  %saveas(gcf, drug_str_pdf{i})
%  close all

end
vectorgen_compiled(bra_sheet, gata_sheet, sox_sheet, dapi_sheet, dist_sheet, num_cells, cm5_master,C2_max, C3_max, C4_max, C2_min, C3_min, C4_min, plate_id)

%% Export Matrix As CSV
%writematrix(information_matrix, append('New_Corr_',plate_id,'.csv'));
%sheet_convert(bra_sheet, gata_sheet, sox_sheet, dist_sheet, plate_id);
