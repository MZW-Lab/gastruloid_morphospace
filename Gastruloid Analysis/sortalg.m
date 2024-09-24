%MAIN
%2021 MZW CQ

%% Initialize
clear all; close all;
plate_id = 'July Test';
% Define file locations
root = append(pwd, "/",plate_id,"/");

%Define Nuclear Marker index (e.g. Dapi channel number or H2B Channel number). 
C1 = 1; %DAPI
C2 = 2; %GATA3
C3 = 3; %BRA
C4 = 4; %SOX2

% Define Channel Names. No spaces or special characters. 
Cnames = {'C2','C3','C4'};

%Approximate the diameter of a nucleus in pixels on the short axis
nuc_diam = 10;

%Set display Flag. 1 to display. 0 to hide. 
flag = 0;
%% Set up files for processing
filenames = dir(root); %Get filenames
filenames = natsortfiles(filenames);
[samplenames, drug_name, drug_str, drug_str_pdf, filenames_mod] = return_drug_str(filenames);
information_matrix = []; gata_sheet = {}; bra_sheet = {}; sox_sheet = {}; dist_sheet = {};
for i = 1:numel(samplenames)
%% Load nuclear image image
    tif_file = loadTiffStack([root+filenames_mod(i).name]);
    im_DAPI = pull_dapi(tif_file,C1);
%Plotting Function
    if flag == 1
        figure(1)
        imshow(-im_DAPI,[])
    end

%% Background Removal
    [colonymask, colonymask2, colonymask3, colonymask4, colonymask5, cm4_CC] = background_rm(nuc_diam, im_DAPI);
%Plotting Function
    if flag == 1
        figure(1)
        imshowpair(colonymask, colonymask5, 'falsecolor')
    end
    
%% Do region-based segmentation
    BW = reg_seg(nuc_diam, im_DAPI, colonymask5);
%Plotting Function
    if flag == 1
        figure(1)
        imshowpair(-im_DAPI, BW, 'blen')
    end
    
%% Get connected components
    [im_props, ri_CC, nuc_info, num_objects] = get_connected(BW, nuc_diam);

%% Calculate distance from edge
    [mask_edges, border, distances, nuc_info] = find_nuc_2_edge(colonymask5, nuc_info, num_objects);

%% Get imaging data
    [col_name, sub_struct, chan_idx] = image_gen(distances, num_objects, tif_file, samplenames, Cnames, drug_name, ri_CC, i, border);
    
%% Generate Information Vector
%Need to put multi-value generation into a for loop
    [info_vec] = vectorgen(num_objects, col_name, sub_struct, colonymask5);
    info_vec = info_vec';
    information_matrix = [information_matrix, info_vec]; 
    
 %% Generate information of each channel
    gata_vec = {sub_struct.(col_name).C2(:).MeanIntensity};
    bra_vec = {sub_struct.(col_name).C3(:).MeanIntensity};
    sox_vec = {sub_struct.(col_name).C4(:).MeanIntensity};
    gata_sheet = [gata_sheet; gata_vec{:}];
    bra_sheet =  [bra_sheet; bra_vec{:}];
    sox_sheet = [sox_sheet; sox_vec{:}];
    dist_vec = {sub_struct.(col_name).C2(:).Dist2Edge};
    dist_sheet = [dist_sheet; {dist_vec}];
    dist_1 = [sub_struct.(col_name).C2(:).Dist2Edge];

end