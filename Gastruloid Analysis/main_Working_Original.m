%MAIN
%2021 MZW CQ

%% Initialize
clear all; close all;
% Define file locations
root = pwd + "/Plate Tifs/Plate A/";

%Define Nuclear Marker index (e.g. Dapi channel number or H2B Channel number). 
C1 = 1;
C2 = 2;
C3 = 3;
C4 = 4;

% Define Channel Names. No spaces or special characters. 
Cnames = {'C2','C3','C4'};

%Approximate the diameter of a nucleus in pixels on the short axis
nuc_diam = 10;

%Set display Flag. 1 to display. 0 to hide. 
flag = 1;
%% Set up files for processing
filenames = dir(root); %Get filenames

filenames = filenames(arrayfun(@(x) x.name(1), filenames) ~= '.'); %Remove '.','..' from array call
%splits the names of the files for access
for ii = 1:numel(filenames)
    tmp = strsplit(filenames(ii).name, {'_','.tif'});
    samplenames{ii} = tmp{end-1};
end
%%
for i = 1:numel(samplenames)
%% Load nuclear image image
    ifile = loadTiffStack([root+filenames(i).name]);
    im = ifile(:,:,C1); %Invert so nuclei are darker
    im = rescale(im);
    %Plotting Function
    if flag == 1
        figure(1)
        imshow(-im,[])
    end

%% Background Removal
    se = strel('disk', nuc_diam)                        %Generates a structured disk element in a disk shape
    colonymask = imbinarize(im, 'global');              %computes mean intensity around a pixel and takes threshold
    colonymask2 = imfill(colonymask, 'holes');          %Fills the holes in image
    colonymask3 = imdilate(colonymask2, se);            %Dilates image using structural element
    colonymask4 = imfill(colonymask3, 'holes');         %further improves image graininess via filling
    colCC = bwconncomp(colonymask4);                    %Finds and counts connected components (i.e. nuclei)
    %Count all the connected elements
    colIDX = 1;
    for ii = 1:(numel(colCC.PixelIdxList)-1)
        if numel(colCC.PixelIdxList{colIDX}) < numel(colCC.PixelIdxList{ii+1})
            colIDX = ii+1;
        end
    end

    colonymask5 = zeros(colCC.ImageSize);               %Create empty mask
    colonymask5(colCC.PixelIdxList{colIDX})=1;          %Fill empty mask in indices = 1
    %Plotting Function
    if flag == 1
        figure(1)
        imshowpair(colonymask, colonymask5, 'falsecolor')
    end

%% Do region-based segmentation
%Perform morphological dilation and erosion with strel "structured element"
%denoted as se. Element size is crucial here and should be approx. equal to
%nuclear radius 
    se = strel('disk', nuc_diam/2);                      %Generate secondary structural element
    Io = imopen(im, se);                                %opening on binary image using structural element
    Io_colony = Io.*colonymask5;                        %Combine with binary colonymask via dot product to further enhance image
    BW = imregionalmax(Io_colony);                      %Return the binary image that identifies the regional maxiam in I

%Plotting Function
    if flag == 1
        figure(1)
        imshowpair(-im, BW, 'blen')
    end

%% Get connected components
    %Generate a list of "nuclei", then measure each and remove outliers
    CC = bwconncomp(BW);                                %Get connected components
    STATS = regionprops(CC, 'Centroid');                %Returns measurements for the set of properties for a connected object in BW
    NumObjects = 0;
    k = 0; to_remove = [];
    lowerbound = pi*(nuc_diam/2).^2/2;                   %Establish the upper and lower limits of nuclei size
    upperbound = pi*(nuc_diam/2).^2*2;
%Filter out large and small objects
    for ii = 1:numel(CC.PixelIdxList)
        nucsize = numel(CC.PixelIdxList{ii});
        if nucsize > lowerbound && nucsize < upperbound %Check Size (weird typo in original)
            NumObjects = NumObjects + 1;                            
            NUCs(NumObjects).PixelIdxList = CC.PixelIdxList{ii};    %Assign Nuclei object to array
            NUCs(NumObjects).Centroid = STATS(ii).Centroid;         %Assign centroid object to array
            idxin(NumObjects) = ii; 
        else
            k = k+1;
            to_remove(k) = ii;
        end
    end

%Update CC with correct number of objects so that nuclei we removed are
%not counted. 
    CC.PixelIdxList(to_remove) = [];
    CC.NumObjects = NumObjects;

%% Calculate distance from edge
    Btmp = bwboundaries(colonymask5);                                %Finds the region for a boundary, returns row and colony coordinates of boundary pixes
    Btmp = Btmp{:}; B = []; D = [];
    B(:,1) = Btmp(:,2);
    B(:,2) = Btmp(:,1);
    %Iterate through the list of nuclei and find distance from the edge
    for ii = 1:NumObjects
        iCentroid = NUCs(ii).Centroid;
        for iii = 1:size(B,1)
            D(iii) = sqrt(sum((B(iii,:) - iCentroid) .^ 2));
        end
        [M,idx] = min(D);
        NUCs(ii).Dist2Edge = M;
        NUCs(ii).ClosestEdge = B(idx,:);
    end
%% Get imaging data
    C = {C2,C3,C4};
    S = struct;
    cidx = cellfun(@(a) a ~= 0, C);
    Csub = C(cidx);

    fprintf('Nuclei Identified')

    for ii = 1:sum(cidx)
        fprintf(['Processing sample: ', samplenames{i}, ', Channel: ', Cnames{ii}])
        im = ifile(:,:,Csub{ii});
        STATS = regionprops(CC,im, 'Centroid', 'MeanIntensity');
        for iii = 1:NumObjects
            disp(sprintf("%d percent completed",round(100*iii/NumObjects)))
            iCentroid = STATS(iii).Centroid;                         %Setting Centroid
            for iiii = 1:size(B,1)
                D(iiii) = sqrt(sum((B(iiii,:) - iCentroid) .^ 2));   %Distance from Centroid
        end
        [M,idx] = min(D);
        STATS(iii).Dist2Edge = M; 
        STATS(iii).ClosestEdge = B(idx,:);
    end
    varname = ['col_' genvarname(samplenames{i})];      %Generates variable name from string
    S.(varname).(Cnames{ii}) = STATS;                   %Setting S varname to regionprops

end
%% Plot Results

L = labelmatrix(CC);
RGB2 = label2rgb(L,'spring','k','noshuffle'); 

figure(2)
subplot(1,2,1)

imshow(RGB2)
hold on
for ii = 1:NumObjects
    edge = NUCs(ii).ClosestEdge;
    cent = NUCs(ii).Centroid;
    hline = line([edge(1), cent(1)], [edge(2), cent(2)]);
    hline.LineStyle = '-';
    hline.Color = 'c';
    plot(B(:,1),B(:,2),'w','LineWidth',3)
end
hold off

ylim([0 1024])
xlim([0 1024])
set(gca,'XTick',[], 'YTick', [])


subplot(1,2,2)
hold on
plot([S.(varname).C2(:).Dist2Edge],[S.(varname).C2(:).MeanIntensity],'o')
plot([S.(varname).C3(:).Dist2Edge],[S.(varname).C3(:).MeanIntensity],'o')
plot([S.(varname).C4(:).Dist2Edge],[S.(varname).C4(:).MeanIntensity],'o')
legend(Cnames(cidx))
hold off

sgtitle(sprintf('Colony %s', varname))
set(gcf, 'Position',  [100, 100, 800, 300])

%savefig(gcf, sprintf('Colony %s.fig', varname))
close all

end
save('data.mat','S')