%% Function Description
%Takes fill path and pulls tiff stack of each channel that has been Z-stack Maxed

%% Load Tiff Stack
function [ TiffStack rows cols slices] = loadTiffStack( ImFilePath )
%mwLoadTiffStack
%   Takes in file path and returns 3D tiff stack
    
    imInfo = imfinfo(ImFilePath);
    TiffStack = zeros(imInfo(1).Width, imInfo(1).Height, numel(imInfo));
    for k = 1:numel(imInfo)
        IM = imread(ImFilePath, k, 'Info', imInfo); %read in the entire tiff
        TiffStack(:,:,k) = IM;
    end
    [rows cols slices] = size(TiffStack);
end

