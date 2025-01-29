clc
clear
close all
%%

global analysisParam;

for PlateNum = 1:analysisParam.NumofPlates
    fprintf(['***********************','\n'])
    fprintf(['Plate:', num2str(PlateNum),'\n'])
    fprintf(['***********************','\n'])
    
for WellNumber = analysisParam.WellsWithData{PlateNum}
    fprintf(['-----------------------','\n'])
    fprintf(['Well:', num2str(WellNumber),'\n'])
    fprintf(['-----------------------','\n'])
for nposition = 1:analysisParam.ImagesperWell
    fprintf(['Position:', num2str(nposition),'\n'])
    fprintf(['-----------------------','\n'])

positionname=['P',num2str(PlateNum),'_','W',num2str(WellNumber,'%03.f'),'_',num2str(nposition),'_MAXProj'];
MAXProj_name = [analysisParam.pathnamesave filesep 'MAXProj/', positionname,'.tif'];
DAPIMask_name = [analysisParam.pathnamesave filesep 'MAXProj_DAPI_ilastik/', positionname,'_Simple Segmentation.h5'];

% substract background
if ~isfield(analysisParam,'bgvalues')
    IP_ComputeBGSubstractionlevelsUsingSegmentation
end 

% Load the 4-channel image
filereader = bfGetReader(MAXProj_name);
nx = filereader.getSizeX; % x
ny = filereader.getSizeY; % y
nt = filereader.getSizeT; % # of time points
nz = filereader.getSizeZ; % # of z slices
nc = filereader.getSizeC; % # of channels
% Separate the channels and background subtraction
for cc = 1:analysisParam.ChannelMaxNum{PlateNum}(WellNumber)
    imauxfluorescencelevels = double(bfGetPlaneAtZCT(filereader,1,cc,1));
    FluorescenceLevelsRawBGnorm(:,:,cc) = imauxfluorescencelevels - analysisParam.bgvalues(analysisParam.MapChannels.ChannelsCoordMatrix{PlateNum,WellNumber}(cc));
end
% Combine channels 1-3 into an RGB image
combinedImage = cat(3, mat2gray(FluorescenceLevelsRawBGnorm(:,:,2)), mat2gray(FluorescenceLevelsRawBGnorm(:,:,3)), mat2gray(FluorescenceLevelsRawBGnorm(:,:,1)));

% Read segmentation data (squeeze removes excessive dimensions) and turn everything that was red into 1 in binary mode
% Segmentation definition, 1: nucleus, 2: background
foregroundLabel =1;
datasetname = '/exported_data';% Path to Ilastiks exported files:
nucleiMask = squeeze(h5read(DAPIMask_name,datasetname)) == foregroundLabel;
nucleiMask=nucleiMask';
nuclei1=nucleiMask;%saved to plot

% Nuclei mask improvement
nucleiMask = bwareaopen(nucleiMask,100);%original(20X)
nucleiMask = imfill(nucleiMask, 'holes');
nuclei2 = nucleiMask; %saved to plot
%imshow(cat(3,nuclei1(:,:,1),nucleiMask,0*nuclei1(:,:,1)),[]) 

% use imopen and imclose with the structuring element se to 
% make edges smooth. It mafikes nuclei round Chose a disc of radius 10.
se = strel('disk',analysisParam.imopendiskradious);
nucleiMask = imopen(nucleiMask, se); 
nuclei3 = nucleiMask; %saved to plot

% Erode and dilate to just divide mask that joins two cells (however, this affects those cells that have moonshape, making them circular) 
seerode = strel('disk',analysisParam.imerodediskradious);
nucleiMask=imerode(nucleiMask,seerode);
sedil = strel('disk',analysisParam.imdilatediskradious);
nucleiMask=imdilate(nucleiMask,sedil);
nuclei4 = nucleiMask; %saved to plot

% Option to use imclose or imclearboarders
se = strel('disk',1);
nucleiMask = imclose(nucleiMask, se);

% Watershed algorithm
nucleiMask = imclearborder(nucleiMask); %Removes nuclei in the border of the image, which can lead to errors in the watershed algorithm
CC=bwconncomp(nucleiMask); %Find the connected components in the segmentation
statsnew = regionprops(CC,'Area','Centroid'); %Computes the area and centroids of the connected components
area = [statsnew.Area];
fusedcandidates = area > mean(area)+std(area); %Fused candidates are such that the area is bigger than the mean+std of the areas
statscell = struct2cell(statsnew);
CentroidsMat = cat(1,statscell{2,:});

sublist = CC.PixelIdxList(fusedcandidates);
sublist = cat(1,sublist{:});
fusedMask = false(size(nucleiMask));
fusedMask(sublist) = 1;

s = round(analysisParam.pwatershed*sqrt(mean(area))/pi);

if any(fusedMask,'all')
nucmin = imerode(fusedMask,strel('disk',s));
%imshow(cat(3,fusedMask,nucmin,0*fusedMask),[]);
outside = ~imdilate(fusedMask,strel('disk',1));
%imshow(outside)
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin,nucmin|outside);
pcolor(basin); shading flat;
L=watershed(basin);
%imshow(L); colormap('jet');caxis([0,20]);
newNuclearMask = L>1 | (nucleiMask - fusedMask);
%imshow(newNuclearMask)

%f=figure('Position',[100 100 1200 800]);

RawIMAGEplot=imadjust(mat2gray(FluorescenceLevelsRawBGnorm(:,:,1)));
masknew=newNuclearMask-imerode(newNuclearMask,strel('disk',2));
mask=nucleiMask(:,:,1)-imerode(nucleiMask(:,:,1)',strel('disk',2));

% subplot(1,3,1)
% imshow(cat(3,mask,RawIMAGEplot,mask),[])
% title('mask')
% 
% subplot(1,3,2)
% imshow(cat(3,masknew,RawIMAGEplot,masknew),[])
% title('masknew')
% 
% subplot(1,3,3)
% imshow(cat(3,newNuclearMask,RawIMAGEplot,masknew),[])
% title('newNuclearMask')
close all
% nuclei mask finished with newNuclearMask

%% Find background for each channel, substracting the brightspots resulting from bad immunostaining  
se = strel('disk',20);

opbackgroundCell = imdilate(newNuclearMask(:,:,1), se);
backgroundCell = not(opbackgroundCell);
IndbackDAPI = find(backgroundCell);

FluorescenceDAPI = FluorescenceLevelsRawBGnorm(:,:,1);
backgroundintDAPI = mean(FluorescenceDAPI(IndbackDAPI));
BackgroundData = struct('DAPI',backgroundintDAPI);
end
% % Display the combined RGB image
% figure;
% imshow(combinedImage);
% hold on;
% 
% % Overlay the nuclear mask in red
% maskOverlay = imoverlay(combinedImage, nucleiMask, [1 0 0]);
% 
% % Display the combined image with mask overlay
% imshow(maskOverlay);
% title('Combined Image with Nuclear Mask Overlay');
% hold off;

%% Radius profiling
% Parameters
numDirections = 6; % Number of sampling directions (4-6)
topPercentage = 0.2; % Top 20% for max intensity calculation
maxDistance = 25; % Maximum distance outward from nucleus

% Get properties of labeled regions (nuclei)
stats = regionprops(newNuclearMask, 'Centroid','MajorAxisLength','Orientation','Area','PixelIdxList');

% Create a meshgrid for distance calculations
[rows, cols] = size(FluorescenceLevelsRawBGnorm(:,:,1));
[x, y] = meshgrid(1:cols, 1:rows);

% Initialize arrays to store final intensity values per cell
finalSignal_mem = zeros(length(stats), 1);
finalSignalN = zeros(length(stats), 1);

% Outcome: 
%   AllDataExperiment.mat: AllDataExperiment{1,platenum}{wellnum} contains
%   [xposition,yposition,rawDAPI,DAPInormalisedChannel2,DAPInormalisedChannel3,DAPInormalisedChannel4,Area,position in well] for every
%   cell in each plate and well


% Loop over each nucleus
MatrixData_perPosition = zeros(length(stats),6);
for k = 1:length(stats)
    MatrixData_perPosition(k,1) = stats(k).Centroid(1);%x coordinate of centroid
    MatrixData_perPosition(k,2) = stats(k).Centroid(2);%y coordinate of centroid
    MatrixData_perPosition(k,4+3) = stats(k).Area;
    
    % Get the centroid of the current nucleus
    centroid = stats(k).Centroid;
    for c = 1:analysisParam.ChannelMaxNum{PlateNum}(WellNumber)

        if matches(analysisParam.MemOrNuc{PlateNum}{WellNumber}{c},'DAPI')
            FluorescenceLevelsRaw = FluorescenceLevelsRawBGnorm(:,:,c);
            MatrixData_perPosition(k,c+2) = mean(FluorescenceLevelsRaw(stats(k).PixelIdxList));%-backgroundintChannel2;
        elseif matches(analysisParam.MemOrNuc{PlateNum}{WellNumber}{c},'Nuclear')
            FluorescenceLevelsRaw = FluorescenceLevelsRawBGnorm(:,:,c);
            MatrixData_perPosition(k,c+2) = mean(FluorescenceLevelsRaw(stats(k).PixelIdxList));%-backgroundintChannel2;
        elseif matches(analysisParam.MemOrNuc{PlateNum}{WellNumber}{c},'Membrane')
        FluorescenceLevelsRaw = FluorescenceLevelsRawBGnorm(:,:,c);
            % Initialize arrays to store top 20% values for each direction
        topValues_mem = zeros(1, numDirections);
        
        % Loop over multiple directions (e.g., 6 directions)
        for d = 1:numDirections
            % Calculate the angle for the current direction (equally spaced angles)
            angle = (d-1) * (2 * pi / numDirections);
            directionVec = [cos(angle), sin(angle)];
            
            % Sample along the current direction from the centroid
            distances = 1:maxDistance;
            sampledPoints = bsxfun(@plus, centroid, distances' * directionVec);
            
            % Make sure sampled points are within image bounds
            validIdx = all(sampledPoints >= 1 & sampledPoints <= [cols, rows], 2);
            sampledPoints = sampledPoints(validIdx, :);
            
            % Check if there are valid sampled points before proceeding
            if isempty(sampledPoints)
                continue; % Skip this direction if no valid points
            end
            
            % Get intensity values along this direction for this Channel 
            sampleIndices = sub2ind([rows, cols], round(sampledPoints(:,2)), round(sampledPoints(:,1)));
            intensities_mem = FluorescenceLevelsRaw(sampleIndices); 
            
            % Check if the intensities array is empty
            if isempty(intensities_mem)
                continue; % Skip this direction if no valid intensity values
            end
            
            % Sort intensity values in descending order and then extract the top 20%
            sortedIntensities_mem = sort(intensities_mem, 'descend');
            
            % Calculate the number of top values to select, ensuring it's at least 1
            numTopValues = max(1, round(topPercentage * length(sortedIntensities_mem)));
            
            % Safely get the top 20% values (ensure there are enough elements)
            topValues_mem(d) = mean(sortedIntensities_mem(1:numTopValues));
        end
        
        % Calculate the average of the top 20% values across all directions
        finalSignal_mem = mean(topValues_mem);
        MatrixData_perPosition(k,c+2) = finalSignal_mem;
        end
    
    end

end

% % Store results in a table or display
% resultsTable = table((1:length(stats))', MatrixData_perPosition(:,2), MatrixData_perPosition(:,3), MatrixData_perPosition(:,4),...
%     'VariableNames', {'NucleusID', 'AvgSignal_E', 'AvgSignal_N', 'AvgSignal_nuclei'});
% disp(resultsTable);

alldata{nposition} = MatrixData_perPosition;
end
name = strcat(analysisParam.savingpathforData,'/Plate',num2str(PlateNum),'_Well',num2str(WellNumber,'%03.f'),'_Data.mat');
save(name,'alldata','analysisParam') 
alldatawells{PlateNum}{WellNumber} = alldata;
end
end

name = strcat(analysisParam.savingpathforData,'/AllDataPlates.mat');
save(name,'alldatawells','analysisParam')
% Outcome: 
%       Platei_Wellj_Data.mat: alldata{position} contains
%       [xposition,yposition,rawDAPI,Channel2,Channel3,Channel4,Area] for every
%       cell
%       AllDataPlates: alldatawells{PlateNum}{WellNumber}{position} contains
%       [xposition,yposition,rawDAPI,Channel2,Channel3,Channel4,Area] for every
%       cell in PlateNum, Wellnumber, position
%% Format data into a data array for each plate
for platenum = 1:analysisParam.NumofPlates
disp(['Plate number:',num2str(platenum)])

AllDataCell = cell(1,1);
AllDataMatrix = cell(1,analysisParam.NumofWells);
CellNumberPosition = cell(1,analysisParam.NumofWells);

for wellnum = analysisParam.WellsWithData{platenum}
    
load([analysisParam.savingpathforData,'/Plate',num2str(platenum),'_Well',num2str(wellnum,'%03.f'),'_Data.mat'],'alldata');
AllDataCell{wellnum} = alldata;

    for positionnumber = 1:size(alldata,2)
        AllDataMatrix{wellnum} = [AllDataMatrix{wellnum};[alldata{positionnumber},positionnumber*ones(size(alldata{positionnumber},1),1)]];

        CellNumberPosition{wellnum}(positionnumber) = size(alldata{positionnumber},1);
    end
      
end

NamesConditions = analysisParam.NamesConditions{platenum};

ChannelMax = 4;
minnormValues = min(AllDataMatrix{analysisParam.WellsWithData{platenum}(1)}(:,3:(ChannelMax+2)));
maxnormValues = max(AllDataMatrix{analysisParam.WellsWithData{platenum}(1)}(:,3:(ChannelMax+2)));

for wellnum = analysisParam.WellsWithData{platenum}
    minnormValues = min([minnormValues;AllDataMatrix{wellnum}(:,3:(ChannelMax+2))]);
    maxnormValues = max([maxnormValues;AllDataMatrix{wellnum}(:,3:(ChannelMax+2))]);
end

analysisParam.minnormValues = minnormValues;
analysisParam.maxnormValues = maxnormValues;

save([analysisParam.savingpathforData,'/Plate',num2str(platenum),'_AllDataMatrix.mat'])
end
%% Save all data in each plate together
% Outcome: 
%       AllDataExperiment.mat: AllDataExperiment{1,platenum}{wellnum} contains
%       [xposition,yposition,rawDAPI,DAPInormalisedChannel2,DAPInormalisedChannel3,DAPInormalisedChannel4,Area,position in well] for every
%       cell in each plate and well
%       Also finds max and min intensities of each channel and saves them
%       in limitschannels (ordered as
%       analysisParam.MapChannels.DifferentChannelsPresent)

for platenum = 1:analysisParam.NumofPlates
    load([analysisParam.savingpathforData,'/Plate',num2str(platenum),'_AllDataMatrix.mat'],'AllDataMatrix');
    AllDataExperiment{1,platenum} = AllDataMatrix;
    clear AllDataMatrix
end

save([analysisParam.savingpathforData,'/AllDataExperiment'])

%%
DA_FindLimitsData_Raw
%% Data Analysis: Clean Data
cd(analysisParam.pathnamesave)
minArea = 60;
maxArea = 1000;

minquantile = 0.0005;
maxquantile = 0.9995;

%If no Background image provided:
options = {'minlimarea',minArea,'maxlimarea',maxArea,'minquantile',minquantile,'maxquantile',maxquantile};

DA_FindLimitsData(options)

%% define parameters for violin plots
close all
cd(analysisParam.pathnamesave)
clear ConditionsSelection

% follow the order of markers (DAPI, ...)
% maxLimBar=[0 0 0 0 1.5 0 1.5 0.8 1.5 1 1 1.5];
maxLimVio=[0 1300 800 1500 0 0 0 0 0 0 0 0 0 0 0];
minLimVio=[0 0 100 0 0 0 0 0 0 0 0 0 0 0 0];
colorOrder=[0 1 2 3 0 0 0 0 0 0 0 0 0 0 0];

plotname_all={'DE Cad','DE SNAIL1'};
channelnums_all = {[2 3],4,};
ConditionsSelection_all = {[3 2 3;159 1 1],[1 2 3;1 1 1]};

for i = 1:length(plotname_all)
    plotname = plotname_all{i};
    channelnums = channelnums_all{i};
    titleplottosave=plotname;title = plotname;
    ConditionsSelection = ConditionsSelection_all{i}; 
    
windowdims = [10 10 length(ConditionsSelection)*50+100 length(channelnums)*200+100];

% Data Analysis: Violin plots
% Optional options
raworclean = 1; %0 if RAW data, 1 if using CLEAN data (0 default)
blackBG = 0; %0 if white background, 1 if black background (0 default)
% channelnums = [5]; %choose channels in analysisParam.MapChannels.DifferentChannelsPresent, (analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,1),ConditionsSelection(2,1)} default)
stdmean = 1; %1 if plotting std of the mean, 0 if plotting std of data;
angleticks = 45; %angle for xticks in barplots

[sidx,ordercolors]=sort(channelnums);
uniformlimits = 0;
options = {'raworclean',raworclean,'blackbg',blackBG,'channels',channelnums,'stdmean',stdmean,'angleticks',angleticks,'title',title,'uniformlimits',uniformlimits,'ordercolors',ordercolors,'ylimdef',maxLimVio,'ylimdeflow',minLimVio,'windowdims',windowdims,'colordef',colorOrder};

DA_ViolinPlotsData_FixedData_ConditionsSelection(ConditionsSelection,titleplottosave,options)

close all

end

%% define parameters for scatter plots
tic
close all
clear ConditionsSelection
%DAPI, E-Cad, N-Cad, SNAIL1, MIXL1，
%FOXA2, SOX17, CDX2, TBX6, BRA，
%ISL1，HAND1, SOX2, OCT4, NANOG
stainingGroup_Channel = {[2 3 4],[5 6 7],[5 8 9],[10 11 12],[10 11 9],[13 14 15]};

plotname_all={'DE','PM','LM'};
stainingGroup_all = {1,1,1};
ConditionsSelection_all = {[2 3;1 1],[2 3;17 27],[2 3;27 57]};


for n = 1:length(plotname_all)
for i = 1:length(stainingGroup_all{n})
    ss=stainingGroup_all{n}(i);
    channelnums = stainingGroup_Channel{ss};
    ConditionsSelection = ConditionsSelection_all{n}; 
    
% Data Analysis: Scatter Plots Separate
close all
cd(analysisParam.pathnamesave)

%Optional options
raworclean = 1; %0 if RAW data, 1 if using CLEAN data (0 default)
blackBG = 0; %0 if white background, 1 if black background (0 default)
% channelnums = [2,3,4]; %choose channels in analysisParam.MapChannels.DifferentChannelsPresent, (analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,1),ConditionsSelection(2,1)} default)
stdmean = 1; %1 if plotting std of the mean, 0 if plotting std of data;
angleticks = 45; %angle for xticks in barplots
[sidx,ordercolors]=sort(channelnums);
uniformlimits = 1;

% options = {'raworclean',raworclean,'blackbg',blackBG,'channels',channelnums,'stdmean',stdmean,'angleticks',angleticks,'title',title,'uniformlimits',uniformlimits,'ordercolors',ordercolors};

% To plot 2 genes
if length(ConditionsSelection)>=5
    windowdims=[10 10 1800 50+240*ceil(length(ConditionsSelection)/5)];
else
    windowdims=[10 10 300+length(ConditionsSelection)*300 50+240*ceil(length(ConditionsSelection)/5)];
end
scatter2order = [3 2]; %G2-BRA FOXA2 SOX17 in the order of SOX17-FOXA2
titleplottosave = [plotname_all{n} ' G' num2str(ss) ' ' analysisParam.MapChannels.DifferentChannelsPresent{channelnums(scatter2order(1))} ' ' analysisParam.MapChannels.DifferentChannelsPresent{channelnums(scatter2order(2))}];

options = {'raworclean',raworclean,'blackbg',blackBG,'channels',channelnums,'stdmean',stdmean,'angleticks',angleticks,'uniformlimits',uniformlimits,'ordercolors',ordercolors,'windowdims',windowdims,'scatter2order',scatter2order};
DA_HeatScatter_FixedData_ConditionsSelection(ConditionsSelection,titleplottosave,options)
titleplottosave = strrep(titleplottosave,'.',' ')
scatter2order = [3 1]; %G2-BRA FOXA2 SOX17 in the order of SOX17-FOXA2
titleplottosave = [plotname_all{n} ' G' num2str(ss) ' ' analysisParam.MapChannels.DifferentChannelsPresent{channelnums(scatter2order(1))} ' ' analysisParam.MapChannels.DifferentChannelsPresent{channelnums(scatter2order(2))}];
titleplottosave = strrep(titleplottosave,'.',' ')

close all
end
end
toc