global analysisParam;

fprintf(1, '%s called  to define params\n',mfilename);

%% Parameters to modify

% Path to your experiment
analysisParam.pathnamesave = '/Users/yezhu/Desktop/Lab file/Data/CadSwitch FIXED/240729 EcadKO 2/nuclearAndMembrane';
analysisParam.pathnamedata = [analysisParam.pathnamesave,'/MAXProj']; %No need to change usually
% Number of plates
analysisParam.NumofPlates = 2;
% Plate Type
analysisParam.NumofWells = 94;
% Number of images per well
analysisParam.ImagesperWell = 8;

% define NamesConditions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = {'E-CadKO 2'};% Define experiment name
treatment = {'DE(a)','none','APS','MPS','PPS',...
             'MPS Ldn','MPS Bmp4 0','MPS Bmp4 80',...
             'APS Iwp2','APS Chir 0','APS Chir 10',...
             'DE(a)-DE(b)','none-PM','none-LM','APS-PM','APS-LM','MPS-PM','MPS-LM','PPS-PM','PPS-LM',...
             'MPS-LM Ldn Ldn','MPS-LM Bmp4 0 0','MPS-LM Ldn Bmp4 30','MPS-LM Bmp4 0 30','MPS-LM Bmp4 80 30','MPS-LM Bmp4 40 Ldn','MPS-LM Bmp4 40 0','MPS-LM Bmp4 40 60',...
             'APS-PM Iwp2 Iwp2','APS-PM Chir 0 0','APS-PM Iwp2 Chir','APS-PM Chir 0 3','APS-PM Chir 10 3','APS-PM Chir Iwp2','APS-PM Chir 4 0','APS-PM Chir 4 10',...
             'mTeSR+','APS Chir 4','MPS Bmp4 40','APS-PM Chir 4 3','MPS-LM Bmp4 40 30'}; % Define treatment names
treatment_plotlabel = treatment; % Define treatment names for plotting

% define immunostaining markers: protein Name, order, membrane or nuclear marker...
staining = {'ECAD NCAD SNAIL1','MIXL1 FOXA2 SOX17','MIXL1 CDX2 TBX6','BRA ISL1 HAND1','BRA ISL1 TBX6','SOX2 OCT4 NANOG'};%in file names
channels = {{'DAPI','E-Cad','N-Cad','SNAIL1'},{'DAPI','MIXL1','FOXA2','SOX17'},...
            {'DAPI','MIXL1','CDX2','TBX6'},{'DAPI','BRA','ISL1','HAND1'},...
            {'DAPI','BRA','ISL1','TBX6'},{'DAPI','SOX2','OCT4','NANOG'}};
OrderChannel = {{1,2,3,4},{1,2,3,4},{1,2,3,4},{1,2,3,4},{1,2,3,4},{1,2,3,4}};
MemOrNuc = {{'DAPI','Membrane','Membrane','Nuclear'},{'DAPI','Nuclear','Nuclear','Nuclear'},{'DAPI','Nuclear','Nuclear','Nuclear'},{'DAPI','Nuclear','Nuclear','Nuclear'},...
            {'DAPI','Nuclear','Nuclear','Nuclear'},{'DAPI','Nuclear','Nuclear','Nuclear'}};

cellLine = {'WT','E-Cad--','E-Cad+-'}; % Cell line
timePoint = {'0hr','24hr'}; % Time point

% define treatment number, marker number and cell line number for each well tested
treatmentNum = {[37 37 37 37 37 37 37 37 37 37 37 37 37],...
                [1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 6 6 6 6 7 7 7 7 39 39 39 39 8 8 8 8 9 9 9 9 10 10 10 10 38 38 38 38 11 11 11 11 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37]};
stainingNum  = {[1 2 3 4 5 1 2 3 4 5 1 2 5],...
                [1 2 1 2 1 2 1 3 4 5 1 3 4 5 1 5 1 3 4 5 1 3 4 5 1 5 1 3 4 5 1 3 4 5 1 5 1 3 4 5 1 3 4 5 1 5 1 4 1 4 1 4 1 4 1 4 1 4 1 4 1 4 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 3 1 2 3 4 5 6 1 2 3 4 5 6 1 2 5 6]};
cellNum      = {[1 1 1 1 1 2 2 2 2 2 3 3 3],...
                [1 1 2 2 3 3 1 1 1 1 2 2 2 2 3 3 1 1 1 1 2 2 2 2 3 3 1 1 1 1 2 2 2 2 3 3 1 1 1 1 2 2 2 2 3 3 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3]};

% Names of conditions in each well and plate
nWell=analysisParam.NumofWells;
% plate_name = cell(1,analysisParam.NumofWells);

NameConditions=cell(1,analysisParam.NumofPlates);
Channelsnames=cell(1,analysisParam.NumofPlates);

for plateNum = 1:analysisParam.NumofPlates
    TIME = timePoint{plateNum};
    disp(['Processing: Plate ' num2str(plateNum) ' ' TIME]);
    
    for wellNum = 1:length(treatmentNum{plateNum})
        disp(['Plate ' num2str(plateNum) ' Well ' num2str(wellNum)]);
        TREATMENT = treatment{treatmentNum{plateNum}(wellNum)};
        STAINING = staining{stainingNum{plateNum}(wellNum)};
        CELL = cellLine{cellNum{plateNum}(wellNum)};
        
        
        name = [TIME ' ' TREATMENT ' ' CELL ' ' STAINING ' ' ];
        plate_name{wellNum} = name;
        
        name_plot = [TREATMENT ' (' TIME ')'];
        %name_plot = [TIME ' ' CELL ' ' TREATMENT];
        plate_name_plot{wellNum} = name_plot;

        PlateChannels{wellNum} = channels{stainingNum{plateNum}(wellNum)};
        PlateOrderChannels{wellNum} = OrderChannel{stainingNum{plateNum}(wellNum)};
        PlateMemOrNuc{wellNum} = MemOrNuc{stainingNum{plateNum}(wellNum)};
    end

    NameConditions{plateNum} = plate_name;
    Channelsnames{plateNum} = PlateChannels;
    OrderChannels{plateNum} = PlateOrderChannels;
    MemOrNuc_{plateNum} = PlateMemOrNuc;

    NameConditions_plot{plateNum} = plate_name_plot;

    clear plate_name
    clear PlateChannels
    clear PlateOrderChannels
    clear PlateMemOrNuc

    clear plate_name_plot

end

analysisParam.NamesConditions = NameConditions;
analysisParam.NamesConditions_plot = NameConditions_plot;
analysisParam.Channelsnames = Channelsnames;
analysisParam.OrderChannels = OrderChannels;
analysisParam.MemOrNuc = MemOrNuc_;

analysisParam.WellsWithData = {1:13, 1:94};
analysisParam.ChannelMaxNum = {};
for ii=1:analysisParam.NumofPlates
    analysisParam.ChannelMaxNum{ii} = cellfun(@length,analysisParam.Channelsnames{ii});
end

%Background images:
%------------------
analysisParam.bgsubstractionopt = 2; %1: use background images given, 2: use min mean background value over images (needs segmentation!), 3: use imopen to substract background, 4: don't substract background
analysisParam.path2BGImages = analysisParam.pathnamedata;
%One image per well and per plate (in this example, the bg image is the
%same for all wells, otherwise you'll need to specify each bg image
bgimagename = 'Background_Exp37_MAXProj.tif';
analysisParam.BGImages = {repmat({bgimagename},1,nWell),repmat({bgimagename},1,nWell)};

%Image processing parameters:
%----------------------------
analysisParam.imopendiskradious = 2;
analysisParam.imerodediskradious = 2;
analysisParam.imdilatediskradious = 2;
analysisParam.pwatershed = 1;


%% Other parameters usually not needed to be modified
analysisParam.savingpathforData = [analysisParam.pathnamesave,'/Matlab_Analysis_Segmentation'];
analysisParam.savingpathforImages = [analysisParam.pathnamesave,'/Visualize_Images_BGSubstracted'];


analysisParam.figDir = 'figures';
mkdir(analysisParam.figDir)
mkdir(analysisParam.savingpathforData)
mkdir(analysisParam.savingpathforImages)

%% Create map channels

IP_CreateMapChannels
