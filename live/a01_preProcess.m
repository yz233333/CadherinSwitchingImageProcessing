clc
clear
close all

% PS 13
experimentName='PS13';
% master folder
masterFolder = '/Users/yezhu/Desktop/Lab file_Ye/ECad-NCad reporter line/221118 PS_bi_13/MATLAB_live';
% raw data folder
rawFolder = [masterFolder filesep 'rawData'];
filelist  = dir(rawFolder);

% output folder
outputDataFolder = '/Users/yezhu/Desktop/Lab file_Ye/ECad-NCad reporter line/ilastik_Pre_MIP_dead/C1C2C3_raw';

%% 0. Specify info: conditions, channels and time points
% conditions
conditionNames = {'01 Pluri','02 PS1','03 MPS','04 MPS_noAct','05 MPS_noBmp','06 MPS-LM','07 MPS-LM_noAct','08 MPS-LM_noBmp','09 blank'}; % name matches the condition folder name in the rawData Folder
nConditions=length(conditionNames);
nColoniesPerCondition = [3,3,3,3,3,3,3,3,1];% # of colonies for each well/condition
nColoniesMax=max(nColoniesPerCondition);
conds   = {'mTeSR PLUS','PS1','MPS','MPS noAct','MPS noBmp','MPS-LM','MPS-LM noAct','MPS-LM noBmp','blank'}; % condition names in plots

% channels
channelNames = {'CAAX:mCFP','E-Cad:mCitrine','N-Cad:mCherry'};%skip the brightfield channel
nChannels = length(channelNames);

% time points
timepoint = [1:1.45:21.3 24.5:1.45:44.8];
nT = length(timepoint); % total time points

%% 1. pre-process - option 1: save multi-t,c maxProj image
backdiskrad = 50;
for ii = 1:nConditions % loop thru each condition
    tic 
    disp(['Condition' int2str(ii) ': ' conds{ii}])

    for jj = 1:nColoniesPerCondition(ii) % loop thru each colony
        disp(['Colony' int2str(jj)])
        colonyName = [conditionNames{1,ii} ' ' num2str(jj,'%02.f')];

        % read raw image
        rawFile = [rawFolder filesep experimentName '_' colonyName '.tif'];
        filereader = bfGetReader(rawFile);
        nt = filereader.getSizeT; % # of time points
        nz = filereader.getSizeZ; % # of z slices
        nc = filereader.getSizeC; % # of channels
        
        tic
        MultiDimImg = zeros(1024-3,1024-2,nc,1,nt,'uint16');

        for tt = 1:nt
            maxZ = MakeMaxZImage(filereader, 1:nc, tt);
            % function maxIntensityImage = MakeMaxZImage(fileReader, channels, timepoint,  z_Range, nSeries)

            for cc = 1:nc
                % pixel shift - raw image (-3x-2 in matlab)
                if cc == 3 % channel of 561 - phase 2 
                    imgnow = maxZ(1:end-3,1:end-2,cc);
                else
                    imgnow = maxZ(4:end,3:end,cc);
                end
                % bg removal
                imgnow = presubBackground_self(imgnow,backdiskrad);

                MultiDimImg(:,:,cc,1,tt)=imgnow;

            end
        end

        
        fiji_descr = ['ImageJ=1.52p' newline ...
                    'images=' num2str(size(MultiDimImg,3)*...
                                      size(MultiDimImg,4)*...
                                      size(MultiDimImg,5)) newline... 
                    'channels=' num2str(size(MultiDimImg,3)) newline...
                    'slices=' num2str(size(MultiDimImg,4)) newline...
                    'frames=' num2str(size(MultiDimImg,5)) newline... 
                    'hyperstack=true' newline...
                    'mode=grayscale' newline...  
                    'loop=false' newline...  
                    'min=0.0' newline...      
                    'max=4095.0'];  % change this to 256 if you use an 8bit image
        tiffName= [outputDataFolder filesep experimentName '_' colonyName '-MAX.tif'];
        t = Tiff(tiffName,'w');
        tagstruct.ImageLength = size(MultiDimImg,1);
        tagstruct.ImageWidth = size(MultiDimImg,2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 16;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.Compression = Tiff.Compression.LZW;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        tagstruct.ImageDescription = fiji_descr;
        for frame = 1:size(MultiDimImg,5)
            for slice = 1:size(MultiDimImg,4)
                for channel = 1:size(MultiDimImg,3)
                    t.setTag(tagstruct)
                    t.write(MultiDimImg(:,:,channel,slice,frame));
                    t.writeDirectory(); % saves a new page in the tiff file
                end
            end
        end
        t.close() 
        toc
    end
end

%% 2. pre-process - option 2: save multi-t,z,c image
backdiskrad = 50;
for ii = 1:nConditions % loop thru each condition
    tic 
    disp(['Condition' int2str(ii) ': ' conds{ii}])

    for jj = 1:nColoniesPerCondition(ii) % loop thru each colony
        disp(['Colony' int2str(jj)])
        colonyName = [conditionNames{1,ii} ' ' num2str(jj,'%02.f')];

        % read raw image
        rawFile = [rawFolder filesep experimentName '_' colonyName '.tif'];
        filereader = bfGetReader(rawFile);
        nt = filereader.getSizeT; % # of time points
        nz = filereader.getSizeZ; % # of z slices
        nc = filereader.getSizeC; % # of channels
        
        tic
        MultiDimImg = zeros(1024-3,1024-2,nc,nz,nt,'uint16');

        for tt = 1:nt
            for zz = 1:nz
                for cc = 1:nc
                    imgnow = bfGetPlaneAtZCT(filereader,zz,cc,tt);
                    % pixel shift - raw image (-3x-2 in matlab)
                    if cc == 3 % channel of 561 - phase 2 
                        imgnow_px = imgnow(1:end-3,1:end-2);
                    else
                        imgnow_px = imgnow(4:end,3:end);
                    end
                    % bg removal
                    imgnow_bgr = presubBackground_self(imgnow_px,backdiskrad);

                    MultiDimImg(:,:,cc,zz,tt)=imgnow_bgr;
   
                end
            end
        end
        
        fiji_descr = ['ImageJ=1.52p' newline ...
                    'images=' num2str(size(MultiDimImg,3)*...
                                      size(MultiDimImg,4)*...
                                      size(MultiDimImg,5)) newline... 
                    'channels=' num2str(size(MultiDimImg,3)) newline...
                    'slices=' num2str(size(MultiDimImg,4)) newline...
                    'frames=' num2str(size(MultiDimImg,5)) newline... 
                    'hyperstack=true' newline...
                    'mode=grayscale' newline...  
                    'loop=false' newline...  
                    'min=0.0' newline...      
                    'max=4095.0'];  % change this to 256 if you use an 8bit image
        tiffName= [outputDataFolder filesep experimentName '_' colonyName '_pre.tif'];
        t = Tiff(tiffName,'w');
        tagstruct.ImageLength = size(MultiDimImg,1);
        tagstruct.ImageWidth = size(MultiDimImg,2);
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 16;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.Compression = Tiff.Compression.LZW;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        tagstruct.ImageDescription = fiji_descr;
        for frame = 1:size(MultiDimImg,5)
            for slice = 1:size(MultiDimImg,4)
                for channel = 1:size(MultiDimImg,3)
                    t.setTag(tagstruct)
                    t.write(MultiDimImg(:,:,channel,slice,frame));
                    t.writeDirectory(); % saves a new page in the tiff file
                end
            end
        end
        t.close() 
        toc
    end
end