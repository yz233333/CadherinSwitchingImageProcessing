%% 0.0 Define directories
clc
clear
close all
%%
% master folder
masterFolder = '/Volumes/research/aw21/Ye/Data/CadSwitch LIVE 8+/Current-IP/ilastik_Pre_MIP_dead';
% raw image folder
rawFolder = [masterFolder filesep 'C1C2C3_raw'];
filelist  = dir(rawFolder);
% ilastik mask folder
maskFolder = [masterFolder filesep 'C1_raw'];
deadCellMaskFolder = [masterFolder filesep 'C2C3_MAX'];

%% 0.1 Specify info: conditions, channels and time points
% PS 13
experimentNum = 13;
experimentName=['PS' num2str(experimentNum)];

% conditions
conditionNames = {'01 Pluri','02 PS1','03 MPS','04 MPS_noAct','05 MPS_noBmp','06 MPS-LM','07 MPS-LM_noAct','08 MPS-LM_noBmp'}; %,'09 blank'}; % name matches the condition folder name in the rawData Folder
nConditions=length(conditionNames);
nColoniesPerCondition = [3,3,3,3,3,3,3,3];%,1];% # of colonies for each well/condition
nColoniesMax=max(nColoniesPerCondition);
conds   = {'mTeSR PLUS','PS1','MPS','MPS noAct','MPS noBmp','MPS-LM','MPS-LM(-ActivinA)','MPS-LM(-Bmp4)'}; %,'blank'}; % condition names in plots
% define colors
colors_plot3 = {'#fdbf6f','#fb9a99','#1f78b4','#a6cee3','#8dd3c7',...
    '#6a3d9a','#cab2d6','#bc80bd','#b15928'};

% channels
channelNames = {'CAAX:mCFP','E-Cad:mCitrine','N-Cad:mCherry'};%skip the brightfield channel
nChannels = length(channelNames);

% time points
timepoint = [1:1.45:21.3 24.5:1.45:44.8];
nT = length(timepoint); % total time points

%% 0.2 preallocate
mem_intensities          = cell(max(nColoniesPerCondition),nConditions);
mem_counter              = cell(max(nColoniesPerCondition),nConditions);
mem_intensities_final    = cell(max(nColoniesPerCondition),nConditions);
cell_intensities         = cell(max(nColoniesPerCondition),nConditions);
cell_counter             = cell(max(nColoniesPerCondition),nConditions);
cell_intensities_final   = cell(max(nColoniesPerCondition),nConditions);
bg_intensities           = cell(max(nColoniesPerCondition),nConditions);
bg_counter               = cell(max(nColoniesPerCondition),nConditions);
bg_intensities_final     = cell(max(nColoniesPerCondition),nConditions);
nonMem_intensities       = cell(max(nColoniesPerCondition),nConditions);
nonMem_counter           = cell(max(nColoniesPerCondition),nConditions);
nonMem_intensities_final = cell(max(nColoniesPerCondition),nConditions);
nonBg_intensities        = cell(max(nColoniesPerCondition),nConditions);
nonBg_counter            = cell(max(nColoniesPerCondition),nConditions);
nonBg_intensities_final  = cell(max(nColoniesPerCondition),nConditions);
memToBg                = cell(max(nColoniesPerCondition),nConditions);
memToCell              = cell(max(nColoniesPerCondition),nConditions);
memToNonmem            = cell(max(nColoniesPerCondition),nConditions);
nonBgToBg              = cell(max(nColoniesPerCondition),nConditions);

%%
parpool(4);
%% 1.1 Analyze intensities and ratios
parfor ii = 1:nConditions % loop thru each condition
    tic 
    disp(['Condition' int2str(ii) ': ' conds{ii}])

%     for jj = 1:nColoniesPerCondition(ii) % loop thru each colony
    for jj = 1:nColoniesMax % loop thru each colony
        disp(['Colony' int2str(jj)])
        colonyName = [experimentName '_' conditionNames{1,ii} ' ' num2str(jj,'%02.f')];
        
        % read raw image
        rawFile = [rawFolder filesep  colonyName '_pre.tif'];
        filereader = bfGetReader(rawFile);
        nx = filereader.getSizeX; % x
        ny = filereader.getSizeY; % y
        nt = filereader.getSizeT; % # of time points
        nz = filereader.getSizeZ; % # of z slices
        nc = filereader.getSizeC; % # of channels

        % read main mask (membrane, background, cell)
        maskFile = [maskFolder filesep colonyName '_pre_C1_Simple Segmentation.h5'];
        [mask_mem,mask_bg,mask_cell] = readIlastikFile3(maskFile,1); 

        % read dead cell mask (dead, background, cell)
        deadcellmaskFile = [deadCellMaskFolder filesep colonyName '-MAX_C2C3_Simple Segmentation.h5'];
        [mask_dead,mask_bg_d,mask_cell_d] = readIlastikFile3(deadcellmaskFile,1);

        if ii == 1
            mask_dead_Final = bwareaopen(mask_dead,200);% remove small objects
            mask_dead_Final = imdilate(mask_dead_Final,strel('disk',15));
        else
            mask_dead_Final = bwareaopen(mask_dead,120);% remove small objects
            mask_dead_Final = imdilate(mask_dead_Final,strel('disk',4)); 
        end
        mask_dead_Final = imfill(mask_dead_Final,'holes');% fill holes
        
        if ii == 1 && experimentNum ~= 13
            % read dead cell mask_pluri_late (dead, background, cell)
            deadcellmaskFile_pluri_late = [deadCellMaskFolder filesep colonyName '-MAX_C2C3_late_Simple Segmentation.h5'];
            [mask_dead_late,mask_bg_d_late,mask_cell_d_late] = readIlastikFile3(deadcellmaskFile_pluri_late,1);
            mask_dead_late_processed = bwareaopen(mask_dead_late,120);% remove small objects
            mask_dead_late_processed = imdilate(mask_dead_late_processed,strel('disk',10));
            mask_dead_late_processed = imfill(mask_dead_late_processed,'holes');% fill holes
        
            mask_dead_Final(:,:,21:nt) = mask_dead_late_processed(:,:,21:nt);
        end

        % remove deadCell pixels from the original masks
        mask_mem=removeDeadCellMask(mask_mem,mask_dead_Final);
        mask_bg=removeDeadCellMask(mask_bg,mask_dead_Final);
        mask_cell=removeDeadCellMask(mask_cell,mask_dead_Final);
        
        % preallocate
        mem_intensities{jj,ii}    = zeros(nt,nc);
        mem_counter{jj,ii}        = zeros(nt,nc);
        cell_intensities{jj,ii} = zeros(nt,nc);
        cell_counter{jj,ii}     = zeros(nt,nc);
        bg_intensities{jj,ii} = zeros(nt,nc);
        bg_counter{jj,ii}     = zeros(nt,nc);
        nonBg_intensities{jj,ii} = zeros(nt,nc);
        nonBg_counter{jj,ii}     = zeros(nt,nc);
        nonMem_intensities{jj,ii} = zeros(nt,nc);
        nonMem_counter{jj,ii}     = zeros(nt,nc);

        for tt = 1:nt
            for zz = 1:nz-1
                masknow_mem = squeeze(mask_mem(:,:,zz,tt));
                masknow_cell = squeeze(mask_cell(:,:,zz,tt));
                masknow_bg = squeeze(mask_bg(:,:,zz,tt));

                for cc = 1:nc
                    imgnow = double(bfGetPlaneAtZCT(filereader,zz,cc,tt));

                    mem_counter{jj,ii}(tt,cc) = mem_counter{jj,ii}(tt,cc) + sum(sum(masknow_mem));
                    mem_intensities{jj,ii}(tt,cc) = mem_intensities{jj,ii}(tt,cc) + sum(sum(masknow_mem.*imgnow));
                    cell_counter{jj,ii}(tt,cc) = cell_counter{jj,ii}(tt,cc) + sum(sum(masknow_cell));
                    cell_intensities{jj,ii}(tt,cc) = cell_intensities{jj,ii}(tt,cc) + sum(sum(masknow_cell.*imgnow));
                    bg_counter{jj,ii}(tt,cc) = bg_counter{jj,ii}(tt,cc) + sum(sum(masknow_bg));
                    bg_intensities{jj,ii}(tt,cc) = bg_intensities{jj,ii}(tt,cc) + sum(sum(masknow_bg.*imgnow));
                    
                    nonBg_counter{jj,ii}(tt,cc) = nonBg_counter{jj,ii}(tt,cc) + sum(sum(masknow_mem)) + sum(sum(masknow_cell));
                    nonBg_intensities{jj,ii}(tt,cc) = nonBg_intensities{jj,ii}(tt,cc) + sum(sum(masknow_mem.*imgnow)) + sum(sum(masknow_cell.*imgnow));
                    nonMem_counter{jj,ii}(tt,cc) = nonMem_counter{jj,ii}(tt,cc) + sum(sum(masknow_cell)) + sum(sum(masknow_bg));
                    nonMem_intensities{jj,ii}(tt,cc) = nonMem_intensities{jj,ii}(tt,cc) + sum(sum(masknow_cell.*imgnow)) + sum(sum(masknow_bg.*imgnow));
                end
            end
        end
        mem_intensities_final{jj,ii} = mem_intensities{jj,ii}./mem_counter{jj,ii};
        cell_intensities_final{jj,ii} = cell_intensities{jj,ii}./cell_counter{jj,ii};
        bg_intensities_final{jj,ii} = bg_intensities{jj,ii}./bg_counter{jj,ii};
        nonBg_intensities_final{jj,ii} = nonBg_intensities{jj,ii}./nonBg_counter{jj,ii};
        nonMem_intensities_final{jj,ii} = nonMem_intensities{jj,ii}./nonMem_counter{jj,ii};

        memToBg{jj,ii} = mem_intensities_final{jj,ii}./bg_intensities_final{jj,ii};
        memToCell{jj,ii} = mem_intensities_final{jj,ii}./cell_intensities_final{jj,ii};
        memToNonmem{jj,ii} = mem_intensities_final{jj,ii}./nonMem_intensities_final{jj,ii};
        nonBgToBg{jj,ii} = nonBg_intensities_final{jj,ii}./bg_intensities_final{jj,ii};

        disp('Colony finished.')
        toc
    end
    disp('Condition finished.')
end
disp('1.1 finished')

% 1.2 Save Result to outputDataFolder
versionNum='V0107';
% output folder
outputDataFolder = [masterFolder filesep 'output_' versionNum];
matName=[outputDataFolder filesep experimentName '_memRatio_z-1_' versionNum '.mat'];

Result = cell(1,9);
Result{1,1} = mem_intensities_final;
Result{1,2} = cell_intensities_final;
Result{1,3} = bg_intensities_final;
Result{1,4} = nonBg_intensities_final;
Result{1,5} = nonMem_intensities_final;
Result{1,6} = memToBg;
Result{1,7} = memToCell;
Result{1,8} = memToNonmem;
Result{1,9} = nonBgToBg;

save(matName, 'Result')
disp(['1.2 Data saved to: ' matName])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.1 Load Result (Run Step 0.0 and 0.1 first)
versionNum='V0107';z='';
% output folder
outputDataFolder = [masterFolder filesep 'output_' versionNum];
matName=[outputDataFolder filesep experimentName '_memRatio_' z versionNum '.mat'];
load(matName)
mem_intensities_final = Result{1,1};
cell_intensities_final = Result{1,2};
bg_intensities_final = Result{1,3};
nonBg_intensities_final = Result{1,4};
nonMem_intensities_final = Result{1,5};
memToBg = Result{1,6};
memToCell = Result{1,7};
memToNonmem = Result{1,8};
nonBgToBg = Result{1,9};
disp(['2.1 Data loaded from: ' matName])
%% 2.2 normalized ratio memToCell by dividing membrane intensity
memToCell_norm = cell(max(nColoniesPerCondition),nConditions);
for ii = 1:nConditions % loop thru each condition
    disp(['Condition' int2str(ii) ': ' conds{ii}])

    for jj = 1:nColoniesPerCondition(ii) % loop thru each colony
%     for jj = 1:nColoniesMax % loop thru each colony
        disp(['Colony' int2str(jj)])

        memToCell_norm{jj,ii}(:,1) = memToCell{jj,ii}(:,1)./memToCell{jj,ii}(:,1);
        memToCell_norm{jj,ii}(:,2) = memToCell{jj,ii}(:,2)./memToCell{jj,ii}(:,1);
        memToCell_norm{jj,ii}(:,3) = memToCell{jj,ii}(:,3)./memToCell{jj,ii}(:,1);
    end
end
disp('2.2 Normalized ratio calculated.')

%% decide if to make optional plot 1-3 
plot1=0;
plot2_1=0;
plot2_2=0;
plot3=1;
%% 2.3 Optional Plot 1: check intensity for each colony (1 image saved per condition)
close all
if plot1==1
% define plot colors
colors ={'#b2df8a','#fdbf6f','#cab2d6','#ffff99'};

% determine ylim_max
ylim_max_mem_intensities_final_ch1=0;
ylim_max_mem_intensities_final_ch2=0;
ylim_max_mem_intensities_final_ch3=0;
for ii=1:nConditions
    for jj=1:nColoniesPerCondition(ii)
        if ylim_max_mem_intensities_final_ch1 < max(mem_intensities_final{jj,ii}(:,1))*1.3
            ylim_max_mem_intensities_final_ch1=max(mem_intensities_final{jj,ii}(:,1))*1.3;
        end
        if ylim_max_mem_intensities_final_ch2 < max(mem_intensities_final{jj,ii}(:,2))*1.2
            ylim_max_mem_intensities_final_ch2=max(mem_intensities_final{jj,ii}(:,2))*1.2;
        end
        if ylim_max_mem_intensities_final_ch3 < max(mem_intensities_final{jj,ii}(:,3))*1.2
            ylim_max_mem_intensities_final_ch3=max(mem_intensities_final{jj,ii}(:,3))*1.2;
        end
    end
end

% membrane intensity plot
for ii = 1 : nConditions
    f=figure; set(gcf,'Position',[0 0 1500 500]);
    hold on
    for cc = 1:nChannels
        subplot(1,3,cc); hold on;
        for jj=1:nColoniesPerCondition(ii)
            plot(timepoint,mem_intensities_final{jj,ii}(:,cc),'.-','LineWidth',3,'MarkerSize',20,'Color',colors{jj});
            
        end
        set(gca,'FontSize',24);title(channelNames{cc});xlabel('Time (hr)');ylabel('Membrane Intensity');xlim([0 48])

        if cc == 1
            legend('Colony 1','Colony 2','Colony 3');legend('Location','best','FontSize',14)
            ylim([0 ylim_max_mem_intensities_final_ch1])
        elseif cc == 2
            ylim([0 ylim_max_mem_intensities_final_ch2])
        elseif cc == 3
            ylim([30 ylim_max_mem_intensities_final_ch3])
        end
    end
    hold off
%SAVE
fDir=[outputDataFolder filesep experimentName '_Plot1_' conditionNames{ii} '_colonyMembraneIntensity_' versionNum '.png'];
saveas(f,fDir,'png')
close
end
end

%% 2.4 Calculate mean, standard deviation, standard error

%[meanRatio,stdRatio,steRatio] = calcMeanStdSte(input,nChannels,nConditions,nColoniesPerCondition,nT)
[mean_memToNonmem,std_memToNonmem,ste_memToNonmem] = calcMeanStdSte(memToNonmem,nChannels,nConditions,nColoniesPerCondition,nT);
[mean_memToCell,std_memToCell,ste_memToCell] = calcMeanStdSte(memToCell,nChannels,nConditions,nColoniesPerCondition,nT);
[mean_memToBg,std_memToBg,ste_memToBg] = calcMeanStdSte(memToBg,nChannels,nConditions,nColoniesPerCondition,nT);
[mean_nonBgToBg,std_nonBgToBg,ste_nonBgToBg] = calcMeanStdSte(nonBgToBg,nChannels,nConditions,nColoniesPerCondition,nT);
[mean_mem,std_mem,ste_mem] = calcMeanStdSte(mem_intensities_final,nChannels,nConditions,nColoniesPerCondition,nT);
[mean_cell,std_cell,ste_cell] = calcMeanStdSte(cell_intensities_final,nChannels,nConditions,nColoniesPerCondition,nT);
[mean_bg,std_bg,ste_bg] = calcMeanStdSte(bg_intensities_final,nChannels,nConditions,nColoniesPerCondition,nT);

[mean_memToCell_norm,std_memToCell_norm,ste_memToCell_norm] = calcMeanStdSte(memToCell_norm,nChannels,nConditions,nColoniesPerCondition,nT);
disp('2.4 Mean, standard deviation, and standard error calculated.')

%% Optional Plot 2.1: Intensity avg plot for each condition (1 image saved per condition)
close all
if plot2_1==1
% define plot color
% membrane = yellow#0072BD, non-mem cell = red#D95319, background = blue#EDB120
colors={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30'}; % MATLAB default color

% determine ylim_max
ylim_max_mean_mem_ch1=0;
ylim_max_mean_mem_ch2=0;
ylim_max_mean_mem_ch3=0;
for ii = 1: nConditions
    if ylim_max_mean_mem_ch1 < max(mean_mem{1,ii})*1.15
        ylim_max_mean_mem_ch1=max(mean_mem{1,ii})*1.15;
    end
    if ylim_max_mean_mem_ch2 < max(mean_mem{2,ii})*1.1
        ylim_max_mean_mem_ch2=max(mean_mem{2,ii})*1.1;
    end
    if ylim_max_mean_mem_ch3 < max(mean_mem{3,ii})*1.1
        ylim_max_mean_mem_ch3=max(mean_mem{3,ii})*1.1;
    end
end

% plot
for ii = 1: nConditions
    f=figure; set(gcf,'Position',[0 0 1500 500]);
    for cc = 1:nChannels
        subplot(1,3,cc); hold on;
        plot(timepoint,mean_mem{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors{3});
        errorbar(timepoint,mean_mem{cc,ii},ste_mem{cc,ii}/2,'.','Color',colors{3},'HandleVisibility','off')  

        plot(timepoint,mean_cell{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors{2});
        errorbar(timepoint,mean_cell{cc,ii},ste_cell{cc,ii}/2,'.','Color',colors{2},'HandleVisibility','off') 

        plot(timepoint,mean_bg{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors{1});
        errorbar(timepoint,mean_bg{cc,ii},ste_bg{cc,ii}/2,'.','Color',colors{1},'HandleVisibility','off') 
        
        set(gca,'FontSize',24);title(channelNames{cc});xlabel('Time (hours)');ylabel('Intensity');xlim([0 48])

        if cc == 1
            legend('Membrane', 'Non-membrane Cell','Background');
            legend('Location','best','FontSize',14)
            ylim([0 ylim_max_mean_mem_ch1])
        elseif cc == 2
            ylim([0 ylim_max_mean_mem_ch2])
        elseif cc == 3
            ylim([30 ylim_max_mean_mem_ch3])
        end
    end

    %SAVE
    fDir=[outputDataFolder filesep experimentName '_Plot2-1_' conditionNames{ii} '_memCellBgIntensity_avg_' versionNum '.png'];
    saveas(f,fDir,'png')
    close
end
end
%% Optional Plot 2.2: Ratio avg plot for each condition (1 image saved per condition)
close all
if plot2_2==1
% define plot color
colors ={'#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99'};

% determine ylim_max
ylim_max_mean_memToBg_ch1=0;
ylim_max_mean_memToBg_ch2=0;
ylim_max_mean_memToBg_ch3=0;
for ii = 1: nConditions
    if ylim_max_mean_memToBg_ch1 < max(mean_memToBg{1,ii})*1.1
        ylim_max_mean_memToBg_ch1=max(mean_memToBg{1,ii})*1.1;
    end
    if ylim_max_mean_memToBg_ch2 < max(mean_memToBg{2,ii})*1.1
        ylim_max_mean_memToBg_ch2=max(mean_memToBg{2,ii})*1.1;
    end
    if ylim_max_mean_memToBg_ch3 < max(mean_memToBg{3,ii})*1.1
        ylim_max_mean_memToBg_ch3=max(mean_memToBg{3,ii})*1.1;
    end
end

% plot
for ii = 1: nConditions
    f=figure; set(gcf,'Position',[0 0 1500 500]);
    for cc = 1:nChannels
        subplot(1,3,cc); hold on;
        plot(timepoint,mean_memToNonmem{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors{1});
        errorbar(timepoint,mean_memToNonmem{cc,ii},ste_memToNonmem{cc,ii}/2,'.','Color',colors{1},'HandleVisibility','off')  

        plot(timepoint,mean_memToBg{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors{2});
        errorbar(timepoint,mean_memToBg{cc,ii},ste_memToBg{cc,ii}/2,'.','Color',colors{2},'HandleVisibility','off')  

        plot(timepoint,mean_nonBgToBg{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors{3});
        errorbar(timepoint,mean_nonBgToBg{cc,ii},ste_nonBgToBg{cc,ii}/2,'.','Color',colors{3},'HandleVisibility','off')  

        plot(timepoint,mean_memToCell{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors{4});
        errorbar(timepoint,mean_memToCell{cc,ii},ste_memToCell{cc,ii}/2,'.','Color',colors{4},'HandleVisibility','off')  
        
        title(channelNames{cc});xlabel('Time (hours)');ylabel('Ratio');xlim([0 48])
       
        if cc == 1
            ylim([0 ylim_max_mean_memToBg_ch1])
        elseif cc == 2
            ylim([0 ylim_max_mean_memToBg_ch2])
        elseif cc == 3
            legend('Mem:Non-mem-px', 'Mem:Bg','Non-bg-px:Bg','Mem:Non-mem-cell');
            legend('Location','best','FontSize',14)
            ylim([0.7 ylim_max_mean_memToBg_ch3])
        end
        set(gca,'FontSize',24);
    end

    %SAVE
    fDir=[outputDataFolder filesep experimentName '_Plot2-2_' conditionNames{ii} '_ratio_avg_' versionNum '.png'];
    saveas(f,fDir,'png')
    close
end
end

%% Optional Plot 3: avg memToCell plot for multiple conditions (1 image saved for all conditions)
close all
% PS 13
plotCondition = {[1 6 7],[1 6 8],[1 6 7 8]};
plotName      = {'MPS-LM_Act','MPS_Bmp','MPS_Act_Bmp'};

for pp=1:length(plotCondition)

% memToCell Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine ylim_max
ylim_max_mean_memToCell_ch1=0;
ylim_max_mean_memToCell_ch2=0;
ylim_max_mean_memToCell_ch3=0;
for ii = 1: nConditions
    if ylim_max_mean_memToCell_ch1 < max(mean_memToCell{1,ii})*1.02
        ylim_max_mean_memToCell_ch1=max(mean_memToCell{1,ii})*1.02;
    end
    if ylim_max_mean_memToCell_ch2 < max(mean_memToCell{2,ii})*1.02
        ylim_max_mean_memToCell_ch2=max(mean_memToCell{2,ii})*1.02;
    end
    if ylim_max_mean_memToCell_ch3 < max(mean_memToCell{3,ii})*1.02
        ylim_max_mean_memToCell_ch3=max(mean_memToCell{3,ii})*1.02;
    end
end

% determine ylim_min
ylim_min_mean_memToCell_ch1=10;
ylim_min_mean_memToCell_ch2=10;
ylim_min_mean_memToCell_ch3=10;
for ii=1:nConditions
    if ylim_min_mean_memToCell_ch1 > min(mean_memToCell{1,ii})*0.7
        ylim_min_mean_memToCell_ch1=min(mean_memToCell{1,ii})*0.7;
    end
    if ylim_min_mean_memToCell_ch2 > min(mean_memToCell{2,ii})*0.98
        ylim_min_mean_memToCell_ch2=min(mean_memToCell{2,ii})*0.98;
    end
    if ylim_min_mean_memToCell_ch3 > min(mean_memToCell{3,ii})*0.98
        ylim_min_mean_memToCell_ch3=min(mean_memToCell{3,ii})*0.98;
    end
end

f=figure; set(gcf,'Position',[0 0 1500 1000]);
% memToCell Plot
for cc = 1:nChannels
    subplot(2,3,cc); hold on;
    for ii = plotCondition{pp}
        plot(timepoint,mean_memToCell{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors_plot3{ii});
        errorbar(timepoint,mean_memToCell{cc,ii},ste_memToCell{cc,ii}/2,'.','Color',colors_plot3{ii},'HandleVisibility','off')  
    end

    set(gca,'FontSize',24);title(channelNames{cc});xlabel('Time (hours)');ylabel('Membrane: Non-Membrane Cell'); xlim([0 48])
    if cc == 1
        legend(conds(plotCondition{pp}),'Location','best','FontSize',14)
        ylim([ylim_min_mean_memToCell_ch1 ylim_max_mean_memToCell_ch1])
    elseif cc == 2
        ylim([ylim_min_mean_memToCell_ch2 ylim_max_mean_memToCell_ch2])
    elseif cc == 3
        ylim([ylim_min_mean_memToCell_ch3 ylim_max_mean_memToCell_ch3])
    end
end

% memToCell Normalized Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine ylim_max
ylim_max_mean_memToCell_norm_ch2=0;
ylim_max_mean_memToCell_norm_ch3=0;
for ii=1:nConditions
    if ylim_max_mean_memToCell_norm_ch2 < max(mean_memToCell_norm{2,ii})*1.02
        ylim_max_mean_memToCell_norm_ch2=max(mean_memToCell_norm{2,ii})*1.02;
    end
    if ylim_max_mean_memToCell_norm_ch3 < max(mean_memToCell_norm{3,ii})*1.02
        ylim_max_mean_memToCell_norm_ch3=max(mean_memToCell_norm{3,ii})*1.02;
    end
end

% determine ylim_min
ylim_min_mean_memToCell_norm_ch2=10;
ylim_min_mean_memToCell_norm_ch3=10;
for ii=1:nConditions
    if ylim_min_mean_memToCell_norm_ch2 > min(mean_memToCell_norm{2,ii})*0.98
        ylim_min_mean_memToCell_norm_ch2=min(mean_memToCell_norm{2,ii})*0.98;
    end
    if ylim_min_mean_memToCell_norm_ch3 > min(mean_memToCell_norm{3,ii})*0.98
        ylim_min_mean_memToCell_norm_ch3=min(mean_memToCell_norm{3,ii})*0.98;
    end
end

% memToCell Normalized Plot
for cc = 1:nChannels
    subplot(2,3,cc+3); hold on;
    for ii = plotCondition{pp}
        plot(timepoint,mean_memToCell_norm{cc,ii},'.-','LineWidth',3,'MarkerSize',20,'Color',colors_plot3{ii});
        errorbar(timepoint,mean_memToCell_norm{cc,ii},ste_memToCell_norm{cc,ii}/2,'.','Color',colors_plot3{ii},'HandleVisibility','off')  
    end

    if cc == 1
        legend(conds(plotCondition{pp}),'Location','best','FontSize',14)
        ylim([0 1.5])
    elseif cc == 2
        ylim([ylim_min_mean_memToCell_norm_ch2 ylim_max_mean_memToCell_norm_ch2])
    elseif cc == 3
        ylim([ylim_min_mean_memToCell_norm_ch3 ylim_max_mean_memToCell_norm_ch3])
    end
    set(gca,'FontSize',24);title(channelNames{cc});xlabel('Time (hours)');ylabel('Membrane: Non-Mem Cell (Normalized)'); xlim([0 48])
end
versionNum='V0107';

%SAVE
fDir=[masterFolder filesep experimentName '_Plot3-2_memToCell_' plotName{pp} '_' z versionNum '.png'];
saveas(f,fDir,'png')
%close
end

%% Main plot: E-Cad and N-Cad of one condition in on a plot
plotCondition = {[1 6 7],[1 6 8],[1 6 7 8]};
plotName      = {'MPS-LM_Act','MPS_Bmp','MPS_Act_Bmp'};

for ii = [1 6 7 8]
    f=figure; set(gcf,'Position',[2200 1800 400 500]);
    left_color = [0 0.62 0.349];%'#009e59');%green
    right_color = [0.835 0.161 0];%'#d52900');%red
    set(f,'defaultAxesColorOrder',[left_color; right_color]);
    hold on
    
    yyaxis left
    plot(timepoint,mean_memToCell_norm{2,ii},'.-','LineWidth',3,'MarkerSize',20,'Color','#009e59');%green
    
    ylim([0.6 1])%ECad
    ylabel('E-Cad:mCitrine (Normalized Mem:Cell)'); 
    
    yyaxis  right
    plot(timepoint,mean_memToCell{3,ii},'.-','LineWidth',3,'MarkerSize',20,'Color','#d52900');%red
    errorbar(timepoint,mean_memToCell{3,ii},ste_memToCell{3,ii}/2,'.','Color','#d52900','HandleVisibility','off')  
    ylim([1 1.7])%NCad
    ylabel('N-Cad:mCherry (Mem:Cell)');
    
    yyaxis left
    errorbar(timepoint,mean_memToCell_norm{2,ii},ste_memToCell_norm{2,ii}/2,'.','Color','#009e59','HandleVisibility','off')  
    title(conds{ii})
    set(gca,'FontSize',20);
    xlabel('Time (hours)');
    xlim([0 45])

%     fDir=[masterFolder filesep experimentName '_' conds{ii} '.png'];
%     saveas(f,fDir,'png')

end

% close all