function DA_HeatScatter_FixedData_ConditionsSelection(ConditionsSelection,titleplottosave,varargin)


global analysisParam;
global threshold;

if ( length(varargin) == 1 )
    varargin = varargin{:};
end

raworclean = 0;
blackBG = 0;
stdmean = 1;
angleticks = 0;
uniformlimits = 1;
distancepoints = 50;%no DAPI norm
% distancepoints = 0.05;% DAPI norm


raworcleantitle = {'RAW','CLEAN'};

Title=[titleplottosave,' ', raworcleantitle{raworclean+1}];


channelnums = analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,1),ConditionsSelection(2,1)};


while ~isempty(varargin)
    switch lower(varargin{1})
          case {'raworclean'}
              raworclean = varargin{2};   
          case {'blackbg'}
              blackBG = varargin{2};  
          case {'channels'}
              channelnums = varargin{2}; 
              ordercolors = 1:length(channelnums);
          case {'stdmean'}
              stdmean = varargin{2}; 
          case {'angleticks'}
              angleticks = varargin{2};  
          case {'title'}
              Title = [varargin{2}];  
          case {'uniformlimits'}
              uniformlimits = [varargin{2}]; 
          case {'ordercolors'}
              ordercolors = [varargin{2}];    
        case {'windowdims'}%Ye
            windowdims = [varargin{2}];
        case{'scatter3order'}%Ye
            scatter3order = [varargin{2}]; scatter2order = [];
        case{'scatter2order'}%Ye
            scatter2order = [varargin{2}]; scatter3order = [];        
    otherwise
        error(['Unexpected option: ' varargin{1}])
    end
      varargin(1:2) = [];
end

analysisParam.figDir = [analysisParam.pathnamesave filesep 'figures'];
mkdir(analysisParam.figDir)


%%

if raworclean
    load([analysisParam.savingpathforData filesep 'AllDataExperimentClean.mat'])
    AnalysisParamScript_IP;%Ye
else
    load([analysisParam.savingpathforData filesep 'AllDataExperiment.mat'])
    AnalysisParamScript_IP;%Ye
end

if blackBG
    
    colorbg = 'k';
    colorfont = 'w';
    colorbgplotname = 'BLACK';
    
else
    
    colorbg = 'w';
    colorfont = 'k';
    colorbgplotname = 'WHITE';
end

%% Check conditions contain the same channels

nCon = size(ConditionsSelection,2);
nChan = length(channelnums);

FindChannelsinConditions = zeros(nCon,nChan);

for ii = 1:nChan
    
    channelinterest = channelnums(ii);
    
    for jj = 1:nCon
        
        auxvar = find(channelinterest == analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,jj),ConditionsSelection(2,jj)});
        
        if isempty(auxvar)
            error(['DA_BarPlotsData_FixedData_ConditionsSelection: Selected condition ',num2str(jj),' does not contain data for channel ',analysisParam.MapChannels.DifferentChannelsPresent{channelinterest}]);
        
        else
            FindChannelsinConditions(jj,ii) = auxvar;
            
        end
        
    end
    
end

%%

figure;
set(gcf,'Position',windowdims)

x=parula;

nrows=ceil(nCon/5);  
s1=scatter2order(1);s2=scatter2order(2);

for condnum = 1:nCon
   
    if nCon<6
        subplot(1,nCon,condnum)

    else
        subplot(nrows,5,condnum)
    end

    DataPlot =AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(:,2+FindChannelsinConditions(condnum,:));
        
    ndimensions = size(DataPlot,2);

    pointsize=5;

    disp('Plotting')

    if ndimensions == 3 && ~isempty(scatter3order)%Ye
        
        scatterLabel = '3genes';%Ye
        s1=scatter3order(1);s2=scatter3order(2);s3=scatter3order(3);%Ye

        scatter3(DataPlot(:,s1),DataPlot(:,s2),DataPlot(:,s3),pointsize,assignnumberneighbours3(DataPlot(:,s1),DataPlot(:,s2),DataPlot(:,s3),distancepoints),'filled','MarkerEdgeAlpha',1)

        colormap(x(50:end,:))

        xlabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(s1)},'Color',colorfont);
        ylabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(s2)},'Color',colorfont);
        zlabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(s3)},'Color',colorfont);
        
        if uniformlimits
            xlim([limitschannels(1,channelnums(s1)),limitschannels(2,channelnums(s1))])
            ylim([limitschannels(1,channelnums(s2)),limitschannels(2,channelnums(s2))])
            zlim([limitschannels(1,channelnums(s3)),limitschannels(2,channelnums(s3))])
        end
        

        set(gca,'LineWidth', 2);
        set(gca,'FontWeight', 'bold')
        set(gca,'FontName','Arial')
        set(gca,'FontSize',18)
        set(gca,'Color',colorbg)

        set(gca,'XColor',colorfont)
        set(gca,'YColor',colorfont) 
        set(gca,'ZColor',colorfont)
            

    
    elseif ndimensions == 3 && ~isempty(scatter2order)%Ye%%%%%
        
        scatterLabel = '2genes';%Ye
        s1=scatter2order(1);s2=scatter2order(2);%Ye

        hold on;
        scatter(DataPlot(:,s1),DataPlot(:,s2),pointsize,assignnumberneighbours2(DataPlot(:,s1),DataPlot(:,s2),distancepoints),'MarkerEdgeAlpha',1)
        xlabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(s1)},'Color',colorfont);
        ylabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(s2)},'Color',colorfont);

        xline(threshold(analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,condnum),ConditionsSelection(2,condnum)}(s1+1),1), '--'); % s1 threshold
        yline(threshold(analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,condnum),ConditionsSelection(2,condnum)}(s2+1),1), '--'); % s2 threshold
        
        if uniformlimits
            xlim([0,limitschannels(2,channelnums(s1))])
            ylim([0,limitschannels(2,channelnums(s2))])
        end
       
        grid on

        set(gca, 'LineWidth', 2);
        set(gca,'FontWeight', 'bold')
        set(gca,'FontName','Arial')
        set(gca,'FontSize',18)
        set(gca,'Color',colorbg)

        set(gca,'XColor',colorfont)
        set(gca,'YColor',colorfont) 
        hold off;

    elseif ndimensions == 2 && ~isempty(scatter2order)%Ye
        
        scatterLabel = '2genes';%Ye
        s1=scatter2order(1);s2=scatter2order(2);%Ye

        hold on;
        scatter(DataPlot(:,s1),DataPlot(:,s2),pointsize,assignnumberneighbours2(DataPlot(:,s1),DataPlot(:,s2),distancepoints),'MarkerEdgeAlpha',1)
        xlabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(s1)},'Color',colorfont);
        ylabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(s2)},'Color',colorfont);
        
        if uniformlimits
            xlim([0,limitschannels(2,channelnums(s1))])
            ylim([0,limitschannels(2,channelnums(s2))])
        end
        
        xline(threshold(analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,condnum),ConditionsSelection(2,condnum)}(s1+1),1), '--'); % s1 threshold
        yline(threshold(analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,condnum),ConditionsSelection(2,condnum)}(s2+1),1), '--'); % s2 threshold
      

        grid on

        set(gca, 'LineWidth', 2);
        set(gca,'FontWeight', 'bold')
        set(gca,'FontName','Arial')
        set(gca,'FontSize',18)
        set(gca,'Color',colorbg)

        set(gca,'XColor',colorfont)
        set(gca,'YColor',colorfont) 

        hold off;
    else
        scatterLabel = '2genes';%Ye

        hold on;
        scatter(DataPlot(:,2),DataPlot(:,2),pointsize,assignnumberneighbours2(DataPlot(:,1),DataPlot(:,2),distancepoints),'MarkerEdgeAlpha',1)
        xlabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(1)},'Color',colorfont);
        ylabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(2)},'Color',colorfont);
        
        if uniformlimits
            xlim([limitschannels(1,channelnums(1)),limitschannels(2,channelnums(1))])
            ylim([limitschannels(1,channelnums(2)),limitschannels(2,channelnums(2))])
        end
        
        grid on

        set(gca, 'LineWidth', 2);
        set(gca,'FontWeight', 'bold')
        set(gca,'FontName','Arial')
        set(gca,'FontSize',18)
        set(gca,'Color',colorbg)

        set(gca,'XColor',colorfont)
        set(gca,'YColor',colorfont) 

        hold off;
    end

        title(analysisParam.NamesConditions_plot{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)},'Color',colorfont)% add _plot by Ye 202040409
end

fig = gcf;
fig.Color = colorbg;
fig.InvertHardcopy = 'off';
set(findall(fig,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)

saveas(fig,strrep([analysisParam.figDir filesep 'DA_HeatScatter-' colorbgplotname '-' scatterLabel '-DAPINorm-' titleplottosave '-' raworcleantitle{raworclean+1}],'.','_'),'svg')% filename ->strrep(filename,'.','_') to replace "." in the file name by Ye 20240425
% saveas(fig,[analysisParam.figDir filesep 'DA_HeatScatter-' colorbgplotname '-' scatterLabel '-DAPINorm-' titleplottosave '-' raworcleantitle{raworclean+1}],'fig')
saveas(fig,strrep([analysisParam.figDir filesep 'DA_HeatScatter-' colorbgplotname '-' scatterLabel '-DAPINorm-' titleplottosave '-' raworcleantitle{raworclean+1}],'.','_'),'png')% filename ->strrep(filename,'.','_') to replace "." in the file name by Ye 20240425

% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % print(fig,'filename','-dpdf','-r0')
% 
% saveas(fig,[analysisParam.figDir filesep 'DA_HeatScatter-' colorbgplotname '-' scatterLabel '-DAPINorm-' titleplottosave '-' raworcleantitle{raworclean+1} '.pdf'],'pdf');


    
    
