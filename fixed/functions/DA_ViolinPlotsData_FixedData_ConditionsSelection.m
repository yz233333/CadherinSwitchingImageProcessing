
function DA_ViolinPlotsData_FixedData_ConditionsSelection(ConditionsSelection,titleplottosave,varargin)

global analysisParam;

if ( length(varargin) == 1 )
    varargin = varargin{:};
end

raworclean = 0;
blackBG = 0;
stdmean = 1;
angleticks = 0;
uniformlimits = 1;

windowdims = [10 10 2000 1000];

raworcleantitle = {'RAW','CLEAN'};

Title=[titleplottosave,' ', raworcleantitle{raworclean+1}];

% channelnums = analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,1),ConditionsSelection(2,1)};
channelnums =[];

colordefs=[];
ylimdefs=[];

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
          case {'ylimdef'}%by Ye 20230316
              ylimdefs=varargin{2};
          case {'ylimdeflow'}%by Ye 20240907
              ylimdefsLow=varargin{2};
          case {'colordef'}%by Ye 20230316
              colordefs=varargin{2};
          case {'windowdims'}%by Ye 20230316
              windowdims = [varargin{2}];

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
times = 1:nCon;

FindChannelsinConditions = zeros(nCon,nChan);

for ii = 1:nChan
    channelinterest = channelnums(ii);
    
    for jj = 1:nCon       
        auxvar = find(channelinterest == analysisParam.MapChannels.ChannelsCoordMatrix{ConditionsSelection(1,jj),ConditionsSelection(2,jj)});
        
        if isempty(auxvar)
            error(['DA_BarPlotsData_FixedData_ConditionsSelection: Selected condition ',num2str(jj),' does not contain data for channel ',analysisParam.MapChannels.DifferentChannelsPresent{ii}]);
        
        else
            FindChannelsinConditions(jj,ii) = auxvar;
            
        end
    end
end

%% Violin plot
figure;
set(gcf,'Position',windowdims)
% colourclusters = {colorconvertorRGB([236,28,36]),colorconvertorRGB([249,236,49]),colorconvertorRGB([64,192,198])};
% colourclusters = {'g','r','b'};

if nChan<4
colourclusters = [colorconvertorRGB([64,192,198]);colorconvertorRGB([185,82,159]);colorconvertorRGB([249,236,49])];
% colourclusters = [[0,0,0];[0,0,0];[0,0,0]];%all black (by Ye20240830
% colourclusters = colourclusters(ordercolors,:);

else
    colourclusters = distinguishable_colors(nChan,{'w','k'});
    colourclusters = colourclusters(ordercolors,:);
    
end

colourclusters = [colorconvertorRGB([64,192,198]);colorconvertorRGB([185,82,159]);[1,0.8,0];...
                  [0.6,0.55,0.55];[0.7,0.4,1];...
                  colorconvertorRGB([185,104,23]);colorconvertorRGB([102,110,21]);colorconvertorRGB([0,130,165])];%Ye

disp('Plotting....')
tiledlayout(nChan,1,'TileSpacing','compact')

alinearrelation = 10;
blinearrelation = 5;

for channelnum = 1:nChan
nexttile

celln = 0;
xticksnum = [];
daystickslabels = {};
plotshandle = [];
legendhandle = {};

    meanData = [];
    stdData = [];  
    for condnum = 1:nCon
       
        %Set the interval in which the violin will spread
        daytickmin = alinearrelation*(times(condnum)-times(1));
        daytickmax = alinearrelation*(times(condnum)-times(1))+blinearrelation;
        daytick = (2*alinearrelation*(times(condnum)-times(1))+blinearrelation)/2;

        DataPlot =AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(:,2+FindChannelsinConditions(condnum,channelnum));
        
%         disp('Computing pdf.....')
        % Obtain the distribution of the data
        [pdf,y] = ksdensity(DataPlot);
        
%         disp('Computed!.....')
        
        % Obtain higher accuracy so that the violin looks nicer:
        ymax = max(DataPlot);%max(y);
        ymin = min(DataPlot);%min(y);
        points = linspace(ymin,ymax,1000);
        
%         disp('Evaluating pdf.....')
        %Repeat ksdensity
        [pdf,y] = ksdensity(DataPlot,points);
%         disp('Evaluated!.....')
        
        %Now, for each point in y, we need to find a change of coordinates
        %between [pdfmin,pdfmax] and [daytickmin,daytickmax]
        
        pdfmax = max(pdf);
        pdfmin = min(pdf);
        
        Aright = (daytick-daytickmax)/(pdfmin-pdfmax);
        bright = daytickmax - Aright*pdfmax;
        
        Aleft = (daytick-daytickmin)/(pdfmin-pdfmax);
        bleft = daytickmin - Aleft*pdfmax;
        
        pdfright = Aright*pdf+bright;
        pdfleft = Aleft*pdf+bleft;
        
        startPoint = 2;

% current violin:       
        if ~isempty(colordefs)%by Ye 20230316
            colorIndex=colordefs(channelnums(channelnum));
            plot(pdfright+startPoint,y,'Color',colourclusters(colorIndex,:),'LineWidth',2);
        else
            plot(pdfright+startPoint,y,'Color',colourclusters(channelnum,:),'LineWidth',2);
        end
        hold on

        if ~isempty(colordefs)%by Ye 20230316
            colorIndex=colordefs(channelnums(channelnum));
            plot(pdfleft+startPoint,y,'Color',colourclusters(colorIndex,:),'LineWidth',2);
        else
            plot(pdfleft+startPoint,y,'Color',colourclusters(channelnum,:),'LineWidth',2);
        end

% original: 
%         plot(pdfright,y,'Color',colourclusters(channelnum,:),'LineWidth',2);
%         hold on
%         plot(pdfleft,y,'Color',colourclusters(channelnum,:),'LineWidth',2);
        
        if stdmean
            meanDataAux = zeros(1,AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(end,end));
            
            for posnum = 1:AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(end,end)
                
                indposnum = find(AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(:,end)==posnum);
                meanDataAux(posnum) = mean(AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(indposnum,2+FindChannelsinConditions(condnum,channelnum)));
            end
            
            meanData = meannonan(meanDataAux);
            stdData = stdnonan(meanDataAux)/sqrt(AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(end,end));
            
        else
                    
            meanData = meannonan(AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(:,2+FindChannelsinConditions(condnum,channelnum)));
            stdData = stdnonan(AllDataExperiment{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)}(:,2+FindChannelsinConditions(condnum,channelnum)));
            
        end

% Current: 
        if ~isempty(colordefs)%by Ye 20230316
            colorIndex=colordefs(channelnums(channelnum));
            errorbar(daytick+startPoint, meanData,stdData/2,'LineWidth',2,'Color',colourclusters(colorIndex,:));
        else
            errorbar(daytick+startPoint, meanData,stdData/2,'LineWidth',2,'Color',colourclusters(channelnum,:));
        end
% original: 
%         errorbar(daytick, meanData,stdData/2,'LineWidth',2,'Color',colourclusters(channelnum,:));

        
        
        % Set appearance of plot:
        xticksnum = [xticksnum, daytick];
        daystickslabels{condnum} = analysisParam.NamesConditions_plot{ConditionsSelection(1,condnum)}{ConditionsSelection(2,condnum)};% add _plot by Ye 202040409
        
        
%             legendhandle{daynum} = ['(',num2str(ids(daynum)),') ',mutregimes{daynum},' D',num2str(days{daynum}),' ', CHIRcond, ' ',LGKcond];
%         if daynum >1
%            plot([daytickaux,daytick], [meanaux,mean(DataPlot)],'Color',colors(conditionnum,:),'LineWidth',1.5) 
%         end
        
        
        
        
    end
    
    
%         xlabel('Hours post treatment')
        ylabel(analysisParam.MapChannels.DifferentChannelsPresent{channelnums(channelnum)},'Color',colorfont);

    if uniformlimits
    ylim([limitschannels(1,channelnums(channelnum)),limitschannels(2,channelnums(channelnum))])
    end


    xlim([0,max(xticksnum)+blinearrelation+startPoint])
    [xticksordered,indicesxticks] = sort(xticksnum);
    xticks(unique(xticksordered)+startPoint)
    if channelnum<nChan
        xticklabels((cell(1,nCon)));
    else
        xticklabels((daystickslabels(1:nCon)));
    end

    if ~isempty(ylimdefs) && isempty(ylimdefsLow)%by Ye 20230316
        ylim([0 ylimdefs(channelnums(channelnum))])
    elseif ~isempty(ylimdefs) && ~isempty(ylimdefsLow)
        ylim([ylimdefsLow(channelnums(channelnum)) ylimdefs(channelnums(channelnum))])
    end
    
    if channelnum==1
    title(Title,'Color',colorfont)
    end
    
    xtickangle(angleticks)
    
    set(gca, 'LineWidth', 2);
    set(gca,'FontWeight', 'bold')
    set(gca,'FontName','Arial')
    set(gca,'FontSize',18)
    set(gca,'Color',colorbg)
    set(gca,'Color',colorbg)
    set(gca,'XColor',colorfont)
    set(gca,'YColor',colorfont)  
    
end



fig = gcf;
fig.Color = colorbg;
fig.InvertHardcopy = 'off';
set(findall(fig,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','LineWidth'),'LineWidth',2)

%saveas(fig,[analysisParam.figDir filesep 'DA_ViolinPlots-',colorbgplotname,'-AllGenes-DAPINorm-', titleplottosave, '-', raworcleantitle{raworclean+1}],'eps')
saveas(fig,strrep([analysisParam.figDir filesep 'DA_ViolinPlots-',colorbgplotname,'-DAPINorm-', titleplottosave, '-', raworcleantitle{raworclean+1}],'.','_'),'svg')% filename ->strrep(filename,'.','_') to replace "." in the file name by Ye 20240425
saveas(fig,strrep([analysisParam.figDir filesep 'DA_ViolinPlots-',colorbgplotname,'-DAPINorm-', titleplottosave, '-', raworcleantitle{raworclean+1}],'.','_'),'png')% filename ->strrep(filename,'.','_') to replace "." in the file name by Ye 20240425
% saveas(fig,[analysisParam.figDir filesep 'DA_ViolinPlots-',colorbgplotname,'-DAPINorm-', titleplottosave, '-', raworcleantitle{raworclean+1}],'fig')

% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % print(fig,'filename','-dpdf','-r0')
% 
% saveas(fig,[analysisParam.figDir filesep 'DA_ViolinPlots-',colorbgplotname,'-DAPINorm-', titleplottosave, '-', raworcleantitle{raworclean+1}, '.pdf'],'pdf');
