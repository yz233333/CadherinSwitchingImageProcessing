function DA_FindLimitsData(varargin)
%% Description
% Function that will identify all the different channels in the experiment,
% and will find the limits for each of them using all the data. It will
% also save the CleanData, that would have gotten rid of outliers

global analysisParam

if ( length(varargin) == 1 )
    varargin = varargin{:};
end

minareaaux = 0;
maxareaaux = 0;
minquantile = 0.005;
maxquantile = 0.995;

while ~isempty(varargin)
    switch lower(varargin{1})
          case {'minlimarea'}
              minareaaux = 1;
              minLimArea = varargin{2};
          case {'maxlimarea'}
              maxareaaux = 1;
              maxLimArea = varargin{2};  
          case {'minquantile'}
              minquantile = varargin{2};
          case {'maxquantile'}
              maxquantile = varargin{2};    

    otherwise
        error(['Unexpected option: ' varargin{1}])
    end
      varargin(1:2) = [];
end

load([analysisParam.savingpathforData,'/AllDataExperiment.mat'],'AllDataExperiment');


%% finding how many channels used in each well

ChannelsPresent = analysisParam.MapChannels.DifferentChannelsPresent;
numchannelspresent = length(ChannelsPresent);
PlateCoordinates = analysisParam.MapChannels.PlateCoordinates;
WellCoordinates = analysisParam.MapChannels.WellCoordinates;
ChannelCoordinates = analysisParam.MapChannels.ChannelCoordinates;
ChannelsMatrix = analysisParam.MapChannels.ChannelsCoordMatrix;


%% generating limits for the panel (per condition)

disp('Generating limits......')

limitschannels = zeros(2,length(ChannelsPresent));

for ChannelNumber = 1:length(ChannelsPresent)
    
    PlateNumAux = PlateCoordinates{ChannelNumber}(1);
    WellNumAux = WellCoordinates{ChannelNumber}(1);
    
    Data = AllDataExperiment{PlateNumAux}{WellNumAux};
    
    % Discriminate by Area if selected (commented out by Ye 20230512)
%     if minareaaux 
%         minareaaux
%         Data = Data(Data(:,end-1)>minLimArea,:);
%     end
%     
%     if maxareaaux 
%         maxareaaux
%         Data = Data(Data(:,end-1)<maxLimArea,:)
%     end
    
    % Take just the data for the specific channel
    DataAll = Data(:,2+ChannelCoordinates{ChannelNumber}(1));
    
    for imageswithchannel = 1:length(PlateCoordinates{ChannelNumber})
        
        for positionnumber = 1:analysisParam.ImagesperWell 
            PlateNumAux = PlateCoordinates{ChannelNumber}(imageswithchannel);
            WellNumAux = WellCoordinates{ChannelNumber}(imageswithchannel);
    
            Data = AllDataExperiment{PlateNumAux}{WellNumAux};
            
            % Discriminate by Area if selected
%             if minareaaux 
%                 Data = Data(Data(:,end-1)>minLimArea,:);
%             end
% 
%             if maxareaaux 
%                 Data = Data(Data(:,end-1)<maxLimArea,:);
%             end

            % Take just the data for the specific channel
            DataAll = [DataAll;Data(:,2+ChannelCoordinates{ChannelNumber}(imageswithchannel))];
            
        end
        
    end
    
    %Compute limits
    
    limitsaux = [quantile(DataAll,minquantile); ...
                quantile(DataAll,maxquantile)];    
            
    limitschannels(:,ChannelNumber) = limitsaux;
     
        
end

%%

disp('Saving Clean Data.....')

AllDataExperimentClean = cell(1,analysisParam.NumofPlates);

for platenumber = 1:analysisParam.NumofPlates
    AllDataExperimentClean{platenumber} = {};
    
for wellnumber = analysisParam.WellsWithData{platenumber}

    AllDataExperimentClean{platenumber}{wellnumber} = [];
    AllDataAux = AllDataExperiment{platenumber}{wellnumber};
    
    nchannel=length(analysisParam.Channelsnames{platenumber}{wellnumber});


% Within the channel loop, the data in AllDataAux is filtered based on
% specified limits for each channel. The line AllDataAux = 
% AllDataAux(AllDataAux(:,2+channelnumber)>limitschannels(1,ChannelsMatrix{platenumber,wellnumber}(channelnumber)),:)
% removes rows where the channel value is below a lower limit, 
% and AllDataAux = AllDataAux(AllDataAux(:,2+channelnumber)<limitschannels(2,ChannelsMatrix{platenumber,wellnumber}(channelnumber)),:); 
% removes rows where the channel value is above an upper limit.
            for channelnumber = 1:nchannel
                AllDataAux = AllDataAux(AllDataAux(:,2+channelnumber)>limitschannels(1,ChannelsMatrix{platenumber,wellnumber}(channelnumber)),:);
                AllDataAux = AllDataAux(AllDataAux(:,2+channelnumber)<limitschannels(2,ChannelsMatrix{platenumber,wellnumber}(channelnumber)),:);
                
            end
            
    AllDataExperimentClean{platenumber}{wellnumber} = AllDataAux;    
end

end
% % test
% AllDataExperiment

quantiles = [minquantile,maxquantile];
limAreas = [minLimArea,maxLimArea];
AllDataExperiment = AllDataExperimentClean;


%%

limitsbars = zeros(2,numchannelspresent);
for channelnum = 1:numchannelspresent

% MatrixToPlot=[];
% grp1=[];
meanData = [];
stdData = [];


    for imageswithchannel = 1:length(PlateCoordinates{channelnum})

            PlateNumAux = PlateCoordinates{channelnum}(imageswithchannel);
            WellNumAux = WellCoordinates{channelnum}(imageswithchannel);
            ChannelNumAux = ChannelCoordinates{channelnum}(imageswithchannel);
    
            Data = AllDataExperiment{PlateNumAux}{WellNumAux};

            meanDataAux = zeros(1,Data(end,end));
            
            for posnum = 1:Data(end,end)
                
                indposnum = find(Data(:,end)==posnum);
                meanDataAux(posnum) = mean(Data(indposnum,2+ChannelNumAux));
            end
            
            meanData = [meanData,meannonan(meanDataAux)];
            stdData = [stdData,stdnonan(meanDataAux)/sqrt(Data(end,end))];
            
  
    end
    
    limitsbars(:,channelnum) = [min((meanData-stdData/2));max((meanData+stdData/2))];
    
end



save([analysisParam.savingpathforData,'/AllDataExperimentClean'],'AllDataExperiment','limitschannels','quantiles','limAreas','limitsbars','analysisParam');

disp('Saved data')
