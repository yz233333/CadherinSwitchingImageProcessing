function [meanRatio,stdRatio,steRatio] = calcMeanStdSte(input,nChannels,nConditions,nColoniesPerCondition,nT)
% Calculate mean, standard deviation, standard error
% Input
% rawData; nChannels; nConditions; nColoniesPerCondition(array); nT

rawData   = zeros(nT,3);
if nT == 1 % for fixed data (using array instead of cell)
    meanRatio = zeros(nChannels,nConditions);
    stdRatio  = zeros(nChannels,nConditions);
    steRatio  = zeros(nChannels,nConditions);
    
    for cc = 1:nChannels
        for ii = 1:nConditions
            if ~isempty(input{1,ii})
                for jj = 1:nColoniesPerCondition(ii)
                    rawData(:,jj)=input{jj,ii}(:,cc);
                end
        
                % mean
                meanRatio(cc,ii)= mean(rawData,2,'omitnan');
        
                % standard deviation
                stdRatio(cc,ii) = std(rawData);
        
                % standard error
                steRatio(cc,ii) = stdRatio(cc,ii)./sqrt(nColoniesPerCondition(ii));
            end
        end
    end

else % for live data with nT>1 (using cell)
    meanRatio = cell(nChannels,nConditions);
    stdRatio  = cell(nChannels,nConditions);
    steRatio  = cell(nChannels,nConditions);
    
    for cc = 1:nChannels
        for ii = 1:nConditions
            if ~isempty(input{1,ii}(:,cc))   
                for jj = 1:nColoniesPerCondition(ii)
                   rawData(:,jj)=input{jj,ii}(:,cc);
                end
        
                % mean
                meanRatio{cc,ii}=mean(rawData,2,'omitnan');
        
                % standard deviation
                for tt = 1:nT
                    stdRatio{cc,ii}(tt,1) = std(rawData(tt,:));
                end
        
                % standard error
                steRatio{cc,ii} = stdRatio{cc,ii}./sqrt(nColoniesPerCondition(ii));
            end
        end
    end
end