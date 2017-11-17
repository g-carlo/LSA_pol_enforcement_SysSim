%%% SON SIMULATOR: Final data analysis function class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%   Last Modification: March, 2015
%

function [minVal, minInd] = getNextEventIdx(eventQueue)

    %%% UE SINR and spectral efficiency 
    timeVector  = [eventQueue.time];
    
    [minVal] = min(timeVector);
    minInd = find(timeVector == minVal);

%     outputStruct = struct( 'eventType', cell(1,numOfEvents),'eventTime', cell(1,numOfEvents),'eventActivity', cell(1,numOfEvents),'eventID', cell(1,numOfEvents));
%     for ii = 1:length(minInd)
%         [eventType, eventTime, eventActivity, eventID] = eventQueue(minInd(ii)).getProperties();
%         outputStruct(ii).eventType = eventType;
%         outputStruct(ii).eventTime = eventTime;
%         outputStruct(ii).eventActivity = eventActivity;
%         outputStruct(ii).eventID = eventID;
%     end
%     eventQueue(minInd).getProperties();
    
    

end






