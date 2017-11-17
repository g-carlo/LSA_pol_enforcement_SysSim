%%% SON SIMULATOR: findUser
%
%   Author: Almir Maric
%
%   Date:   August 2013
%
%   Description: return vector of cells which the user belongs to at different frequencies, and return the user's numbers at those frequencies

function [usersV,cellsV]=findUser_v2(freq,users,numOfCarriers,mtxUserVsFreq,...
    carrierArray)


    numUes = length(users);
    numCols = length(mtxUserVsFreq(freq,:));
    
    %Rows: number of carriers
    %Columns: number of users that are attached to the current cell
    
    usersV=zeros(numOfCarriers,numUes);
    cellsV=zeros(numOfCarriers,numUes);
    
    userAtFirsFreq=zeros(size(users));
    for k=1:numUes
        ueInd = users(k);
        for i=1:numCols
            if sum(mtxUserVsFreq(freq,1:i)) == ueInd
                userAtFirsFreq(k)=i;
                break;
            end
        end
    end

    
    for i=1:numOfCarriers
        for k=1:numUes
            tmpUser=0;
            tmpCell=0;
            if mtxUserVsFreq(i,userAtFirsFreq(k))~=0
                tmpUser=sum(mtxUserVsFreq(i,1:userAtFirsFreq(k)));
            end;
            if tmpUser~=0
                varCellsInCarrier = carrierArray(1).cellsInCarrier;
                varTmp = dataVar.userVsCell_allocationMtx(varCellsInCarrier,tmpUser);
                [~, bb]=sort(varTmp);
                
                %[~, bb]=sort(clsV{1,i}(1,1).dataVar.userVsCell_allocationMtx(:,tmpUser));
                tmpCell=bb(end,1);
            end;
            usersV(i,k)=tmpUser;
            cellsV(i,k)=tmpCell;
        end
    end