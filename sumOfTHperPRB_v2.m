function [sum_]=sumOfTHperPRB_v2(obj, cls,carrierNum,numOfCarriers,...
    mtxUserVsFreq, carrierArray)
        
    [users, cells]=findUser_v2(carrierNum,obj.usersList,numOfCarriers,...
        mtxUserVsFreq,carrierArray);
    
    for j=1:size(users,2)%number of users
        sum_temp=0;
        for i=1:size(users,1)%number of carriers
            if cells(i,j)~=0
                varInd = cells(i,j);
                varCellInd = carrierArray(i).cellsInCarrier(varInd);
                varCell = cls(varCellInd);
                [r, c]=find(varCell.usersList==users(i,j));
                if ~isempty(r) && ~isempty(c)
                    
                    
                    cellUEindex=varCell.usersList;
                    
                    av_sp = varCell.spectrumAvailable(i);
                    av_sp = av_sp(av_sp>0);
                    
                    
                    estimatedTHperPRB_mtx = carrierArray(i).C_mtx_extim_fsu(cellUEindex,av_sp);
                    tmp_est=estimatedTHperPRB_mtx(r,:);
                    
                    if i==carrierNum
                        if size(sum_temp,1)>size(tmp_est,1)
                            sum_temp=sum_temp+[tmp_est;(zeros(size(sum_temp,1)-size(tmp_est,1),1))];
                        else
                            sum_temp=[sum_temp;(zeros(abs(size(sum_temp,1)-size(tmp_est,1)),1))]+tmp_est;
                        end
                    end                    
                end
            end
        end
        sum_(j,:)=sum_temp;
    end