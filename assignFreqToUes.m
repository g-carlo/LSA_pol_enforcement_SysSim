function mtxUserVsFreq = assignFreqToUes(totNumUEs,carrierArray,confUeRandom,seedNumber )

% in order to select users that can use a given carrier
    numOfCarriers = length(carrierArray);
    seed1 = seedNumber .* rand(1, numOfCarriers); %seed: selection of users on different frequencies 

%%  randomly assign frequencies to users:
%   for example: 3 frequncy, [7,4,3] users of frequencies
%
%   mtxUserVsFreq                                   mat
%   users:  |  1   2   3   4   5   6   7            users:  |  1   2   3   4   5   6   7
%   -------------------------------------           -------------------------------------
%   f1:     |  1   1   1   1   1   1   1            f1:     |  1   2   3   4   5   6   7
%   f2:     |  0   1   1   1   1   0   0            f2:     |  2   3   4   5   0   0   0
%   f3:     |  1   0   0   1   0   0   1            f3:     |  1   4   7   0   0   0   0

    maxLength = totNumUEs;
    
    noOfCarriers = length(carrierArray);
    mtxUserVsFreq=zeros(noOfCarriers,maxLength);
    mtxUserVsFreq(1,:)=1;
    
    mat=zeros(noOfCarriers,maxLength);
    mat(1,1:maxLength)=1:maxLength;
    
    carrierArray(1).usersInCarrier = 1:totNumUEs;
    
    %switch configV(1).general.usersRandomization
    switch confUeRandom
        case 1
            for ii=2:noOfCarriers
                % to achive quasi stationary positions of users there is
                % developt an algorithm.
                
                % if the number of users on f2 is smaller than a half of number
                % of users on f1, then:
                % for ex.: 
                % nOfU on f1 = 7
                % nOfU on f2 = 3
                %                                 ¡ nOfU on f2                     ¡nOfU on f1
                % A:     |  0.7  |  0.3   |  0.1  ||       |       |       |       ||
                %                 
                %                    ¡nOfU on f2      ¡nOfU on f1
                % B:     |_0_|_1_|_1_||_0_|_1_|_0_|___||
                %                                   ^ from this place to nOfU on f1 write zeros
                
                
                %else:        
                % for ex.: 
                % nOfU on f1 = 7
                % nOfU on f2 = 4
                %                                          ¡nOfU on f2             ¡nOfU on f1 on f1
                % A:     |  0.7  |  0.3  |   0.1  |   0.6  ||      |       |       ||
                %                 
                %                        ¡nOfU on f2  ¡nOfU on f1
                % B:     |_0_|_1_|_1_|_0_||_1_|_0_|_0_||_1_||
                %                                ^ - this one, or in other
                %                                cases, more of them, needs to
                %                                be translated to the zeros
                %                                places in front of ¡nOfU
                %                                on f1 mark. This has to
                %                                be done randomly.
                
                
                RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed1(ii)));
                % tempLength=length(uesV{i});
                
                tempLength = round(maxLength * carrierArray(ii).uePercentage);
                               
                A=(rand(1,tempLength));
                B=zeros(1,2*tempLength);
                for j=1:tempLength
                    if A(j)>=0.5
                        B(j*2-1)=0;
                    else
                        B(j*2-1)=1;
                    end;
                    B(j*2)=~B(j*2-1);
                end;
                mtxUserVsFreq(ii,1:maxLength)= B(1:maxLength);
                if maxLength<tempLength*2
                    for j=maxLength+1:tempLength*2
                        if B(j)==1
                            flag=0;
                            while flag==0
                                k=randi([1,maxLength],1);
                                if mtxUserVsFreq(ii,k)==0
                                    mtxUserVsFreq(ii,k)=1;
                                    flag=1;
                                end
                            end
                        end
                    end
                end
                [C,D]=sort(mtxUserVsFreq(ii,:),2,'descend');
                mat(ii,:)=C.*D;
                carrierArray(ii).usersInCarrier = mat(ii,1:tempLength); 
            end
        case 2
            for ii=2:noOfCarriers
                %RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed1(i)));
                tempLength = round(maxLength * carrierArray(ii).uePercentage);
                [~, IX]=sort(rand(1,maxLength));
                [B, ~]=sort(IX(1:tempLength));
                mtxUserVsFreq(ii,B)=1;
                mat(ii,1:tempLength)=B;
                carrierArray(ii).usersInCarrier = B; 
            end
    end
end
