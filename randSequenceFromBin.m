%%% SON SIMULATOR: random draw of a sequence of K numbered marbles among a
%%% set of N numbered murbles
%
%   Author: Carlo Galiotto
%
%   Date:   September 2010
%
%   Description: starting from a bin with N marbles numbered from 1 to N,
%   this function returns a sequence of the indices of K marbles randomly
%   drawn from the bin.


function y = randSequenceFromBin(k,n)

    if k >n 
       
        error('k must be greater than 0 and not greater than n');
        
    end
    
    if k <= 0 
        
        error('k must be greater than 0 and not greater than n');
        
    elseif k == 1;
        
        y = randi(1,1,[1 n]);
        
    else
        
        bin = rand(1,n);
    
        [~, bin_sorted_ind] = sort(bin);
    
        y = bin_sorted_ind(1:k);
                
    end




end 


