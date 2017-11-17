function seedNumber = getSeedNumber(varGeneratorMode, varIsEiRP, ...
    varFileName)
    %%% random number generator
    switch varGeneratorMode
        case  1
            % clock returns a six element date vector containing the current
            % time and date in decimal form:  
            % [year month day hour minute seconds]

            seedNumber = 100*sum(clock);

            % RandStream is a class that represents a random number stream. In
            % the following stament a new object of this class is instantiated 
            % by the constructor RandStream('mt19937ar','seed',seedNumber). 
            % The first argument indicates the algorithm 'mt19937ar' that will 
            % be used to generate the stream, the second argument indicates the 
            % input parameter name that will be passed to the algorithm, and the 
            % third one indicates the parameter value. 

            % once the constructor creates the new stream, this is passed to
            % the function setGlobalStream to specify this new stream as the
            % global stream that will be used by functions such as: rand, randi,
            % randn.         

            if ~varIsEiRP
                % creates the file with the name specified as first argument
                % and saves the values of variables that are specified as the
                % following arguments. In this case, value of variable 
                % 'seedNumber' is saved. The file format is MATLAB specific, so
                % that the variable values are recovered using the MATLAB
                % function load.            
                save(varFileName,'seedNumber');
            end;
        case  2
            load(varFileName);        
        otherwise
            seedNumber = 100*sum(clock);        
    end
end