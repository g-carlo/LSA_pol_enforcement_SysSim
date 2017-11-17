function objArray = instantiateObjectFromInputParam(structName,...
    varClassName,dim1,dim2)
    
    numOfElem = dim1*dim2;
    objArray = feval(varClassName);    
    names=fieldnames(structName)';
    noOfFields=length(names);
    
    values = cell(noOfFields,1);
    for m = 1:noOfFields
        if ~isstruct(structName.(names{m}))
            values{m} = getfield(structName,names{m});
            if iscell(values{m})
                vartmp = values{m};
            else
                vartmp = values{m}(1:dim1,1:dim2);
            end
            values{m} = reshape(vartmp,1,numOfElem);          
        end
    end
    
    
    for n = 1:numOfElem
        objArray(n) = feval(varClassName);
        for m = 1:noOfFields
            if ~isstruct(structName.(names{m}))
                objArray(n) = setfield(objArray(n),names{m},...
                    values{m}(n));
            end
        end        
    end    
end
