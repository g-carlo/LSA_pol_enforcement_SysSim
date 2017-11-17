function dispSimTime( n,N ) 
    
    if fix(n/N*100/10) > fix((n-1)/N*100/10) 
        
        disp([ num2str(fix(n/N*100/10)*10) '% completed'])
        
    end


end