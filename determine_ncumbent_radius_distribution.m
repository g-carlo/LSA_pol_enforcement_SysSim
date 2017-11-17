% determine incumbent radius distribution

indexNonZero = min(find(tmpSave(1,:)<=0))-1; 
plot(tmpSave(1,1:indexNonZero),tmpSave(2,1:indexNonZero)/57,'.')