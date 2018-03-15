function out = runSnapshotSimulationQuick(conf,cellsPerCarrier,netOperatorEventGen)

global macroInd;
global picoInd;
global smallInd;

% throughput B_1

cellCapacity    = conf.operator.cellSpectralEfficiency;
cellArea        = sqrt(3)*2*conf.deployment.intersiteDistance^2; 

throughput      = zeros(length(conf.carrier.frequency),conf.deployment.numOfNetworkOperators);

for band_idx = 1:length(conf.carrier.frequency)

    for n =   1:conf.deployment.numOfNetworkOperators

        numOfActiveCells        = sum(cellsPerCarrier(band_idx,netOperatorEventGen.cell_idx_per_NO(n).list));
        throughput(band_idx,n)  = numOfActiveCells*cellCapacity*conf.carrier.availableBW(band_idx);        

    end

end

%%% save results here
out = throughput; 

end