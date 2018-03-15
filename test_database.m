% test database
clear classes;
clear all;
close all;

RandStream.setGlobalStream(RandStream('mt19937ar','seed',0));

conf = getInputParameters();
num_of_inc = conf.deployment.numOfIncumbents;

% generate database and incuments 
paramVar.deployment = conf.deployment;
paramVar.database = conf.database;
database_var = dataBase_Cls(conf.database.numOfLSAChannels,conf.deployment.numOfIncumbents,paramVar);      

n_steps = 10000;
inc_var(1) = incumbent_Cls(1,conf.database.lambdaIdle,conf.database.lambdaActive);
inc_var(2) = incumbent_Cls(2,conf.database.lambdaIdle,conf.database.lambdaActive);

% initialize database and incuments

database_var.initializeDataBase( inc_var);

inc_var(1).initialize( 100 , database_var);
inc_var(2).initialize( 200 , database_var);

in_state = zeros(1,n_steps) ;

for n= 1:n_steps
    
    for i_inc = 1:num_of_inc
        inc_var(i_inc).update();
    end
    in_state(1,n) = inc_var(1).activityStatus;
    in_state(2,n) = inc_var(2).activityStatus;
    if any( n == [170 753 947 1783] )
        disp('ciao');
    end
    
end