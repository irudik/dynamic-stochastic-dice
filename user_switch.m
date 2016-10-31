if user == 'ivan'
    warning('off','all');
    rmpath(genpath('C:\Users\irudik\Desktop\Git\robust_control'));
    rmpath(genpath('C:\Users\irudik\Desktop\Git\robust-control'));
    rmpath(genpath('C:\Users\irudik\Desktop\Git\robust-control-mat'));
    rmpath(genpath('C:\Users\irudik\Desktop\Git\lrr_mcmc_dice'));
    rmpath(genpath('C:\Users\irudik\Desktop\Git\lrr_mcmc_dice'));
    rmpath(genpath('C:\Users\irudik\Desktop\Git\lrr_mcmc_dice'));
    addpath('C:\Users\irudik\Desktop\Git\dynamic_stochastic_dice\');
    addpath(genpath('C:\Users\irudik\Desktop\Git\dynamic_stochastic_dice\smolyak'));
    addpath(genpath('C:\Users\irudik\Desktop\Git\dynamic_stochastic_dice\compecon'));
    addpath(genpath('C:\Users\irudik\Desktop\Git\dynamic_stochastic_dice\results'));
    
    warning('on','all');
   
    %%% Define directories
    dir_root = 'C:\Users\irudik\Desktop\Git\dynamic_stochastic_dice';
elseif user == 'derek'
    
elseif user == 'max'
    
end