set(0,'DefaultFigureWindowStyle','docked')
cd /store/users/skukies/Projects/2026UnfoldRIDE/
addpath("/store/users/skukies/StudentProjects/2024_MSc_RIDE/scripts/matlab_RIDE/RIDE_call")

data_path = "./hdf-data/";

%import data
import_data_path = data_path + "2026-07-23-Simulations-S_amp_0C_amp_10R_amp_0-NL_5-TN_200" + ".h5";
data = h5read((import_data_path),"/dataset_data");
rt = h5read((import_data_path),"/dataset_rt");
%chanlocs = load("chanlocs.mat");


%% configuration of parameters and run RIDE

for section = 1:1
    cfg = [];%initialization
    cfg.samp_interval = 10;
    cfg.epoch_twd = [-300,1600];%time window for the epoched data (relative to stimulus)
    cfg.comp.name = {'s','c','r'};%component names
    cfg.comp.twd = {[-100,300],[300,800],[-100,400]}; %time windows for extracting components, for 's' and 'c' it is raltive to stimulus, for 'r' it is relative to RT
    cfg.comp.latency = {0,'unknown',rt*10};%latency for each RIDE component
    cfg.bl =  20;
    cfg.prg = 0;
    % cfg.rwd = 20;
    % cfg.high_cutoff = 100;

    cfg = RIDE_cfg(cfg);%standardize

    results = RIDE_call(data,cfg); %run RIDE
end

%% Plot the time courses of RIDE results

%chan_index  = find(strcmpi({chanlocs.labels},'Pz'));%select which channel to plot
chan_index = 1;

%plot erp and RIDE components superimposed together
figure;RIDE_plot(results,{'erp','s','c','r'},chan_index);
