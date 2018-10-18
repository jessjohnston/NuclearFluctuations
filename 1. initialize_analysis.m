%% refresh matlab
close all;
clear all;
clc;

% set location of dv moviesstrain_all=dir(rootpath);
rootpath='C:\nuclei\data\';
strain_all=dir(rootpath);
strain_names={strain_all(3:end).name};
strain_selected=strain_names(cellfun(@(x)~isempty(x),...
    strfind(strain_names,'cdc25')));
%%

for ii=1:length(strain_selected)
    strainpath=strain_selected{ii};
   % movienames=dir(fullfile(rootpath,strainpath,'*.dv'));
    movienames=dir(fullfile(rootpath,strainpath,'cdc25_2.5_10min_05_R3D.dv'));%% edit here for 1 movie
    for jj=1:length(movienames)
        nm=nucmem3(fullfile(rootpath,strainpath,movienames(jj).name));
        % load movie
        nm.endframe=101;
        nm.loadmovie(101*10);
        % select the threshold for analysis
        nm.get_centroid_firstframe;
        nm.remove_badcentroid(1);
        nm.remove_badcentroid(50);
        nm.remove_badcentroid(101);
        % initialze
        nm.initialize();
        %
        nm.pickorientation;
        %
        nm.save_contour(1);
    end
end

