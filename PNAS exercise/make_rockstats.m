% Generate structure containing experimental particle characteristics shown
% in SI Appendix, Table S2. 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear

% Read in data
stats{1}=csvread('Data/rockstats.small.csv');
stats{2}=csvread('Data/rockstats.medium.csv');
stats{3}=csvread('Data/rockstats.large.csv');

% Create cell arrays used to populate fields
rocksize={'small','med','large'};
fields={'Mass','Volume','IntermediateAxis','Density'};
originalunits={'g','cm^3','cm','g/cm^3'};
conversion={'./1e3','./1e6','./1e2','.*1e3'};
newunits={'kg','m^3','m','kg/m^3'};

% Populate fields of rockstats structure
for i=1:3
    stats{i}(stats{i}==0)=NaN;
    for f=1:4
        rockstats(i).RockSize=rocksize{i};
        eval(['rockstats(i).',fields{f},'.mean=nanmean(stats{i}(:,f)',conversion{f},');']); %mean values (converted to SI units)
        eval(['rockstats(i).',fields{f},'.stdev=nanstd(stats{i}(:,f)',conversion{f},');']); %standard deviations (converted to SI units)
        eval(['rockstats(i).',fields{f},'.units=newunits{f};']); %SI units for reported values
        eval(['rockstats(i).',fields{f},'.values=stats{i}(:,f)',conversion{f},';']); %raw measured values (converted to SI units)
    end
end

save('Data/rockstats.mat','rockstats')

disp('Experimental particle data saved to Data/rockstats.mat');