% Generate structure containing empirical particle travel distance
% distributions, raw experimental data and metadata. 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

clear

% Read in data
load('Data/rockstats.mat');
B17 = csvread('Data/B17.rockdist.csv');
B20 = csvread('Data/B20.rockdist.csv');
B28 = csvread('Data/B28.rockdist.csv');
V14 = csvread('Data/V14.rockdist.csv');
V20 = csvread('Data/V20.rockdist.csv');
V24 = csvread('Data/V24.rockdist.csv');
V39 = csvread('Data/V39.rockdist.csv');

sites = who; %generate list of site names
rocksize = {'small' 'med' 'large'};

theta=0.1; %left truncation point (m)
trunc=@(data) (data.x(data.x>theta)-theta); %anonymous function for left-truncating and re-zeroing data
CensorPoint = csvread('Data/censorpoints.csv',1,1)-theta; %right-censor points shifted by theta to match re-zeroed data

for s=1:7
    for r=1:3 
        % Store raw experimental data and metadata in rockdrop.data structure
        rockdrop{s,r}.data.SiteName = [sites{s},'_',rocksize{r}];
        rockdrop{s,r}.data.SlopeDeg=str2num(sites{s}(2:3));
        rockdrop{s,r}.data.Slope=tand(rockdrop{s,r}.data.SlopeDeg);
        rockdrop{s,r}.data.RockSize=rocksize{r};
        rockdrop{s,r}.data.IntermediateAxis=rockstats(r).IntermediateAxis.mean;
        rockdrop{s,r}.data.Units='meters';
        rockdrop{s,r}.data.x=eval([sites{s},'(:,r)']); %Raw experimental travel distances

        % Make and store empirical distributions in rockdrop.distributions structure
        [pd.X,pd.PDF,pd.CDF,pd.R,pd.P,pd.Xemp,pd.CDFemp,pd.CensorPoint] = distributions(trunc(rockdrop{s,r}.data),'dX',0.001,'CensorPoint',CensorPoint(s,r));
        pd.DistName='Empirical';
        pd.Method='distributions';
        rockdrop{s,r}.distributions=pd;
        rockdrop{s,r}.distributions.theta=theta;
    end
end

save('Data/rockdrop.mat','rockdrop');

disp('Particle travel distance data and empirical distributions saved to Data/rockdrop.mat');