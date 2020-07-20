% Calculate bootstrapped standard error of optimized lomax parameters A, B,
% and B/|A| for all sites and particle sizes. 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
%
% Written by Danica Roth, Colorado School of Mines, May 2019.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Check whether bootstrap.mat file exists already to avoid rerunning on accident
files=dir('Data/*mat');
filenames={files.name};

% If bootstrap.mat does exist, ask for confirmation before rerunning
if ismember('bootstrap.mat',filenames)
    wmsg='Bootstrap error algorithm has already been run. Rerunning it may take significant time.';
    warning(wmsg);
    
    disp('Rerun bootstrap algorithm? (Y/N) *required to proceed');
    drawnow;
   
    prompt=['\n',wmsg,'\n\nRerun anyway? (Y/N) \n\n'];
    runboot = input(prompt,'s');
    
    if ismember(runboot,{'N' 'n'}) || isempty(runboot) %if no input or N/n is entered
        disp('Bootstrap algorithm canceled. Previously saved error values will be used instead.');

    elseif ~ismember(runboot,{'N' 'n' 'Y' 'y' ''}) %if input other than Y/y or N/n is entered
        disp('Response not recognized. Bootstrap algorithm canceled. Rerun section if this was a mistake.');
    end
end  
% If bootstrap.mat does not exist, or if Y/y was entered above, proceed with calculation
if ~ismember('bootstrap.mat',filenames) || [ismember('bootstrap.mat',filenames) && ismember(runboot,{'Y' 'y'})]

    % Read in data
    load('Data/results.mat','X','A_opt','B_opt','theta','CensorPoint','fitmethod','conditions','site');
      
    for s=1:7 %slope identifier    
        for r=1:3 %rock size identifier
            
            disp(['[',datestr(datetime),'] Site ',site{s},' starting...']);

            % Get precision (number of decimal places) of travel distance X
            % measurements
            [~,ca]=find(num2str(X{s,r})=='.');
            prec=max(size(num2str(X{s,r}),2)-ca);

            % Run bootstrap error calculation algorithm
            [Ase(s,r),Bse(s,r),xrse(s,r),A_se{s,r},B_se{s,r},xr_se{s,r}] = lomaxbootstrap(A_opt(s,r),B_opt(s,r),theta,CensorPoint(s,r),prec,'fitmethod',fitmethod, 'conditions',conditions);%,'MaxIter',10,'MaxFunEvals',10,'iboot',10);  %uncomment to test
         end
    end 

    save('Data/bootstrap.mat','A_se','B_se','xr_se','Ase','Bse','xrse');

    disp('Bootstrapped standard error arrays saved to Data/bootstrap.mat');

end