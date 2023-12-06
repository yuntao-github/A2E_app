% the e0 and sigma influence the inversion results. how about using them as
% variables to be inverted using a inversion method.

% update history, 
% 2023-06-13, changes: (1) changed the number-index for different age
% systems, so as to make use it for sort the age, (2) changed the sorting
% methods, taking into account the age methods and elevations,
% 2023-12-01, fixed a bug in calculating the size of the length along
% longitude, xl. 

clear
close all

%% setting up the input and output folder
codeFolder = pwd;
idcs = strfind(codeFolder,filesep);
currentFolder = codeFolder(1:idcs(end)-1);
input_data_folder = strcat(currentFolder,filesep,"data_A2Epaper");
output_folder = strcat(currentFolder,filesep,"000_output4code_version_no_rounding20231121");
if exist(output_folder,"dir")
else
    mkdir(output_folder)
end
%% below are the datasets of Tibetan transects, and constants and parameter setup
which_one2process_input = "Fitzgerald1995"; % if input "all" - run all inversion; giving a project name, the code will work on the specific project 
G_present_authority = 38.9; % if G_present_authority is defined, it will be used in the following calculations
start_time_authority = 25.5; % if start_time_authority is defined, it will be used in the following calculations and the iterations for this value will not be run
dt2plot_authority = 1.5; % if dt2plot_authority is defined, it will be used in the following calculations
% % constants and parameter setup
h_mean_authority = 4.02; % mean elevation
zmax = 80;  % maximum depth of crustal block
Nz = 161;  % 
ite = 5;   % iteration
lapse_rate = 6; % C/km
kappa = 30; % thermal diffusivity, km2/Myr 
age0 = 110; % the maximum age used for inversion 
plottingOrNot = 1; %set as 1 if plotting, 0 otherwise
removing_previous_results = 0;
misfit_mode = 0; %1 for using age misfit, 0 for (misfit_age + misfit_G/sqrt(n)), and 2 for using (misfit_age + misfit_G/n), and 3 for (misfit_age + misfit_G)
T0_authority=-12;
n_model4steps2N3 = 20; % the number of models to run for optimizing G0 and e0 
error_type = 0; outputfoldername = "error_type0";
sigma_pr = 0.5; % the priori variance of the mean exhumation rate
consider_basal_flux = 1; % 1 and 0 means considering or not 

%% below are the datasets of Tibetan transects, and constants and parameter setup
% which_one2process_input = "DhauladarRange"; % if input "all" - run all inversion; giving a project name, the code will work on the specific project 
% G_present_authority = 45; % if G_present_authority is defined, it will be used in the following calculations
% G_present_error_authority = 8;
% start_time_authority = 19.5; % if start_time_authority is defined, it will be used in the following calculations and the iterations for this value will not be run
% dt2plot_authority = 1.5; % if dt2plot_authority is defined, it will be used in the following calculations
% % % constants and parameter setup
% zmax = 60;  % maximum depth of crustal block
% Nz = 121;  % 
% ite = 5;   % iteration
% lapse_rate = 6; % C/km
% kappa = 30; % thermal diffusivity
% age0 = 19; % the maximum age used for inversion 
% plottingOrNot = 1; %set as 1 if plotting, 0 otherwise
% removing_previous_results = 0;
% misfit_mode = 0; %1 for using age misfit, 0 for (misfit_age + misfit_G/sqrt(n)), and 2 for using (misfit_age + misfit_G/n)
% T0_authority = 5;
% n_model4steps2N3 = 20; % the number of models to run for optimizing G0 and e0 
% error_type = 0; outputfoldername = "error_type0"; 
% sigma_pr = 0.5; % the priori variance of the mean exhumation rate
% consider_basal_flux = 1;
%% for testing based on the K1TB, 0-1000 m
% which_one2process_input = "K1TB";  % if input "all" - run all inversion; giving a project name, the code will work on the specific project 
% G_present_authority = 27.5; % if G_present_authority is defined, it will be used in the following calculations
% dt2plot_authority = 10; % if dt2plot_authority is defined, it will be used in the following calculations
% start_time_authority = 120; % should be n*dt2plot_authority. if Tgrd_initial_authority is defined, it will be used in the following calculations and the iterations for this value will not be run
% T0_authority = 0; % degree C, if T0_authority is defined, it will be used in the following calculations 
% zmax = 60;  % maximum depth of crustal block
% Nz = 121;  % 
% ite = 5;  % iteration 
% lapse_rate = 6; % C/km
% kappa = 30; % thermal diffusivity
% age0 = 100; % the maximum age used for inversion 
% plottingOrNot = 1;
% n_model4steps2N3 = 20; % the number of models to run for optimizing G0 and e0 
% misfit_mode = 0; %1 for using age misfit, 0 for (misfit_age + misfit_G/sqrt(n)), and 2 for using (misfit_age + misfit_G/n)
% h_mean_authority = 0; % mean elevation
% error_type = 0; outputfoldername = "error_type0";
% sigma_pr = 0.5; % the priori variance of the mean exhumation rate
% consider_basal_flux = 1;
%% 
cd (input_data_folder)
files = dir();
dirFlags = [files.isdir];
subFolders = files(dirFlags);
cd ..
%% 
    n1=3;
    n2=length(subFolders(:,1));
    which_one2process = table('Size',[n2,4],...
        'VariableNames',{'name','long','lat','G_present'}, ...
        'VariableTypes', {'string','double','double','double'});
    for ii=n1:n2  
        which_one2process{ii-2,1} = {subFolders(ii,1).name};
        datafile=strcat(input_data_folder,'/',which_one2process.name(ii-2),'/','DATA_',which_one2process.name(ii-2),'.csv');
        data = readtable(datafile);
        data = data(~isnan(data{:,5}),:);
        which_one2process{ii-2,2}=mean(data.longi);
        which_one2process{ii-2,3}=mean(data.lati);
    end
    if which_one2process_input=="all"
    else
        [Lia,Locb] = ismember(which_one2process_input, which_one2process{:,1});
        if sum(Lia)== length(which_one2process_input)
            which_one2process = which_one2process(Locb,:);
        else
            disp ("some of the transects are not in the 'which_one2process' folder generated by the 'master' code" )
        end
    end


% G_present_summary=table('Size',[n2,4],'VariableTypes', {'double','double','double','string'}, ...
%     'VariableNames',{'long','lat','G_present','ID'});
% %% get the heat flow data
% for ii=n1:n2
%     folder_ii=subFolders(ii,1).name; 
%     datafile=strcat(datafolder,'/',folder_ii,'/','DATA_',folder_ii,'.csv');
%     data = readtable(datafile);
%     G_present_summary{ii,1}=mean(data{~isnan(data{:,3}),3});
%     G_present_summary{ii,2}=mean(data{~isnan(data{:,2}),2});
%     G_present_summary{ii,4}=strcat("transect", num2str(ii-2));
% end

n_transects2process = size(which_one2process,1);

nstart = 1;
for ii=nstart:n_transects2process
    folder_ii=which_one2process.name(ii); 
    % remove previous run
    datafolder_ii = strcat(input_data_folder,'/',folder_ii);
    files_i = dir(datafolder_ii);
    for k = 1:length(files_i)
        if contains(files_i(k).name, "RUN_") && removing_previous_results==1
            rmdir(strcat(datafolder_ii,'/',files_i(k).name), 's');
        else
        end
    end
    fprintf(strcat("Running inversion for the project: ", folder_ii,'\n'));
    datafile=strcat(datafolder_ii,'/','DATA_',folder_ii,'.csv');
    demfile=strcat(datafolder_ii,'/','DEM_',folder_ii,'.csv');
    nxfile=strcat(datafolder_ii,'/','nx',folder_ii,'.csv');
    nyfile=strcat(datafolder_ii,'/','ny',folder_ii,'.csv');
    data = readtable(datafile);
    data = data(:,1:7);
%     data(data.elevation<0,:) = []; % remove borehole samples
%     data(contains(data.sys,'ZHe'),:) = []; % remove borehole samples
    demdata = load(demfile);
    data = data(~isnan(data{:,5}),:);
    data = data(data{:,5}<age0,:);
    N = length(data.lati);
    nx = load(nxfile);
    ny = load(nyfile);
    system = zeros(N,1);
    n4each_sys = zeros(10,1);
    for i=1:N
        if strcmpi(data.sys(i),'aft')
            system(i) = 2;n4each_sys(2)=n4each_sys(2)+1;
        elseif strcmpi(data.sys(i),'zft')
            system(i) = 4;n4each_sys(4)=n4each_sys(4)+1;
        elseif strcmpi(data.sys(i),'ahe')
            system(i) = 1;n4each_sys(1)=n4each_sys(1)+1;
        elseif strcmpi(data.sys(i),'zhe')
            system(i) = 3;n4each_sys(3)=n4each_sys(3)+1;
        elseif strcmpi(data.sys(i),'HAr')
            system(i) = 10;n4each_sys(10)=n4each_sys(10)+1;
        elseif strcmpi(data.sys(i),'MAr')
            system(i) = 9;n4each_sys(9)=n4each_sys(9)+1;
        elseif strcmpi(data.sys(i),'BAr')
            system(i) = 8;n4each_sys(8)=n4each_sys(8)+1;
        elseif strcmpi(data.sys(i),'KAr')
            system(i) = 7;n4each_sys(7)=n4each_sys(7)+1; 
        elseif strcmpi(data.sys(i),'THe')
            system(i) = 5;n4each_sys(5)=n4each_sys(5)+1;  
        elseif strcmpi(data.sys(i),'TFT') || strcmpi(data.sys(i),'SFT')
            system(i) = 6;n4each_sys(6)=n4each_sys(6)+1; 
        end
    end
    data.('system') = system;
    data2use = sortrows(data,[5,8,4]); %important for calculating the p
    rng(3)
    data2use.age(mod(data2use.age,dt2plot_authority)==0) = data2use.age(mod(data2use.age,dt2plot_authority)==0)+0.01*(2*randi(2)-3);
    system = data2use.system;
    n_systems = length(unique(system));
    N = length(data2use.lati);
    Age = data2use{:,5};
    Error = data2use{:,6};
    Elevation = data2use{:,4}/1000;
    sys = data2use{:,7};
    %% set up initial parameters
    if exist('h_mean_authority', 'var')
        h_mean = h_mean_authority;
    else
        h_mean = mean(demdata(:,3))/1000;
%         h_mean = mean(data2use.elevation)/1000; 
    end
    h = Elevation-h_mean;
%     h_data = Elevation-h_mean_data;
    z_mean = 0;

    if (exist('G_present_authority', 'var'))
        G_present = G_present_authority; % if G_present_authority is defined, it will be used in the following calculations
    else
        G_present = which_one2process.G_present(ii);
    end
    if (exist('G_present_error_authority', 'var'))
        G_present_error = G_present_error_authority;
    else
        G_present_error=G_present*0.1;   % error of present geothermal gradient [K/km] 
    end
    
    if exist ('T0_authority', 'var')
        T0 = T0_authority;
    else
        if h_mean < 0 % using depth
%             T0 = 10;
            disp("Please set up the T0 parameter")
        else
            T0 = 30-lapse_rate*h_mean;  % surface temperature using a temperature - elevation model 
        end
    end
    % set up e0
        n_Age=length(Age);
        zz=zeros(n_Age,1);
    for i=1:n_Age
        if strcmpi(sys(i),'aft')
            zz(i)=110/G_present;
        elseif strcmpi(sys(i),'zft')
            zz(i)=220/G_present;
        elseif strcmpi(sys(i),'ahe')
            zz(i)=70/G_present;
        elseif strcmpi(sys(i),'zhe')
            zz(i)=180/G_present;
        elseif strcmpi(sys(i),'HAr')
            zz(i)=450/G_present;
        elseif strcmpi(sys(i),'MAr')
            zz(i)=400/G_present;
        elseif strcmpi(sys(i),'BAr')
            zz(i)=350/G_present;
        elseif strcmpi(sys(i),'KAr')
            zz(i)=350/G_present;  
        elseif strcmpi(sys(i),'THe')
            zz(i)=200/G_present;
%         elseif strcmpi(sys(i),'TFT') || strcmpi(sys(i),'SFT')
%             zz(i)=300/G_present;            
        end
    end
    if exist('e0_pr_authority', 'var')
        e0_pr = e0_pr_authority;
    else
%         e0_pr = mean((zz+mean(Elevation))./Age);  % assumed initial erosion rate [km/my]
        e0_pr = mean(zz./Age);  % assumed initial erosion rate [km/my]
    end

    tmax_age=max(Age);
    tmin_age=min(Age);
%     if tmax_age <= 10 
%         tmax = tmax_age+10; % the starting time of calculation, which must be older than the maximum age of the sample
%     else
%         tmax = tmax_age*1.4; % the starting time of calculation, which must be older than the maximum age of the sample
%     end
    if exist('dt2plot_authority', 'var')
        dt2plot = dt2plot_authority;
    else
        dt2plot = ceil(median(Age)*0.2);
%         if tmax_age < 4  && tmin_age < 1 && length(Age)>5
%             dt2plot = .5; %Ma 
%             nt4each_step = 5;
%         elseif tmax_age>=4 && tmax_age<15 || tmin_age < 2
%             dt2plot = 1; %Ma
%             nt4each_step = 5;
%         elseif tmax_age>=15 && tmax_age<30 || tmin_age < 4
%             dt2plot = 2; %Ma
%             nt4each_step = 10;
%         elseif tmax_age>=30 && tmax_age<50 || tmin_age< 6
%             dt2plot = 3; %Ma  
%             nt4each_step = 30;
%         elseif tmax_age>=50 && tmax_age<80 || tmin_age < 20
%             dt2plot = 5; %Ma
%             nt4each_step = 20;
%         else 
%             dt2plot = 10; %Ma  
%             nt4each_step = 20;
%         end
    end
    if exist('start_time_authority', 'var')
        tmax = start_time_authority;
    else
        tmax = tmax_age+10;
        tmax = (fix(tmax/dt2plot)+1)*dt2plot;
    end
    
%     dt2plot = dt2plot*2; % can we used a larger time step to avoid
%     negative rates?
%     dt = dt2plot/nt4each_step;  % step size of age. A large dt may result in large misfit.   
%     dz = zmax*1000/Nz;
%     kappa=kappa*(1000^2.);
%     dt = (dz^2.)/kappa/0.02;
    Nt = fix(tmax/dt2plot);  % number of steps of calculation   

    %% test the effect of geothermal gradient on the inversion results and find the best G0
    fprintf('Running inversion with different G0 \n');
    e0 = e0_pr;  % assumed initial erosion rate [km/my]
    sigma = e0*sigma_pr;  % error of the assumed initial erosion rate [km/my]

    misfit1=zeros(n_model4steps2N3,4); % misfit of age, misfit of thermal gradient, thermal gradient input and e0 input
    if exist('Tgrd_initial_authority', 'var')
        Tgrd = Tgrd_initial_authority;
    else
        G_fraction = linspace(0.5, 1.2, n_model4steps2N3);
        for irun = 1:n_model4steps2N3
            Tgrd=G_present*G_fraction(irun); 
            fprintf('  G0 = %d\n', Tgrd);
            cd (codeFolder)
            inversion_clc2;
            cd (output_folder)
            if exist(folder_ii,'dir')
            else
                mkdir(folder_ii)
            end
            cd (folder_ii)
            if exist(outputfoldername,'dir')
            else
                mkdir(outputfoldername)
            end
            cd (outputfoldername)
            dirname=strcat('masterRUN_G', num2str(irun));
            mkdir(dirname);
            cd (dirname)
            misfit1(irun,1)=misfit_age_postior;
            misfit1(irun,2)=misfit_G_postior;
            misfit1(irun,3)=Tgrd;
            misfit1(irun,4)=e0;
%             clear plottingOrNot;
            save run.mat
            cd (currentFolder)
        end
%         if misfit_mode == 1
%             [M,Index] = min(misfit1(:,1)); 
%         elseif misfit_mode == 0
%             [M,Index] = min(sum(misfit1(:,1)+misfit1(:,2)/n,2));
%         elseif misfit_mode == 2
%             [M,Index] = min(sum(misfit1(:,1)+misfit1(:,2)/sqrt(n),2));
%         else
%             disp("please defining a misfit mode")
%             return
%         end
        [M,Index] = min(misfit1(:,2)); 

        Tgrd=misfit1(Index,3);
    end

    %% test the effect of prior e and error on the inversion results, it seems the effect is minor

    fprintf('Running inversion with different e0 \n');
    misfit2=zeros(n_model4steps2N3,4); % misfit of age, misfit of thermal gradient, thermal gradient input and e0 input
    e_fraction = linspace(0.1, 2.0, n_model4steps2N3);
    for irun = 1:n_model4steps2N3
        if exist('G_present_authority', 'var') && exist('e0_pr_authority', 'var')
            e0 = e0_pr_authority;
        else
            e0=e0_pr*e_fraction(irun);   % change between 0.2*e0_pr and 2*e0_pr
        end
        sigma=e0*sigma_pr;
        formatSpec = '  e0 = %4.2f, with an error = %4.2f \n';
        fprintf(formatSpec,e0,sigma)
        cd (codeFolder)
        inversion_clc2;
        cd (output_folder)
        if exist(folder_ii,'dir')
        else
            mkdir(folder_ii)
        end
        cd (folder_ii)
        if exist(outputfoldername,'dir')
        else
            mkdir(outputfoldername)
        end
        cd (outputfoldername)
        dirname=strcat('masterRUN_e', num2str(irun));
        mkdir(dirname);
        cd (dirname)
        misfit2(irun,1)=misfit_age_postior;
        misfit2(irun,4)=e0;
        misfit2(irun,2)=misfit_G_postior;
        misfit2(irun,3)=Tgrd;
%         clear plottingOrNot;
        save run.mat
        cd (currentFolder)
    end

    if plottingOrNot == 0
    else
        cd(codeFolder)
        plot_figures_during_batchRUN4master_no_rounding20224paper1
    end

    fprintf('------------------\n');
end
