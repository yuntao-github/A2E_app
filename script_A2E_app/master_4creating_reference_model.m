% the e0 and sigma influence the inversion results. how about using them as
% variables to be inverted using a inversion method.

% update history, 
% 2023-06-13, changes: (1) changed the number-index for different age
% systems, so as to make use it for sort the age, (2) changed the sorting
% methods, taking into account the age methods and elevations

clear
close all

%% setting up the input and output folder
codeFolder = pwd;
idcs = strfind(codeFolder,filesep);
currentFolder = codeFolder(1:idcs(end)-1);
input_data_folder = strcat(currentFolder,filesep,"data_A2Epaper");
output_folder = strcat(currentFolder,filesep,"000_output4code_version_no_rounding2023");
if exist(output_folder,"dir")
else
    mkdir(output_folder)
end

%% for testing based on the Fitzgerald1995 data
which_one2process_input = "Fitzgerald1995"; % if input "all" - run all inversion; giving a project name, the code will work on the specific project 
G_present_authority = 38.9; % if G_present_authority is defined, it will be used in the following calculations
dt2plot_authority = 2.5; % if dt2plot_authority is defined, it will be used in the following calculations
e0_pr_authority = 0.5; % if e0_pr_authority is defined, it will be used in the following calculations
e0_pr_err_authority = 0.15; % if e0_pr_err_authority is defined, it will be used in the following calculations
Tgrd_initial_authority = 24; % if Tgrd_initial_authority is defined, it will be used in the following calculations and the iterations for this value will not be run
start_time_authority = 25; % should be n*dt2plot_authority. if Tgrd_initial_authority is defined, it will be used in the following calculations and the iterations for this value will not be run
T0_authority = -12; % degree C, if T0_authority is defined, it will be used in the following calculations 
h_mean_authority = 4.02; % mean elevation
zmax = 80;  % maximum depth of crustal block
Nz = 161;  % 
ite = 5;  % iteration
lapse_rate = 6; % C/km
kappa = 30; % thermal diffusivity, km2/Myr 
age0 = 110; % the maximum age used for inversion 
plottingOrNot = 1;
outputfolder_prefix = "RefModel"; 
misfit_mode = 0; % 1 for using age misfit and 0 for using age and Geothermal gradient (misfit_age + misfit_G/n)
error_type = 0;
consider_basal_flux = 1;
%% 
cd (input_data_folder)
files = dir();
dirFlags = [files.isdir];
subFolders = files(dirFlags);
cd ..    
    %% find which ones to process and build the which_one2process matrix
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
    demdata = load(demfile);
    data = data(~isnan(data{:,5}),:);
    data = data(data{:,5}<age0,:);
    nx = load(nxfile);
    ny = load(nyfile);
    N = length(data.lati);
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
%         h_mean = mean(demdata(:,3))/1000;
        h_mean = mean(Elevation); 
    end
    h = Elevation-h_mean;
%     h_data = Elevation-h_mean_data;
    z_mean = 0;
%%
    if (exist('G_present_authority', 'var'))
        G_present = G_present_authority; % if G_present_authority is defined, it will be used in the following calculations
    else
        G_present = which_one2process.G_present(ii);
    end
    G_present_error=G_present*0.1;   % error of present geothermal gradient [K/km]    
    
    if exist ('T0_authority', 'var')
        T0 = T0_authority;
    else
        if h_mean < 0 % using depth
            T0 = 10;
        else
            T0 = 30-lapse_rate*h_mean;  % surface temperature using a temperature - elevation model 
        end
    end
    %% model set up
        n_Age=length(Age);
        zz=zeros(n_Age,1);
    for i=1:n_Age
        if strcmpi(sys(i),'aft')
            zz(i)=110/G_present;
        elseif strcmpi(sys(i),'zft')
            zz(i)=200/G_present;
        elseif strcmpi(sys(i),'ahe')
            zz(i)=50/G_present;
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
        end
    end
    if exist('e0_pr_authority', 'var')
        e0_pr = e0_pr_authority;
    else
        e0_pr = mean(zz./Age);  % assumed initial erosion rate [km/my]
    end
    
%     sigma = e0_pr*0.5;  % error of the assumed initial erosion rate [km/my]
    tmax_age=max(Age);
    tmin_age=min(Age);
    if exist('start_time_authority', 'var')
        tmax = start_time_authority;
    else
        tmax = tmax_age+10;
    end
%     if tmax_age <= 10 
%         tmax = tmax_age+10; % the starting time of calculation, which must be older than the maximum age of the sample
%     else
%         tmax = tmax_age*1.4; % the starting time of calculation, which must be older than the maximum age of the sample
%     end
    if exist('dt2plot_authority', 'var')
        dt2plot = dt2plot_authority;
    else
        if tmax_age < 4  && tmin_age < 1 && length(Age)>5
            dt2plot = .5; %Ma 
            nt4each_step = 5;
        elseif tmax_age>=4 && tmax_age<15 || tmin_age < 2
            dt2plot = 1; %Ma
            nt4each_step = 5;
        elseif tmax_age>=15 && tmax_age<30 || tmin_age < 4
            dt2plot = 2; %Ma
            nt4each_step = 10;
        elseif tmax_age>=30 && tmax_age<50 || tmin_age< 6
            dt2plot = 3; %Ma  
            nt4each_step = 30;
        elseif tmax_age>=50 && tmax_age<80 || tmin_age < 20
            dt2plot = 5; %Ma
            nt4each_step = 20;
        else 
            dt2plot = 10; %Ma  
            nt4each_step = 20;
        end
    end
%     dt2plot = dt2plot*2; % can we used a larger time step to avoid
%     negative rates?
%     dt = dt2plot/nt4each_step;  % step size of age. A large dt may result in large misfit.   
%     dz = zmax*1000/Nz;
%     kappa=kappa*(1000^2.);
%     dt = (dz^2.)/kappa/0.02;

    Nt = fix(tmax/dt2plot);  % number of steps of calculation  


    %% test the effect of prior e and error on the inversion results, it seems the effect is minor

    fprintf('Running inversion with different e0 \n');
    misfit2=zeros(1,4); % misfit of age, misfit of thermal gradient, thermal gradient input and e0 input

        e0 = e0_pr_authority;
        Tgrd = Tgrd_initial_authority;
        if exist('e0_pr_err_authority', 'var')
            sigma = e0_pr_err_authority;
        else
            sigma=e0*0.5;
        end
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
        output_dirname=strcat(outputfolder_prefix,'_masterRUN');
        mkdir(output_dirname);
        cd (output_dirname)
        misfit2(1,1)=misfit_age_postior;
        misfit2(1,4)=e0;
        misfit2(1,2)=misfit_G_postior;
        misfit2(:,3)=Tgrd;
%         clear plottingOrNot;
        save run.mat

    if plottingOrNot == 0
    else
        cd(codeFolder)
        plot_figures_during_batchRUN4reference_model
    end

    fprintf('------------------\n');
end
