clear
close all

%% setting up the input and output folder
codeFolder2 = pwd;
idcs = strfind(codeFolder2,filesep);
currentFolder = codeFolder2(1:idcs(end)-1);
input_data_folder = strcat(currentFolder,filesep,"data_A2Epaper");
output_folder2 = strcat(currentFolder,filesep,"000_output4code_version_no_rounding2023");
addpath(strcat(codeFolder2,filesep,"othercolor",filesep,"othercolor"))

%% plotting
plottingOrNot = 1;
which_one2process_input = "Fitzgerald1995";
panel_numb = {'(a)','(b)','(c)','(d)','(e)','(f)'};

dirname=strcat(input_data_folder,'/',which_one2process_input);
cd(dirname)

wid = 600;
hei = 350;
%% Modeling using different e0 values and a constant G0
ne_run=6;
f=figure(1);
f.Position(3:4) = [wid hei];
for pp=1:ne_run
    dirname=strcat('RUN_e', num2str(pp));
    cd (output_folder2)
    cd (which_one2process_input)
    cd (dirname)
    load run
    e=flipud(e(:, ite+1));
    e=[e;e(end)];
    de=flipud(sqrt(diag(Cpo)));
    de=[de;de(end)];
    time=0:dt2plot:tmax;
    
    subplot(2,3,pp)
    yyaxis left
    hold on
    stairs(time, e-de,'-','LineWidth',0.5,'Color','b');  
    stairs(time, e+de,'-','LineWidth',0.5,'Color','b'); 
    stairs(time, e,'-','LineWidth',1.5,'Color','b'); 
%     scatter(Age, zeros(length(Age),1))
    set(gca,'xdir','reverse')
    xlim([0 tmax])
    ylim([-.2 1.2])
    box on
    hold off
    if pp ==1 || pp ==4
        ylabel('Erosion rate (km/my)')
    else
    end
    if  pp>3
        xlabel('Time (Ma)')
    else
    end

    text2print1 = strcat(panel_numb(pp), " e0 = ", num2str(round(e0,2)), "\pm", num2str(round(sigma,2))); 
    text2print2 = strcat("\it\phi_{\tau, pr}","\rm = ", num2str(round(misfit_age_prior,2)),", \it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
%     text2print3 = strcat("\it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
    text2print4 = strcat("\it\phi_{\gamma, pr}","\rm = ", num2str(round(misfit_G_prior,2)),", \it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
%     text2print5 = strcat("\it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
%     text(15, 1.0, {text2print1;text2print2;text2print3;text2print4;text2print5},'fontsize', 7.5);
    text(.07, .82, {text2print1;text2print2;text2print4},'fontsize', 7.5,'Units','normalized');
    
    yyaxis right
    hold on
    plot(dTdz_1km_po(:,1), dTdz_1km_po(:,2),'r-','LineWidth',1)
    plot(dTdz_1km_pr(:,1), dTdz_1km_pr(:,2),'r--','LineWidth',1)
    set(gca,'FontSize',8)
    if pp ==3 || pp ==6
        ylabel('Geothermal gradient (^oC/km)')
    else
    end

end
sgtitle(strcat("Modeling using different e0 values and a G0 = ",num2str(Tgrd)),'fontsize', 10)
cd(codeFolder2)
print(gcf,'testing_e0','-dpdf') % then print it

%% Modeling using different sigma0 values and constant e0 and G0 values
ne_run=6;
f=figure(2);
f.Position(3:4) = [wid hei];
for pp=1:ne_run
    dirname=strcat('RUN_sigma', num2str(pp));
    cd (output_folder2)
    cd (which_one2process_input)
    cd (dirname)
    load run
    e=flipud(e(:, ite+1));
    e=[e;e(end)];
    de=flipud(sqrt(diag(Cpo)));
    de=[de;de(end)];
    time=0:dt2plot:tmax;
    
    subplot(2,3,pp)
    yyaxis left
    hold on
    stairs(time, e-de,'-','LineWidth',0.5,'Color','b');  
    stairs(time, e+de,'-','LineWidth',0.5,'Color','b'); 
    stairs(time, e,'-','LineWidth',1.5,'Color','b'); 
%     scatter(Age, zeros(length(Age),1))
    set(gca,'xdir','reverse')
    xlim([0 tmax])
    ylim([-.2 1.2])
    box on
    hold off
    if pp ==1 || pp ==4
        ylabel('Erosion rate (km/my)')
    else
    end
    if  pp>3
        xlabel('Time (Ma)')
    else
    end

    text2print1 = strcat(panel_numb(pp), " e0 = ", num2str(round(e0,2)), "\pm", num2str(round(sigma,2))); 
    text2print2 = strcat("\it\phi_{\tau, pr}","\rm = ", num2str(round(misfit_age_prior,2)),", \it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
%     text2print3 = strcat("\it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
    text2print4 = strcat("\it\phi_{\gamma, pr}","\rm = ", num2str(round(misfit_G_prior,2)),", \it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
%     text2print5 = strcat("\it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
%     text(15, 1.0, {text2print1;text2print2;text2print3;text2print4;text2print5},'fontsize', 7.5);
    text(.07, .82, {text2print1;text2print2;text2print4},'fontsize', 7.5,'Units','normalized');

    yyaxis right
    hold on
    plot(dTdz_1km_po(:,1), dTdz_1km_po(:,2),'r-','LineWidth',1)
    plot(dTdz_1km_pr(:,1), dTdz_1km_pr(:,2),'r--','LineWidth',1)
    set(gca,'FontSize',8)
    if pp ==3 || pp ==6
        ylabel('Geothermal gradient (^oC/km)')
    else
    end

end
sgtitle(strcat("Modeling using different sigma0 values and a e0 = ",num2str(e0)),'fontsize', 10)
cd(codeFolder2)
print(gcf,'testing_sigma0','-dpdf') % then print it

%% Modeling using different G0 values and constant e0
nG_run=6;
f=figure(3);
f.Position(3:4) = [wid hei];
for pp=1:nG_run
    dirname=strcat('RUN_G', num2str(pp));
    cd (output_folder2)
    cd (which_one2process_input)
    cd (dirname)
    load run
    e=flipud(e(:, ite+1));
    e=[e;e(end)];
    de=flipud(sqrt(diag(Cpo)));
    de=[de;de(end)];
    time=0:dt2plot:tmax;
    
    subplot(2,3,pp)
    yyaxis left
    hold on
    stairs(time, e-de,'-','LineWidth',0.5,'Color','b');  
    stairs(time, e+de,'-','LineWidth',0.5,'Color','b'); 
    stairs(time, e,'-','LineWidth',1.5,'Color','b'); 
%     scatter(Age, zeros(length(Age),1))
    set(gca,'xdir','reverse')
    xlim([0 tmax])
    ylim([-.2 1.2])
    box on
    hold off
    if pp ==1 || pp ==4
        ylabel('Erosion rate (km/my)')
    else
    end
    if  pp>3
        xlabel('Time (Ma)')
    else
    end


    text2print1 = strcat(panel_numb(pp), " G0 = ", num2str(round(Tgrd,2))); 
    text2print2 = strcat("\it\phi_{\tau, pr}","\rm = ", num2str(round(misfit_age_prior,2)),", \it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
%     text2print3 = strcat("\it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
    text2print4 = strcat("\it\phi_{\gamma, pr}","\rm = ", num2str(round(misfit_G_prior,2)),", \it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
%     text2print5 = strcat("\it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
%     text(20, 1.0, {text2print1;text2print2;text2print3;text2print4;text2print5},'fontsize', 7.5);
    text(.07, .82, {text2print1;text2print2;text2print4},'fontsize', 7.5,'Units','normalized');

    yyaxis right
    hold on
    plot(dTdz_1km_po(:,1), dTdz_1km_po(:,2),'r-','LineWidth',1)
    plot(dTdz_1km_pr(:,1), dTdz_1km_pr(:,2),'r--','LineWidth',1)

    set(gca,'FontSize',8)
    if pp ==3 || pp ==6
        ylabel('Geothermal gradient (^oC/km)')
    else
    end

end
sgtitle(strcat("Modeling using different G0 values and a e0 = ",num2str(e0)),'fontsize', 10)
cd(codeFolder2)
print(gcf,'testing_G0','-dpdf') % then print it

%% Modeling using different dt values
nG_run=6;
f=figure(4);
f.Position(3:4) = [wid hei];
for pp=1:nG_run
    dirname=strcat('RUN_deltat', num2str(pp));
    cd (output_folder2)
    cd (which_one2process_input)
    cd (dirname)
    load run
    e=flipud(e(:, ite+1));
    e=[e;e(end)];
    de=flipud(sqrt(diag(Cpo)));
    de=[de;de(end)];
    time=0:dt2plot:tmax;
    
    subplot(2,3,pp)
    yyaxis left  
    hold on
    stairs(time, e-de,'-','LineWidth',0.5,'Color','b');  
    stairs(time, e+de,'-','LineWidth',0.5,'Color','b'); 
    stairs(time, e,'-','LineWidth',1.5,'Color','b'); 
%     scatter(Age, zeros(length(Age),1))
    set(gca,'xdir','reverse')
    xlim([0 tmax])
    ylim([-.2 1.2])
    box on
    hold off
    if pp ==1 || pp ==4
        ylabel('Erosion rate (km/my)')
    else
    end
    if  pp>3
        xlabel('Time (Ma)')
    else
    end


    text2print1 = strcat(panel_numb(pp), " \Delta{t} = ", num2str(round(dt2plot,1))); 
    text2print2 = strcat("\it\phi_{\tau, pr}","\rm = ", num2str(round(misfit_age_prior,2)),", \it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
%     text2print3 = strcat("\it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
    text2print4 = strcat("\it\phi_{\gamma, pr}","\rm = ", num2str(round(misfit_G_prior,2)),", \it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
%     text2print5 = strcat("\it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
%     text(20, 1.0, {text2print1;text2print2;text2print3;text2print4;text2print5},'fontsize', 7.5);
    text(.07, .82, {text2print1;text2print2;text2print4},'fontsize', 7.5,'Units','normalized');

    yyaxis right
    hold on
    plot(dTdz_1km_po(:,1), dTdz_1km_po(:,2),'r-','LineWidth',1)
    plot(dTdz_1km_pr(:,1), dTdz_1km_pr(:,2),'r--','LineWidth',1)

    set(gca,'FontSize',8)
    if pp ==3 || pp ==6
        ylabel('Geothermal gradient (^oC/km)')
    else
    end
end
sgtitle(strcat("Modeling using different dt values and a e0 = ",num2str(e0)),'fontsize', 10)
cd(codeFolder2)
print(gcf,'testing_dt','-dpdf') % then print it


