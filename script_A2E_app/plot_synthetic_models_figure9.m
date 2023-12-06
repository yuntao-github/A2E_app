clear
close all
%% setting up the input and output folder
codeFolder2 = pwd;
idcs = strfind(codeFolder2,filesep);
currentFolder = codeFolder2(1:idcs(end)-1);
output_folder2 = strcat(currentFolder,filesep,"outputs_synthetic_models",filesep,"with_boreholesamples");
addpath(strcat(codeFolder2,filesep,"othercolor",filesep,"othercolor"))

% if ispc
%     addpath("D:\MACandPC\matlab_tools\othercolor\othercolor")
%     output_folder2 = 'D:\MACandPC\software\modeling_vertical_profiles\000_output4code_version_no_rounding2023';
% elseif ismac
%     addpath("/Users/yuntao/Documents/MACandPC/matlab_tools/othercolor/othercolor")
%     output_folder2 = '/Users/yuntao/Documents/MACandPC/software/modeling_vertical_profiles/000_output4code_version_no_rounding20231121';
% end

outputfoldername_all = ["error_type0";"error_type0";"error_type0";"error_type0"];
age_symboles = {'d', 'p', 's', 'o','^','<','>','v','*','+'};
age_types = {'AHe','AFT','ZHe','ZFT','THe','TFT','KAr','BAr','MAr','HAr'};
the3transects=["synthetic_model_a","synthetic_model_b","synthetic_model_c","synthetic_model_d"];
the3transectsNames=["Synthetic dataset-a","Synthetic dataset-b","Synthetic dataset-c","Synthetic dataset-d"];
synthetic_history = [[1 0.1];[1 0.3]; [0.3 1]; [0.1 1]];
synthetic_history = [synthetic_history,synthetic_history(:,end)];
synthetic_time = [20 5 0];
the3transectsNames=["Synthetic dataset-a","Synthetic dataset-b","Synthetic dataset-c","Synthetic dataset-d"];
ylim_value = [1.5, 1.5, 1.5, 1.5];
abc=["a", "b", "c","d"];
figure('Renderer', 'painters', 'Position', [50 50 750 1100])
n_cols_subplot = 4;
for iii=1:length(the3transects)
    cd (output_folder2)
    cd (the3transects(iii))
    cd (outputfoldername_all(iii))
    lastrun = strcat('masterRUN_e',num2str(20));
    cd (lastrun)
    load run
    misfit_mode = 0;
    % [Min2,Index2] = min(sum(misfit1(:,2),2));
    % dirname=strcat('RUN_G', num2str(Index2));
        if misfit_mode == 1
            [M,Index] = min(misfit2(:,1)); 
        elseif misfit_mode == 2
            [M,Index] = min(sum(misfit2(:,1)+misfit2(:,2)/n,2));
        elseif misfit_mode == 0
            [M,Index] = min(sum(misfit2(:,1)+misfit2(:,2)/sqrt(n),2));
        elseif misfit_mode == 3
            [M,Index] = min(sum(misfit2(:,1)+misfit2(:,2),2));
        else
            disp("please defining a misfit mode")
            return
        end
        dirname=strcat('masterRUN_e', num2str(Index));
    cd ..
    % close all
    % cd (datafolder_ii)
    % [Min1,Index1] = min(sum(misfit1,2));
    % [Min2,Index2] = min(sum(misfit2,2));
    % if Min1<Min2
    %     dirname=strcat('RUN_G', num2str(Index1));
    % else
    %     dirname=strcat('RUN_e', num2str(Index2));
    % end
    cd (dirname)
    load run
    cd ..
    y=flipud(e(:, ite+1));
    y=[y;y(end)];
    dy=flipud(sqrt(diag(Cpo)));
    dy=[dy;dy(end)];
    
    %%
    subplot(5,n_cols_subplot,iii)
    hold on
    for i=1:length(age_symboles)
        errorbar(Age(system==i),Elevation(system==i),Error(system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','black','MarkerFaceColor','white','Color','black')
        plot(age_po(system==i),Elevation(system==i),age_symboles{i},'MarkerEdgeColor','blue','MarkerFaceColor','white')
    end
    box on
    if iii==3
        ylim([(floor(min(Elevation)/.500)-.25)*.5 (ceil(max(Elevation)/.500)+0.25)*.500])
    else
        ylim([(floor(min(Elevation)/.500)-.75)*.5 (ceil(max(Elevation)/.500)+.75)*.500])
    end
    xlim([0 tmax])
    xlabel('Age (Ma)')
    
    if iii==1
        legend('Obs (black)','Prd (blue)','Location','southeast')
        ylabel('Elevation (km)')
    else
    end
    set(gca,'FontSize',8)
    title(strcat("(",abc(iii),") ", the3transectsNames(iii)),'FontSize',10,'Interpreter','none');
    
    subplot(5,n_cols_subplot,iii+n_cols_subplot)
    hold on
    for i=1:length(age_symboles)
        hh(i) = errorbar(Age(system==i),age_po(system==i),Error(system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','k','MarkerFaceColor','white','Color','black');
    end
    box on
    xlabel('Observed age (Ma)')
    if iii ==1
        ylabel('Predicted age (Ma)')
    else
    end
    line([0 tmax], [0 tmax])
    legend(hh(sort(unique(system))), age_types{sort(unique(system))}, 'Location','southeast')
    xlim([0 tmax])
    ylim([0 tmax])
    set(gca,'FontSize',8)
    
    subplot(5,n_cols_subplot,iii+2*n_cols_subplot)
    yyaxis right
    hold on
    fill([0 tmax tmax 0],[G_present-G_present_error,G_present-G_present_error, G_present+G_present_error,G_present+G_present_error],'c','EdgeColor','none','FaceAlpha',.2)
    plot([0 tmax],[G_present,G_present],'-c', 'LineWidth',1.5)
    plot(dTdz_1km_po(:,1), dTdz_1km_po(:,2),'r-','LineWidth',1)
    plot(dTdz_1km_pr(:,1), dTdz_1km_pr(:,2),'r--','LineWidth',1)
    set(gca,'FontSize',8)
    if iii == n_cols_subplot
        ylabel('Geothermal gradient (^oC/km)')
    else
    end
    % title(strcat("Erosion history, transect ", folder_ii),'Interpreter','none')
    hold off
    
    yyaxis left
    hold on
    stairs(synthetic_time, synthetic_history(iii,:),'-','LineWidth',1.5,'Color','k'); 
    stairs(0:dt2plot:tmax, y-dy,'-','LineWidth',0.5,'Color','b');  
    stairs(0:dt2plot:tmax, y+dy,'-','LineWidth',0.5,'Color','b'); 
    stairs(0:dt2plot:tmax, y,'-','LineWidth',1.5,'Color','b'); 
%     scatter(Age, zeros(length(Age),1))
    text2print1 = strcat("e0 = ", num2str(round(e0,2)), "\pm", num2str(round(sigma,2)),", G0 = ", num2str(round(Tgrd,2))); 
    text2print2 = strcat("\it\phi_{\tau, pr}","\rm = ", num2str(round(misfit_age_prior,2)),", \it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
    % text2print3 = strcat("\phi_{po, age}"," = ", num2str(round(misfit_age_postior,2))); 
    text2print4 = strcat("\it\phi_{\gamma, pr}","\rm = ", num2str(round(misfit_G_prior,2)), ", \it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
    % text2print5 = strcat("\phi_{po, G}"," = ", num2str(round(misfit_G_postior,2))); 
    % text(.2, .8, {text2print1;text2print2;text2print3;text2print4;text2print5},'fontsize', 8,'Units','normalized');
    text(.16, .83, {text2print1;text2print2;text2print4},'fontsize', 7,'Units','normalized','Interpreter','tex');
    
%     set(gca,'xdir','reverse')
    xlim([0 tmax])
    ylim([0 ylim_value(iii)])
    box on
    xlabel('Time (Ma)')
    if iii ==1
        ylabel('Erosion rate (km/my)')
    else
    end


    
    R2 = R;
    timeXY = [0 tmax]-dt2plot/2;
    
    cm = othercolor('RdBu9', 40);
    thinning = 1;
    % cm(cm>0.99) = 1;
    % colors = [cm(1:thinning:fix(length(cm)/2),:);[1 1 1]; cm(fix(length(cm)/2)-thinning:thinning:end,:)];
    % cm = acc_colormap('cmo_diff');
    % R2(R2<0.01) = NaN;
    
    subplot(5,n_cols_subplot,iii+3*n_cols_subplot)
    Rplot = imagesc(fliplr(timeXY),fliplr(timeXY),R2);
    %     set(gca,'Ydir','reverse')
%     set(gca,'Xdir','reverse')
    xlabel('Time (Ma)')
    if iii == 1
        ylabel('Time (Ma)')   
    else
    end
     
    % title('Resolution matrix')
    % colormap(cm)
    colormap( cm(1:thinning:end,:))
    caxis([-0.75 0.75])
    if iii==1
        c = colorbar('location', 'eastoutside','Orientation','horizontal');
        c.Position(1) = 0.19;
        c.Position(2) = 0.3;
        c.Position(3) = 0.08;
        c.Position(4) = 0.005;
        c.Label.String = "Resolution";
        c.Label.FontSize = 10;
        c.Label.Position = [0 4.2 0];
    else
    end
    xlabel('Time (Ma)')
    if iii == 1
        ylabel('Time (Ma)')   
    else
    end
    set(gca,'FontSize',8)
    % title('Correlation matrix')
    
    subplot(5,n_cols_subplot,iii+4*n_cols_subplot)
    Cpo_scaled = zeros(Nt, Nt);
    for i=1:Nt
        for j=1:Nt
            Cpo_scaled(i,j) = Cpo(i,j)/sqrt(Cpo(i,i))/sqrt(Cpo(j,j));
        end
    end
    imagesc(fliplr(timeXY),fliplr(timeXY),Cpo_scaled);
%     set(gca, 'xdir', 'reverse')
    colormap( cm(1:thinning:end,:))
    caxis([-1 1])
    if iii==1
        c = colorbar('location', 'eastoutside','Orientation','horizontal');
        c.Position(1) = 0.15;
        c.Position(2) = 0.127;
        c.Position(3) = 0.08;
        c.Position(4) = 0.005;
        c.Label.String = "Correlation";
        c.Label.FontSize = 10;
        c.Label.Position = [0 4.2 0];
    else
    end
    xlabel('Time (Ma)')
    ylabel('Time (Ma)') 
    % title('Correlation matrix')
    set(gca,'FontSize',8)

end
cd(codeFolder2)
% set(gcf,'PaperSize',[25 20]); %set the paper size to what you want  
set(gcf,'Position', [50 50 730 810]);
% saveas(gcf, 'plots_figures_4the3transects.pdf')
% exportgraphics(gcf,'plots_figures_synthetic_transects.eps') % then print it
print(gcf,'plots_figures_synthetic_transects','-dpdf') % then print it



