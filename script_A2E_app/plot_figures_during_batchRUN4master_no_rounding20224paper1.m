
close all

%% setting up the input and output folder
codeFolder2 = pwd;
addpath(strcat(codeFolder2,filesep,"othercolor",filesep,"othercolor"))

% if ispc
%     addpath("D:/MACandPC/matlab_tools/othercolor/othercolor")
% elseif ismac
%     addpath("/Users/yuntao/Documents/MACandPC/matlab_tools/othercolor/othercolor")
% elseif isunix
%     addpath("/home/thermo/matlab_tools/othercolor/othercolor")    
% end
age_symboles = {'p', 'o', 'd', 's','^','<','>','v','*','+'};
age_types = {'AHe','AFT','ZHe','ZFT','THe','TFT','KAr','BAr','MAr','HAr'};
cd (output_folder)
cd (folder_ii)
% [Min2,Index2] = min(sum(misfit1(:,2),2));
% dirname=strcat('RUN_G', num2str(Index2));

        if misfit_mode == 1
            [M,Index] = min(misfit2(:,1)); 
        elseif misfit_mode == 2
            [M,Index] = min(sum(misfit2(:,1)+misfit2(:,2)/N,2));
        elseif misfit_mode == 0
            [M,Index] = min(sum(misfit2(:,1)+misfit2(:,2)/sqrt(N),2));
        elseif misfit_mode == 3
            [M,Index] = min(sum(misfit2(:,1)+misfit2(:,2),2));
        else
            disp("please defining a misfit mode")
            return
        end
        dirname=strcat('masterRUN_e', num2str(Index));

% close all
% cd (datafolder_ii)
% [Min1,Index1] = min(sum(misfit1,2));
% [Min2,Index2] = min(sum(misfit2,2));
% if Min1<Min2
%     dirname=strcat('RUN_G', num2str(Index1));
% else
%     dirname=strcat('RUN_e', num2str(Index2));
% end
cd (outputfoldername)
cd (dirname)
load run
cd ..

y=flipud(e(:, ite+1));
y=[y;y(end)];
dy=flipud(sqrt(diag(Cpo)));
dy=[dy;dy(end)];
% y = [y; y(end)];
% dy = [dy; dy(end)];
% calculate the binned erosion hisotry
% the linear inversion method requires a small dt to make sure the
% calculation converges. But the data often does not have the resolution.
% Therefore, the results, derived using a small dt, are binned using a
% reasonable time interval.
% 
% [bined_e,bined_err,time_bin] = binning(Nt,dt2plot,dt,y,dy) ;

%%
figure(2)
plot(1:ite,zc(1,1:end-1),'o-',1:ite,zc(end,1:end-1),'o-');
xlabel('iteration')
ylabel('mean zc')
legend('Lowest sample','Highest sample')
% saveas(gcf,'fig2.png')
% 
% figure(3)
% plot(1:ite,Tc(3,:),'o-',1:ite,Tc(11,:),'o-');
% xlabel('iteration')
% ylabel('Tc (C)')
% legend('5.2Ma','9.1Ma')
% 
% figure(4)
% plot(1:M, e(:, 1),'o-',1:M, e(:,ite+1),'o-');

%%
figure(7)
yyaxis left
plot(h*1000,'b-o')
hold on
plot(p(:,end),'r-o')
yyaxis right
plot(zc(:,ite+1)*1000,'--')
legend('h','p','zmi','location', 'best')
% saveas(gcf,'fig7.png')
    
%%    
figure('Renderer', 'painters', 'Position', [200 200 600 240])
for i=1:length(age_symboles)
    subplot(1,2,1)
    hold on
    errorbar(Age(data2use.system==i),Elevation(data2use.system==i),Error(data2use.system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','black','Color','black')
    plot(age_po(data2use.system==i),Elevation(data2use.system==i),age_symboles{i},'MarkerEdgeColor','blue')
    
    subplot(1,2,2)
    hold on
    errorbar(Age(data2use.system==i),age_po(data2use.system==i),Error(data2use.system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','k','Color','black');
%     hh(i) = errorbar(Age(data2use.system==i),age_po(data2use.system==i),Error(data2use.system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','k','Color','black');
end
subplot(1,2,1)
xlabel('Age (Ma)')
ylabel('Elevation (km)')
legend('Obs (black)','Prd (blue)')
% ylim([1.5 6.5])
% xlim([0 20])
hold off
box on
subplot(1,2,2)
box on
xlabel('Observed age (Ma)')
ylabel('Predicted age (Ma)')  
line([0 1.1*max([Age,age_po])], [0 1.1*max([Age,age_po])])
% legend(hh(sort(unique(system))), age_types{sort(unique(system))}, 'Location','Best')
hold off
sgtitle(strcat("The model yilelds a \phi_{po, age} = ",num2str(round(misfit_age_postior,2)), " and a \phi_{po, G} = ",num2str(round(misfit_G_postior,1))), 'fontsize',10)

%%    
figure('Renderer', 'painters', 'Position', [50 50 250 1300])
subplot(5,1,1)
hold on
for i=1:length(age_symboles)
    errorbar(Age(data2use.system==i),Elevation(data2use.system==i),Error(data2use.system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','black','Color','black')
    plot(age_po(data2use.system==i),Elevation(data2use.system==i),age_symboles{i},'MarkerEdgeColor','blue')
end
box on
% ylim([2 6])
xlabel('Age (Ma)')
ylabel('Elevation (km)')
legend('Obs (black)','Prd (blue)','Location','southeast')

subplot(5,1,2)
hold on
for i=1:length(age_symboles)
    errorbar(Age(data2use.system==i),age_po(data2use.system==i),Error(data2use.system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','k','Color','black');
    hh(i) = errorbar(Age(data2use.system==i),age_po(data2use.system==i),Error(data2use.system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','k','Color','black');
end
box on
xlabel('Observed age (Ma)')
ylabel('Predicted age (Ma)')
line([0 tmax], [0 tmax])
legend(hh(sort(unique(system))), age_types{sort(unique(system))}, 'Location','southeast')
% xlim([0 12])
% ylim([0 12])

subplot(5,1,3)
yyaxis left
hold on
stairs(0:dt2plot:tmax, y-dy,'-','LineWidth',0.5,'Color','b');  
stairs(0:dt2plot:tmax, y+dy,'-','LineWidth',0.5,'Color','b'); 
stairs(0:dt2plot:tmax, y,'-','LineWidth',2,'Color','b'); 
scatter(Age, zeros(length(Age),1))
text2print1 = strcat("e0 = ", num2str(round(e0,2)), "\pm", num2str(round(sigma,2)),", G0 = ", num2str(round(Tgrd,2))); 
text2print2 = strcat("\phi_{pr, age}"," = ", num2str(round(misfit_age_prior,2)),", \phi_{po, age}"," = ", num2str(round(misfit_age_postior,2))); 
% text2print3 = strcat("\phi_{po, age}"," = ", num2str(round(misfit_age_postior,2))); 
text2print4 = strcat("\phi_{pr, G}"," = ", num2str(round(misfit_G_prior,2)), ", \phi_{po, G}"," = ", num2str(round(misfit_G_postior,2))); 
% text2print5 = strcat("\phi_{po, G}"," = ", num2str(round(misfit_G_postior,2))); 
% text(.2, .8, {text2print1;text2print2;text2print3;text2print4;text2print5},'fontsize', 8,'Units','normalized');
text(.1, .8, {text2print1;text2print2;text2print4},'fontsize', 8,'Units','normalized');

set(gca,'xdir','reverse')
% xlim([0 12])
box on
xlabel('Time (Ma)')
ylabel('Erosion rate (km/my)')
% title(strcat("Erosion history, transect ", folder_ii),'Interpreter','none')
hold off

yyaxis right
hold on
fill([0 tmax tmax 0],[G_present-G_present_error,G_present-G_present_error, G_present+G_present_error,G_present+G_present_error],'c','EdgeColor','none','FaceAlpha',.2)
plot([0 tmax],[G_present,G_present],'-c', 'LineWidth',1)
plot(dTdz_1km_po(:,1), dTdz_1km_po(:,2),'r-','LineWidth',1)
plot(dTdz_1km_pr(:,1), dTdz_1km_pr(:,2),'r--','LineWidth',1)
ylabel('Geothermal gradient (oC/km)')


R2 = R;
timeXY = [0 tmax]+dt2plot/2;

cm = othercolor('RdBu9', 40);
thinning = 1;
% cm(cm>0.99) = 1;
% colors = [cm(1:thinning:fix(length(cm)/2),:);[1 1 1]; cm(fix(length(cm)/2)-thinning:thinning:end,:)];
% cm = acc_colormap('cmo_diff');
% R2(R2<0.01) = NaN;

subplot(5,1,4)
Rplot = imagesc(fliplr(timeXY),fliplr(timeXY),R2);
%     set(gca,'Ydir','reverse')
set(gca,'Xdir','reverse')
xlabel('Time (Ma)')
ylabel('Time (Ma)')    
% title('Resolution matrix')
% colormap(cm)
colormap( cm(1:thinning:end,:))
caxis([-0.7 0.7])
c = colorbar('location', 'eastoutside','Orientation','horizontal');
c.Position(1) = 0.54;
c.Position(2) = 0.3;
c.Position(3) = 0.27;
c.Position(4) = 0.005;
c.Label.String = "Resolution";
c.Label.FontSize = 12;
c.Label.Position = [0 3.5 0];
xlabel('Time (Ma)')
ylabel('Time (Ma)') 
% title('Correlation matrix')

subplot(5,1,5)
Cpo_scaled = zeros(Nt, Nt);
for i=1:Nt
    for j=1:Nt
        Cpo_scaled(i,j) = Cpo(i,j)/sqrt(Cpo(i,i))/sqrt(Cpo(j,j));
    end
end
imagesc(fliplr(timeXY),fliplr(timeXY),Cpo_scaled);
set(gca, 'xdir', 'reverse')
colormap( cm(1:thinning:end,:))
caxis([-1 1])
c = colorbar('location', 'eastoutside','Orientation','horizontal');
c.Position(1) = 0.54;
c.Position(2) = 0.127;
c.Position(3) = 0.27;
c.Position(4) = 0.005;
c.Label.String = "Correlation";
c.Label.FontSize = 12;
c.Label.Position = [0 3.5 0];
xlabel('Time (Ma)')
ylabel('Time (Ma)') 
% title('Correlation matrix')


