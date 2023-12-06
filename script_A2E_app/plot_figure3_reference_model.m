clear
close all
%% setting up the input and output folder
codeFolder2 = pwd;
idcs = strfind(codeFolder2,filesep);
currentFolder = codeFolder2(1:idcs(end)-1);
output_folder = strcat(currentFolder,filesep,"000_output4code_version_no_rounding2023");
addpath(strcat(codeFolder2,filesep,"othercolor",filesep,"othercolor"))

% 
% codeFolder2 = pwd;
% if ispc
%     output_folder = 'D:\MACandPC\software\modeling_vertical_profiles\000_output4code_version_no_rounding2023';
%     addpath("D:\MACandPC\matlab_tools\othercolor\othercolor")
% elseif ismac
%     output_folder = '/Users/yuntao/Documents/MACandPC/software/modeling_vertical_profiles/000_output4code_version_no_rounding2023';
%     addpath("/Users/yuntao/Documents/MACandPC/matlab_tools/othercolor/othercolor")
% end

age_symboles = {'d', 'p', 's', 'o','^','<','>','v','*','+'};
age_types = {'AHe','AFT','ZHe','ZFT','THe','TFT','KAr','BAr','MAr','HAr'};

% [Min2,Index2] = min(sum(misfit2(:,1),2));
% dirname=strcat('RUN_e', num2str(Index2));

% close all
% cd (datafolder_ii)
% [Min1,Index1] = min(sum(misfit1,2));
% [Min2,Index2] = min(sum(misfit2,2));
% if Min1<Min2
%     dirname=strcat('RUN_G', num2str(Index1));
% else
%     dirname=strcat('RUN_e', num2str(Index2));
% end
cd (output_folder)
cd ("Fitzgerald1995")
% cd (which_one2process_input)
cd ('RefModel_masterRUN')
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
yyaxis left
plot(1:ite,zc(3,1:end-1),'o-',1:ite,zc(10,1:end-1),'o--');
xlabel('iteration')
ylabel('mean zc')
legend('Lowest sample','Highest sample')
ylim([3.5 4.5])
yyaxis right
plot(1:ite,Tc(3,1:end),'o-',1:ite,Tc(10,1:end),'o--');
xlabel('iteration')
ylabel('mean Tc')
legend('Sample of 5.2 Ma','Sample of 9.2 Ma')
ylim([105 125])
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
figure('Renderer', 'painters', 'Position', [50 50 350 300])
subplot(2,2,1)
hold on

for i=2:2
    errorbar(Age(system==i),Elevation(system==i),Error(system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','black','Color','black','MarkerFaceColor','w','MarkerSize',9);
    plot(age_po(system==i),Elevation(system==i),age_symboles{i},'MarkerEdgeColor','blue','MarkerSize',9);
end
box on
ylim([floor(min(Elevation)/0.5-1)*0.5 ceil(max(Elevation)/0.5+1)*0.5])
xlim([0 floor(max(Age)/5+1)*5])
xlabel('Age (Ma)')
ylabel('Elevation (km)')
legend({'Observation','Prediction'},'Location','southeast')
% axes('Position',[.6 .82 .3 .06])
% hold on
% box on
% for i=1:length(age_symboles)
%     hh(i) = errorbar(Age(system==i),age_po(system==i),Error(system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','k','Color','black');
% end
% xlabel('Observation (Ma)')
% ylabel('Prediction (Ma)')
% line([3 18], [3 18])
% % legend(hh(sort(unique(system))), age_types{sort(unique(system))}, 'Location','southeast')
% xlim([3 18])
% ylim([3 18])

% subplot(5,1,2)
% hold on
% for i=1:length(age_symboles)
%     hh(i) = errorbar(Age(system==i),age_po(system==i),Error(system==i),age_symboles{i},'horizontal','CapSize',0,'MarkerEdgeColor','k','Color','black');
% end
% box on
% xlabel('Observed age (Ma)')
% ylabel('Predicted age (Ma)')
% line([3 18], [3 18])
% legend(hh(sort(unique(system))), age_types{sort(unique(system))}, 'Location','southeast')
% xlim([3 18])
% ylim([3 18])

subplot(2,2,2)
yyaxis left
hold on
stairs(0:dt2plot:tmax, y-dy,'-','LineWidth',0.5,'Color','b');  
stairs(0:dt2plot:tmax, y+dy,'-','LineWidth',0.5,'Color','b'); 
stairs(0:dt2plot:tmax, y,'-','LineWidth',2,'Color','b'); 
scatter(Age, zeros(length(Age),1))
set(gca,'xdir','reverse')
xlim([0 tmax])
box on
text2print2 = strcat("\it\phi_{\tau, pr}","\rm = ", num2str(round(misfit_age_prior,2)),", \it\phi_{\tau, po}","\rm = ", num2str(round(misfit_age_postior,2))); 
text2print4 = strcat("\it\phi_{\gamma, pr}","\rm = ", num2str(round(misfit_G_prior,2)),", \it\phi_{\gamma, po}","\rm = ", num2str(round(misfit_G_postior,2))); 
text(23, .68, {text2print2;text2print4},'fontsize', 9);

xlabel('Time (Ma)')
ylabel('Erosion rate (km/my)')
% title(strcat("Erosion history, transect ", folder_ii),'Interpreter','none')
hold off

yyaxis right
hold on
plot(dTdz_1km_po(:,1), dTdz_1km_po(:,2),'r-','LineWidth',1)
plot(dTdz_1km_pr(:,1), dTdz_1km_pr(:,2),'r--','LineWidth',1)
ylabel('Geothermal gradient (^oC/km)')


R2 = R;
timeXY = [0 tmax]+dt2plot/2;
% cm = othercolor('PiYG11', 40);
% cm = flipud(cm);
cm = othercolor('RdBu7', 40);
thinning = 1;
% cm(cm>0.99) = 1;
% colors = [cm(1:thinning:fix(length(cm)/2),:);[1 1 1]; cm(fix(length(cm)/2)-thinning:thinning:end,:)];
% cm = acc_colormap('cmo_diff');
% R2(R2<0.01) = NaN;

subplot(2,2,3)
Rplot = imagesc(fliplr(timeXY),fliplr(timeXY),R2);
%     set(gca,'Ydir','reverse')
set(gca,'Xdir','reverse')
xlim([0 tmax])
ylim([0 tmax])
xlabel('Time (Ma)')
ylabel('Time (Ma)')    
% title('Resolution matrix')
% colormap(cm)
colormap( cm(1:thinning:end,:))
caxis([-0.7 0.7])
c = colorbar('location', 'eastoutside','Orientation','horizontal');
c.Position(1) = 0.3;
c.Position(2) = 0.15;
c.Position(3) = 0.15;
c.Position(4) = 0.015;
c.Label.String = "Resolution";
c.Label.FontSize = 10;
c.Label.Position = [0 3.2 0];
%     saveas(gcf,'fig1.png')

subplot(2,2,4)
%%
Cpo_scaled = zeros(Nt, Nt);
for i=1:Nt
    for j=1:Nt
        Cpo_scaled(i,j) = Cpo(i,j)/sqrt(Cpo(i,i))/sqrt(Cpo(j,j));
    end
end
% timeXY = [0 tmax-dt2plot]+dt2plot/2;
imagesc(fliplr(timeXY),fliplr(timeXY),Cpo_scaled);
set(gca, 'xdir', 'reverse')
xlim([0 tmax])
ylim([0 tmax])
colormap( cm(1:thinning:end,:))
caxis([-1 1])
c = colorbar('location', 'eastoutside','Orientation','horizontal');
c.Position(1) = 0.74;
c.Position(2) = 0.15;
c.Position(3) = 0.15;
c.Position(4) = 0.015;
c.Label.String = "Correlation";
c.Label.FontSize = 10;
c.Label.Position = [0 3.2 0];
xlabel('Time (Ma)')
ylabel('Time (Ma)') 
cd(codeFolder2)
print(gcf,strcat(which_one2process_input,'ref_model_output'),'-dpdf') % then print it
