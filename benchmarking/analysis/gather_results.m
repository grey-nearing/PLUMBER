clear all
close all
clc
iFig = 0;
addpath('../../matlab_tools')

% get colors
figure(100);
h = plot(randn(20));
colors = get(h,'Color');
close(100);

% information resolutions
resolutions = [0.01,0.02,0.05,0.10];
resNames = [{'1% info resolution'},{'2% info resolution'},{'5% info resolution'},{'10% info resolution'}];
Nres = length(resolutiont);

% site names
sites = [{'Amplero'},{'Blodgett'},{'Bugac'   },{'ElSaler2'},{'ElSaler' },...
         {'Espirra'},{'FortPeck'},{'Harvard' },{'Hesse'   },{'Howlandm'},...
         {'Howard' },{'Hyytiala'},{'Kruger'  },{'Loobos'  },{'Merbleue'},...
         {'Mopane' },{'Palang'  },{'Sylvania'},{'Tumba'   },{'UniMich' }];
Nsites = length(sites);

% Models
modelNames = [ ...
              {'CABLE 2.0'}
              {'CABLE 2.0 SLI'}
              {'CHTESSEL'}
              {'COLASSiB 2.0'}
              {'ISBA SURFEX 7.3 3l'}
              {'ISBA SURFEX 7.3 dif'}
              {'JULES 3.1'}
              {'JULES 3.1 (altP)'}
              {'Manabe Bucket'}
              {'Mosaic'}
              {'NOAH 2.7.1'}
              {'Noah 3.2'}
              {'NOAH 3.3'}
%          {'ORCHIDEE.trunk_r1401'}
%          {'SUMMA.1.0.exp.01.000'}
%          {'SUMMA.1.0.exp.01.001'}
%          {'SUMMA.1.0.exp.01.002'}
%          {'SUMMA.1.0.exp.01.003'}
%          {'SUMMA.1.0.exp.01.004'}
%          {'SUMMA.1.0.exp.01.005'}
%          {'SUMMA.1.0.exp.01.006'}
%          {'SUMMA.1.0.exp.01.007'}
%          {'SUMMA.1.0.exp.01.008'}
%          {'SUMMA.1.0.exp.01.009'}
          {'Penman Monteith'}];
Nmodels = length(modelNames);

% stuff we need for the file name
Nepochs = 10000;
trnfctn = 'trainscg';

% experiment targets
targNames = [{'Qe'},{'Qh'}];%,{'NEE'}];
Ntargs = length(targNames);

% number of boodstraps
Nboots = 2;

% loop through experiments
for iTarg = 1:Ntargs
 for iBoot = 1:Nboots

  % screen report 
  tic;

  % init storage
  Yloo = []; 
  Yobs = [];
  Ysit = []; 
  Nsite = zeros(Nsites,1);

  for iSite = 1:Nsites

%% --------------------------
   % file name 
   LOOfname = strcat('../results/LOO_all_',targNames{iTarg},'_',sites{iSite},'_',num2str(iBoot),'_',num2str(Nepochs),'_',trnfctn,'.mat');
   SITEfname = strcat('../results/Site_all_',targNames{iTarg},'_',sites{iSite},'_',num2str(iBoot),'_',num2str(Nepochs),'_',trnfctn,'.mat');

   % load file
   try
    load(LOOfname);  LOO  = results; clear results;
    load(SITEfname); SITE = results; clear results;
   catch
    continue
   end 

   % make sure there was no stupid mistake
   assert(max(abs(LOO.test.Yobs(:)-SITE.Yobs(:)))==0);

   % store all data
   Yloo = [Yloo(:);LOO.test.Yhat(:)];;
   Yobs = [Yobs(:);LOO.test.Yobs(:)];
   assert(length(Yloo) == length(Yobs));
   Nsite(iSite) = length(LOO.test.Yobs(:));

   % store all data
   Ysit = [Ysit(:);SITE.Yhat(:)];;
   assert(length(Ysit) == length(Yobs));

  end % iSite

%% -------------------------- 
  % load models
  Ymod = zeros(length(Ysit),Nmodels);
  
  edex = 0;
  for iSite = 1:Nsites

    % file name
    fname = strcat('../../data/lagged_data/models_',sites{iSite},'.mat');
    load(fname);
    if Nsite(iSite) > 0
     sdex = edex+1; edex = sdex+Nsite(iSite)-1;
     Ymod(sdex:edex,:) = squeeze(model(:,3+iTarg,:)); 
    end

  end % iSite

%% --------------------------

  % calculate loo statistics
  if isempty(Yloo)
   INFO(1,iTarg,iBoot,1:Nres) = 0./0;
  else
   for r = 1:Nres; tic;
    Bw = (max(Yobs)-min(Yobs))*resolutions(r);
    Bhat = (min(Yloo)-Bw):Bw:max(Yloo+Bw);
    Bobs = (min(Yobs)-Bw):Bw:max(Yobs+Bw);
    [INFO(1,iTarg,iBoot,r),H(iTarg,iBoot,r)] = info(Yloo,Yobs,Bhat,Bobs,2);
    t = toc; fprintf('Finished loo stats - iBoot = %d/%d - iRes = %d/%d - time = %f\n',iBoot,Nboots,r,Nres,t);
   end
  end

  % calculate site statistics
  if isempty(Ysit)
   INFO(2,iTarg,iBoot,1:Nres) = 0./0;
  else
   for r = 1:Nres; tic;
    Bw = (max(Yobs)-min(Yobs))*resolutions(r);
    Bhat = (min(Ysit)-Bw):Bw:max(Ysit+Bw);
    Bobs = (min(Yobs)-Bw):Bw:max(Yobs+Bw);
    INFO(2,iTarg,iBoot,r) = info(Ysit,Yobs,Bhat,Bobs,2);
    t = toc; fprintf('Finished site stats - iBoot = %d/%d - iRes = %d/%d - time = %f\n',iBoot,Nboots,r,Nres,t);
   end
  end

  % calculate model statistics
  if iBoot == 1
   for iMod = 1:Nmodels
    if isempty(Ysit)
     INFO(2+iMod,iTarg,iBoot,1:Nres) = 0./0;
    else
     for r = 1:Nres; tic
      Bw = (max(Yobs)-min(Yobs))*resolutions(r);
      Bhat = (min(Ymod(:,iMod))-Bw):Bw:max(Ymod(:,iMod)+Bw);
      Bobs = (min(Yobs)-Bw):Bw:max(Yobs+Bw);
      INFO(2+iMod,iTarg,iBoot,r) = info(Ymod(:,iMod),Yobs,Bhat,Bobs,2);
      t = toc; fprintf('Finished model stats - iMod = %d/%d - iRes = %d/%d - time = %f\n',iMod,Nmodels,r,Nres,t);
     end
    end
   end
  end

  % screen report
  t = toc; fprintf('Finished: iBoot = %d - iTarg = %d - time = %f\n',iBoot,iTarg,t);
  % screen report
  t = toc; fprintf('Finished: iBoot = %d - iTarg = %d - time = %f\n',iBoot,iTarg,t);

 end % iBoot
end % iTarg

% save results
save('loo_site_results.mat','INFO');


%% ---------------------------------------


% plot
iFig = iFig+1;
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');

mu  = squeeze(nanmean(INFO,3));
sig = squeeze(nanstd(INFO,[],3));

[numgroups,numbars] = size(squeeze(mu(1,:,:)));
groupwidth = min(0.8,numbars/(numbars+1.5));
xx = [];
for i = 1:numbars
 xx = [xx;(1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars)]; 
end
xx = xx(:);

% plot
h = bar(squeeze(mu(1,:,:))); hold on;
%h = bar(squeeze(mu(1,:,:)),'facecolor',0.6*ones(3,1)); hold on;
for iRes = 1:Nres
 h(iRes).FaceColor = 1-iRes/Nres * ones(3,1);
end
for iMod = 1:Nmodels
 mdat = squeeze(mu(2+iMod,:,:))'; mdat = mdat(:);
 hm(iMod) = plot(xx,mdat,'o-','linewidth',1,'color',colors{iMod});  
end

% labels
set(gca,'xticklabel',targNames);%,'fontsize',16);
ylabel('Info Ratio','fontsize',16);
grid on;
title('Out-of-Sample Benchmark','fontsize',18)

% legend
leg = legend([h,hm],[resNames';modelNames]);
%set(leg,'fontsize',8)


%% ---------------------------------------


iFig = iFig+1;
figure(iFig); close(iFig); figure(iFig);
set(gcf,'color','w');

for iTarg = 1:Ntargs
 subplot(1,Ntargs,iTarg)

 mu  = squeeze(nanmean(INFO(1:2,iTarg,:,:),3));
 sig = squeeze(nanstd(INFO(1:2,iTarg,:,:),[],3));

 loss = nanmean(squeeze((INFO(2,iTarg,:,:)-INFO(1,iTarg,:,:))./INFO(2,iTarg,:,:)),1);
% mu = cat(1,mu,loss);

 [numgroups,numbars] = size(squeeze(mu));
 groupwidth = min(0.8,numbars/(numbars+1.5));
 xx = [];
 for i = 1:numbars
  xx = [xx;(1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars)]; 
 end
 xx = xx';

 % plot
 [ax,h1,h2] = plotyy(xx(:),mu(:),1+xx(2,:),loss,'bar','bar'); hold on;
 h1.FaceColor = 'none';
 h2.FaceColor = 'none';
 h1.EdgeColor = 'none';
 h2.EdgeColor = 'none';
 h3 = bar(cat(1,mu,loss));
 h4 = bar(mu);
 for iRes = 1:Nres
  frac = (iRes-1)/Nres;
  h3(iRes).FaceColor = (1-frac)*[0.850,0.325,0.098];
  frac = (iRes)/Nres;
  h4(iRes).FaceColor = (1-frac) * [1,1,1];
 end
 ax(2).YColor = colors{2};

 if iTarg == 1; leg = legend(h4,resNames,'location','nw'); end;
 title(targNames{iTarg},'fontsize',18);
 set(gca,'xtick',[1,2,3],'xticklabel',[{'LOO'},{'in-Sample'},{'fractional loss'}]);
 ax(1).YLabel.String = 'Info Ratio';
 ax(1).YLabel.FontSize = 18;
 ax(2).YLabel.String = 'Info Loss Fraction';
 ax(2).YLabel.FontSize = 18;
 set(gca, 'xlim',[0.5,3.5]);

end


%% ---------------------------------------

for iTarg = 1:Ntargs
 iFig = iFig+1;
 figure(iFig); close(iFig); figure(iFig);
 set(gcf,'color','w');

 for iMod = 1:Nmodels
  subplot(2,Nmodels/2,iMod)

  Hmu = squeeze(nanmean(H(          iTarg,:,4),2));
  Umu = squeeze(nanmean(INFO(2     ,iTarg,:,4),3));
  Mmu = squeeze(nanmean(INFO(2+iMod,iTarg,:,4),3)); 
  
  data = [Hmu-Umu,Umu-Mmu,Mmu]/Hmu; sum(data)
  labels = [{'Info Mising from Inputs'},{'Info Lost from Model'},{'Info Provided by Model'}];
  explode = [1,1,1];
  p = pie3(data,explode);
%  p(1).FaceColor = colors{1}; 
%  p(3).FaceColor = colors{7}; 
%  p(5).FaceColor = colors{3}; 
%  p(2).BackgroundColor = 0.8*ones(3,1); 
%  p(4).BackgroundColor = 0.8*ones(3,1); 
%  p(6).BackgroundColor = 0.8*ones(3,1); 

  if iMod == 1; legend(labels); end; 

  title(modelNames{iMod},'fontsize',12);
  frac_loss(iTarg,iMod) = (Umu-Mmu)/Mmu; 
  frac_total(iTarg,iMod) = Mmu/Umu; 
  frac_error(iTarg,iMod) = Mmu/(Umu-Mmu); 

 end

 if iTarg == 1; mtit('Latent Heat Flux','fontsize',20); end;
 if iTarg == 2; mtit('Sensible Heat Flux','fontsize',20); end;

end















