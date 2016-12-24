clear all
close all
clc

modelNames = [ ...
              {'CABLE.2.0'}
              {'CABLE_2.0_SLI.vxh599_r553'}
              {'CHTESSEL'}
              {'COLASSiB.2.0'}
              {'ISBA_SURFEX_3l.SURFEX7.3'}
              {'ISBA_SURFEX_dif.SURFEX7.3'}
              {'JULES.3.1'}
              {'JULES3.1_altP'}
              {'Mosaic.1'}
              {'NOAH.2.7.1'}
              {'Noah.3.2'}
              {'NOAH.3.3'}
%              {'ORCHIDEE.trunk_r1401'}
%              {'SUMMA.1.0.exp.01.000'}
%              {'SUMMA.1.0.exp.01.001'}
%              {'SUMMA.1.0.exp.01.002'}
%              {'SUMMA.1.0.exp.01.003'}
%              {'SUMMA.1.0.exp.01.004'}
%              {'SUMMA.1.0.exp.01.005'}
%              {'SUMMA.1.0.exp.01.006'}
%              {'SUMMA.1.0.exp.01.007'}
%              {'SUMMA.1.0.exp.01.008'}
%              {'SUMMA.1.0.exp.01.009'}
              {'Manabe_Bucket.2 '}
              {'Penman_Monteith.1'}];

siteNames =  [{'Amplero'}
              {'Blodgett'}
              {'Bugac'}
              {'ElSaler2'}
              {'ElSaler'}
              {'Espirra'}
              {'FortPeck'}
              {'Harvard'}
              {'Hesse'}
              {'Howlandm'}
              {'Howard'}
              {'Hyytiala'}
              {'Kruger'}
              {'Loobos'}
              {'Merbleue'}
              {'Mopane'}
              {'Palang'}
              {'Sylvania'}
              {'Tumba'}
              {'UniMich'}];

% dimensions
Nsites = length(siteNames);
Nmodels = length(modelNames);

% loop through sites
for s = 1:Nsites   

 % get FluxNet data
 fname = strcat('pals_data/extracted/',siteNames{s},'.txt'); 
 odata = load(fname);
 To = size(odata,1);

 % extract dates
 odates = odata(:,1:3);

 % init storage
 model = zeros(To,6,Nmodels)-9999;

 % get model data
 for m = 1:Nmodels 

  % screen report
  fprintf('Loading data at site %s from model %s ... ',siteNames{s},modelNames{m}); tic;

  % get model data
  fname = strcat('model_data/',modelNames{m},'/',modelNames{m},'_',siteNames{s},'Fluxnet.1.4.nc');

  % extract data
  try
   time = squeeze(ncread(fname,'time')); 
   Qe = squeeze(ncread(fname,'Qle'));  
   Qh = squeeze(ncread(fname,'Qh'));  
   SM = squeeze(ncread(fname,'SoilMoist'));  

size(SM)
keyboard

   SM1 = SM(:,1); SM2 = SM(:,2);
  catch
   fprintf('Failed !!!! \n')
   error('Did not work !!!')
%   continue
  end

  try
   NEE = squeeze(ncread(fname,'NEE'));  
  catch
   NEE = zeros(size(Qe))./0;
  end

  % extract dates
  time = squeeze(ncread(fname,'time'));  
  Tm = length(time);
  clear yr dy hr
  for t = 1:Tm
   [yr(t),dy(t),hr(t)] = convert_time(time(t));
  end
  yr = yr'; dy = dy'; hr = hr';

  % find matching dates in observation files
  for t = 1:To
   Id = find(dy         == odates(t,2));
   Ih = find(hr(Id)     == odates(t,3));
   Iy = find(yr(Id(Ih)) == odates(t,1));
   if ~isempty(Iy) 
    if length(Iy>1); Iy = Iy(1); end; % wtf, chtessl?
    model(t,:,m) = [yr(Id(Ih(Iy))),dy(Id(Ih(Iy))),hr(Id(Ih(Iy))),Qe(Id(Ih(Iy))),Qh(Id(Ih(Iy))),NEE(Id(Ih(Iy)))]; 
   end
  end

  % number of points found
  Npts = length(find(model(:,1,m)>0));

  % fix CHTESSEL
  if strcmpi(modelNames{m},'CHTESSEL')
   model(:,4:5,m) = -model(:,4:5,m);
  end

  % screen report
  t = toc; fprintf('finished - # points = %d - time = %f \n',Npts,t);

 end % models

 % file name
 fname = strcat('model_data/extracted/',siteNames{s},'.mat');

 % save
 save(fname,'model'); 

end % sites












