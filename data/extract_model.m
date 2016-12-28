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

% fixed dimensions
Dd = 3;
Du = 5;
Dz = 5;
D = Dd+Du+Dz;

% loop through sites
for s = 1:Nsites   

 % get FluxNet data
 fname = strcat('pals_data/extracted/',siteNames{s},'.txt'); 
 pals = load(fname);
 Npals = size(pals,1);

 % init storage
 model = zeros([size(pals),Nmodels])./0;

 % get model data
 for m = 1:Nmodels 

  % screen report
  fprintf('%s - %s ... ',siteNames{s},modelNames{m}); tic;

  % get model data
  fname = strcat('model_data/',modelNames{m},'/',modelNames{m},'_',siteNames{s},'Fluxnet.1.4.nc');

  % extract data
  try
   time = squeeze(ncread(fname,'time')); 
   Qe = squeeze(ncread(fname,'Qle'));  
   Qh = squeeze(ncread(fname,'Qh'));  
  catch
   error('Failed to find Qe and Qh variables !!!')
  end

  try
   SM = squeeze(ncread(fname,'SoilMoist'));  
   if(size(SM,2)>size(SM,1)); SM = SM'; end;
   assert(size(SM,1)==length(Qe));
   SM1 = SM(:,1); SM2 = SM(:,2);
  catch
   SM1 = zeros(size(Qe))./0;
   SM2 = zeros(size(Qe))./0;
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
  Npts = 0;
  for t = 1:size(pals,1)
   Id = find(dy         == pals(t,2));
   Ih = find(hr(Id)     == pals(t,3));
   Iy = find(yr(Id(Ih)) == pals(t,1));
   if ~isempty(Iy) 
    if length(Iy>1); Iy = Iy(1); end; % wtf, chtessel?
    pals(t,Dd+Du+1) = Qe(Id(Ih(Iy)));
    pals(t,Dd+Du+2) = Qh(Id(Ih(Iy)));
    pals(t,Dd+Du+3) = NEE(Id(Ih(Iy)));
    pals(t,Dd+Du+4) = SM1(Id(Ih(Iy)));
    pals(t,Dd+Du+5) = SM2(Id(Ih(Iy)));
    Npts = Npts + 1;
   else
    pals(t,Dd+Du+1) = 0/0;
    pals(t,Dd+Du+2) = 0/0;
    pals(t,Dd+Du+3) = 0/0;
    pals(t,Dd+Du+4) = 0/0;
    pals(t,Dd+Du+5) = 0/0;
   end
   model(t,:,m) = pals(t,:);
  end

  % fix CHTESSEL
  if strcmpi(modelNames{m},'CHTESSEL')
   model(:,Dd+Du+(1:2),m) = -model(:,Dd+Du+(1:2),m);
  end

  % screen report
  t = toc; fprintf('Npoints = %d - time = %f \n',Npts,t);

 end % models

 % file name
 fname = strcat('model_data/extracted/',siteNames{s},'.mat');

 % save
 save(fname,'model'); 

end % sites












