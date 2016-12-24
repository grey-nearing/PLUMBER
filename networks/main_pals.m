clear all
close all
clc
restoredefaultpath
addpath('../matlab_tools');

% site names
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

% number of sites
Ns = length(siteNames);

% variable names
varNames = [{'Ws'},{'Ta'},{'Rn'},{'Rh'},{'PP'},{'Qe'},{'Qh'},{'NEE'},{'SM1'},{'SM2'}];
Du = 5; % number of forcing data

% number of histogram bins
Nb = 50;

% temporal lags
lags = 48*[1/48,1/2,1,7,30,120];
Nl = length(lags);

% init storage
D = 10;
T = zeros(Ns,Nl,D,D)-0/0;
H = zeros(Ns,Nl,D,D)-0/0;
S = zeros(Ns,Nl,D,D)-0/0;

% loop through sites
for s = 1:20

 % load data
 fname = strcat('../data/pals_data/extracted/',siteNames(s),'.txt');
 data  = load(strcat(fname{1}));
 data(:,1:3) = [];  % remove dates

 % data dimensions
 [Nt,Dz] = size(data);
 assert(D==Dz);

 % histogram bin sizes
 for d = 1:D
  Bmin = min(data(:,d))-1e-6;
  Bmax = max(data(:,d))+1e-6;
  B(:,d) = linspace(Bmin,Bmax,Nb); 
 end

 % process networks at the different scales
 for l = 1:Nl
  tic;

  for x = 1:D
   if ~isempty(find(data(:,x)<-9990,1,'first')); continue; end;
   Xw = window_average(data(:,x),lags(l));

   for y = Du+1:D
    if x==y; continue; end;
    if ~isempty(find(data(:,y)<-9990,1,'first')); continue; end;

    Yw = window_average(data(:,y),lags(l));
    Xt = Xw(1:end-1,:); 
    Yt = Yw(1:end-1,:); 
    Ys = Yw(2:end,:); 

    [Isxt,Isx,Ist,Ixt,Hs,Hx,Ht] = mutual_info_3(Ys,Xt,Yt,B(:,y),B(:,x),B(:,y));
    T(s,l,x,y) = Isx-Isxt;
    H(s,l,x,y) = Hs;
    S(s,l,x,y) = Isxt;
   end
  end

  % save in Circos format
  TT = squeeze(T(s,l,:,:)); 
  fname = strcat('results/pals/PALS_',siteNames(s),'_',num2str(l),'.txt');
  saveCircos(TT,varNames,Du,D,fname{1});

  % save in Circos format
  TT = squeeze(T(s,l,:,:)./H(s,l,:,:)); 
  fname = strcat('results/pals/PALS_Normalized_',siteNames(s),'_',num2str(l),'.txt');
  saveCircos(TT,varNames,Du,D,fname{1});

  % screen report
  time = toc;
  [s/Ns,l/Nl,time]

 end

 % save progress
 save('results/pals/savePALSprogress.mat');

end




