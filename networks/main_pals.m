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
lags = round(logspace(log10(1),log10(48*180),15));
%lags = 48*[1/48,1/2,1,7,30,120];
Nl = length(lags);

% init storage
D = 10;
T = zeros(Ns,Nl,D,D)-0/0;
H = zeros(Ns,Nl,D,D)-0/0;
S = zeros(Ns,Nl,D,D)-0/0;

% loop through sites
for s = 1:Ns

 % screen report
 fprintf('Working on site %d of %d (%s) ...',s,Ns,siteNames{s}); tsite = tic;

 % load data
 fname = strcat('../data/pals_data/extracted/',siteNames(s),'.txt');
 data  = load(strcat(fname{1}));
 dates = data(:,1:3);
 data(:,1:3) = [];  % remove dates

 % find start and end dates of actual data
 If = find(~isnan(dates(:,1)),1,'first');
 Il = find(~isnan(dates(:,1)),1,'last');
 dates = dates(If:Il,:); data = data(If:Il,:);

 % find discontinuous missing data 
 imiss = isempty(find(isnan(dates(:)),1,'first'));
 count = 0;
 while ~imiss
  If = find(isnan(dates(:,1)),1,'first');
  if If>0.5*length(dates(:,1));
   dates(If:end,:) = []; data(If:end,:) = []; 
  else
   dates(1:If,:) = []; data(1:If,:) = []; 
  end

  count = count+1;
  if count > 1; error('strange patterns of missing data: (%d,%d)',If,Il); end;
  imiss = isempty(find(isnan(dates(:)),1,'first'));
 end

 % data dimensions
 [Nt,Dz] = size(data);
 assert(D==Dz);

 % screen report
 if count == 0 
  fprintf('. Ndata = %d -- Any Missing? No. ',Nt);
 else
  fprintf('. Ndata = %d -- Any Missing? Yes!! ',Nt);
 end

 % histogram bin sizes
 for d = 1:D
  Bmin = min(data(:,d))-1e-6;
  Bmax = max(data(:,d))+1e-6;
  B(:,d) = linspace(Bmin,Bmax,Nb); 
 end

 % process networks at the different scales
 for l = 1:Nl

  for x = 1:D
%   if ~isempty(find(isnan(data(:,x)),1,'first')); continue; end;
   XX = grandma_smoothing(data(:,x),lags(l));
   if isempty(XX); continue; end;
   Xw = window_average(XX,lags(l));

   for y = Du+1:D
    if x==y; continue; end;
%    if ~isempty(find(isnan(data(:,y)),1,'first')); continue; end;
    YY = grandma_smoothing(data(:,y),lags(l));
    assert(length(XX) == length(YY));
    Yw = window_average(YY,lags(l));

    % time shift
    Xt = Xw(1:end-1,:); 
    Yt = Yw(1:end-1,:); 
    Ys = Yw(2:end,:); 

    % make sure we have dealt with missing values
    assert(isempty(find(isnan(Xt),1,'first')));
    assert(isempty(find(isnan(Yt),1,'first')));
    assert(isempty(find(isnan(Ys),1,'first')));

    % do the actaul calculations
    [T(s,l,x,y),H(s,l,x,y),S(s,l,x,y)] = transfer_entropy_window_average(Ys,Xt,Yt,B(:,x),B(:,y));

    % screen report
    fprintf('\n    Site: %s - Lag = %d - %s -> %s = %f - Ndata = %d',siteNames{s},lags(l),varNames{x},varNames{y},T(s,l,x,y)./H(s,l,x,y),length(YY));

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

 end

 % save progress
 save('results/pals/savePALSprogress.mat');

 % screen report
 t = toc; fprintf('\n');
 fprintf('Finished Site %s - time = %f\n',siteNames{s},t);

end




