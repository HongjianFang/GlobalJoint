% Matlab script
% For constructing G matrix for surface waves.
% Written by Hongjian Fang @ MIT, Feb 2 2018. hfang@mit.edu

%% First Gridding
% parameterization, refer to makegrid.m
% single precision is enough
%showwaitbar = 0;
addpath '/scratch1/hfang/Tomo_clean/mksensm/surfdata'
latmin = -90;
latmax = 90;
nlat = 513; % 2^n+1
longmin = 0;
longmax = 360;
nlon = 1025;
depmin = 0;
depmax = 2891.5;
ndep = 65;
latrange = linspace(latmin,latmax,nlat);
nlat = nlat - 1;
dlat = latrange(2)-latrange(1);
longrange = linspace(longmin,longmax,nlon)-180;
nlon = nlon - 1;
dlon = longrange(2)-longrange(1);
deprange = linspace(depmin,depmax,ndep);
ddep = deprange(2)-deprange(1);
ndep = ndep - 1;

% data input (sources and stations distribution & )
lnvs1 = load('../../depthkernel/lnVs.dat'); % fundamental mode
disper = load('../../depthkernel/disper.dat'); % fundamental mode
periods = load('../periods.dat');
%meas = load('test13.dat');
%[nmeas,~] = size(meas);
nmeas = 500000;
nnzero = floor(0.0001*nlat*nlon*ndep*nmeas);
disp(['max number of nnzeros:',num2str(nnzero)]);
S = single(zeros(nnzero,1));
C = single(zeros(nnzero,1));
L = single(zeros(nmeas,1));
N = single(zeros(nmeas,1));
D = single(zeros(nmeas,1));
lin = single(zeros(nmeas,1));
linsize = [0,0];

step = 30000; %minimu step for ray path, unit:feet
disp(['begining ...'])
% if showwaitbar
%     h = waitbar(0,'Please wait...');
% end
nzall = 0;
% ray tracing, loop over for all pairs
Block = 1;
fid = fopen('SL2013_3_2.dat','r');
iter = 0;
while (~feof(fid))                               % For each block:
    fprintf('Block: %s\n', num2str(Block))           % Print block number to the screen
    InputText = textscan(fid,'%s',2,'delimiter','\n');  % Read 2 header lines
    HeaderLines = InputText{1};
    %disp(HeaderLines);                        % Display header lines
    
    InputText = textscan(fid,'%f');     % Read the numeric value
    Data = InputText{1};                          % Specify that this is the
    
    srcinf = str2num(HeaderLines{2});
    nm = srcinf(1);
    Data = reshape(Data,3,nm);
    Data = Data';
    lat_s = srcinf(4);
    lon_s = srcinf(5);
    lat_r = srcinf(6);
    lon_r = srcinf(7);
    
    type = 1;
    iuse = 1;
    
    %if (type == 1)
    [lats, lons]=generate_great_circle_path(lat_s,lon_s,lat_r,lon_r,step);
    idx = find(lons<0);
    lons(idx) = lons(idx) + 360;
    %else
    %    [lats, lons]=generate_great_circle_path_l2(lat_s,lon_s,lat_r,lon_r,step);
    %    idx = find(lons<0);
    %    lons(idx) = lons(idx) + 360;
    %end
    
    idxlat = ceil((lats+90)/dlat);
    idxlon = ceil(lons+eps/dlon);
    
    row = zeros(nlat*nlon,1);
    npts = length(lats);
    for ii = 1:npts-1
        %idx = (idxlon(ii+1)-1)*nlat+idxlat(ii+1);
        idx = (idxlon(ii)-1)*nlat+idxlat(ii);
        row(idx) = row(idx)+distance(lats(ii),lons(ii),lats(ii+1),lons(ii+1));
    end
    row = row/180*pi*6371;
    
    % for checking the results
    dis = distance(lat_s,lon_s,lat_r,lon_r)/180*pi*6371;
    sdis = sum(row);
    if (abs(dis-sdis)>1e0)
        warning('this path is not accurate');
        iuse = 0;
        % go out the loop and continue?
    end
    
    % incorporating depth kernel, read from
    % from fundamental to overtones
    % lat 1 S-->N
    % lon 2 W-->E
    % dep 3 surface-->deep
    roundf = Data(:,3);
    for ii = 1:length(roundf)
        [~,midx] = min(abs(roundf(ii)-periods));
        roundf(ii) = periods(midx);
    end
    [roundf,ia,ic] = unique(roundf,'stable');
    nm = length(roundf);
    data_tmp = zeros(nm,1);
    for ii = 1:nm-1
        data_tmp(ii) = mean(Data(ia(ii):ia(ii+1)-1,2));
    end
    data_tmp(end) = mean(Data(ia(nm):end,2));
    
    for ifreq = 1:nm
        iter = iter + 1;
        %perd = find(periods == round(Data(ifreq,3)));
        perd = find(periods == roundf(ifreq));
        data = data_tmp(ifreq);
        lrow = zeros(ndep*nlat*nlon,1);
        nidx = find(row>1e0);
        for ii = 1:length(nidx)
            lrow(nidx(ii):nlat*nlon:(length(lnvs1(:,1))-2)*nlat*nlon+nidx(ii)) = row(nidx(ii))*lnvs1(1:end-1,perd);
        end
        
        nidx = find(lrow>1e0);
        nz = length(nidx);
        S(nzall+1:nzall+nz) = lrow(nidx);
        C(nzall+1:nzall+nz) = nidx;
        L(iter) = nzall+nz;
        D(iter) = dis/data*1000-dis/disper(perd);
        N(iter) = iter;
        lin(iter) = iuse;
        nzall = nzall + nz;
%         if showwaitbar
%             h = waitbar(iter/nmeas,h,...
%                 ['remaining percent =',num2str((nmeas-iter)/nmeas*100), 'percent']);
%         end
    end
    Block = Block + 1;

end
%%
linsize(1) = iter;
linsize(2) = nzall;
S = S(1:nzall);
C = C(1:nzall);
L = L(1:iter);
N = N(1:iter);
D = D(1:iter);
lin = lin(1:iter);


% write into binary files: C,S,L,N,lin...
addpath '/scratch1/hfang/Tomo_clean/kmbin'
writeb(S,'S.bin','float');
writeb(C,'C.bin','int');
writeb(L,'L.bin','int');
writeb(N,'N.bin','int');
writeb(lin,'lin.bin','int');
writeb(D,'D.bin','float');
writeb(linsize,'linsize.bin','int');
disp('Finished')
