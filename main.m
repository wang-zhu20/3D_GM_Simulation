

clear all;
clc;

% Rrup:     Rupture distance (km)
% Rhyp:     Hypocentral distance (km)
% inVs30:   Average shear wave velocity with 30m in surface
% Mw:       Minimum Moment Magnitude
% nsmpl:    # of samples to generate



offset=0;

%% conditions for each case
ind=1;

% MM=[ 6 7 8];
% RrupM = [ 10 20 50];
% Vs30M = [270 450 760];
% ncase = length(MM)*length(RrupM)*length(Vs30M);
% number = 150;
% nsmpl=ones(1,ncase).*number;

MM=[ 6];
RrupM = [ 10];
Vs30M = [270];
ncase = length(MM)*length(RrupM)*length(Vs30M);
number = 1;
nsmpl=ones(1,ncase).*number;


for i =1:length(MM)
    for j = 1:length(RrupM)
        for k = 1:length(Vs30M)
            N = (i-1)*9+(j-1)*3+k;
            Mw(N)=MM(i);Rrup(N)=RrupM(j);Vs30(N) = Vs30M(k);
        end
    end
end
Rhyp = Rrup;




%%
% initializing random number generater
rs.TotalSign = RandStream('mt19937ar');
rs.MinorRand = RandStream('mt19937ar');
rs.MajorLoc  = RandStream('mt19937ar');
rs.MajorAmp  = RandStream('mt19937ar');
inrand=ceil(rand(1)*100);
for i=1:1:inrand
    rand(rs.TotalSign);
    rand(rs.MinorRand);
    rand(rs.MajorLoc);
    rand(rs.MajorAmp);
end

% loop for each case
for j=1:1:ncase
    % making directory for output files
    filename=sprintf('SimDatabase\\M%03.1f_Rrup%06.2f_Vs%05.1f',Mw(j),Rrup(j),Vs30(j));
    system(sprintf('mkdir %s',filename));
    filename1=sprintf('M%03.1f_Rrup%06.2f_Vs%05.1fH1',Mw(j),Rrup(j),Vs30(j));
    filename2=sprintf('M%03.1f_Rrup%06.2f_Vs%05.1fH2',Mw(j),Rrup(j),Vs30(j));
    filename3=sprintf('M%03.1f_Rrup%06.2f_Vs%05.1fV',Mw(j),Rrup(j),Vs30(j));
    % 初始化计算出的数据
    PGA1=[];PGV1=[];PGD1=[];PSA1=[];
    IA1=[];CAV1=[];
    PGA2=[];PGV2=[];PGD2=[];PSA2=[];
    IA2=[];CAV2=[];
    PGA3=[];PGV3=[];PGD3=[];PSA3=[];
    IA3=[];CAV3=[];
    % loop for each sample
    for i=1:1:nsmpl(j)
        for iflg = 1:3
            disp(['computing... smpl ' num2str(i,'% 6d') '/' num2str(nsmpl(j),'% 6d') '| case' num2str(j,'% 6d') '/' num2str(ncase,'% 6d')]);
            [prmcoef0]=fn_PredictWaveletPara(Mw(j),Rrup(j),Rhyp(j),Vs30(j));
            prmcoef=prmcoef0(iflg);
            [th, dt] = fn_getsim33_wz(Mw(j),Rrup(j),Vs30(j),rs,Rhyp(j),prmcoef);
            if iflg==1
                fn_writeACCinNGA(sprintf('.\\%s\\%s_#%05d.AT2',filename,filename1,(i+offset)),length(th),dt,th, Mw(j), Rrup(j), Rhyp(j), Vs30(j))
                th1=th;
            elseif iflg==2
                fn_writeACCinNGA(sprintf('.\\%s\\%s_#%05d.AT2',filename,filename2,(i+offset)),length(th),dt,th, Mw(j), Rrup(j), Rhyp(j), Vs30(j))
                th2=th;
            else
                fn_writeACCinNGA(sprintf('.\\%s\\%s_#%05d.AT2',filename,filename3,(i+offset)),length(th),dt,th, Mw(j), Rrup(j), Rhyp(j), Vs30(j))
                th3=th;
            end
        end
    end
end



