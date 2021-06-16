clear all; close all; clc;
load lakegpp.mat; % load data -- also available in .xlsx form 
lp = log10(p); % gross primary production, gC/d
la = log10(a); % area, m2
lv = log10(v); % volume, m3

%% estimate light and temperature by best-fit Fourier series and correct area accordingly

I =  293 + 131.2*cos(y*0.03251) -4.08*sin(y*0.03251);
Io = 293 + 131.2*cos(58.5*0.03251) -4.08*sin(58.5*0.03251);
T = -137.2 + 209.1*cos(y*0.01607) + 1.002*sin(y*0.01607)-44.86*cos(2*y*0.01607)-0.2613*sin(2*y*0.01607);
atic = a.*exp(-.32./(8.6e-5.*(273+T)))./exp(-.32./(8.6e-5.*(273))).*I./Io;
latic = log10(atic);

%% bootstrap to get scaling estimates, uncertainties, and probabilities

for j = 1:10000;
    i = randi(73,1,73);
    vb = lv(i);
    ab = la(i);
    pb = lp(i);
    tb = latic(i);
    [rv(j),svb(j),~] = regression(vb',pb');
    [ra(j),sab(j),~] = regression(ab',pb');
    [rav(j),savb(j),~] = regression(ab',vb');
    [rti(j),stib(j),~] = regression(tb',pb');
end

ae = mean(svb) %.744
au = std(svb) %.048
be = mean(sab) %1.109
bu = std(sab) %.048
ce = mean(savb) %1.329
cu = std(savb) %.073
de = mean(stib) %1.002
du = std(stib) %.037

p34vs23 = sum(svb>2/3)./length(svb) %94.9%
psuper = sum(sab>1)./length(sab) %98.4%

%% get r^2's and RMSEs for scaling relationships

mdl = fitlm(lv,lp);
rmsea = mdl.RMSE % 0.635
r2a = mdl.Rsquared.Ordinary % 0.726
mdl = fitlm(la,lp);
rmseb = mdl.RMSE % 0.504
r2b = mdl.Rsquared.Ordinary % 0.828
mdl = fitlm(la,lv);
rmsec = mdl.RMSE % 0.423
r2c = mdl.Rsquared.Ordinary % 0.908
mdl = fitlm(latic,lp);
rmsed = mdl.RMSE % 0.435
r2d = mdl.Rsquared.Ordinary % 0.872

%% check for correlation of other variables with residuals

[~,m,b] = regression(la',lp');
resid1 = lp-(m.*la+b); 
[~,m,b] = regression(latic',lp');
resid2 = lp-(m.*latic+b); 
clear m b;
[~,presid] = corr([resid1 resid2 c k n y],'type','Kendall','rows','complete');
presid = presid(1:2,3:end);

%% multivariate regression with volume and other variables

mdl = stepwiselm([log10(v) c log10(c) k log10(k) n log10(n) y log10(y)],log10(p));
mdl.Rsquared.Ordinary % 0.804

%% global productivity

clear all; close all; clc;
load hydrolakes.mat; % area (a) and latitude (y) from HydroLAKES saved as .mat file
a = 1e6.*a; % km2-->m2
beta = 0.24;
beta_u = 0.03;

I =  293 + 131.2*cos(y*0.03251) -4.08*sin(y*0.03251);
Io = 293 + 131.2*cos(58.5*0.03251) -4.08*sin(58.5*0.03251);
T = -137.2 + 209.1*cos(y*0.01607) + 1.002*sin(y*0.01607)-44.86*cos(2*y*0.01607)-0.2613*sin(2*y*0.01607);
D = 365.25.*exp(-((T-27)/34.52).^2); % ice free days per year, as a Gaussian function of T
atic = a.*exp(-.32./(8.6e-5.*(273+T)))./exp(-.32./(8.6e-5.*(273))).*I./Io;
global_P = nansum(beta.*atic.*D)./1e12 % global annual lake GPP in Tg:  523
global_P_unc = nansum(beta_u.*atic.*D)./1e12 %65
GCm2y = 1e12.*global_P./nansum(a) %179
GCm2y_unc = 1e12.*global_P_unc./nansum(a) %22
