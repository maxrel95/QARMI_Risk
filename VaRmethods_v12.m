%==========================================================================
% Quantitative Asset and Risk Management
% Risk Management Project : Comparison of VaR estimation techniques
%
% Maxime Borel, Joe Habib, Yoann Joye & Axel Bandelier
% Date: 29 May 2020
%==========================================================================
%the last section of the code "graph" are the same graphs as among the code
%it is juste for grouping them for the report. indivudal graph are still
%available in the result folder. 
%please before running the code put the Kevin sheppard toolbok in the
%same folder with the name as in line 18. 
clear; 
clc; 
close all;

restoredefaultpath      % restore search path to default
addpath(genpath('Kevin Sheppard Toolbox'))
clear RESTOREDEFAULTPATH_EXECUTED

data=importdata('classeur5.xlsx');
date.all=datetime(data.data(:,1),'ConvertFrom','excel');
date.all=date.all(2:end,1);
price=data.data(:,6);
rt=diff(log(price)); %compute log return
T=size(rt,1);
    
evalwindow=T-3343+1:T; %defin evalwindow keep 1000 data to backtest out of sample
date.evalwindow=date.all(evalwindow); %corresponding date
estimwindow=1:T-3343; %sample to estimate parameters
r.evalwindow=rt(evalwindow,1); %define the return in each window
r.estimwindow=rt(estimwindow,1);
t=size(r.evalwindow,1);
varlvl=[0.05,0.01,0.005,0.001,0.0001]; %different varlvl 
cv=ones(1,size(varlvl,2)).*chi2inv(0.95,1); %critical values for test
cv2=ones(1,size(varlvl,2)).*chi2inv(0.95,2);
VarNames1={'VaR estimate','Failure rate','LR tStat','LRCCI','CV5','CC','CV25'}; %names for table
VarNames2={'Average VaR estimate','Failure rate','LR tStat','LRCCI','CV5','CC','Chi225'};

%descriptive statistics on whole sample
r.stat.mu=mean(rt);
r.stat.median=median(rt);
r.stat.max=max(rt);
r.stat.min=min(rt);
r.stat.sd=std(rt);
r.stat.skew=skewness(rt);
r.stat.ekurt=kurtosis(rt)-3;
r.stat.table=mat2dataset([T;r.stat.mu;r.stat.median;r.stat.max;r.stat.min;...
    r.stat.sd;r.stat.skew;r.stat.ekurt],'ObsNames',...
    {'Number of obs.','Mean','Median','Max','Min','Standard dev.','Skewness','Kurtosis'});
export(r.stat.table,'file','results/sumstatreturn.xls')

%check autocorrelation on eps^2 for arch effect
centered_rt=rt-mean(rt);
r.acf=sacf(centered_rt.^2,10,0,0);
r.pacf=spacf(centered_rt.^2,10,0,0);
[r.Ljung.stat,r.Ljung.pv]=ljungbox(centered_rt.^2,10);
r.temporaltable=[r.acf r.pacf r.Ljung.stat r.Ljung.pv];
r.temporaltable=mat2dataset(r.temporaltable,'ObsNames',...
    {'1','2','3','4','5','6','7','8','9','10'},'VarNames',{'ACF','PACF','LJungBox','pValue'});
export(r.temporaltable,'file','results/sumstatautocorr.xls')

g1=figure(); %return overall period
plot(date.all,rt)
title('log Return of SMI index overall period')
ylabel('Return')
xlabel('Time')
saveas(g1,'results/logreturnSMIall','png');

g2=figure(); %in term of histogram to see distribution
histfit(rt)
title('log Return of SMI index overall period')
ylabel('# Occurence')
xlabel('Return')
saveas(g2,'results/distrreturn','png');
%% STATIC
mdl=fitdist(r.estimwindow,'normal'); %fit a normal distribution
vari.normal.para=[mdl.mu,mdl.sigma]; %get the parameter mu and sigma for the norma
vari.normal.cov=sqrt(diag(mdl.ParameterCovariance)); %covariance matrix of the paramters
vari.normal.tstat=vari.normal.para./vari.normal.cov'; %perform t test
stat.normal=[vari.normal.para',vari.normal.tstat']; %table of result
stat.normal=mat2dataset(stat.normal,'varnames',{'Estimate','tStat'},...
    'obsnames',{'mu','sigma'});
export(stat.normal,'file','results/tablenormal.xls');

mdl2=fitdist(r.estimwindow,'tLocationScale'); %do the same for the student t
varmdl2=mdl2.sigma.^2*mdl2.nu/(mdl2.nu-2);
vari.student.para=[mdl2.mu,mdl2.sigma,mdl2.nu];
vari.student.cov=sqrt(diag(mdl2.ParameterCovariance));
vari.student.tstat=vari.student.para./vari.student.cov';
stat.student=[vari.student.para',vari.student.tstat'];
stat.student=mat2dataset(stat.student,'varnames',{'Estimate','tStat'},...
    'obsnames',{'mu','sigma','nu'});
export(stat.student,'file','results/tablestudent.xls');

%normal
FR.static.normal=zeros(size(varlvl)); %initilize vector
LR.static.normal=zeros(size(varlvl));
Rlt.static.normal=zeros(size(varlvl));
LRCCI.static.normal=zeros(size(varlvl));
CC.static.normal=zeros(size(varlvl));

for i=varlvl %using own function failure rate estime different stats we report for the normal
    VaR.static.normal(varlvl==i)=-(exp(vari.normal.para(1)+vari.normal.para(2)*norminv(i))-1);
    vecvar=ones(t,1).*VaR.static.normal(varlvl==i); %faillure rate needs a vector of VaR 
    [FR.static.normal(varlvl==i),LR.static.normal(varlvl==i),Rlt.static.normal(varlvl==i),...
        LRCCI.static.normal(varlvl==i),CC.static.normal(varlvl==i)]...
        =failurerate(r.evalwindow,vecvar,1-i,0.95);
end

f1=figure(); %plot of the return eval window vs the var
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-ones(t,1).*VaR.static.normal(1,1),'linewidth',1.5,'color','r')
hold off
title('Static VaR 5% with normal distribution')
xlim([date.evalwindow(1),date.evalwindow(end)])
legend('return in evaluation window','Normal VaR 5% ')
ylabel('Return')
xlabel('Time')
saveas(f1,'results/staticvar5normal','png');

%table of result 
static.result.normal=[VaR.static.normal.*100;FR.static.normal.*t;LR.static.normal;LRCCI.static.normal;cv;...
    CC.static.normal;cv2]';
static.result.normal=mat2dataset(static.result.normal,'ObsNames',string(varlvl),'VarNames',VarNames1);
export(static.result.normal,'file','results/staticnormal.xls')

%student
FR.static.student=zeros(size(varlvl)); %same as for normal 
LR.static.student=zeros(size(varlvl));
Rlt.static.student=zeros(size(varlvl));
LRCCI.static.student=zeros(size(varlvl));
CC.static.student=zeros(size(varlvl));

for i=varlvl %same as for normal
    VaR.static.student(varlvl==i)=-(exp(vari.student.para(1)+vari.student.para(2)*tinv(i,vari.student.para(3)))-1);
    vecvar=ones(t,1).*VaR.static.student(varlvl==i); % compute the VaR for all level
    [FR.static.student(varlvl==i),LR.static.student(varlvl==i),Rlt.static.student(varlvl==i),...
        LRCCI.static.student(varlvl==i),CC.static.student(varlvl==i)]...
        =failurerate(r.evalwindow,vecvar,1-i,0.95);
end

f2=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-ones(t,1).*VaR.static.student(1,1),'linewidth',1.5,'color','r')
hold off
title('Static VaR 5% with student distribution')
legend('return in evaluation window','Student VaR 5%')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f2,'results/staticvar5student','png');

static.result.student=[VaR.static.student.*100;FR.static.student.*t;LR.static.student;LRCCI.static.student;cv;...
    CC.static.student;cv2]';
static.result.student=mat2dataset(static.result.student,'ObsNames',string(varlvl),'VarNames',VarNames1);
export(static.result.student,'file','results/staticstudent.xls')

%Skewed T
options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'MaxIterations',10000,'ConstraintTolerance',1.0000e-6, ...
    'OptimalityTolerance',1.0000e-6,'MaxFunctionEvaluations',...
    100000,'Display','none');

x0=[mean(r.estimwindow),var(r.estimwindow),2.5,0]; %starting values
lb=[-Inf,0.000001,2.001,-0.999999]; %lowerbond for the paramters
ub=[Inf,Inf,Inf,0.999999]; %upper bonds


fun=@(para)LLskewt(r.estimwindow,para); %loglik of the skewed t, function return sigma2  

[vari.skewt.para,fval,~,~,~,~,vari.skewt.hessian]=fmincon(fun,x0,[],[],[],[],lb,ub,[],options); %run optim to find parameter
%take the hessian to have the covariance matrix 

vari.skewt.cov=inv(vari.skewt.hessian); %covariance matrix from the hessian
vari.skewt.para=vari.skewt.para'; %param for skewed t
vari.skewt.tstat=vari.skewt.para./sqrt(diag(vari.skewt.cov)); %perfom t test
stat.skewt=[vari.skewt.para,vari.skewt.tstat]; %result
stat.skewt=mat2dataset(stat.skewt,'varnames',{'Estimate','tStat'},...
    'obsnames',{'mu','sigma','nu','lambda'});
export(stat.skewt,'file','results/tableskewt.xls');

FR.static.skewt=zeros(size(varlvl));%same as normal
LR.static.skewt=zeros(size(varlvl));
Rlt.static.skewt=zeros(size(varlvl));
LRCCI.static.skewt=zeros(size(varlvl));
CC.static.skewt=zeros(size(varlvl));

for i=varlvl %we take the sqrt of para 2 as it is sigma2 
    VaR.static.skewt(varlvl==i)=-(exp(vari.skewt.para(1)+sqrt(vari.skewt.para(2))*...
        skewtinv(i,vari.skewt.para(3),vari.skewt.para(4)))-1); %compute the the Value at risk 
    
    vecvar=ones(t,1).*VaR.static.skewt(varlvl==i); 
    [FR.static.skewt(varlvl==i),LR.static.skewt(varlvl==i),Rlt.static.skewt(varlvl==i),...
        LRCCI.static.skewt(varlvl==i),CC.static.skewt(varlvl==i)]...
        =failurerate(r.evalwindow,vecvar,1-i,0.95);
end

f17=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-ones(t,1).*VaR.static.skewt(1,1),'linewidth',1.5,'color','r')
hold off
title('Static VaR 5% with Skewed T distribution')
legend('return in evaluation window','Skewed T VaR 5% ')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f17,'results/staticvar5skewt','png');

static.result.skewt=[VaR.static.skewt.*100;FR.static.skewt.*t;LR.static.skewt;LRCCI.static.skewt;cv;...
    CC.static.skewt;cv2]';
static.result.skewt=mat2dataset(static.result.skewt,'ObsNames',string(varlvl),'VarNames',VarNames1);
export(static.result.skewt,'file','results/staticskewt.xls')
%% GARCH
mu_garch.estimwindow=mean(r.estimwindow); %uncond. mean for have unexpected return
std_mugarch=sqrt((var(r.estimwindow)/length(r.estimwindow)));
mu_garch_tstat=mu_garch.estimwindow/std_mugarch;
eps.estimwindow=r.estimwindow-mu_garch.estimwindow; %get the eps in estim window
eps.garch=rt(T-t+1:T,1)-mu_garch.estimwindow; %get eps in eval window
figure()
plot(eps.estimwindow.^2)

%normal 
[garch_para.normal,~,garch_ht.normal,gev.rvcv,garch_VCV.normal]=tarch(eps.estimwindow,1,0,1,'NORMAL'); %estim garch normal
garch_tstat.normal=garch_para.normal./sqrt(diag(gev.rvcv)); %perfom ttest
stat.garch_normal=[garch_para.normal,garch_tstat.normal]; %table of result 
stat.garch_normal=mat2dataset(stat.garch_normal,'varnames',{'Estimate','tStat'},...
    'obsnames',{'Omega','Alpha','Beta',});  %rajouter la somme alpha+beta latex
export(stat.garch_normal,'file','results/tablestatnormalgarch.xls')
garch.normal.stdab=sqrt((gev.rvcv(2,2)+gev.rvcv(3,3)+2*gev.rvcv(2,3)));
garch.normal.ab_tstat=(garch_para.normal(2)+garch_para.normal(3)-1)/garch.normal.stdab;

sigma_garch.normal=zeros(t,1); %initialize vector
sigma_garch.normal(1,1)=garch_para.normal(1)+garch_para.normal(2)*eps.estimwindow(end).^2+...
    garch_para.normal(3)*garch_ht.normal(end);%compute first value for estimated sigma^2 in eval window

for i=2:t %compute the rest of sigma^2 for eval window
  sigma_garch.normal(i,1)=garch_para.normal(1,1)+garch_para.normal(2,1)*eps.garch(i-1,1).^2+...
      garch_para.normal(3,1)*sigma_garch.normal(i-1,1);
end

FR.garch.normal=zeros(size(varlvl));%initiliaze vector for faillure rate function
LR.garch.normal=zeros(size(varlvl)); %kupiec test 
Rlt.garch.normal=zeros(size(varlvl));
VaR.garch.normal=zeros(t,size(varlvl,2));
LRCCI.garch.normal=zeros(size(varlvl)); %independance test
CC.garch.normal=zeros(size(varlvl)); %christofersen test

sigma_garch.normal=[garch_ht.normal(end);sigma_garch.normal(1:end-1)]; %in t we compute VaRt+1 we need the last sigma in
%estimation window and all the rest but not the last one
for i=varlvl
    VaR.garch.normal(:,varlvl==i)=-(exp(mu_garch.estimwindow+sqrt(sigma_garch.normal(1:end)).*norminv(i))-1); 
    [FR.garch.normal(varlvl==i),LR.garch.normal(varlvl==i),Rlt.garch.normal(varlvl==i),...
        LRCCI.garch.normal(varlvl==i),CC.garch.normal(varlvl==i)]...
        =failurerate(r.evalwindow,VaR.garch.normal(:,varlvl==i),1-i,0.95);
end

garch.avgvar.normal=mean(VaR.garch.normal); %avg VaR for report 

f4=figure(); %GARCH normal VaR vs return in eval window graphical repres. of faillure rate
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-VaR.garch.normal(:,1),'linewidth',1.2,'color','r')
hold off
legend('return in evaluation window','GARCH normal VaR 5%')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
title('VaR 5% with GARCH Normal')
saveas(f4,'results/garchvar5normal','png');

%export table of result 
garch.result.normal=[garch.avgvar.normal.*100;FR.garch.normal.*t;LR.garch.normal;LRCCI.garch.normal;cv;...
    CC.garch.normal;cv2]';
garch.result.normal=mat2dataset(garch.result.normal,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(garch.result.normal,'file','results/garchnormal.xls')

%Student
%we use robust covariance as we are using non normal distribution from the
%QMLE estimaton methods
[garch_para.student,~,garch_ht.student,garch_rVCV.student]=tarch(eps.estimwindow,1,0,1,'STUDENTST');
garch_tstat.student=garch_para.student./sqrt(diag(garch_rVCV.student));
stat.garch_student=[garch_para.student,garch_tstat.student];
stat.garch_student=mat2dataset(stat.garch_student,'varnames',{'Estimate','tStat'},...
    'obsnames',{'Omega','Alpha','Beta','nu'}); %rajouter la somme alpha+beta latex
export(stat.garch_student,'file','results/tablestatstudentgarch.xls')
garch.student.stdab=sqrt((garch_rVCV.student(2,2)+garch_rVCV.student(3,3)+2*garch_rVCV.student(2,3)));
garch.student.ab_tstat=(garch_para.student(2)+garch_para.student(3)-1)/garch.student.stdab;

sigma_garch.student=zeros(t,1);
sigma_garch.student(1)=garch_para.student(1)+garch_para.student(2).*eps.estimwindow(end).^2+...
    garch_para.student(3).*garch_ht.student(end);

for i=2:t
  sigma_garch.student(i,1)=garch_para.student(1,1)+garch_para.student(2,1)*eps.garch(i-1,1).^2+...
      garch_para.student(3,1)*sigma_garch.student(i-1,1);
end

FR.garch.student=zeros(size(varlvl));
LR.garch.student=zeros(size(varlvl));
Rlt.garch.student=zeros(size(varlvl));
VaR.garch.student=zeros(t,size(varlvl,2));
LRCCI.garch.student=zeros(size(varlvl));
CC.garch.student=zeros(size(varlvl));

sigma_garch.student=[garch_ht.student(end);sigma_garch.student(1:end-1)];
for i=varlvl
    VaR.garch.student(:,varlvl==i)=-(exp(mu_garch.estimwindow+sqrt(sigma_garch.student(1:end)).*...
        tinv(i,garch_para.student(4,1)))-1);
    [FR.garch.student(varlvl==i),LR.garch.student(varlvl==i),Rlt.garch.student(varlvl==i),...
        LRCCI.garch.student(varlvl==i),CC.garch.student(varlvl==i)]...
        =failurerate(r.evalwindow,VaR.garch.student(:,varlvl==i),1-i,0.95);
end

garch.avgvar.student=mean(VaR.garch.student);

f6=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-VaR.garch.student(:,1),'linewidth',1.2,'color','r')
hold off
legend('return in evaluation window','GARCH Student VaR 5%')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
title('VaR 5% with GARCH Student t')
saveas(f6,'results/garchvar5student','png');

garch.result.student=[garch.avgvar.student.*100;FR.garch.student.*t;LR.garch.student;LRCCI.garch.student;cv;...
    CC.garch.student;cv2]';
garch.result.student=mat2dataset(garch.result.student,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(garch.result.student,'file','results/garchstudent.xls')

%SkwedT
%we use robust covariance as we are using non normal distribution
[garch_para.skewt,~,garch_ht.skewt,garch_rVCV.skewt]=tarch(eps.estimwindow,1,0,1,'SKEWT');
garch_tstat.skewt=garch_para.skewt./sqrt(diag(garch_rVCV.skewt));
stat.garch_skewdt=[garch_para.skewt,garch_tstat.skewt];
stat.garch_skewdt=mat2dataset(stat.garch_skewdt,'varnames',{'Estimate','tStat'},...
    'obsnames',{'Omega','Alpha','Beta','Nu','Lambda'}); %rajouter la somme alpha+beta latex
export(stat.garch_skewdt,'file','results/tablestatskewtgarch.xls')
garch.skewt.stdab=sqrt((garch_rVCV.skewt(2,2)+garch_rVCV.skewt(3,3)+2*garch_rVCV.skewt(2,3)));
garch.skewt.ab_tstat=(garch_para.skewt(2)+garch_para.skewt(3)-1)/garch.skewt.stdab;

sigma_garch.skewt=zeros(t,1);
sigma_garch.skewt(1)=garch_para.skewt(1)+garch_para.skewt(2).*eps.estimwindow(end).^2+...
    garch_para.skewt(3).*garch_ht.skewt(end);

for i=2:t
  sigma_garch.skewt(i,1)=garch_para.skewt(1,1)+garch_para.skewt(2,1)*eps.garch(i-1,1).^2+...
      garch_para.skewt(3,1)*sigma_garch.skewt(i-1,1);
end

FR.garch.skewt=zeros(size(varlvl));
LR.garch.skewt=zeros(size(varlvl));
Rlt.garch.skewt=zeros(size(varlvl));
VaR.garch.skewt=zeros(t,size(varlvl,2));
LRCCI.garch.skewt=zeros(size(varlvl));
CC.garch.skewt=zeros(size(varlvl));

sigma_garch.skewt=[garch_ht.skewt(end);sigma_garch.skewt(1:end-1)];
for i=varlvl
    VaR.garch.skewt(:,varlvl==i)=-(exp(mu_garch.estimwindow+sqrt(sigma_garch.skewt(1:end)).*...
        skewtinv(i,garch_para.skewt(4,1),garch_para.skewt(5,1)))-1); %changer le quantile
    [FR.garch.skewt(varlvl==i),LR.garch.skewt(varlvl==i),Rlt.garch.skewt(varlvl==i),...
        LRCCI.garch.skewt(varlvl==i),CC.garch.skewt(varlvl==i)]...
        =failurerate(r.evalwindow,VaR.garch.skewt(:,varlvl==i),1-i,0.95);
end

garch.avgvar.skewt=mean(VaR.garch.skewt);

j2=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-VaR.garch.skewt(:,1),'linewidth',1.2,'color','r')
hold off
title('GARCH Skew T Conditional volatility wiht error terms')
legend('return in evaluation window','GARCH Skew T VaR 5%')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
title('VaR 5% with GARCH Skewed t')
saveas(j2,'results/garchvar5skewt','png');

garch.result.skewt=[garch.avgvar.skewt.*100;FR.garch.skewt.*t;LR.garch.skewt;LRCCI.garch.skewt;cv;...
    CC.garch.skewt;cv2]';
garch.result.skewt=mat2dataset(garch.result.skewt,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(garch.result.skewt,'file','results/garchskewt.xls')
%% Historical simulation baisc
%define window size
hs.window=[250,500,1000,T-t];

hs.var.s=zeros(size(T-t:T-1,2),size(varlvl,2)); %initialize vectors
hs.var.m=zeros(size(T-t:T-1,2),size(varlvl,2));
hs.var.b=zeros(size(T-t:T-1,2),size(varlvl,2));
hs.var.t=zeros(size(T-t:T-1,2),size(varlvl,2));
interval=T-t:T-1;

for i=interval %compute the VaR in each rolling window for different size of rolling window and 
    hs.var.s(interval==i,:)=-quantile(rt(i-hs.window(1)+1:i),varlvl);
    hs.var.m(interval==i,:)=-quantile(rt(i-hs.window(2)+1:i),varlvl);
    hs.var.b(interval==i,:)=-quantile(rt(i-hs.window(3)+1:i),varlvl);
    hs.var.t(interval==i,:)=-quantile(rt(i-hs.window(4)+1:i),varlvl);
end

f7=figure(); %HS VaR for different size of rolling window
plot(date.evalwindow,hs.var.s(:,1),'linewidth',1.2,'color','r')
hold on
plot(date.evalwindow,hs.var.b(:,1),'linewidth',1.2,'color','b')
legend(' VaR 5% n=250','VaR 5% n=1000')
title('HS VaR 5% estimation 250 vs 1000 windows')
hold off
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('VaR')
xlabel('Time')
saveas(f7,'results/hsvar5250vs1000','png');

hs.avgvar.s=mean(hs.var.s); %mean of VaR for report purpose
hs.avgvar.m=mean(hs.var.m);
hs.avgvar.b=mean(hs.var.b);
hs.avgvar.t=mean(hs.var.t);

%250
FR.hs.s=zeros(size(varlvl)); %initialize vector for faillure rate
LR.hs.s=zeros(size(varlvl)); %kupiece
Rlt.hs.s=zeros(size(varlvl));
LRCCI.hs.s=zeros(size(varlvl)); %independence
CC.hs.s=zeros(size(varlvl)); %christoffersen

for i=varlvl
    [FR.hs.s(varlvl==i),LR.hs.s(varlvl==i),Rlt.hs.s(varlvl==i),...
        LRCCI.hs.s(varlvl==i),CC.hs.s(varlvl==i)]...
        =failurerate(r.evalwindow,hs.var.s(:,varlvl==i),1-i,0.95);
end

f8=figure(); %VaR 5% with rolling window of size 250
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.var.s(:,1),'linewidth',1.2,'color','r')
title('Return and VaR 5% window 250')
hold off
legend('return in evaluation window','HS VaR 5% n=250')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f8,'results/hsvar5w250','png');

hs.result.s=[hs.avgvar.s.*100;FR.hs.s.*t;LR.hs.s;LRCCI.hs.s;cv;CC.hs.s;cv2]';
hs.result.s=mat2dataset(hs.result.s,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(hs.result.s,'file','results/hs250.xls')

%500
FR.hs.m=zeros(size(varlvl)); %same as for 250
LR.hs.m=zeros(size(varlvl));
Rlt.hs.m=zeros(size(varlvl));
LRCCI.hs.m=zeros(size(varlvl));
CC.hs.m=zeros(size(varlvl));

for i=varlvl
    [FR.hs.m(varlvl==i),LR.hs.m(varlvl==i),Rlt.hs.m(varlvl==i),...
        LRCCI.hs.m(varlvl==i),CC.hs.m(varlvl==i)]...
        =failurerate(r.evalwindow,hs.var.m(:,varlvl==i),1-i,0.95);
end

f9=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.var.m(:,1),'linewidth',1.2,'color','r')
title('Return and VaR 5% window 500')
hold off
legend('return in evaluation window','HS VaR 5% n=500')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f9,'results/hsvar5w500','png');

hs.result.m=[hs.avgvar.m.*100;FR.hs.m.*t;LR.hs.m;LRCCI.hs.m;cv;CC.hs.m;cv2]';
hs.result.m=mat2dataset(hs.result.m,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(hs.result.m,'file','results/hs500.xls')

%1000
FR.hs.b=zeros(size(varlvl)); %same as for 250
LR.hs.b=zeros(size(varlvl));
Rlt.hs.b=zeros(size(varlvl));
LRCCI.hs.b=zeros(size(varlvl));
CC.hs.b=zeros(size(varlvl));

for i=varlvl
    [FR.hs.b(varlvl==i),LR.hs.b(varlvl==i),Rlt.hs.b(varlvl==i),...
        LRCCI.hs.b(varlvl==i),CC.hs.b(varlvl==i)]...
        =failurerate(r.evalwindow,hs.var.b(:,varlvl==i),1-i,0.95);
end

f10=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.var.b(:,1),'linewidth',1.2,'color','r')
title('Return and VaR 5% window 1000')
hold off
legend('return in evaluation window','HS VaR 5% n=1000')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f10,'results/hsvar5w1000','png');

hs.result.b=[hs.avgvar.b.*100;FR.hs.b.*t;LR.hs.b;LRCCI.hs.b;cv;CC.hs.b;cv2]';
hs.result.b=mat2dataset(hs.result.b,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(hs.result.b,'file','results/hs1000.xls')

%T-t=4037
FR.hs.t=zeros(size(varlvl)); %same as for 250
LR.hs.t=zeros(size(varlvl));
Rlt.hs.t=zeros(size(varlvl));
LRCCI.hs.t=zeros(size(varlvl));
CC.hs.t=zeros(size(varlvl));

for i=varlvl
    [FR.hs.t(varlvl==i),LR.hs.t(varlvl==i),Rlt.hs.t(varlvl==i),...
        LRCCI.hs.t(varlvl==i),CC.hs.t(varlvl==i)]...
        =failurerate(r.evalwindow,hs.var.t(:,varlvl==i),1-i,0.95);
end

f11=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.var.t(:,1),'linewidth',1.2,'color','r')
title('Return and VaR 5% window 3025')
hold off
legend('return in evaluation window','HS VaR 5% n=4025')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f11,'results/hsvar5w3025','png');

hs.result.t=[hs.avgvar.t.*100;FR.hs.t.*t;LR.hs.t;LRCCI.hs.t;cv;CC.hs.t;cv2]';
hs.result.t=mat2dataset(hs.result.t,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(hs.result.t,'file','results/hs3025.xls')
%% age-weighted HS
hs.age.lambda=0.98;%set value for lambda for the weight value from FRM
vec250=flip(1:hs.window(1))';
vec500=flip(1:hs.window(2))';
vec1000=flip(1:hs.window(3))';
veclong=flip(1:hs.window(4))';

for i=interval
    mat250=[vec250,rt(i-hs.window(1)+1:i),zeros(hs.window(1),2)]; %table for the own function varageweighted
    mat500=[vec500,rt(i-hs.window(2)+1:i),zeros(hs.window(2),2)]; %first vector is 1:length of the rolling window
    mat1000=[vec1000,rt(i-hs.window(3)+1:i),zeros(hs.window(3),2)]; %it determine the age the return, second column
    mat3025=[veclong,rt(i-hs.window(4)+1:i),zeros(hs.window(4),2)]; %are the return, 2 other columns to put the 
    hs.age.var.s(interval==i,:)=varageweighted(mat250,hs.age.lambda,varlvl); %weight 3rd and cum sum weight 4th 
    hs.age.var.m(interval==i,:)=varageweighted(mat500,hs.age.lambda,varlvl); %determine the var with the cumsum of weight
    hs.age.var.b(interval==i,:)=varageweighted(mat1000,hs.age.lambda,varlvl);
    hs.age.var.t(interval==i,:)=varageweighted(mat3025,hs.age.lambda,varlvl);
end

hs.age.avgvar.s=mean(hs.age.var.s); %for report purpose 
hs.age.avgvar.m=mean(hs.age.var.m);
hs.age.avgvar.b=mean(hs.age.var.b);
hs.age.avgvar.t=mean(hs.age.var.t);

f12=figure();
plot(date.evalwindow,hs.age.var.b(:,1),'linewidth',0.9,'color','r') %report al the VaR for different varlvl n=500
hold on
plot(date.evalwindow,hs.age.var.b(:,2),'linewidth',0.9,'color','b')
hold on
plot(date.evalwindow,hs.age.var.b(:,3),'linewidth',0.9,'color','k')
hold on 
plot(date.evalwindow,hs.age.var.b(:,4),'linewidth',0.9,'color','#D95319')
hold off
legend('5%','1%','0.5%','0.1%')
title('Age-Weighted VaR with n=1000')
ylabel('VaR')
xlabel('Time')
saveas(f12,'results/AgeWeightedVaR500','png');

%250
FR.hs.age.s=zeros(size(varlvl)); %initialize vector for faillure rate function
LR.hs.age.s=zeros(size(varlvl));
Rlt.hs.age.s=zeros(size(varlvl));
LRCCI.hs.age.s=zeros(size(varlvl));
CC.hs.age.s=zeros(size(varlvl));

for i=varlvl
    [FR.hs.age.s(varlvl==i),LR.hs.age.s(varlvl==i),Rlt.age.hs.s(varlvl==i),...
        LRCCI.hs.age.s(varlvl==i),CC.hs.age.s(varlvl==i)]...
        =failurerate(r.evalwindow,hs.age.var.s(:,varlvl==i),1-i,0.95);
end

f13=figure(); %VaR ageweighted 5% vs return eval window
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.age.var.s(:,1),'linewidth',1.2,'color','r')
title('Return and  HS Age-WeightedVaR 5% window 250')
hold off
legend('return in evaluation window','HS Age-Weighted VaR 5%')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f13,'results/hsagevar5w250','png');

%table of result
hs.result.age.s=[hs.age.avgvar.s.*100;FR.hs.age.s.*t;LR.hs.age.s;LRCCI.hs.age.s;cv;CC.hs.age.s;cv2]';
hs.result.age.s=mat2dataset(hs.result.age.s,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(hs.result.age.s,'file','results/hsage250.xls')

%500
FR.hs.age.m=zeros(size(varlvl));%same as for 250
LR.hs.age.m=zeros(size(varlvl));
Rlt.hs.age.m=zeros(size(varlvl));
LRCCI.hs.age.m=zeros(size(varlvl));
CC.hs.age.m=zeros(size(varlvl));

for i=varlvl
    [FR.hs.age.m(varlvl==i),LR.hs.age.m(varlvl==i),Rlt.hs.age.m(varlvl==i),...
        LRCCI.hs.age.m(varlvl==i),CC.hs.age.m(varlvl==i)]...
        =failurerate(r.evalwindow,hs.age.var.m(:,varlvl==i),1-i,0.95);
end

f14=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.age.var.m(:,1),'linewidth',1.2,'color','r')
title('Return and Age-Weighted VaR 5% window 500')
hold off
legend('return in evaluation window','HS Age-Weighted VaR 5%')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f14,'results/hsagevar5w500','png');

hs.result.age.m=[hs.age.avgvar.m.*100;FR.hs.age.m.*t;LR.hs.age.m;LRCCI.hs.age.m;cv;CC.hs.age.m;cv2]';
hs.result.age.m=mat2dataset(hs.result.age.m,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(hs.result.age.m,'file','results/hsage500.xls')

%1000
FR.hs.age.b=zeros(size(varlvl)); %same as for 250
LR.hs.age.b=zeros(size(varlvl));
Rlt.hs.age.b=zeros(size(varlvl));
LRCCI.hs.age.b=zeros(size(varlvl));
CC.hs.age.b=zeros(size(varlvl));

for i=varlvl
    [FR.hs.age.b(varlvl==i),LR.hs.age.b(varlvl==i),Rlt.hs.age.b(varlvl==i),...
        LRCCI.hs.age.b(varlvl==i),CC.hs.age.b(varlvl==i)]...
        =failurerate(r.evalwindow,hs.age.var.b(:,varlvl==i),1-i,0.95);
end

f15=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.age.var.b(:,1),'linewidth',1.2,'color','r')
title('Return and Age-Weighted VaR 5% window 1000')
hold off
legend('return in evaluation window','HS Age-Weighted VaR 5%')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f15,'results/hsagevar5w1000','png');

hs.result.age.b=[hs.age.avgvar.b.*100;FR.hs.age.b.*t;LR.hs.age.b;LRCCI.hs.age.b;cv;CC.hs.age.b;cv2]';
hs.result.age.b=mat2dataset(hs.result.age.b,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(hs.result.age.b,'file','results/hsage1000.xls')

%T-t=4037
FR.hs.age.t=zeros(size(varlvl)); %same as for 250
LR.hs.age.t=zeros(size(varlvl));
Rlt.hs.age.t=zeros(size(varlvl));
LRCCI.hs.age.t=zeros(size(varlvl));
CC.hs.age.t=zeros(size(varlvl));

for i=varlvl
    [FR.hs.age.t(varlvl==i),LR.hs.age.t(varlvl==i),Rlt.hs.age.t(varlvl==i),...
        LRCCI.hs.age.t(varlvl==i),CC.hs.age.t(varlvl==i)]...
        =failurerate(r.evalwindow,hs.age.var.t(:,varlvl==i),1-i,0.95);
end

f16=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.age.var.t(:,1),'linewidth',1.2,'color','r')
title('Return and Age-Weighted VaR 5% n=4037')
hold off
legend('return in evaluation window','HS Age Weighted VaR 5%')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f16,'results/hsagevar5w3025','png');

hs.result.age.t=[hs.age.avgvar.t.*100;FR.hs.age.t.*t;LR.hs.age.t;LRCCI.hs.age.t;cv;CC.hs.age.t;cv2]';
hs.result.age.t=mat2dataset(hs.result.age.t,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(hs.result.age.t,'file','results/hsage3025.xls')
%% Extreme value theory 
%GEV 
gev.r=-r.estimwindow-mean(-r.estimwindow); %centralized loss function 
gev.z=gev.r./sqrt(garch_ht.normal); %get the z process using garch from above
[gev.zautocorr, gev.pvalue]=ljungbox(gev.z,10); %check autocorrelation innovation process
gev.g=60; %choose subsample size of 1 quarter
gev.l=length(gev.z); 
gev.end=gev.g:gev.g:gev.l-gev.g-1; %end value of all subsample but not the last one
gev.sample=zeros(gev.g,length(gev.l-gev.g)); 
gev.ljunrest=mat2dataset([gev.zautocorr, gev.pvalue],'ObsNames',...
    {'1','2','3','4','5','6','7','8','9','10'},'VarNames',{'LJungBox','pValue'});
export(gev.ljunrest,'file','results/zautocorr.xls')


for i=gev.end %get return of all subsample in column vector
    gev.sample(:,gev.end==i)=gev.z(i-gev.g+1:i,1);
end
gev.lastsample=gev.z((gev.l-gev.g):end); %compute the last subsample in order to have all value included

gev.max=max(gev.sample,[],1)'; %compute the max of all subsample 
gev.max(end+1,1)=max(gev.lastsample); %max of th last sumbsample which is different in size 
gev.max=sort(gev.max); 
gev.ecdf=homemade_ecdf(gev.max); %compute the ecdf to have qqplot
gev.gumbel=-log(-log(gev.ecdf));
f19=figure();
plot(gev.gumbel(1:end-1),gev.max)
title('QQplot')
xlabel('Gumbel quantile')
ylabel('Empirical quantile')
saveas(f19,'results/QQplotGev','png');

%test the garch process using QMLE
gev.garch.tstat=garch_para.normal ./sqrt(diag(gev.rvcv));
gev.garch.result=[garch_para.normal,gev.garch.tstat];
gev.garch.result=mat2dataset(gev.garch.result,'varnames',{'Estimate','tStat'},...
    'obsnames',{'Omega','Alpha','Beta',});  %rajouter la somme alpha+beta latex
export(gev.garch.result,'file','results/tablestatgpdgarch.xls')

gev.distr=fitdist(gev.max,'GeneralizedExtremeValue'); %fite the GEV 
gev.ksi=gev.distr.k; gev.psi=gev.distr.sigma; gev.mu=gev.distr.mu; %get all para and test it 
gev.tstat=gev.distr.ParameterValues'./sqrt(diag(gev.distr.ParameterCovariance)); %compute the tstat
gev.result=[gev.distr.ParameterValues',gev.tstat];
gev.result=mat2dataset(gev.result,'varnames',{'Estimate','tStat'},...
    'obsnames',{'Ksi','Psi','Mu'}); 
export(gev.result,'file','results/tablestatGEV21.xls')

gev.zq=(gev.mu+(gev.psi/gev.ksi).*(((-gev.g.*log(1-varlvl)).^(-gev.ksi))-1)); %compute the qz quantile 
FR.tail.gev=zeros(size(varlvl));
LR.tail.gev=zeros(size(varlvl));
Rlt.tail.gev=zeros(size(varlvl));
LRCCI.tail.gev=zeros(size(varlvl));
CC.tail.gev=zeros(size(varlvl));

for i=varlvl %compute the var and construct the faillure rate
    gev.var(:,varlvl==i)=(exp(mean(r.estimwindow)+sqrt(sigma_garch.normal).*gev.zq(varlvl==i))-1); 
    [FR.tail.gev(varlvl==i),LR.tail.gev(varlvl==i),Rlt.tail.gev(varlvl==i),...
        LRCCI.tail.gev(varlvl==i),CC.tail.gev(varlvl==i)]...
        =failurerate(r.evalwindow,gev.var(:,varlvl==i),1-i,0.95);
end

f19=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-gev.var(:,1),'linewidth',1.2,'color','r')
hold off
title('Tail VaR 5% with GEV n=60')
legend('return in evaluation window','tail GEV VaR 5% n=60')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f19,'results/tailvar5gev','png');

gev.avgvar=mean(gev.var);
tail.result.gev=[gev.avgvar.*100;FR.tail.gev.*t;LR.tail.gev;LRCCI.tail.gev;cv;...
    CC.tail.gev;cv2]';
tail.result.gev=mat2dataset(tail.result.gev,'ObsNames',string(varlvl),'VarNames',VarNames1);
export(tail.result.gev,'file','results/taillgev.xls')

%GPD
gpd.r=-(r.estimwindow); %take the loss function 
gpd.mu=mean(gpd.r); %compute the mean assume time invarient
gpd.eps=gpd.r-gpd.mu; %centralized  to get eps
[gpd.garch.para,~,gpd.garch.ht,gpd.garch.VCV]=tarch(gpd.eps,1,0,1,'NORMAL'); % estim garch to find the zt
gpd.garch.tstat=gpd.garch.para./sqrt(diag(gpd.garch.VCV));
gpd.garch.result=[gpd.garch.para,gpd.garch.tstat];
gpd.garch.result=mat2dataset(gpd.garch.result,'varnames',{'Estimate','tStat'},...
    'obsnames',{'Omega','Alpha','Beta',});  %rajouter la somme alpha+beta latex
export(gpd.garch.result,'file','results/tablestatgpdgarch.xls')

gpd.z=gpd.eps./sqrt(gpd.garch.ht); %using the garch ht we find the zt
gpd.srtz=sort(gpd.z);
gpd.q=round(size(gpd.srtz,1)*0.1,-2); %defin thershold, 10% of the data in estimwindow 
gpd.upper=gpd.srtz(end-gpd.q+1:end);
gpd.u=gpd.srtz(end-gpd.q);
gpd.excessr=gpd.upper-gpd.u; %get the return above u for zt


gpd.dist.para=mle(gpd.excessr,'distribution','gp'); %estime by mle the parameter
gpd.dist.cov=mlecov(gpd.dist.para,gpd.excessr,'pdf',@gppdf); %retrieve the covariance with the above parameter
gpd.ksi=gpd.dist.para(1);gpd.psi=gpd.dist.para(2);
gpd.dist.tstat=gpd.dist.para'./sqrt(diag(gpd.dist.cov));%perform tstat
gpd.dist.result=[gpd.dist.para',gpd.dist.tstat];
gpd.dist.result=mat2dataset(gpd.dist.result,'varnames',{'Estimate','tStat'},...
    'obsnames',{'Ksi','Psi'});  %rajouter la somme alpha+beta latex 
export(gpd.dist.result,'file','results/tablestatgpddist.xls')

gpd.zq=gpd.u+(gpd.psi/gpd.ksi).*((((T-t)/gpd.q).*varlvl).^(-gpd.ksi)-1); %qz quantile

gpd.feps=rt(T-t+1:T,1)-gpd.mu; %get the central return in evalwindow

gpd.garch.esigma=zeros(t,1); %get expected sigma 
gpd.garch.esigma(1,1)=gpd.garch.para(1)+gpd.garch.para(2)*gpd.eps(end).^2+... %
    gpd.garch.para(2)*garch_ht.normal(end); %estime first sigmat from t-1 t is first date in evalwindow

for i=2:t %find estim of sigmat+1 for all eval window
  gpd.garch.esigma(i,1)=gpd.garch.para(1,1)+gpd.garch.para(2,1)*gpd.feps(i-1,1).^2+...
      gpd.garch.para(3,1)*gpd.garch.esigma(i-1,1);
end

FR.gpd=zeros(size(varlvl));
LR.gpd=zeros(size(varlvl));
Rlt.gpd=zeros(size(varlvl));
LRCCI.gpd=zeros(size(varlvl));
CC.gpd=zeros(size(varlvl));

gpd.garch.esigma=[gpd.garch.ht(end);gpd.garch.esigma(1:end-1)]; %get all the the corresponding sigma as prevision
for i=varlvl %compute the var
    gpd.var(:,varlvl==i)=(exp(gpd.mu+sqrt(gpd.garch.esigma(1:end)).*gpd.zq(varlvl==i))-1); 
    [FR.gpd(varlvl==i),LR.gpd(varlvl==i),Rlt.gpd(varlvl==i),...
        LRCCI.gpd(varlvl==i),CC.gpd(varlvl==i)]...
        =failurerate(r.evalwindow,gpd.var(:,varlvl==i),1-i,0.95);
end

gpd.avgvar=mean(gpd.var);

f21=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-gpd.var(:,1),'r','Linewidth',1.5)
hold on 
plot(date.evalwindow,-gev.var(:,1),'linewidth',1.2,'color','b')
hold off
legend('return in evaluation window','Tail GPD VaR 5%','GEV')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f21,'results/gpdvar5','png');

gpd.result=[gpd.avgvar.*100;FR.gpd.*t;LR.gpd;LRCCI.gpd;cv;...
    CC.gpd;cv2]';
gpd.result=mat2dataset(gpd.result,'ObsNames',string(varlvl),'VarNames',VarNames2);
export(gpd.result,'file','results/gpd.xls')

%% graph
%Static
f28=figure();
subplot(2,2,1); %plot of the return eval window vs the var
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-ones(t,1).*VaR.static.normal(1,1),'linewidth',1.5,'color','r')
hold off
title('Static VaR 5% with normal distribution')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
legend('return in evaluation window','Normal VaR 5% ','location','southoutside')

subplot(2,2,2);
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-ones(t,1).*VaR.static.student(1,1),'linewidth',1.5,'color','r')
hold off
title('Static VaR 5% with student distribution')
legend('Return in evaluation window','Student VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,3);
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-ones(t,1).*VaR.static.skewt(1,1),'linewidth',1.5,'color','r')
hold off
title('Static VaR 5% with Skewed T distribution')
legend('Return in evaluation window','Skewed T VaR 5% ','location','southoutside')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f28,'results/resumestatic','png');

% GARCH
f29=figure(); %GARCH normal VaR vs return in eval window graphical repres. of faillure rate
subplot(2,2,1)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-VaR.garch.normal(:,1),'linewidth',0.7,'color','r')
hold off
title('Dynamic VaR 5% with Normal distribution')
legend('Return in evaluation window','GARCH normal VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,2)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-VaR.garch.student(:,1),'linewidth',0.7,'color','r')
hold off
title('Dynamic VaR 5% with Student T distribution')
legend('return in evaluation window','GARCH Student VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,3)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-VaR.garch.skewt(:,1),'linewidth',0.7,'color','r')
hold off
title('Dynamic VaR 5% with Skewed T distribution')
legend('return in evaluation window','GARCH Skew T VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])
saveas(f29,'results/resumegarch','png');

%HS
f30=figure(); %VaR 5% with rolling window of size 250
subplot(2,2,1)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.var.s(:,1),'linewidth',0.7,'color','r')
title('Return and VaR 5% window 250')
hold off
legend('return in evaluation window','HS VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,2)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.var.m(:,1),'linewidth',0.7,'color','r')
title('Return and VaR 5% window 500')
hold off
legend('return in evaluation window','HS VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,3)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.var.b(:,1),'linewidth',0.7,'color','r')
title('Return and VaR 5% window 1000')
hold off
legend('return in evaluation window','HS VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,4)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.var.t(:,1),'linewidth',0.7,'color','r')
title('Return and VaR 5% window 3025')
hold off
legend('return in evaluation window','HS VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])
saveas(f30,'results/resumehs','png');

%age-weighted
f31=figure(); %VaR ageweighted 5% vs return eval window
subplot(2,2,1)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.age.var.s(:,1),'linewidth',0.7,'color','r')
title('Return and  HS Age-WeightedVaR 5% window 250')
hold off
legend('return in evaluation window','HS Age-Weighted VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,2)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.age.var.m(:,1),'linewidth',0.7,'color','r')
title('Return and Age-Weighted VaR 5% window 500')
hold off
legend('return in evaluation window','HS Age-Weighted VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,3)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.age.var.b(:,1),'linewidth',0.7,'color','r')
title('Return and Age-Weighted VaR 5% window 1000')
hold off
legend('return in evaluation window','HS Age-Weighted VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(2,2,4)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-hs.age.var.t(:,1),'linewidth',0.7,'color','r')
title('Return and Age-Weighted VaR 5% window 3025')
hold off
legend('return in evaluation window','HS Age Weighted VaR 5%','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])
saveas(f31,'results/resumeageweighted','png');

%tail index
f32=figure();
subplot(1,2,1)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-ones(t,1).*gev.var(:,1),'linewidth',0.7,'color','r')
hold off
title('Tail VaR 5% with GEV n=60')
legend('return in evaluation window','tail GEV VaR 5% n=60','location','southoutside')
ylabel('Return')
xlabel('Time')
xlim([date.evalwindow(1),date.evalwindow(end)])

subplot(1,2,2)
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-gpd.var(:,1),'linewidth',0.7,'color','r')
hold off
legend('return in evaluation window','Tail GPD VaR 5% N=400','location','southoutside')
ylabel('Return')
xlabel('Time')
title('Tail VaR 5% GPD')
xlim([date.evalwindow(1),date.evalwindow(end)])
saveas(f32,'results/resumetail','png');

f35=figure;
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-ones(t,1).*gev.var(:,1),'linewidth',0.9,'color','r')
hold on
plot(date.evalwindow,-gpd.var(:,1),'linewidth',0.9,'color','b')
hold off
legend('return in evaluation window','Tail GEV VaR 5%','Tail GPD VaR 5% N=400','location','southoutside')
ylabel('Return')
xlabel('Time')
title('Tail VaR 5% GPD vs Tail VaR 5% GEV')
xlim([date.evalwindow(1),date.evalwindow(end)])
saveas(f35,'results/compgevgpd5','png');

f33=figure();
plot(-5:0.01:5,skewtpdf(-5:0.01:5,garch_para.skewt(4),garch_para.skewt(5)))
hold on
plot(-5:0.01:5,tpdf(-5:0.01:5,garch_para.student(4)))
title('GARCH Student t vs Skewed t')
legend('Skewed t','Student t')
ylabel('pdf')
xlabel('Domain')
saveas(f33,'results/garchskewedtvsstudent','png');

f23=figure();
plot(date.evalwindow,r.evalwindow)
hold on
plot(date.evalwindow,-gpd.var(:,2), 'LineWidth',1.5,'Color','black')
hold on
plot(date.evalwindow,-VaR.garch.student(:,2),'.r','MarkerSize',0.2)
hold off
legend('return in evaluation window','Tail GPD VaR 1%','GARCH Skewed T VaR 1%')
title('Comparison tail approach and GARCH Skewed T')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f23,'results/gpdvsgarchnormal','png');

f34=figure;
plot(0:0.1:7,gevpdf(0:0.1:7,gev.ksi,gev.psi,gev.mu),'LineWidth',1.5,'Color','black')
hold on 
plot(0:0.1:7,gppdf(0:0.1:7,gpd.ksi,gpd.psi,gpd.u),'LineWidth',1.5,'Color','r')
legend('GEV','GPD')
ylabel('f(x)')
xlabel('Domain')
title('Comparaison of the GEV and GPD distribution')
saveas(f34,'results/distrgpdvsgev','png')

%illustration of violation of the GEV
test=r.evalwindow;
test(test<-gev.var(:,1))=1;
test(test~=1)=0;
f36=figure();
plot(date.evalwindow,r.evalwindow,'k')
hold on
plot(date.evalwindow,-gev.var(:,1),'linewidth',1.2,'color','r')
hold on
plot(date.evalwindow,test.*0.01,'m')
hold off
title('Tail VaR 5% with GEV n=60')
legend('return in evaluation window','tail GEV VaR 5% n=60','Violations')
xlim([date.evalwindow(1),date.evalwindow(end)])
ylabel('Return')
xlabel('Time')
saveas(f36,'results/violationillustration','png');


