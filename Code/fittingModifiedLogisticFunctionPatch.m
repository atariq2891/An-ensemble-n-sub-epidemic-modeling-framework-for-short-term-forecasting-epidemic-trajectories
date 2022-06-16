function [Phatss,npatches,onset_thr,curves,fittedcurves,bestfit,data1]=fittingModifiedLogisticFunction(datafilename1,DT,epidemic_period,M,scalar1,flagX,numstartpoints)

global flag1 timevect ydata yfit

global I0 npatches onset_thr

flag1=flagX;

global invasions
global timeinvasions
global Cinvasions

global npatches_fixed

global onset_fixed

global method1
global dist1
global factor1

global smoothfactor1

global LBe UBe


close all

% ABC fitting first to obtain estimates for npatches and onset_thr

[npatches,onset_thr,P0]=fittingModifiedLogisticFunctionPatchABC(datafilename1,DT,epidemic_period,M,scalar1,flag1,numstartpoints)

%pause

close all

% fitting

Phats=[];

goodnessfit=[];

modelfits=[];

SSstats1=[];

SSstats2=[];

data1=load(datafilename1)

data1=data1(epidemic_period,:);

data1(:,2)=data1(:,2)/scalar1;

I0=data1(1,2); % initial condition

if I0==0
    
    data1=data1(2:end,:);
    
end

data=data1(:,2);

[max1,index1]=max(data);

timevect=(data1(:,1))*DT;

rs1=P0(1:npatches);

ps1=P0(npatches+1:2*npatches);

as1=P0(2*npatches+1:3*npatches);

Ks1=P0(3*npatches+1:4*npatches);

alpha=P0(end-1);

d=P0(end);


z=[rs1 ps1 as1 Ks1 alpha d];

[LB1,UB1]=getbounds(npatches);

LB=[LB1 LBe]; %r p a K alpha d

UB=[UB1 UBe];



I0=data(1); % initial condition

hold on

%options=optimset('MaxFunEvals',3000,'MaxIter',3000,'Algorithm','trust-region-reflective','TolFun',1e-6,'TolX',1e-6);

%options=[];

%options = optimoptions('lsqcurvefit','UseParallel',true,...
%    'TolX',10^(-5),'TolFun',10^(-5));

%[P,resnorm,residual,exitflag,output,lambda,J]=lsqcurvefit(@plotModifiedLogisticGrowthPatch1,z,timevect,smooth(data,smoothfactor1),LB,UB,options,I0,npatches,onset_thr,flag1);

% resnorm is the SSE which is given by sum(residual.^2)
% P is the vector with the estimated parameters


%method1=3; %LSQ=0, MLE Poisson=1, Pearson chi-squared=2, MLE (Neg Binomial)=3


ydata=smooth(data,smoothfactor1);

options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-9,'MaxFunEvals',20000,'MaxIter',20000);

%options=optimoptions('fmincon','Algorithm','sqp','tolfun',10^-6,'TolX',10^-6,'MaxFunEvals',3200,'MaxIter',3200);

f=@plotModifiedLogisticGrowthPatchMethods1;

problem = createOptimProblem('fmincon','objective',f,'x0',z,'lb',LB,'ub',UB,'options',options);

%ms = MultiStart('PlotFcns',@gsplotbestf);
ms = MultiStart('Display','final');
%ms=MultiStart;

%rpoints = RandomStartPointSet('NumStartPoints',numstartpoints); % start with a few random starting sets in addition to the guess supplied by the user (z)

pts = z;
tpoints = CustomStartPointSet(z);
%allpts = {tpoints,rpoints};
allpts = {tpoints};

[P,fval,flagg,outpt,allmins] = run(ms,problem,allpts);



rs_hat=P(1,1:npatches);
ps_hat=P(1,npatches+1:2*npatches);
as_hat=P(1,2*npatches+1:3*npatches);
Ks_hat=P(1,3*npatches+1:4*npatches);

alpha_hat=P(end-1);
d_hat=P(end);

IC=zeros(npatches,1);

if onset_fixed==0
    IC(1,1)=I0;
    IC(2:end,1)=1;
    
    invasions=zeros(npatches,1);
    timeinvasions=zeros(npatches,1);
    Cinvasions=zeros(npatches,1);
    
    invasions(1)=1;
    timeinvasions(1)=0;
    Cinvasions(1)=0;
else
    IC(1:end,1)=I0./length(IC(1:end,1));
    
    invasions=zeros(npatches,1);
    timeinvasions=zeros(npatches,1);
    Cinvasions=zeros(npatches,1);
    
    invasions(1:end)=1;
    timeinvasions(1:end)=0;
    Cinvasions(1:end)=0;
end


[~,x]=ode15s(@modifiedLogisticGrowthPatch,timevect,IC,[],rs_hat,ps_hat,as_hat,Ks_hat,npatches,onset_thr,flag1);
%x=ode5(@modifiedLogisticGrowthPatch,timevect,IC,r_hat,p_hat,a_hat,K_hat,npatches,onset_thr,q_hat,flag1);


figure(10)

for j=1:npatches
    
    incidence1=[x(1,j);diff(x(:,j))];
    
    plot(timevect,incidence1)
    hold on
    
end

y=sum(x,2);

totinc=[y(1,1);diff(y(:,1))];

if onset_fixed==0
    totinc(1)=totinc(1)-(npatches-1);
end

bestfit=totinc;

plot(timevect,totinc,'r')
hold on
plot(timevect,data,'bo')
xlabel('Time (days)');
ylabel('Data and best SIR fit')

Ptrue=[rs_hat ps_hat as_hat Ks_hat alpha_hat d_hat]

'generate simulation study to derive parameter uncertainty..'


% generate simulation study to derive parameter uncertainty

yi=cumsum(totinc);

z=Ptrue; %r p a Ks alpha(MLE-Neg binomial)

Phatss=zeros(M,4*npatches+2);
SSEs=zeros(M,1);
curves=[];

fittedcurves=zeros(length(yi),M);


for real=1:M
    
   yirData=AddPoissonError(yi,1,dist1,factor1,d_hat);

    %
    %     yirData=zeros(length(yi),1);
    %
    %     yirData(1)=yi(1);
    %
    %     for t=2:length(yi)
    %         lambda=abs(yi(t)-yi(t-1));
    %         yirData(t,1)=poissrnd(lambda,1,1);
    %     end
    %
    curves=[curves yirData];
    
    rs1=P0(1:npatches);
    
    ps1=P0(npatches+1:2*npatches);
    
    as1=P0(2*npatches+1:3*npatches);
    
    Ks1=P0(3*npatches+1:4*npatches);
    
    alpha=P0(end-1);
    
    d=P0(end);
    
    z=[rs1 ps1 as1 Ks1 alpha d];
    
    
    [LB1,UB1]=getbounds(npatches);
    
    LB=[LB1 LBe]; %r p a K alpha d
    
    UB=[UB1 UBe];

    
    I0=data(1); % initial condition
    
    hold on
    
    %options=optimset('MaxFunEvals',3000,'MaxIter',3000,'Algorithm','trust-region-reflective','TolFun',1e-5,'TolX',1e-5);
    
    %options=[];
    
    %options = optimoptions('lsqcurvefit','UseParallel',true,...
    %    'TolX',10^(-5),'TolFun',10^(-5));
    
    
    ydata=yirData;
    
    options=optimoptions('fmincon','Algorithm','sqp','StepTolerance',1.0000e-9,'MaxFunEvals',20000,'MaxIter',20000);
    
    f=@plotModifiedLogisticGrowthPatchMethods1;
    
    problem = createOptimProblem('fmincon','objective',f,'x0',z,'lb',LB,'ub',UB,'options',options);
    
    %ms = MultiStart('PlotFcns',@gsplotbestf);
    %ms = MultiStart('Display','final');
    ms = MultiStart('Display','off');

    %ms=MultiStart;
    
    %one random start guess
    rpoints = RandomStartPointSet('NumStartPoints',1); % start with a few random starting sets in addition to the guess supplied by the user (z)
    
    pts = z;
    tpoints = CustomStartPointSet(z);
    allpts = {tpoints,rpoints};
    %allpts = {tpoints};
    
    [P,fval,flagg,outpt,allmins] = run(ms,problem,allpts);
    
    fittedCurve1=yfit;
    
    fittedcurves(:,real)=fittedCurve1;
    
    Phatss(real,:)=P;
    
end

