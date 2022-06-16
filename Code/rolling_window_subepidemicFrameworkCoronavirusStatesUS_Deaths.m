function rolling(outbreakx,caddate1)

%%
% <============================================================================>
% <=================== Declare global variables ===============================>
% <============================================================================>

global method1
global npatches_fixed

global onset_fixed

global dist1
global factor1
global smoothfactor1

global calibrationperiod1

% <============================================================================>
% <================================ Datasets properties =======================>
% <============================================================================>
 
%outbreakx=2;
%caddate1='04-20-20';
%caddate1='06-22-20';
%caddate1='06-29-20';
%caddate1='11-15-21';

        
cadregion='USA'; % string indicating the region of the time series (USA, Chile, Mexico, Nepal, etc)

caddisease='coronavirus'; % string indicating the name of the disease

datatype='deaths'; % string indicating the nature of the data (cases, deaths, hospitalizations, etc)

DT=1; % temporal resolution in days (1=daily data, 7=weekly data).

if DT==1
    cadtemporal='daily';
elseif DT==7
    cadtemporal='weekly';
end

% Name of the file containing the cumulative time series data (rows=time, cols=regions)
datafilename1=strcat('cumulative-',cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-',caddate1,'.txt'); %data file with all time series across areas/regions

% Name of the file for the adjusted incidence data file for a specific region and
% after removing early zeros.
datafilename2=strcat(cadtemporal,'-',caddisease,'-',datatype,'-',cadregion,'-state-',num2str(outbreakx),'-',caddate1,'.txt'); % specific adjusted data file for a given region/area

datevecfirst1=[2020 02 27]; % date corresponding to the first data point in time series data


% <============================================================================>
% <============================Adjustments to data ============================>
% <============================================================================>

smoothfactor1=7; % <smoothfactor1>-day rolling average smoothing of the case series

calibrationperiod1=90; % calibrates model using the most recent <calibrationperiod1> days  where calibrationperiod1<length(data1)


% <=============================================================================>
% <=========================== Statistical method ==============================>
% <=============================================================================>

method1=0; % Type of estimation method: 0 = LSQ

% LSQ=0,
% MLE Poisson=1,
% Pearson chi-squared=2,
% MLE (Neg Binomial)=3, with VAR=mean+alpha*mean;
% MLE (Neg Binomial)=4, with VAR=mean+alpha*mean^2;
% MLE (Neg Binomial)=5, with VAR=mean+alpha*mean^d;

% Define dist1 which is the type of error structure:
switch method1
    
    case 0
        
        % Normnal distribution to model error structure
        dist1=0;
       
        %dist1=1; % error structure type (Poisson=1; NB=2)
        
        %factor1=1; % scaling factor for VAR=factor1*mean
        
        
    case 3
        dist1=3; % VAR=mean+alpha*mean;
        
    case 4
        dist1=4; % VAR=mean+alpha*mean^2;
        
    case 5
        dist1=5; % VAR=mean+alpha*mean^d;
        
end

numstartpoints=10; % Number of initial guesses for optimization procedure using MultiStart

topmodelsx=4; % number of best fitting models (based on AICc) that will be generated to derive an ensemble model

M=300; % number of bootstrap realizations to characterize parameter uncertainty


% <==============================================================================>
% <========================= Mathematical model =================================>
% <==============================================================================>

npatches_fixed=2; % maximum number of subepidemics considered in epidemic wave model fit

if npatches_fixed==1
    topmodelsx=1;
end

%block1=1; %type of growth model used to describe a subepidemic

GGM=0;  % 0 = GGM
GLM=1;  % 1 = GLM
GRM=2;  % 2 = GRM
LM=3;   % 3 = LM
RICH=4; % 4 = Richards

flag1=[RICH GLM]; % Sequence of subepidemic growth models considered in epidemic wave


onset_fixed=0; % flag to indicate if the onset timing of subepidemics fixed at time 0 (onset_fixed=1) or not (onset_fixed=0).


% <==============================================================================>
% <============ Load data and proceed to parameter estimation ===================>
% <==============================================================================>

%pause(10*rand)
data=load(strcat('./input/',datafilename1)); % load time series dataset (rows=time, cols=regions/areas)

dataprov=data';

timelags=0; % keeps track of the epidemic onset timing (nonzero case)

for outbreak1=outbreakx
    
    close all
    
    data1=dataprov(outbreak1,:)'; % Cumulative curve
    
    data1=[data1(1);diff(data1)]; % Incidence curve
    
    clear dataprov
    clear data
    
    % <==================================================================================>
    % <============ select length of the calibration period =============================>
    % <==================================================================================>
    
    % select length of the calibration period
    if calibrationperiod1<length(data1)
        
        data1(1:1:end-calibrationperiod1)=0;
        
    end
    
    %plot(data1)
    %hold on
    
    % <=============================================================================================>
    % <========== remove early zeros until a non-zero data point is found in the time series========>
    % <=============================================================================================>
    
    index1=find(data1>0);
    
    data1=data1(index1(1):end,1);
    
    timelags=index1(1)-1;
    
    
    %     if 0
    %         % remove early low testing period and start from the monotonic increasing trend
    %
    %         [max2, index2]=max(data1);
    %
    %         [min1, index1]=min(data1(1:index2));
    %         [max2, index2]=max(data1(index1:end));
    %
    %         index2=index1+index2-1;
    %
    %         [index1, min1]=find(data1(1:index2)==min1);
    %         index1=index1(length(index1));
    %
    %
    %         if index1<index2
    %             data1=data1(index1:end);
    %
    %             timelags=timelags+index1(1)-1;
    %
    %         else
    %             data1=[];
    %             length(data1)
    %             outbreak1
    %             'small number of data points..1'
    %
    %             continue
    %         end
    %
    %         if data1(1)==0
    %             data1=data1(2:end,1);
    %             timelags=timelags+1;
    %         end
    %     end
    
    outbreakx
    
    caddate1
    
    data1
    
    % <=============================================================================================>
    % <================ If time series is too short, do not proceed further ========================>
    % <=============================================================================================>
    
    if length(data1)<10
        length(data1)
        'small number of data points..'
        continue
    end
    
    
    if (method1==0 & dist1==2)  % estimate <factor1> which is the variance to mean ratio to scale the NB error structure
        
        binsize1=7;
        
        [ratios,~]=getMeanVarianceRatio(data1,binsize1,2);
        
        index1=find(ratios(:,1)>0);
        
        factor1=mean(ratios(index1,1));
        
        factor1
    end
    
    
    % <==============================================================================================>
    % <================ Save adjusted incidence data file for further analysis ======================>
    % <==============================================================================================>
    
    data1=[(0:1:length(data1)-1)' data1];
    
    %save(strcat('./output/',datafilename2),'data1','-ascii')
    
    data=data1;
    
    t_window=1:length(data(:,1));
    
    data1=data(t_window',:)
    
    plot(data1(:,1),data1(:,2),'ko')
    
    % <=======================================================================================>
    % <===================== Derive the best fitting sub-epidemic models ====================>
    % <=======================================================================================>
    
    [RMSESx,PS,npatches,onset_thr,~]=fittingModifiedLogisticFunctionPatchABC(datafilename2,data1,DT,t_window,M,flag1,numstartpoints);

    % <===========================================================================================>
    % <======= Derive uncertainty for the <topmodelsx> best fitting models and save results ======>
    % <===========================================================================================>
    
        
    for rank1=1:topmodelsx
        
        [Phatss,npatches,onset_thr,curves,bestfit,data1,P0,AICc_best,factor1,d]=fittingModifiedLogisticFunctionPatchMultiple(RMSESx,PS,data1,DT,t_window,M,flag1,numstartpoints,rank1);
        
        cadfilename1=strcat('./output/modifiedLogisticPatch-ensem-npatchesfixed-',num2str(npatches_fixed),'-onsetfixed-',num2str(onset_fixed),'-smoothing-',num2str(smoothfactor1),'-',datafilename2(1:end-4),'-flag1-',num2str(flag1(1)),'-flag1-',num2str(flag1(2)),'-method-',num2str(method1),'-dist-',num2str(dist1),'-calibrationperiod-',num2str(calibrationperiod1),'-rank-',num2str(rank1),'.mat');
   
        save(cadfilename1,'-mat')
        
        strcat('rank-',num2str(rank1),'-complete!')
    end
    
end

clear RMSESx PS P
clear Phatss curves bestfit data1 P0 AICc_best d

