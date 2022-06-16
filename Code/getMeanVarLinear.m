function [coefs,ns]=getMeanVarLinear(data1,fitcurve_model1d,bins)


residuals1=(data1-fitcurve_model1d);

incidence1=fitcurve_model1d;

predictedranges=linspace(min(incidence1),max(incidence1)+1,bins)';

variances1=[];

ns=[];

for j=1:length(predictedranges(:,1))-1
    
    index1=find(incidence1>=predictedranges(j,1) & incidence1<predictedranges(j+1,1));
    
    
    variances1=[variances1;[mean(incidence1(index1)) var(residuals1(index1))]];

    %variances1=[variances1;[(predictedranges(j,1)+predictedranges(j+1,1))/2 var(residuals1(index1))]];

    ns=[ns;length(index1)];
end

index1=find(ns>3);

length(variances1)

variances1=variances1(index1,:);

length(variances1)


figure(102)

subplot(1,2,1)
line1=plot(incidence1,residuals1,'rx')
set(line1,'MarkerSize',9,'LineWidth',2)
hold on

[rho,t,P]=spear(incidence1,residuals1)

line1=[min(incidence1) 0;max(incidence1)+1 0];
line2=plot(line1(:,1),line1(:,2),'k--')
set(line2,'LineWidth',2)


for j=1:length(predictedranges(:,1))
    
    line1=[predictedranges(j,1) min(residuals1);predictedranges(j,1) max(residuals1)];
    
    line1=plot(line1(:,1),line1(:,2),'k--')
    set(line1,'LineWidth',2,'MarkerSize',9)

    
    hold on
    
end

    
xlabel('Predicted values');
ylabel('Residuals')

set(gca,'FontSize', 24);
set(gcf,'color','white')

subplot(1,2,2)
line1=plot(variances1(:,1),variances1(:,2),'bs')
set(line1,'LineWidth',2,'MarkerSize',9)
hold on

[rho,t,P]=spear(variances1(:,1),variances1(:,2))


X=[ones(length(variances1(:,1)),1) variances1(:,1)];

%X=[variances1(:,1)];


[B,BINT,R,RINT,STATS]=regress(variances1(:,2),X)

line1=plot(X(:,2),X*B,'k-');
set(line1,'LineWidth',2)
hold on

%vector STATS containing, in
%the following order, the R-square statistic, the F statistic and p value
%for the full model, and an estimate of the error variance.

xlabel('Mean')
ylabel('Variance')

title(strcat('R^2=',num2str(STATS(1)*100),'; P=',num2str(STATS(3))))

coefs=B(2);
BINT

set(gca,'FontSize', 24);
set(gcf,'color','white')
