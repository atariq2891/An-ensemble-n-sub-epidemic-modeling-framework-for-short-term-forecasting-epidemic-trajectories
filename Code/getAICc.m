function [AICc,part1,part2,numparams]=getAICc(method1,npatches,flag1,fixI0,fval,n)

% last revised: 01 June 2022

%flag 0=GGM 1=GLM 2=GRM 3=Logistic 4=RIC 5=Gompertz

numparams=0;

for j=1:npatches
    
    switch flag1(j) % model indicator
        
        case 0
            numparams=numparams+2;
            
        case 1
            numparams=numparams+3;
            
        case 2
            numparams=numparams+4;
            
        case 3
            numparams=numparams+2;
            
        case 4
            numparams=numparams+3;
            
        case 5
            numparams=numparams+2;
            
    end
end


if npatches>1
    
    numparams=numparams+1; %Cthr
    
end

if fixI0==0 %fix initial datum or estimated
    
    numparams=numparams+1;
    
end

    
if method1==3 | method1==4 %Neg. Binomial requires one more parameter (alpha)
    
    numparams=numparams+1;
    
elseif method1==5
    
    numparams=numparams+2;  %Neg. Binomial requires 2 more parameters (alpha,d)
    
end


switch method1
    
    case 0
        
        AICc= n*log(fval) + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1= n*log(fval);
        
        part2= 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
    otherwise
        
        AICc=-2*(-fval) + 2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);
        
        part1=-2*(-fval);
        
        part2=2*numparams + (2*numparams*(numparams+1))/(n-numparams-1);

end


%method1

%fval

%n

%numparams
