function error = SEIRCDworlderror(AllCases,AllDeaths,p,P,M,N)
    p=abs(p);
    ic= zeros(6*N,1);
    for i=1:N 
        %ic(1)=p(1), ic(2)=0,.... ic(7)=p(2)...ic(13)=p(3)
        %ic(2)=1....ic(8)=1....ic(14)=1
        ic(6*(i-1)+1)=P(i);
        ic(6*(i-1)+2)=p(3);
        %ic(6*(i-1)+2)=1; % Suscepticle is population, Exposed is 1, rest of 6 variables are 0
    end
    %T=500;
    T=size(AllCases,2);
    inter = [0,T-1];
    n = 10*T; %multiply by 10, then skip 10 to get all t values (1:10:end-1) smaller time step, more accurate
    %p = [0.015, 0.1, 10]; 
    %n = 10*T;
    %inter = [0,T-1];
    %ic = [330e6 0 0 0 0]; % S I R C D leave out E b/c to we don't know whos Exposed 
    [x,t] = SEIRCDworld2(inter,ic,n,p,M);
    error = 0;
    for i=1:N
        C=x(:,6*(i-1)+5); %matlab index starts at 1 (time values,Cases)
        D=x(:,6*(i-1)+6); %collect all cases and deaths for all countries
        error = error + mean((C(1:10:end-1)-AllCases(i,:)').^2+5*(D(1:10:end-1)-AllDeaths(i,:)').^2); % Cases And DEaths real while C and D are experimental
        %(ith country, all time values), AllCases and AllDeaths are row
        %vectors so add ' to transpose, C and D are column vectors
    end
% compare this model to 114days normally   
end