loadfulldata;
ic= zeros(6*N,1);
for i=1:N 
    %ic(1)=p(1), ic(2)=0,.... ic(7)=p(2)...ic(13)=p(3)
    %ic(2)=1....ic(8)=1....ic(14)=1
    ic(6*(i-1)+1)=P(i);% Suscepticle is population,
    ic(6*(i-1)+2)=1; % Exposed is 1, resot of 6 variables are 0
end
T=500;
inter = [0,T-1];
n = 10*T;
p = [0.015, 0.1, 10]; 
[x,t] = SEIRCDworld(inter,ic,n,p);

%problem 3 start

C = zeros(N,N);
M = zeros(N,N);
for i=1:N
    for j=1:N
        if Continent(i) ~= Continent(j)
            C(i,j) = 1;
        end 
    end
end
theta1 = 1/sum(P)^2;
theta2 = 5;
for i=1:N
    for j=1:N
        M(i,j) = (theta1)*(P(i)*P(j))/(1+theta2*C(i,j));
    end
end

%problem 3 end 

SEIRCDworld2(inter,ic,n,p,M); %problem 4 graph

%problem 5 start

f = @(p) SEIRCDworlderror(AllCases,AllDeaths,p,P,M,N);
options = optimset('MaxIter',20)
p0 = fminsearch(f,p,options)
p0 = fminsearch(f,p0,options)
p0 = fminsearch(f,p0,options)
p0 = fminsearch(f,p0,options)
p0=abs(p0)
SEIRCDworld2(inter,ic,n,p0,M);

%problem 5 end