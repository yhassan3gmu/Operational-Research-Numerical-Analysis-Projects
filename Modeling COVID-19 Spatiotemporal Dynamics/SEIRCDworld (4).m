% Program 6.3 Animation program for pendulum
% Inputs: time interval inter,
% initial values ic = [y(1,1) y(1,2)], number of steps n
% Calls a one-step method such as trapstep.m
% Example usage: pend2([0 10],[pi/2 0],200)
function [x, t] = SEIRCDworld(inter,ic,n,p)
    %ic 6*N
    %%% Setup the initial conditions
    
    h=(inter(2)-inter(1))/n;
    t = zeros(n,1);
    x = zeros(n,length(ic));
    x(1,:)=ic; 
    t(1)=inter(1);
    E0 = p(3);
    x(1,2) = E0;
    x(1,3) = ic(2);
    x(1,4) = ic(3);
    x(1,5) = ic(4);
    x(1,6) = ic(5);
    
    %%% Solve the differential equation
    for k=1:n
        t(k+1)=t(k)+h;  
        x(k+1,:)=trapstep(t(k),x(k,:),h,p);
    end
    
    %%% Plot if no outputs are asked for
    if (nargout==0)
        figure;
        plot(t,x,'linewidth',2);
        l=legend('S','E','I','R','C','D');
        set(l,'fontsize',18);
        set(gca,'fontsize',18);
        xlabel('Time','fontsize',22);
        xlim(inter);
    end
end

function y = trapstep(t,x,h,p)
    %one step of the Trapezoid Method
    z1=ydot(t,x,p);
    g=x+h*z1;
    z2=ydot(t+h,g,p);
    y=x+h*(z1+z2)/2;
end

function y=eulerstep(t,y,h)
    %one step of the Euler Method
    %Input: current time t, current vector y, step size h
    %Output: the approximate solution vector at time t+h
    y=y+h*ydot(t,y);
end

function xdot=ydot(t,x,p)
    %g=9.81;length=1;
    %step 1
    %incubation=3;
    %R0=3
    N =length(x)/6;
    incubation=5;
    mu=p(1);tau=p(2);E0=p(3);
    sigma=1/incubation;gamma=1/4;R0=2.2;B=R0*gamma;
    for i=1:N
        S=x(6*(i-1)+1);E=x(6*(i-1)+2);I=x(6*(i-1)+3);R=x(6*(i-1)+4);C=x(6*(i-1)+5);D=x(6*(i-1)+6);
        xdot(6*(i-1)+1) = -(B*S*I)/(S+E+I+R);
        xdot(6*(i-1)+2) = ((B*S*I)/(S+E+I+R))-(sigma*E);
        xdot(6*(i-1)+3) = (sigma*E)-(gamma*I)-(gamma*tau*I)-(gamma*mu*I);
        xdot(6*(i-1)+4) = (gamma*I);
        xdot(6*(i-1)+5) = (gamma*tau*I);
        xdot(6*(i-1)+6) = (gamma*mu*I);
    end
end
