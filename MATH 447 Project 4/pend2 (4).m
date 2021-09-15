% Program 6.3 Animation program for pendulum
% Inputs: time interval inter,
% initial values ic = [y(1,1) y(1,2)], number of steps n
% Calls a one-step method such as trapstep.m
% Example usage: pend2([0 10],[pi/2 0],200)
function [x, t] = pend2(inter,ic,n)

    %%% Setup the initial conditions
    h=(inter(2)-inter(1))/n;
    t = zeros(n,1);
    x = zeros(n,length(ic));
    x(1,:)=ic; 
    t(1)=inter(1);
    
    %%% Solve the differential equation
    for k=1:n
        t(k+1)=t(k)+h;  
        x(k+1,:)=trapstep(t(k),x(k,:),h);
    end
    
    %%% Plot if no outputs are asked for
    if (nargout==0)
        figure;
        plot(t,x,'linewidth',2);
        l=legend('\theta','d\theta/dt');
        set(l,'fontsize',18);
        set(gca,'fontsize',18);
        xlabel('Time','fontsize',22);
        xlim(inter);
    end
end

function y = trapstep(t,x,h)
    %one step of the Trapezoid Method
    z1=ydot(t,x);
    g=x+h*z1;
    z2=ydot(t+h,g);
    y=x+h*(z1+z2)/2;
end

function y=eulerstep(t,y,h)
    %one step of the Euler Method
    %Input: current time t, current vector y, step size h
    %Output: the approximate solution vector at time t+h
    y=y+h*ydot(t,y);
end

function z=ydot(t,y)
    g=9.81;length=1;
    z(1) = y(2);
    z(2) = -(g/length)*sin(y(1));
end