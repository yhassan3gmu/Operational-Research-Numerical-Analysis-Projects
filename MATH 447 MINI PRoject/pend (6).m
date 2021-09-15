% Program 6.3 Animation program for pendulum
% Inputs: time interval inter,
% initial values ic = [y(1,1) y(1,2)], number of steps n
% Calls a one-step method such as trapstep.m
% Example usage: pend([0 10],[pi/2 0],200)
function pend(inter,ic,n)
    h=(inter(2)-inter(1))/n; % plot n points in total
    y(1,:)=ic; % enter initial conds in y
    y2(1,:)=ic;
    t(1)=inter(1);
    fig=figure(1);clf(fig);
    set(gca,'xlim',[-1.2 1.2],'ylim',[-1.2 1.2], ...
    'XTick',[-1 0 1],'YTick',[-1 0 1])
    bob=animatedline('color','r','Marker','.','markersize',40);
    rod=animatedline('color','b','LineStyle','-','LineWidth',3);
    bob2=animatedline('color','g','Marker','.','markersize',40);
    rod2=animatedline('color','b','LineStyle','-','LineWidth',3);
    axis square % make aspect ratio 1 - 1
    for k=1:n
        tic;
        t(k+1)=t(k)+h;
        y(k+1,:)=eulerstep(t(k),y(k,:),h);
        xbob = sin(y(k+1,1)); ybob = -cos(y(k+1,1));
        xrod = [0 xbob]; yrod = [0 ybob];
        clearpoints(bob);addpoints(bob,xbob,ybob);
        clearpoints(rod);addpoints(rod,xrod,yrod);
        
        y2(k+1,:)=trapstep(t(k),y2(k,:),h);
        xbob = sin(y2(k+1,1)); ybob = -cos(y2(k+1,1));
        xrod = [0 xbob]; yrod = [0 ybob];
        clearpoints(bob2);addpoints(bob2,xbob,ybob);
        clearpoints(rod2);addpoints(rod2,xrod,yrod);
        drawnow;
        z=toc;
        pause(max(0,h-z))
    end
    fig=figure(2);clf(fig);
    plot(t,y2(:,1),'b','linewidth',2);
    hold on;
    plot(t,y2(:,2),'r','linewidth',2);
    plot(t,y(:,1),'b--','linewidth',1.5)
    plot(t,y(:,2),'r--','linewidth',1.5)
    l=legend('Trapezoid, \theta','Trapezoid, d\theta/dt','Euler, \theta','Euler, d\theta/dt');
    set(l,'fontsize',18);
    set(gca,'fontsize',18);
    xlabel('Time','fontsize',22);
    xlim(inter);
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
    g=9.81;length=1;d=1;A=15;
    z(1) = y(2);
    z(2) = -(g/length)*sin(y(1))-d*y(2)-A*sin(t);
end
% function z=ydot(t,y)
%     g=9.81;length=1;
%     z(1) = y(2);
%     z(2) = -(g/length)*sin(y(1));
% end
