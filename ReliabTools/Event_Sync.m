% Calculates event synchronization and event delay from two spike trains "x_spikes" and "y_spikes"
%
% Information can be found online under "http://www.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/Event-Sync.html" and/or in
% Quian Quiroga R, Kreuz T, Grassberger P:
% Event Synchronization: A simple and fast method to measure synchronicity and time delay patterns.
% Phys.Rev. E, 66, 041904 (2002).
% 
% Example call:
%
% x=[1:100];                   % first spike train (completely periodic)
% y=x+(rand(1,length(x))-0.5); % second spike train (first spike train with some random jitter) 
% [es,ed]=event_sync(x,y)      % es and ed give you the event synchronization and the event delay between spike train x and spike train y
% 
%
% For questions and/or comments please contact me at "tkreuz (at) ucsd.edu".  

function [es,ed]=Event_Sync(x_spikes,y_spikes)

plotmode=0;             % plotting figure?                      (0-no,1-yes)
numfig=2;               % number of this figure

printmode=0;            % saving to postscript file?            (0-no,1-yes)
filename='Event-Sync-Example2.ps';  % name of this postscript file

xy_name={'X','Y'};      % name of spike trains
xy_col='br';            % colors for the respective spike trains

precision = 6;            % - logarithm_10 ( sampling interval )

taumode  = 2;             % how is the allowed time lag defined?            (1-fixed interval (tau),2-adaptive)
timemode = 2;            % cumulative, windowed or cumulative normalized?  (1-cumulative,2-windowed,3:cumulative normalized)

tau = 2;                  % maximum time delay for which events are still considered to be synchronous (for taumode=1)
window = 10;                % window size (for timemode=2)

if size(x_spikes,1)>size(x_spikes,2) x_spikes=x_spikes'; end
if size(y_spikes,1)>size(y_spikes,2) y_spikes=y_spikes'; end

x_spikes=sort(round(x_spikes*exp(precision*log(10)))/exp(precision*log(10)));
y_spikes=sort(round(y_spikes*exp(precision*log(10)))/exp(precision*log(10)));

dt=10^-precision;

tedge=[0.98*min([x_spikes y_spikes]) max([x_spikes y_spikes])*1.02];
tlim=[fix(tedge(1)/dt)*dt ceil(tedge(2)/dt)*dt];
time=[tlim(1):dt:tlim(2)];
len=length(time);

events=zeros(2,len);
events(1,round((x_spikes-tlim(1))/dt))=1;
events(2,round((y_spikes-tlim(1))/dt))=1;
eventsynchro=zeros(1,len);
if taumode==1                                       % maximum time lag fixed
    [r,c]=find(events==1);
    for efc=1:length(r)-1
        efc2=efc+find(c(efc+1:end)==c(efc));
        for efc3=1:length(efc2)
            eventsynchro(c(efc))=0.5;
        end   
        efc2=efc+find(c(efc+1:end)<=c(efc)+tau & c(efc+1:end)>c(efc) & r(efc+1:end)~=r(efc));
        for efc3=1:length(efc2)
            eventsynchro(c(efc2(efc3)))=sign(r(efc2(efc3))-r(efc));
        end
    end
else                                                % adaptive
    [r,c]=find(events==1);
    for efc=1:length(r)-1
        efc2=efc+find(c(efc+1:end)==c(efc));
        for efc3=1:length(efc2)
            eventsynchro(c(efc))=0.5;
        end
        xdelay=abs(c(efc)-c(find(r==r(efc))));
        xmin=min(xdelay(xdelay>0))/2;
        if isempty(xmin)
            xmin=len;
        end
        efc2=efc+find(c(efc+1:end)<=c(efc)+xmin & c(efc+1:end)>c(efc) & r(efc+1:end)~=r(efc));
        for efc3=1:length(efc2)
            ydelay=abs(c(efc2(efc3))-c(find(r==r(efc2(efc3)))));
            ymin=min(ydelay(ydelay>0))/2;
            if isempty(ymin)
                ymin=len;
            end
            tau=min([xmin ymin]);
            if ismember(c(efc2(efc3)),c(efc)+1:min([c(efc)+tau,size(events,2)]))>0
                eventsynchro(c(efc2(efc3)))=sign(r(efc2(efc3))-r(efc));
            end
        end
    end
end

cumevents=cumsum(events~=0,2);   
cumevents_x=cumevents(1,:);
cumevents_y=cumevents(2,:);
windowsupport=window+1:len-window;
windowlen=length(windowsupport);
winevents_x=cumevents_x(windowsupport+window)-cumevents_x(windowsupport-window);
winevents_y=cumevents_y(windowsupport+window)-cumevents_y(windowsupport-window);

cumeventsynchro=cumsum(abs(ceil(eventsynchro)),2);
cumeventdelay=cumsum(fix(eventsynchro),2);

normcumeventsynchro=cumeventsynchro./(sqrt(cumevents_x.*cumevents_y)+((sqrt(cumevents_x.*cumevents_y))==0));
normcumeventdelay=cumeventdelay./(sqrt(cumevents_x.*cumevents_y)+((sqrt(cumevents_x.*cumevents_y))==0));

wineventsynchro=(cumeventsynchro(windowsupport+window)-cumeventsynchro(windowsupport-window))./(sqrt(winevents_x.*winevents_y)+...
    ((sqrt(winevents_x.*winevents_y))==0));
wineventdelay=(cumeventdelay(windowsupport+window)-cumeventdelay(windowsupport-window))./(sqrt(winevents_x.*winevents_y)...
    +((sqrt(winevents_x.*winevents_y))==0));
wineventsynchro(wineventsynchro>1)=1;
wineventdelay(wineventdelay>1)=1;
wineventdelay(wineventdelay<-1)=-1;

xevc=cumevents(1,end);
yevc=cumevents(2,end);  
delayc=cumeventdelay(end);
synchroc=cumeventsynchro(end);

es=(synchroc./(sqrt(xevc.*yevc)+((sqrt(xevc.*yevc))==0)))';
ed=(delayc./(sqrt(xevc.*yevc)+((sqrt(xevc.*yevc))==0)))';

if plotmode==1
    figure(numfig); clf; hold on
    set(gcf,'Position',[1 31 1024 662])
    set(gcf,'Name',['Event-Sync: ',xy_name{1},'  -  ',xy_name{2}])
    
    % here the data could be plotted (x_data and y_data have to be assigned first)
    % plot(time,1.32+(x_data-min(x_data))./(max(x_data)-min(x_data))*0.26,'k');
    % plot(time,1.02+(y_data-min(y_data))./(max(y_data)-min(y_data))*0.26,'k');

    ylim([-0.05 2.65])
    
    for sc=1:length(x_spikes)
        line(x_spikes(sc)*ones(1,2),[1.33 1.57],'Color',xy_col(1));
    end
    for sc=1:length(y_spikes)
        line(y_spikes(sc)*ones(1,2),[1.03 1.27],'Color',xy_col(2));
    end
    
    if timemode==1
        maxval=max(cumeventsynchro);
        maxval2=max(abs(cumeventdelay));
        plot(time,1.6+cumeventsynchro/(maxval+(maxval==0)),'k')
        plot(time,0.5+cumeventdelay/(maxval2+(maxval2==0))/2,'k')
    elseif timemode==2
        maxval=1; maxval2=1;
        plot(time(windowsupport),1.6+wineventsynchro,'k')
        plot(time(windowsupport),0.5+wineventdelay/2,'k')
    else
        maxval=1; maxval2=1;
        plot(time,1.6+normcumeventsynchro,'k')
        plot(time,0.5+normcumeventdelay/2,'k')
    end
    
    xlim(tlim)
    xl=xlim; yl=ylim;
    line(xl(1):xl(2)-xl(1):xl(2),2.6*ones(1,2),'Color','k','LineStyle',':')
    line(xl(1):xl(2)-xl(1):xl(2),1.6*ones(1,2),'Color','k','LineStyle',':')
    line(xl(1):xl(2)-xl(1):xl(2),1.3*ones(1,2),'Color','k','LineStyle',':')
    line(xl(1):xl(2)-xl(1):xl(2),ones(1,2),'Color','k','LineStyle',':')
    line(xl(1):xl(2)-xl(1):xl(2),0.5*ones(1,2),'Color','k','LineStyle',':')
    line(xl(1):xl(2)-xl(1):xl(2),zeros(1,2),'Color','k','LineStyle',':')
    Ylabelstr{1}=num2str(-maxval2,3); Ylabelstr{2}=num2str(0); Ylabelstr{3}=num2str(maxval2,3);
    Ylabelstr{4}=num2str(0); Ylabelstr{5}=num2str(maxval,3);
    set(gca,'YTick',[0 0.5 1 1.6 2.6],'YTickLabel',Ylabelstr,'FontSize',11);
    
    text(xl(1)-0.06*(xl(2)-xl(1)),1.15,[xy_name{2}],'Color','k','FontSize',12,'FontWeight','bold')
    text(xl(1)-0.06*(xl(2)-xl(1)),1.45,[xy_name{1}],'Color','k','FontSize',12,'FontWeight','bold')
    
    text(xl(1)-0.07*(xl(2)-xl(1)),yl(1)+0.9*(yl(2)-yl(1)),'Sync','Color','k','FontSize',12,'FontWeight','bold')
    text(xl(1)-0.08*(xl(2)-xl(1)),yl(1)+0.1*(yl(2)-yl(1)),'Delay','Color','k','FontSize',12,'FontWeight','bold')
    title(['Event Synchronization: ',num2str(es),'   ---    Event-Delay: ',num2str(ed)],'Color','k','FontSize',14,'FontWeight','bold')
    xlabel('Time','Color','k','FontSize',11,'FontWeight','bold')
    set(gcf,'Color','w'); set(gca,'Color','w');
    box off
    if printmode==1                                                                     % Create postscript file
        set(gcf,'PaperOrientation','Landscape'); set(gcf,'PaperType', 'A4');
        set(gcf,'PaperUnits','Normalized','PaperPosition', [0 0 1.0 1.0]);
        print(gcf,'-dpsc',filename);
    end
end