function [result, smoothTrace] = saccadeDetector(time,trace)
% [result, smoothTrace] = saccadeDetector(time,trace)

%% parameters to tweak
pos_filter_length=5; % short filter to smooth for noisy eye traces
trace=filtfilt(ones(pos_filter_length,1)/pos_filter_length, 1, trace')';
filter_length=20; %number of amples to average for baseline velocity
detect_thresh=100; % threshold in deg/s to detect a saccade
start_thresh=5; % threshold in deg/s to determine the start and end of a saccade
min_isi=50; %minimum number of samples between any two saccades
minDur=5; %minimum duration of a saccade
blink_isi=50; %ignore saccades this many samples around a blink
reproject_onto_saccade_direction=true; %work on signed velocity instead of speed

%% calculate veocities
tps=median(diff(time));
vel=[[0;0],(trace(:,3:end)-trace(:,1:end-2))/(tps*2),[0;0]];
speed=sqrt(sum(vel.^2));

a = 1;
b = ones(1,filter_length/2)/filter_length/2;
speedf = filter(b,a,speed);


speedspeedf=speed-speedf;

potential_saccades=diff([0 (speedspeedf)>detect_thresh])==1;%check data start condition, consider loader more data before...
potential_saccades=find(potential_saccades);
if numel(potential_saccades)>1
    potential_saccades([min_isi diff(potential_saccades)]<min_isi)=[];
end


for iSaccade=1:length(potential_saccades)
    
    if iSaccade>1 && potential_saccades(iSaccade)<endindex(iSaccade-1)
        startindex(iSaccade) = potential_saccades(iSaccade);
        endindex(iSaccade) = potential_saccades(iSaccade);
        baseline_speed(iSaccade) = speedf(startindex(iSaccade));
        continue;
    end
    
    %something like this:
    %1. goto beginning of speed-speedf segment, either in a loop, or using find
    startindex(iSaccade)=potential_saccades(iSaccade);
    while speedspeedf(startindex(iSaccade))> start_thresh
        startindex(iSaccade) = startindex(iSaccade)-1;
    end
    %2. determine pursuit velocity
    baseline_speed(iSaccade) = speedf(startindex(iSaccade));
    endindex(iSaccade)=potential_saccades(iSaccade);
    %would probabbly be better to project onto motion direction and work on
    %signed speedocity. if this becomes necessary simply use this as the
    %starting point, extract samples between these, then project onto
    %sifference in this short bit only.
    while speed(endindex(iSaccade)) - baseline_speed(iSaccade) > start_thresh
        endindex(iSaccade) = endindex(iSaccade)+1;
    end
    
%     rem=endindex+blink_isi > size(speed,2);
%     endindex(rem)=[];
%     startindex(rem)=[];
    if endindex(iSaccade)+blink_isi > size(speed,2) || isnan(sum(speed(endindex(iSaccade)+(1:blink_isi)))) || startindex(iSaccade) <= blink_isi || isnan(sum(speed(startindex(iSaccade)-(1:blink_isi))))
        startindex(iSaccade) = potential_saccades(iSaccade);
        endindex(iSaccade) = potential_saccades(iSaccade);
        continue;
    end
    
	if reproject_onto_saccade_direction
    %3.now with these estimators of the actual saccade, we will rebase the
    %velocity along the saccade trajectory to allow using a velocity
    %instead of speed
        saccade_vector=trace(:,endindex(iSaccade))-trace(:,startindex(iSaccade));
        expected_basis = saccade_vector/sqrt(sum(saccade_vector.^2));
        vxy=expected_basis'*vel(:,(startindex(iSaccade)-filter_length):(endindex(iSaccade)+filter_length));
        baseline_speed(iSaccade) = sum(vxy(1:filter_length))/filter_length;
        baseline_endspeed(iSaccade)= sum(vxy(end-filter_length+1:end))/filter_length;

        old_startindex=startindex(iSaccade);
        startindex(iSaccade)=potential_saccades(iSaccade);

        while vxy(1,startindex(iSaccade)-old_startindex+filter_length+1) - baseline_speed(iSaccade)> start_thresh
            startindex(iSaccade) = startindex(iSaccade)-1;
        end

        endindex(iSaccade)=potential_saccades(iSaccade);

        while vxy(1,endindex(iSaccade)-old_startindex+filter_length+1) - baseline_endspeed(iSaccade)> start_thresh
            endindex(iSaccade) = endindex(iSaccade)+1;
        end
	end
     
end

sacdur=endindex-startindex;
sacsize=sqrt(sum((trace(:,endindex)-trace(:,startindex)).^2));
%do a hist of saccade durarion, hopefullt a clear cutoff
%  hist(sacdur, 0:350)
%  hist(sacsize,0:60)
 
 startindex(sacdur<minDur) = [];
 endindex(sacdur<minDur) = [];
 baseline_speed(sacdur<minDur) = [];
 potential_saccades(sacdur<minDur) = [];
 
 sacdur=endindex-startindex;
 sacsize=sqrt(sum((trace(:,endindex)-trace(:,startindex)).^2));

%  toc
%  [a,b]=max(sacdur);

% figure;
% for iSaccade= 1:length(potential_saccades)
%     this_trace=trace(:,startindex(iSaccade)+(-100 : 2000));
%     this_speed=speed(startindex(iSaccade)+(-100 : 2000));
%     this_speedspeedf=speedspeedf(startindex(iSaccade)+(-100 : 2000));
%     this_velocity = vel(:,startindex(iSaccade)+(-100 : 2000));
% %     this_pupil = pupil(startindex(iSaccade)+(-100 : 2000));
%     subplot(2,1,1)
%     hold off;
%     plot(-100:2000,this_trace');
%     hold on;
%     plot([0 0], [min(min(this_trace)) max(max(this_trace))], '-k');
%     plot([endindex(iSaccade) endindex(iSaccade)]-startindex(iSaccade), [min(min(this_trace)) max(max(this_trace))], '--k');
%     subplot(2,1,2)
%     hold off;
%     plot(-100:2000,this_speed);
%     hold all;
% %     plot(-100:2000,this_speedspeedf);
%     plot(-100:2000,this_velocity');
%     plot([0 0], [min(this_speed) max(this_speed)], '-k');
%     plot([endindex(iSaccade) endindex(iSaccade)]-startindex(iSaccade), [min(this_speed) max(this_speed)], '--k');
% %     subplot(3,1,3)
% %     hold off;
% %     plot(-100:2000,this_pupil);
%     
%     waitforbuttonpress
% end

result=[time(startindex); time(endindex); sacdur; sacsize; trace(:,startindex); trace(:,endindex); startindex; endindex];
if nargout>1
    smoothTrace=trace;
end
