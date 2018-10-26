classdef spikes
    
    properties
        
        folder
        timestamp
        nchannels
        waveforms
        peak
        peak1
        peak2
        trough
        width
        num_trials
        trial_id
        
    end % properties
    
    methods
        
        function obj = spikes(varargin)
            
            if nargin < 1
                obj.timestamp = [ ];
                obj.nchannels = [ ];
                obj.waveforms = [ ];
                obj.peak = [ ];
                obj.peak1 = [ ];
                obj.peak2 = [ ];
                obj.trough = [ ];
                obj.width = [ ];
                obj.folder = [ ];
                obj.trial_id = [];
                obj.num_trials = [];
                
            else
                session_info = varargin{1};
                
                
                n = varargin{2}
                tetrode = varargin{3}
                
                obj.nchannels = session_info.nchannels{n}(tetrode);
                
                
                data_folder = session_info.folder{n};
               
                
                switch session_info.filetype{n}
                    
                    case 'laser event file jakob'
                        
                        if obj.nchannels == 4 % tetrode case
                        
                            fname = [data_folder filesep 'Sc' int2str(tetrode) '.ntt'];
                        
                        else % stereotrode case
                            fname = [data_folder filesep 'Sc' int2str(tetrode) '.nst'];
                        end
                        
                        data = read_cheetah_data(fname);
                        
                        obj.timestamp = data.ts;
                        obj.waveforms = data.waveforms;
                        
                         obj.folder = [data_folder filesep 'Sc' int2str(tetrode) '.cut'];
                        
                         
                    case 'Neuralynx'
                        
                        if obj.nchannels == 4 % tetrode case
                        
                            fname = [data_folder filesep 'Sc' int2str(tetrode) '.ntt'];
                        
                        else % stereotrode case
                            fname = [data_folder filesep 'Sc' int2str(tetrode) '.nst'];
                        end
                        
                        data = read_cheetah_data(fname);
                        
                        obj.timestamp = data.ts;
                        obj.waveforms = data.waveforms;
                        
                         obj.folder = [data_folder filesep 'Sc' int2str(tetrode) '.cut'];
                        
                    case 'Neuralynx G'
                        
                        if obj.nchannels == 4 % tetrode case
                        
                            fname = [data_folder filesep 'TT' int2str(tetrode) '.ntt']
                        
                        else % stereotrode case
                            fname = [data_folder filesep 'TT' int2str(tetrode) '.nst'];
                        end
                        
                        data = read_cheetah_data(fname);
                        
                        obj.timestamp = data.ts;
                        obj.waveforms = data.waveforms;
                        
                         obj.folder = [data_folder filesep 'TT' int2str(tetrode) '.cut'];
                         
                         
                   case 'Neuralynx saturated'
                        
                        fname = [data_folder filesep 'Sc' int2str(tetrode) '.ntt']
                        
                        data = read_cheetah_data(fname);
                        
                        obj.timestamp = data.ts;
                        obj.waveforms = data.waveforms;
                        
                        obj.folder = [data_folder filesep 'Sc' int2str(tetrode) '.cut'];
                        
                        
                   case 'TRN Data'
                        
                         if obj.nchannels == 4 % tetrode case
                        
                            fname = [data_folder filesep 'Sc' int2str(tetrode) '.Ntt'];
                        
                        else % stereotrode case
                            fname = [data_folder filesep 'Sc' int2str(tetrode) '.Nst'];
                        end
                        
                        data = read_cheetah_data(fname);
                        
                        obj.timestamp = data.ts;
                        obj.waveforms = data.waveforms;
                        
                         obj.folder = [data_folder filesep 'Sc' int2str(tetrode) '.cut'];
                        
                        
                    case 'Digital Lynx'
                        
                        fname = [data_folder filesep 'TT' int2str(tetrode) '.ntt'];
                        
                        data = read_cheetah_data(fname);
                        
                        obj.timestamp = data.ts;
                        obj.waveforms = data.waveforms;
                        
                        obj.folder = [data_folder filesep 'TT' int2str(tetrode) '.cut'];
                        
                          case 'Mouse 1 Neuralynx'
                        
                        fname = [data_folder filesep 'TT' int2str(tetrode) '.ntt'];
                        
                        data = read_cheetah_data(fname);
                        
                        obj.timestamp = data.ts;
                        obj.waveforms = data.waveforms;
                        
                       obj.folder = [data_folder filesep 'TT' int2str(tetrode)];
                        
                    case 'MWL'
                        
                        switch tetrode
                            case 1
                                file = 'a107';
                            case 2
                                file = 'a207';
                            case 3
                                file = 'b107';
                            case 4
                                file = 'b207';
                            case 5
                                file = 'c107';
                            case 6
                                file = 'c207';
                        end
                        
                        
                        fname = [data_folder filesep file filesep file '.tt'];
                        
                        obj.folder = [obj.folder filesep file];
                        
                        f = mwlopen(fname);
                        
                        data = load(f);
                        
                        obj.timestamp = double(data.timestamp)'./10000;
                        obj.waveforms = double(data.waveform);
                        
                        
                end
                
                [mx,ind1] = max(obj.waveforms,[],2);
                [mn,ind2] = min(obj.waveforms,[],2);
                
                %size(obj.waveforms)
                
                [mx1] = max(obj.waveforms(:,1:16,:),[],2);
                [mx2] = max(obj.waveforms(:,17:32,:),[],2);
                
                obj.peak = squeeze(mx);
                obj.peak1 = squeeze(mx1);
                obj.peak2 = squeeze(mx2);
                obj.trough = squeeze(mn);
                obj.width = squeeze(ind2) - squeeze(ind1);
                obj.num_trials = 1;
                
            end
            
            
        end %
        
        function obj2 = crop(obj,t1,t2)
            
            % t1, t2 = time in seconds
            
            t = obj.timestamp;
            
            ok_times = find(t > t1 & t < t2);
            
            obj2 = spikes();
            
            obj2.timestamp = obj.timestamp(ok_times);
            obj2.waveforms = obj.waveforms(:,:,ok_times);
            obj2.peak = obj.peak(:,ok_times);
            obj2.peak1 = obj.peak1(:,ok_times);
            obj2.peak2 = obj.peak2(:,ok_times);
            obj2.trough = obj.trough(:,ok_times);
            obj2.width = obj.width(:,ok_times);
            obj2.num_trials = obj.num_trials;
            obj2.trial_id = obj.trial_id;
            
        end
        
        function obj3 = plus(obj1, obj2)
            
            obj3 = spikes();
            
            obj3.timestamp = cat(1,obj1.timestamp,obj2.timestamp);
            obj3.waveforms = cat(3, obj1.waveforms, obj2.waveforms);
            obj3.peak = cat(2, obj1.peak, obj2.peak);
            obj3.peak1 = cat(2, obj1.peak1, obj2.peak1);
            obj3.peak2 = cat(2, obj1.peak2, obj2.peak2);
            obj3.trough = cat(2, obj1.trough, obj2.trough);
            obj3.width = cat(2, obj1.width, obj2.width);
            
            obj3.trial_id = cat(1, obj1.trial_id, obj2.trial_id);
            
            obj3.num_trials = numel(unique(obj3.trial_id));
            
        end
        
        function obj2 = align(obj,laser_info,indices,pre_stim,post_stim)
            
            obj2 = spikes();
            t = obj.timestamp;
            
            for n = 1:length(indices)
                
                start_time = laser_info.start_time(indices(n));
                
                t1 = start_time - pre_stim;
                t2 = start_time + post_stim;
                
                ok_times = find(t > t1 & t < t2);
                
                trial_id = ones(1,length(ok_times)).*n;
                
                obj2.timestamp = cat(1,obj2.timestamp,obj.timestamp(ok_times)-start_time);
                obj2.waveforms = cat(3,obj2.waveforms,obj.waveforms(:,:,ok_times));
                obj2.peak = cat(2,obj2.peak,obj.peak(:,ok_times));
                obj2.peak1 = cat(2,obj2.peak1,obj.peak1(:,ok_times));
                obj2.peak2 = cat(2,obj2.peak2,obj.peak2(:,ok_times));
                obj2.trough = cat(2,obj2.trough,obj.trough(:,ok_times));
                obj2.width = cat(2,obj2.width,obj.width(:,ok_times));
                obj2.trial_id = [obj2.trial_id trial_id];
                
            end
            
            obj2.num_trials = n;
            
        end
        
        function obj2 = filter(obj,varargin)
            
            %
            % filter spikes for various parameters
            %
            
            obj2 = spikes();
            
            for n = 1:2:nargin-1
                
                %disp(n)
                
                param = varargin{n};
                value_range = varargin{1+n};
                
                if strcmp(param,'width')
                    ok_spikes = find(min(obj.(param)) > value_range(1) & min(obj.(param)) < value_range(2));
                elseif strcmp(param,'timestamp')
                    ok_spikes = find(obj.(param) > value_range(1) & obj.(param) < value_range(2));
                    
                    
                else
                    
                    
                    
                    ok_spikes = find(max(obj.(param)) > value_range(1) & max(obj.(param)) < value_range(2));
                    
                end
                
                obj2.timestamp = obj.timestamp(ok_spikes);
                obj2.waveforms = obj.waveforms(:,:,ok_spikes);
                obj2.peak = obj.peak(:,ok_spikes);
                obj2.peak1 = obj.peak1(:,ok_spikes);
                obj2.peak2 = obj.peak2(:,ok_spikes);
                obj2.trough = obj.trough(:,ok_spikes);
                obj2.width = obj.width(:,ok_spikes);
                obj2.num_trials = obj.num_trials;
                
            end
            
        end
        
        function obj2 = cluster(obj,cluster_n,varargin)
            
            obj2 = spikes();
            
            if nargin < 3
                cluster_mode = 'isolate';
            else
                cluster_mode = 'remove';
            end
            
           % if strcmp(obj.folder(end-3:end),'.cut')
                
                cluster_info = load([obj.folder '']);
                
                cl.id = find(cluster_info == cluster_n)-1;
                
           % else
                
           %     fname = [obj.folder filesep 'cl-' int2str(cluster_n)];
            
           %     cl = load(mwlopen(fname));
           % end
            
            switch cluster_mode
                case 'isolate'
                    obj2.timestamp = obj.timestamp(cl.id+1);
                    obj2.waveforms = obj.waveforms(:,:,cl.id+1);
                    obj2.peak = obj.peak(:,cl.id+1);
                    obj2.peak1 = obj.peak1(:,cl.id+1);
                    obj2.peak2 = obj.peak2(:,cl.id+1);
                    obj2.trough = obj.trough(:,cl.id+1);
                    obj2.width = obj.width(:,cl.id+1);
                    obj2.num_trials = obj.num_trials;
                case 'remove'
                    
                    obj2 = obj;
                    
                    obj2.timestamp(cl.id+1) = [ ];
                    obj2.waveforms(:,:,cl.id+1) = [ ];
                    obj2.peak(:,cl.id+1) = [ ];
                    obj2.peak1(:,cl.id+1) = [ ];
                    obj2.peak2(:,cl.id+1) = [ ];
                    obj2.trough(:,cl.id+1) = [ ];
                    obj2.width(:,cl.id+1) = [ ];
                    obj2.num_trials = obj.num_trials;
                    
            end
            
        end
        
        function output = waveshape(obj)
            
            output.peak = [mean(obj.peak,2) std(obj.peak,[],2)];
            output.trough = [mean(obj.trough,2) std(obj.trough,[],2)];
            output.pt_ratio = [mean(obj.peak./obj.trough,2) std(obj.peak./obj.trough,[],2)];
            output.width = [mean(obj.width,2) std(obj.width,[],2)]./32;
            
        end
        
        function plot(obj,varargin)
            
            if nargin < 2
                plottype = 'scatter';
                c = 'k';
                n_bins = 200;
                hist_type = 'count';
            elseif nargin < 3
                plottype = varargin{1};
                c = 'k';
                n_bins = 200;
                hist_type = 'count';
            elseif nargin < 4
                plottype = varargin{1};
                c = varargin{2};
                n_bins = 200;
                hist_type = 'count';
            elseif nargin < 5
                plottype = varargin{1};
                c = varargin{2};
                n_bins = varargin{3};
                hist_type = 'count';
            else
                plottype = varargin{1};
                c = varargin{2};
                n_bins = varargin{3};
                hist_type = varargin{4};
            end
            
            
            switch plottype
                
                case 'ISIhist'
                    hist(diff(obj.timestamp),100);
                    axis tight
                    
                case 'hist'
                    [h,x] = hist(obj.timestamp,n_bins);
                    
                    switch hist_type
                        case 'count'
                            b = bar(x,h,c);
                        case 'percent'
                            b = bar(x,h./obj.num_trials,c);
                    end
                    
                    
                    
                    set(b,'LineStyle','none');
                    set(b,'BarWidth',1.0);
                    axis tight
                    xlabel('time (s)')
                    switch hist_type
                        case 'count'
                            ylabel('spike count')
                        case 'percent'
                            ylabel('spike probability')
                    end
                    
                case 'raster'
                    baseline = 0;
                    for n = 1:length(obj.timestamp)
                        t = obj.timestamp(n);
                        %if n > 1 && t < obj.timestamp(n-1)
                         %   baseline = baseline + 1;
                        %end
                        trial = obj.trial_id(n);
                        plot([t t],[trial trial+1],'color',c)
                        hold on
                    end
                    
                    hold off
                    xlabel('time (s)')
                    ylabel('trial number')
                    
                case 'scatter'
                    R = [1 2;...
                        1 3;...
                        1 4;...
                        2 3;...
                        2 4;...
                        3 4];
                    for plt = 1:6
                        subplot(2,3,plt)
                        plot(obj.peak(R(plt,1),:),obj.peak(R(plt,2),:),'.','color',c,'MarkerSize',0.5);
                    end
                    
                case 'waveforms'
                    
                    for subplt = 1:4
                        
                        subplot(1,4,subplt)
                        
                        plot(squeeze(obj.waveforms(subplt,:,:)),'color',[0.8 0.8 0.8])
                        hold on
                        plot(mean(squeeze(obj.waveforms(subplt,:,:)),2),'color','k','LineWidth',3.0)
                        hold off
                    end
                    
            end
            
        end
        
        % alternate plotting functions
        function hist(obj)
            plot(obj,'hist')
        end
        
        function raster(obj)
            plot(obj,'raster')
        end
        
        
    end % methods
    
    
    
end % classdef




