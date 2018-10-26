classdef sessions

    properties

        ID
        name
        time
        folder
        filetype
        down_sample_rate
        description
        nchannels
        tetrode_depths
        

    end % properties

    methods

        function obj = sessions(varargin)

            % creates a sessions object
            
        end

        function obj = add(obj,varargin)

            %
            % varargin = mouse name, experiment time,
            %            folder path, data type, and description
            %

            n = length(obj.ID) + 1;

            obj.ID(n) = sum(datevec(now));

            obj.name{n} = varargin{1};
            obj.time(n,:) = varargin{2};
            obj.folder{n} = varargin{3};
            obj.filetype{n} = varargin{4};
            obj.down_sample_rate(n) = varargin{5};
            obj.description{n} = varargin{6};
            obj.nchannels{n} = varargin{7};
            obj.tetrode_depths{n} = varargin{8};

        end

        function obj = remove(obj,session_numbers)

            %
            % n = vector of sessions to remove
            %

            for m = 1:length(session_numbers);
                
                n = session_numbers(m) - (m-1);

                obj.ID(n) = [ ];

                obj.name(n) = [ ];
                obj.time(n,:) = [ ];
                obj.folder(n) = [ ];
                obj.filetype(n) = [ ];
                obj.down_sample_rate(n) = [ ];
                obj.description(n) = [ ];
                obj.nchannels(n) = [ ]
                obj.tetrode_depths(n) = [ ];

            end

        end
        
        function view(obj,varargin)
            
            if nargin == 1
                session_numbers = 1:numel(obj.ID);
            else
                session_numbers = varargin{1};
            end
            
            for n = session_numbers
                fprintf(['\t\n Session number: ' int2str(n)])
                fprintf(['\t\n Mouse name: ' obj.name{n}])
                fprintf(['\t\n Date and time: ' datestr(obj.time(n,:))])
                fprintf(['\t\n Description: ' obj.description{n} '\n\n'])
            end
                
        end
        
        function out = filter(obj,varargin)
            
            prop = varargin{1};
            name = varargin{2};
            
            out = [ ];
            
            for n = 1:numel(obj.ID)
                
                val = obj.(prop){n};
                
                if strcmp(val,name)
                    
                    out = [out n];
                    
                end
            end
        end
        
        function save(obj)
            
            save '/home/jsiegle/MATLAB/data/sessions.mat' obj
            
        end

    end % methods



end % classdef
