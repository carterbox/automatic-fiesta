classdef statsdata2
%STATSDATA is an object for organizing statistical variables
% This class stores the data in a 3D cell where:
% column 1 : header
% column 2 : data
% column 3 : units

% EXAMPLE:
% s.variables = {'size', [100x1 double], 'm'};
% s.variables = {'color', [100x3 double], 'RBG'};

% METHODS:
% addData(s,header,data,unit)
%% -----------------------------------------------------------------------
    properties
        variables = cell(20,3);
        vcount = 0; % number of categories
        ocount = 0; % number of observations
        other = {};
    end
    
    methods
        function s = reduce(s, subset)
            % Reduce the total number of observations
            for i = 1:s.vcount
                s.variables{i,2} = s.variables{i,2}(subset);
            end
            s.ocount = numel(s.variables{1,2});

        end
        function s = convert(s, oldformat)
           % Convert old format of stats data to new format of statsdata.
           
           s.variables = [oldformat.predictors; oldformat.responses];
           s.vcount = oldformat.pcount + oldformat.rcount;
           s.ocount = oldformat.ocount;
        end

        function s = addData(s,header,data,unit)
            % Add data to the struct.
            
            % See how big the matrix has to be.
            if s.ocount > 0 && s.ocount ~= numel(data) 
                error('New data is not the same length!');    
            else
                s.ocount = numel(data);
            end
            
            % Check if already exists.
            index = s.findData(header);
            if index < 1
               s.vcount = s.vcount + 1;
               index = s.vcount;
            else
                warning('%s already exists and will be replaced.', header);
            end
            
            % Add the new data to the struct
            fprintf(1, 'ADD DATA: %s, %s\n', header, unit);
            s.variables{index,1} = header;
            s.variables{index,2} = data;
            s.variables{index,3} = unit;   
        end
        
        function [r,headers,units] = getData(s, names)
            % Return variables by index or by string header. 0 returns all
            % of the available options.
            % Returns one lxm matrix of the responsess in the vector n.
            
            n = s.findData(names);
            
            m = numel(n); % number of requested variables
            
            % Check for valid indices
            for i = 1:m
                if n(i) < 1 || n(i) > s.vcount
                    disp(names);
                    error('Requested variable number %i does not exist!', i);
                end
            end
            
            r = zeros(s.ocount,m);
            headers = cell(1,m);
            units = cell(1,m);
            for i = 1:m
                headers{i} = s.variables{n(i),1};
                r(:,i) = s.variables{n(i),2};
                units{i} = s.variables{n(i),3};
            end
        end
        
        function [name] = findData(s, name)
            % Return the index of the variable given by name.
            if isa(name, 'numeric')
                if numel(name) == 0 && name == 0
                    name = 1:s.vcount;
                end
                return
            end
            
            if isa(name,'char'), name = {name}; end
            
            m = numel(name);
            new_n = zeros(1,m);
            for i = 1:m
                for j = 1:s.vcount
                    if(strcmp(s.variables{j,1},name{i}))
                        new_n(i) = j;
                        break;
                    end
                end
                if new_n(i) == 0
                    new_n(i) = -1;
                end
            end
            name = new_n;
        end
            
        function s = scale(s,names,factor,newunit)
            % Multiplies a data by a factor
            
            n = s.findData(names);
            
            for i = n
                if i > 0
                    fprintf(1, 'SCALE: %s, %s x %f = %s\n',...
                            s.variables{i,1}, s.variables{i,3},...
                            factor, newunit);
                    s.variables{i,2} = s.variables{i,2}*factor;
                    s.variables{i,3} = newunit;
                else
                    warning('Variable does not exist!');
                end
            end
        end
        
        function s = addLabel(s,header,label,unit)
            % Adds a uniform label to the data in the object.
            data = repmat(label,s.ocount,1);
            s = s.addData(header,data,unit);
        end
        
        function s = addResponse(s,header,data,unit)
            % Adds data to the responses category
            s = s.addData(header,data,unit);
        end
        
        function s = addPredictor(s,header,data,unit)
            % Adds data to the predictor category
            s = s.addData(header,data,unit);
        end
        
        function [r,headers,units] = getResponse(s,n)
            % Gets data to the responses category
            [r,headers,units] = getData(s,n);
        end
        
        function [r,headers,units] = getPredictor(s,n)
            % Gets data to the predictor category
            [r,headers,units] = getData(s,n);
        end
        
        function newobj = combine(s,manyobjs)
            % Combines many statdata objects into one object by combining
            % all of the predictors and responses that exist in all of the
            % datasets.
            
            if isa(manyobjs, 'statsdata2')
                newobj = manyobjs;
                return
            end

            newobj = statsdata2();
            new_ocount = 0;
            for stats = manyobjs
               new_ocount = new_ocount + stats{1}.ocount; 
            end
            
            % determine which variables are shared between all objects
            for var = 1:manyobjs{1}.vcount
                
                header = manyobjs{1}.variables{var, 1};
                data = zeros(new_ocount, 1, 'like', manyobjs{1}.variables{var, 2});
                unit = manyobjs{1}.variables{var, 3};
                
                % concatenate the shared variables
                lo = 1; hi = 0; fail = false;
                for i = 1:numel(manyobjs)
                   hi = hi + manyobjs{i}.ocount;
                   try
                       this_data = manyobjs{i}.getData(header);
                   catch
                       fail = true;
                       break;
                   end
                   data(lo:hi) = this_data;
                   lo = lo + manyobjs{i}.ocount;
                end
                if fail
                    continue;                 
                else
                    newobj = newobj.addData(header,data,unit);
                end
            end
        end
    end
    
end

