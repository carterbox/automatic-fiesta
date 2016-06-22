classdef statsdata
    %STATSDATA is an object for organizing responses and predictors
    % This class stores the data in two 3D cells where:
    % column 1 : header
    % column 2 : data
    % column 3 : units
    
    % EXAMPLE:
    % s.responses = {'size', [100x1 double], 'm'};
    % s.predictors = {'color', [100x3 double], 'RBG'};

%% -----------------------------------------------------------------------
    properties
        responses = cell(10,3);
        rcount = 0; % number of responses categories
        predictors = cell(10,3);
        pcount = 0; % number of predictor categories
        ocount = 0; % number of observations
    end
    
    methods
        function s = addData(s,RorP,header,data,unit)
            % See how big the matrix has to be.
            if s.ocount ~= numel(data) && s.ocount > 0
                error('New data is not the same length!');    
            else
                s.ocount = numel(data);
            end
            
            if RorP
                s.rcount = s.rcount+1;
                s.responses{s.rcount,2} = data;
                s.responses{s.rcount,1} = header;
                s.responses{s.rcount,3} = unit;
            else
                s.pcount = s.pcount+1;
                s.predictors{s.pcount,2} = data;
                s.predictors{s.pcount,1} = header;
                s.predictors{s.pcount,3} = unit;
            end          
        end
        
        function [r,headers,units] = getData(s,RorP,n)
            % Returns responses by index or by string header. 0 returns all
            % of the available options.
            
            % Parse inputs
            if RorP
                data = s.responses;
                dcount = s.rcount;
            else
                data = s.predictors;
                dcount = s.pcount;
            end
            if isa(n,'char'), n = {n}; end
            
            % convert strings to indices
            if isa(n,'cell')
                for j = 1:numel(n)
                    name = n{j};
                    n{j} = -1;
                    for i = 1:dcount
                        if(strcmp(data{i,1},name))
                            n{j} = i;
                            break;
                        end
                    end
                    if n{j} < 1
                        error('Variable "%s" does not exist!', name);
                    end
                end
                n = cell2mat(n);
            end
            
            % Returns one lxm matrix of the responsess in the vector n.
            if isa(n, 'numeric') & n == 0, n = 1:dcount; end
            m = numel(n);
            % Check for valid indices
            for i = 1:m
            if n(i) < 1 || n(i) > dcount
                error('Variable number %i does not exist!', n(i));
            end
            end
            
            r = zeros(s.ocount,m);
            headers = cell(1,m);
            units = cell(1,m);
            for i = 1:m
                headers{i} = data{n(i),1};
                r(:,i) = data{n(i),2};
                units{i} = data{n(i),3};
            end
        end
        
        function s = scale(s,n,factor,newunit)
            % Multiplies a data by a factor
            if(isa(n,'char'))
                name = n;
                n = -1;
                for i = 1:s.pcount
                    if(strcmp(s.predictors{i,1},name))
                        n = i;
                        s.predictors{i,2} = s.predictors{i,2}*factor;
                        s.predictors{i,3} = newunit;
                        break;
                    end
                end
                if n < 1
                    error('Response "%s" does not exist!', name);
                end
            end
        end
        
        function s = addLabel(s,header,label,unit)
            % Adds a uniform label to the data in the object.
            data = repmat(label,s.ocount,1);
            s = s.addData(false,header,data,unit);
        end
        
        function s = addResponse(s,header,data,unit)
            % Adds data to the responses category
            s = s.addData(true,header,data,unit);
        end
        
        function s = addPredictor(s,header,data,unit)
            % Adds data to the predictor category
            s = s.addData(false,header,data,unit);
        end
        
        function [r,headers,units] = getResponse(s,n)
            % Gets data to the responses category
            [r,headers,units] = getData(s,true,n);
        end
        
        function [r,headers,units] = getPredictor(s,n)
            % Gets data to the predictor category
            [r,headers,units] = getData(s,false,n);
        end
        
        function s = combine(s,manyobjs)
            % Combines many statdata objects into one object by combining
            % all of the predictors and responsess that exist in all of the
            % datasets. #methodstub
            
            manyobjs{end + 1} = s; 
            
            % Determine the end length of the data
            m = numel(manyobjs);
            s.ocount = 0;
            for i = 1:m
                s.ocount = s.ocount + manyobjs{i}.ocount;
            end
            
            % Preallocate the arrays
            for i = 1:manyobjs{1}.rcount
                s.responses{i,1} = manyobjs{1}.responses{i,1};
                s.responses{i,2} = zeros(0,1,'like',manyobjs{1}.responses{i,2});
                s.responses{i,3} = manyobjs{1}.responses{i,3};
            end
            s.rcount = manyobjs{1}.rcount;
            for i = 1:manyobjs{1}.pcount
                s.predictors{i,1} = manyobjs{1}.predictors{i,1};
                s.predictors{i,2} = zeros(0,1,'like',manyobjs{1}.predictors{i,2});
                s.predictors{i,3} = manyobjs{1}.predictors{i,3};
            end
            s.pcount = manyobjs{1}.pcount;
            
            %strjoin(s.predictors(1:3))
            
            % Concatenate the arrays
            for i = 1:m
                for j = 1:manyobjs{i}.rcount
                    s.responses{j,2} = cat(1,s.responses{j,2},manyobjs{i}.responses{j,2});
                end
                for j = 1:manyobjs{i}.pcount
                    s.predictors{j,2} = cat(1,s.predictors{j,2},manyobjs{i}.predictors{j,2});
                end
            end
            
        end
    end
    
end

