function opt = output2file(opt, state, pop, type, varargin)
% Function: opt = output2file(opt, state, pop, type, varargin)
% Description: Output the population 'pop' to file. The file name is
% specified by 'opt.outputfile' field.
% Parameters: 
%   type : output type.  -1 = the last call, close the opened file. 
%          others(or no exist) = normal output
%   varargin : any parameter define in the options.outputfuns cell array.
%
%         LSSSSWC, NWPU
%    Revision: 1.2  Data: 2011-07-13
%*************************************************************************


if(isempty(opt.outputfile))
    return;  % the output file name is not specified, return directly
end

if( isfield(opt, 'outputfileFID') )
    fid = opt.outputfileFID;
else
    fid = [];
end

%*************************************************************************
% 1.Open the output file and output some population info
%*************************************************************************
if( isempty(fid) )
    fid = fopen(opt.outputfile, 'w');
    if( fid == 0)
        error('NSGA2:OutputFileError', 'Can not open output file!! file name:%s', opt.outputfile);
    end
    opt.outputfileFID = fid;
    
    % Output some infomation
    fprintf(fid, '#NSGA2\r\n');

    fprintf(fid, 'popsize %d\r\n', opt.popsize);
    fprintf(fid, 'maxGen %d\r\n', opt.maxGen);
    fprintf(fid, 'numVar %d\r\n', opt.numVar);
    fprintf(fid, 'numObj %d\r\n', opt.numObj);
    fprintf(fid, 'numCons %d\r\n', opt.numCons);
    
    % Output state field names
    fprintf(fid, 'stateFieldNames\t');
    names = fieldnames(state);
    for i = 1:length(names)
        fprintf(fid, '%s\t', names{i});
    end
    fprintf(fid, '\r\n');
    
    fprintf(fid, '#end\r\n\r\n\r\n');
end

%*************************************************************************
% 2. If this is the last call, close the output file
%*************************************************************************
if(type == -1)
    fclose(fid);
    rmfield(opt, 'outputfileFID');
    return
end

%*************************************************************************
% 3. Output population to file
%*************************************************************************
fprintf(fid, '#Generation %d / %d\r\n', state.currentGen, opt.maxGen);

% output each state field
names = fieldnames(state);
for i = 1:length(names)
    fprintf(fid, '%s\t%g\r\n', names{i}, getfield(state, names{i}));
end
fprintf(fid, '#end\r\n');

for i = 1:opt.numVar
    fprintf(fid, '%s\t', opt.nameVar{i});
end
for i = 1:opt.numObj
    fprintf(fid, 'Obj%d\t', i);
end
for i = 1:opt.numCons
    fprintf(fid, 'Cons%d\t', i);
end
for i = 1:opt.numProp
     fprintf(fid, '%s\t', opt.nameProp{i});
end
fprintf(fid, '\r\n');

for p = 1 : opt.popsize
    for i = 1:opt.numVar
        fprintf(fid, '%g\t', pop(p).var(i) );
    end
    for i = 1:opt.numObj
        fprintf(fid, '%g\t', pop(p).obj(i) );
    end
    for i = 1:opt.numCons
        fprintf(fid, '%g\t', pop(p).cons(i));
    end
    for i = 1:opt.numProp
        fprintf(fid, '%g\t', pop(p).properties(i));
    end
    fprintf(fid, '\r\n');
    
    %% visualise package of pareto
    vehicleplot=01;
    if vehicleplot==1 & pop(p).nViol==0
        vehicle=vehicleconcept(pop(p).var);
        h2=figure;
        h=axes;
        axis([-vehicle.Body.rearoverhang-vehicle.Body.wheelbase vehicle.Body.wheelbase+vehicle.Body.frontoverhang -1500 1500 -150 3000]);
        set(h,'Units','normalized')
        set(h,'Position',[0 0 1 1])
        % 2) make figure fill screen
        set(h2,'Units','normalized')
        set(h2,'Position',[0 0 1 1])
        %axis([-6000 6000 -1500 1500 -150 3000])
        h.Visible = 'off';
        vehicle.package(h)
        frame = getframe(h2);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); %
        filename=strcat('pareto.gif');
        % Write to the GIF File
        n=p; % n==1 if file does not exist
        if n == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        delete(h)
        delete(h2)
    end
    
end

fprintf(fid, '\r\n\r\n\r\n');





