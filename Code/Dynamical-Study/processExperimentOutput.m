function [dataOut,ME] = processExperimentOutput(dataCell,param,varargin)
%This function collects the data from a cell array compiled by
%readParasiteExperimentData.m.  It computes a statistic on the data for
%each model.  Obviously this function is going to be extremely specific
%to that particular script and this particular experiment.  I am
%attempting to make it as general as possible by taking an input
%parameter structure that will account for everything I need.


%param is super complicated.  also varargin is probabling going to be a
%mess.

%varargin is any extra cell arrays needed to compute the statistic or
%errorbars.

%Requires function handles for calculating statistics and error bars.

%checking out what is inside the param structure:

ME = struct();
%Check for x input values
try
    xVals = param.xVals;
    nXVal = numel(xVals);
catch errXVals
    ME.errXVals = errXVals;
    return
end

%Check if we are showing a plot
try
    makePlot = param.makePlot;
catch errMakePlot
    ME.makePlot = errMakePlot;
    makePlot = false;
end

%Check if we are saving a plot
try
    savePlot = param.savePlot;
catch errSavePlot
    ME.errSavePlot = errSavePlot;
    savePlot = false;
end

%Check if we are going to calculate a regression
try 
    if param.regressionParam.doRegression
        doReg = true;
        regFcn = param.regressionParam.cmd;
        regOpts = param.regressionParam.regOpts;
        saveReg = param.regressionParam.saveReg;
    end
catch
    doReg = false;
    saveReg = false;
end

%Check if we are going to save the data to make our plots elsewhere.
try
    saveData = param.saveData;
    filepath = param.dataParam.filePath;
    
catch errSaveData
    ME.errSaveData = errSaveData;
    saveData = false;
end

%Make sure we can calculate statistics
if ~isa(param.statFcn,'function_handle')
    errorMessage = 'You must input function handle in param.statFun\n';
    error(errorMessage);
end

%Make sure we can calculate error bars
try
    
    errBarsFcn = param.errBars.errBarsFcn;
    errBarsParam = param.errBars.param;
    
    lNull = zeros(2,1);
    uNull = zeros(2,1);
    yNull = zeros(2,1); 
catch errErrBars
    ME.errBars = errErrBars;
    errBars = false;
end

%Make sure we have plotting parameters
try
    plotParam = param.plotParam;
    marksMatrix = plotParam.marksMatrix;
    legendEntries = plotParam.legendEntries;
    plotRange = plotParam.plotRange;
catch errPlotParam
    ME.errPlotParam = errPlotParam;
    makePlot = false;
    savePlot = false;
    warning('No plot parameters supplied (or something was missing); no plots will be made or saved.\n')
end

% Make sure we have the model codes input
try
    modelCodes = param.modelCodes;
    [nModels,nTreatments] = size(modelCodes);
catch errModelCodes
    ME.errModelCodes = errModelCodes;
    return
end

%Check if we are saving the regressions (Save it as a table of 16
%regressions with some statistics.
if saveReg
    filename = sprintf('%s%sRegressions.textab',filepath,param.dataAnalyzed);
    regFidTable = fopen(filename,'w');
    fprintf(regFidTable,'\\begin{tabular}{c |c c| c c |c c |c}\n\\hline\n');
    fprintf(regFidTable,'$\\beta_0$&$p$-value&\\$beta_1$&$p$-value&$\\beta_2$&$p$-value&Adj. $R^2$\\\\\n \\hline\n');
    %Going to include a column to include quadratics just to see if they
    %are worth including or not.  I choose quadratic vs. linear model based
    %on Adjusted r-squared; if linear is chosen, the column for the
    %quadratic term is left blank.
end


if makePlot
    fig = figure('Position', [0,0,1440,900]);
elseif savePlot&&(~makePlot)
    fig = figure('Position', [0,0,1440,900],'Visible','off');
end

dataOut = struct('center',zeros(nXVal,5,4),...
    'upper',zeros(nXVal,5,4),...
    'lower',zeros(nXVal,5,4));

%Here we end up hard-coding the factors anyway.  I'll elave this; maybe
%I'll come back to this/figure out a wayto generalize it.
factor1 = false(nModels,2);
factor2 = false(nModels,2);
factor3 = false(nModels,2);
factor4 = false(nModels,2);

%defining the factors.. did some code obfuscation here... 0.o
factor1(:,1) = modelCodes(:,1) == 1; %Small Free
factor1(:,2) = modelCodes(:,1) == 2; %Big Free
factor2(:,1) = modelCodes(:,2) == 1; %Big Para
factor2(:,2) = modelCodes(:,2) == 2; %Small Para
factor3(:,1) = modelCodes(:,3) == 1; %Don't include free LIving
factor3(:,2) = modelCodes(:,3) == 2; %Include Free Living
factor4(:,1) = modelCodes(:,4) == 1; %Don't include concomittant
factor4(:,2) = modelCodes(:,4) == 2; %Include concomittant

AXW = fullfact([2,2]);

%This subplot structure probably needs to be hard-coded.  that's how I'm
%going to want it, anyway.
yNull = zeros(2,1);
yNullVec = cell(2,1);

nullChoice = factor1(:,2)+1;
ymax = 0;
ymin = inf;
y = zeros(size(xVals));
yData = zeros(size(dataCell{1}));
for subplotChoice = 1:4
    
    if makePlot||savePlot
        subplot(2,2,subplotChoice);
        hold on
    end
    count = 1;
    dataOut.center(:,1,subplotChoice) = xVals;
    dataOut.upper(:,1,subplotChoice) = xVals;
    dataOut.lower(:,1,subplotChoice) = xVals;
    
    for factor1Choice = 1:2
        for factor2Choice = 1:2
            count = count+1;
            if makePlot||savePlot
                plotStyles = marksMatrix{factor2Choice,factor1Choice};
            end
            
            pickYs = factor2(:,factor2Choice)&...
                factor3(:,AXW(subplotChoice,1))&...
                factor4(:,AXW(subplotChoice,2))&...
                factor1(:,factor1Choice);
            if nargin >=3
                y0 = param.statFcn(dataCell{pickYs},varargin{1}(:),param.statFcnOpt);
            else
                y0 = param.statFcn(dataCell{pickYs},param.statFcnOpt);
            end
            
            if nargin ==4
                [ql0, qu0] = errBarsFcn(dataCell{pickYs},errBarsParam,varargin{2}{pickYs});
            else
                [ql0, qu0] = errBarsFcn(dataCell{pickYs},errBarsParam);
            end
            
            if numel(y0) == nXVal
                
                lNull(nullChoice(pickYs)) = ql0(1);
                uNull(nullChoice(pickYs)) = qu0(1);
                yNull(nullChoice(pickYs)) = y0(1);
                y = y0;
                ql = ql0;
                qu = qu0;
            else
                
                ql = [lNull(nullChoice(pickYs)) ql0];
                qu = [uNull(nullChoice(pickYs)) qu0];
                y = [yNull(nullChoice(pickYs)) y0];
                
            end
            
            dataOut.center(:,count,subplotChoice) = y;
            dataOut.upper(:,count,subplotChoice) = qu;
            dataOut.lower(:,count,subplotChoice) = ql;
            
            if makePlot||savePlot
                errorbar(xVals,y,ql,qu,plotStyles{:})
                ymax = max([ymax,max(y+qu)]);
                ymin = min([ymin,min(y-ql)]);
            end
            %Are we doing the regression?
            if doReg
                %Need to get the f**king null model (0% parasites) with y.
                yData0 = dataCell{pickYs};
                nRow = length(yData0);
                %Picking the right null model as before
                if numel(yData0) == nXVal*nRow

                    yNullVec{nullChoice(pickYs)} = yData0(:,1);
                    yData(:) = yData0(:);
                else
                
                    yData(:) = [yNullVec{nullChoice(pickYs)} yData0];
                    
                end
                
                %Repeat the xvalues as necessary
                nRow = length(yData);
                X = repmat(xVals,nRow,1);
                %Do the regressions.  If we don't specify options just do
                %this.  if we do, not going to compare the two options.
                if numel(fieldnames(regOpts))==0
                    linRegMdl = regFcn(X(:),yData(:));
                    quadRegMdl = regFcn(X(:),yData(:),'PureQuadratic');
                else
                    regMdl = regFcn(X(:),ydata(:),regOpts);
                end
                
                %If we are going to save the regressions...
                if saveReg
                    %pick the better regression
                    
                    %Multi row for each subplot
                    if (factor1Choice*factor2Choice) == 1
                        fprintf(regFidTable,'\\\\multirow{4}{*}{subplot %c}',64+subplotChoice);
                    end
                    if linRegMdl.Rsquared.Adjusted > quadRegMdl.Rsquared.Adjusted
                        %We choose linear here
                        b0 = linRegMdl.Coefficients{1,1};
                        b1 = linRegMdl.Coefficients{2,1};
                        pVal0 = linRegMdl.Coefficients{1,4};
                        pVal1 = linRegMdl.Coefficients{2,4};
                        fprintf(regFidTable,'&%.3f&%.3e&%.3f&%.3e&&&%.3f\\\\\n'...
                            ,b0,pVal0,b1,pVal1,linRegMdl.Rsquared.Adjusted);
                    else
                        b0 = quadRegMdl.Coefficients{1,1};
                        b1 = quadRegMdl.Coefficients{2,1};
                        b2 = quadRegMdl.Coefficients{3,1};
                        
                        pVal0 = quadRegMdl.Coefficients{1,4};
                        pVal1 = quadRegMdl.Coefficients{2,4};
                        pVal2 = quadRegMdl.Coefficients{3,4};
                        fprintf(regFidTable,'&%.3f&%.3e&%.3f&%.3e&%.3f&%.3e&%.3f\\\\\n'...
                            ,b0,pVal0,b1,pVal1,b2,pVal2,...
                            quadRegMdl.Rsquared.Adjusted);
                    end
                        
                end
            end
            grid on
        end
    end
    if saveReg
        fprintf(regFidTable,'\\hline\n');
    end
end

%Should probably automate all of this too..
if makePlot||savePlot
    
    ymax = 1.2*ymax;
    plotRange(4) = ymax;
    plotRange(3) = ymin;
    hold off
    
    
    subplot(2,2,1);
    title(param.plotParam.title1);
    ylabel(param.plotParam.ylabel)
    legend(legendEntries,'Location','Best')
    axis(plotRange);
    
    subplot(2,2,2);
    title(param.plotParam.title2);
    axesPosition = get(gca,'Position');
    ylabel(param.plotParam.ylabel)
    axis(plotRange);
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
    'Box','off');
    ylabel(hNewAxes,param.plotParam.y2label1,'FontWeight','bold')
    
    
    subplot(2,2,3);
    ylabel(param.plotParam.ylabel)
    xlabel(param.plotParam.xlabel)
    axis(plotRange);
    
    subplot(2,2,4);
    xlabel(param.plotParam.xlabel)
    ylabel(param.plotParam.ylabel)
    axesPosition = get(gca,'Position');
    axis(plotRange);
    hNewAxes = axes('Position',axesPosition,...
        'Color','none',...
        'YAxisLocation','right',...
        'XTick',[],...
        'YTick',[],...
        'Box','off');
    ylabel(hNewAxes,param.plotParam.y2label2,'FontWeight','bold')
end

if savePlot
    figFilename = param.plotParam.filename;
    format = param.plotParam.format;
    saveas(fig,figFilename,format)
end

if saveData
    headers = param.dataParam.headers;
    for condition = 1:4
        
        filename = sprintf('%s%s%s',filepath,param.dataAnalyzed);
        
        fid = fopen(sprintf('%s-small-small-subplot-%u',filename,condition),'w');
        fprintf(fid,'%s\n',headers);
        fprintf(fid,'%.5f,%.5f,%.5f,%.5f\n',...
            [dataOut.center(:,1,condition)';
             dataOut.center(:,2,condition)';
             dataOut.upper(:,2,condition)';
             dataOut.lower(:,2,condition)']);
        fclose(fid);
        
        fid = fopen(sprintf('%s-small-big-subplot-%u',filename,condition),'w');
        fprintf(fid,'%s\n',headers);
        fprintf(fid,'%.5f,%.5f,%.5f,%.5f\n',...
            [dataOut.center(:,1,condition)';
             dataOut.center(:,3,condition)';
             dataOut.upper(:,3,condition)';
             dataOut.lower(:,3,condition)']);
        fclose(fid);
        
        fid = fopen(sprintf('%s-big-small-subplot-%u',filename,condition),'w');
        fprintf(fid,'%s\n',headers);
        fprintf(fid,'%.5f,%.5f,%.5f,%.5f\n',...
            [dataOut.center(:,1,condition)';
             dataOut.center(:,4,condition)';
             dataOut.upper(:,4,condition)';
             dataOut.lower(:,4,condition)']);
        fclose(fid);
        
        fid = fopen(sprintf('%s-big-big-subplot-%u',filename,condition),'w');
        fprintf(fid,'%s\n',headers);
        fprintf(fid,'%.5f,%.5f,%.5f,%.5f\n',...
            [dataOut.center(:,1,condition)';
             dataOut.center(:,5,condition)';
             dataOut.upper(:,5,condition)';
             dataOut.lower(:,5,condition)']);
        fclose(fid);
        
    end
    fprintf(regFidTable,'\\hline\n\\end{tabular}');
end