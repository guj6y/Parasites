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

%figured out why it's name-value pairs.  varargin (2n-1) is the field
%name of the parameter strucutre and varargin(2n) is the value.  duh!

%checking out what is inside the param structure:

ME = struct();
try
    xVals = param.xVals;
    nXVal = numel(xVals);
catch errXVals
    ME.errXVals = errXVals;
    return
end

try
    makePlot = param.makePlot;
catch errMakePlot
    ME.makePlot = errMakePlot;
    makePlot = false;
end

try
    savePlot = param.savePlot;
catch errSavePlot
    ME.errSavePlot = errSavePlot;
    savePlot = false;
end

try
    saveData = param.saveData;
catch errSaveData
    ME.errSaveData = errSaveData;
    saveData = false;
end

if ~isa(param.statFcn,'function_handle')
    errorMessage = 'You must save input function handle in param.statFun\n';
    error(errorMessage);
end

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

try
    plotParam = param.plotParam;
    marksMatrix = plotParam.marksMatrix;
    legendEntries = plotParam.legendEntries;
    plotRange = plotParam.plotRange;
catch errPlotParam
    ME.errPlotParam = errPlotParam;
    makePlot = false;
    savePlot = false;
    warning('No plot parameters supplied; no plots will be made or saved.\n')
end

try
    modelCodes = param.modelCodes;
    [nModels,nTreatments] = size(modelCodes);
catch errModelCodes
    ME.errModelCodes = errModelCodes;
    return
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
nullChoice = factor1(:,2)+1;
ymax = 0;
ymin = inf;
for subplotChoice = 1:4;
    
    if makePlot||savePlot
        subplot(2,2,subplotChoice);
        hold on
    end
    count = 1;
    dataOut.center(:,1,subplotChoice) = xVals;
    dataOut.upper(:,1,subplotChoice) = xVals;
    dataOut.lower(:,1,subplotChoice) = xVals;
    
    for factor1Choice = 1:2;
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
                y = param.statFcn(dataCell{pickYs},varargin{1}(:),param.statFcnOpt);
            else
                y = param.statFcn(dataCell{pickYs},param.statFcnOpt);
            end
            
            if nargin ==4
                [ql, qu] = errBarsFcn(dataCell{pickYs},errBarsParam,varargin{2}{pickYs});
            else
                [ql, qu] = errBarsFcn(dataCell{pickYs},errBarsParam);
            end
            
            if numel(y) == nXVal;
                
                lNull(nullChoice(pickYs)) = ql(1);
                uNull(nullChoice(pickYs)) = qu(1);
                yNull(nullChoice(pickYs)) = y(1);
                
            else
                
                ql = [lNull(nullChoice(pickYs)) ql];
                qu = [uNull(nullChoice(pickYs)) qu];
                y = [yNull(nullChoice(pickYs)) y];
                
            end
            
            dataOut.center(:,count,subplotChoice) = y;
            dataOut.upper(:,count,subplotChoice) = qu;
            dataOut.lower(:,count,subplotChoice) = ql;
            
            if makePlot||savePlot
                errorbar(xVals,y,ql,qu,plotStyles{:})
                ymax = max([ymax,max(y+qu)]);
                ymin = min([ymin,min(y-ql)]);
            end
            grid on
        end
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
    filename = param.plotParam.filename;
    format = param.plotParam.format;
    
    saveas(fig,filename,format)
end

if saveData
    headers = param.dataParam.headers;
    for condition = 1:4
        filepath = param.dataParam.filePath;

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
    
end