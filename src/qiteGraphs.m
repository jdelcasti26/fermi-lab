function qiteGraphs(arrStuMTQ, arrStuQ, xAxis, yAxis)
    % Define palette (okabeIto)
    palette = [ ...
    0.000 0.000 0.000;  % Black
    0.902 0.624 0.000;  % Orange        (230,159,0)
    0.337 0.706 0.914;  % Sky blue      (86,180,233)
    0.000 0.620 0.451;  % Bluish green  (0,158,115)
    0.941 0.894 0.259;  % Yellow        (240,228,66)
    0.000 0.447 0.698;  % Blue          (0,114,178)
    0.835 0.369 0.000;  % Vermilion     (213,94,0)
    0.800 0.475 0.655]; % Reddish purple(204,121,167)

    % Retrieve size arrays of structs
    noSamples = size(arrStuMTQ, 2);
    % Number of Trotter steps
    noTS = size(arrStuMTQ(1).energy, 1);

    % Allocate arrays for data processing (axis, log, average and standard deviation)
    graphData = struct('X', zeros(noTS, 1), 'Y', zeros(noTS, 1));
    % Preallocation
    axisDataMTQ = repmat(graphData, noSamples, 1); 
    axisDataQ = repmat(graphData, noSamples, 1);
    % Average and standard deviation of the above
    avgMTQ = repmat(graphData, 1, 1); 
    avgQ = repmat(graphData, 1, 1);
    stdMTQ = repmat(graphData, 1, 1); 
    stdQ = repmat(graphData, 1, 1);

    for s = 1 : noSamples
        % Define X axis according to specified option in string xAxis
        if(xAxis == "trotterSteps")
            axisDataMTQ(s).X = (1 : noTS)';
            axisDataQ(s).X = (1 : noTS)';
            xAxisLabel = "Trotter step";
        elseif(xAxis == "costTotal")
            axisDataMTQ(s).X = arrStuMTQ(s).costLS + arrStuMTQ(s).costEE;
            axisDataQ(s).X = arrStuQ(s).costLS + arrStuQ(s).costEE;
            xAxisLabel = "Pauli meas.";
        elseif(xAxis == "costLS")
            axisDataMTQ(s).X = arrStuMTQ(s).costLS;
            axisDataQ(s).X = arrStuQ(s).costLS;
            xAxisLabel = "Pauli meas. (LS)";
        elseif(xAxis == "opCount")
            axisDataMTQ(s).X = arrStuMTQ(s).opCount;
            axisDataQ(s).X = arrStuQ(s).opCount;
            xAxisLabel = "Depth (CC operators)";
        elseif(xAxis == "PSCount")
            axisDataMTQ(s).X = arrStuMTQ(s).PSCount;
            axisDataQ(s).X = arrStuQ(s).PSCount;
            xAxisLabel = "Depth (Pauli rotations)";
        else
            error("Requested X axis unknown");
        end
        % Define Y axis according to specified option in string yAxis
        if(yAxis == "energy")
            axisDataMTQ(s).Y = arrStuMTQ(s).energy;
            axisDataQ(s).Y = arrStuQ(s).energy;
            yAxisLabel = "Energy (Ha)";
        elseif(yAxis == "fidelity")
            axisDataMTQ(s).Y = log10(1 - arrStuMTQ(s).fidelity);
            axisDataQ(s).Y = log10(1 - arrStuQ(s).fidelity);
            yAxisLabel = "1 - Fidelity";
        elseif(yAxis == "accuracy")
            axisDataMTQ(s).Y = log10(arrStuMTQ(s).accuracy);
            axisDataQ(s).Y = log10(arrStuQ(s).accuracy);
            yAxisLabel = "Accuracy";
        else
            error("Requested Y axis unknown");
        end
    end % for s = 1 : noSamples

    % Average and std calculation
    mtqX  = horzcat(axisDataMTQ.X);   
    mtqY  = horzcat(axisDataMTQ.Y); 
    qX  = horzcat(axisDataQ.X);   
    qY  = horzcat(axisDataQ.Y); 
    avgMTQ.X = mean(mtqX, 2, 'omitnan');
    avgMTQ.Y = mean(mtqY, 2, 'omitnan');
    avgQ.X = mean(qX, 2, 'omitnan');
    avgQ.Y = mean(qY, 2, 'omitnan');
    stdMTQ.Y = std(mtqY, 0, 2, 'omitnan');
    stdQ.Y = std(qY, 0, 2, 'omitnan');

    figure
    hold on
    % Samples are not plotted on computational cost graphs
    if(xAxis ~= "costTotal" && xAxis ~= "costLS")
        for s = 1 : noSamples
            plot(axisDataQ(s).X, axisDataQ(s).Y, '.', 'Color', palette(2,:));
            plot(axisDataMTQ(s).X, axisDataMTQ(s).Y, '.', 'Color', palette(4,:));
        end
        hQ = errorbar(avgQ.X, avgQ.Y, stdQ.Y, stdQ.Y, ':', 'Color', palette(2,:), 'LineWidth', 1.8, 'Marker','s');
        hMTQ = errorbar(avgMTQ.X, avgMTQ.Y, stdMTQ.Y, stdMTQ.Y, ':', 'Color', palette(4,:),'LineWidth', 1.8, 'Marker','s');
    else
        hQ = plot(avgQ.X, avgQ.Y, ':', 'Color', palette(2,:), 'LineWidth', 1.8, 'Marker','s');
        hMTQ = plot(avgMTQ.X, avgMTQ.Y, ':', 'Color', palette(4,:),'LineWidth', 1.8, 'Marker','s');
    end
    if(yAxis == "fidelity" || yAxis == "accuracy")
            yt = yticks; % Works
            yt = unique(floor(yt));
            yticks(yt); % Works for particular case
        yticklabels(num2str(10 .^ yt(:)));
    %else
    %    yticklabels(num2str(yt(:)));
    end
    if(xAxis ~= "costTotal" && xAxis ~= "costLS")
        hQ.Bar.LineWidth = 0.8;   % error bar stems
        hQ.Cap.LineWidth = 0.8;   % caps
        hMTQ.Bar.LineWidth = 0.8;   % error bar stems
        hMTQ.Cap.LineWidth = 0.8;   % caps
    end
    hold off

    xl = xlabel(xAxisLabel);
    yl = ylabel(yAxisLabel);
    legend([hQ, hMTQ], {"QITE", "MT-QITE"});
    xl.FontSize = 22;
    yl.FontSize = 22;
    lgd.FontSize = 22;
    ax = gca;
    ax.FontSize = 22;
end