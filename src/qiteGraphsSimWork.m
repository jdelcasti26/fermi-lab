function qiteGraphs(stuMTQ3P, stuQ3P, stuMTQ2P, stuQ2P, xAxis, yAxis)
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
    noSamples = size(stuMTQ3P, 2);
    if(noSamples > 1)
        error("This graph does not expect arrays of samples");
    end
    % Number of Trotter steps
    noTS = size(stuMTQ3P(1).energy, 1);

    % Allocate arrays for data processing (axis, log, average and standard deviation)
    graphData = struct('X', zeros(noTS, 1), 'Y', zeros(noTS, 1));
    % Preallocation
    axisDataMTQ3P = repmat(graphData, noSamples, 1); 
    axisDataQ3P = repmat(graphData, noSamples, 1);
    axisDataMTQ2P = repmat(graphData, noSamples, 1); 
    axisDataQ2P = repmat(graphData, noSamples, 1);

    for s = 1 : noSamples
        % Define X axis according to specified option in string xAxis
        if(xAxis == "trotterSteps")
            axisDataMTQ3P.X = (1 : noTS)';
            axisDataQ3P.X = (1 : noTS)';
            axisDataMTQ2P.X = (1 : noTS)';
            axisDataQ2P.X = (1 : noTS)';            
            xAxisLabel = "Trotter step";
        elseif(xAxis == "costTotal")
            axisDataMTQ3P.X = stuMTQ3P.costLS + stuMTQ3P.costEE;
            axisDataQ3P.X = stuQ3P.costLS + stuQ3P.costEE;
            axisDataMTQ2P.X = stuMTQ2P.costLS + stuMTQ2P.costEE;
            axisDataQ2P.X = stuQ2P.costLS + stuQ2P.costEE;            
            xAxisLabel = "Pauli meas.";
        elseif(xAxis == "costLS")
            axisDataMTQ3P.X = stuMTQ3P.costLS;
            axisDataQ3P.X = stuQ3P.costLS;
            axisDataMTQ2P.X = stuMTQ2P.costLS;
            axisDataQ2P.X = stuQ2P.costLS;            
            xAxisLabel = "Pauli meas. (LS)";
        elseif(xAxis == "opCount")
            axisDataMTQ3P.X = stuMTQ3P.opCount;
            axisDataQ3P.X = stuQ3P.opCount;
            axisDataMTQ2P.X = stuMTQ2P.opCount;
            axisDataQ2P.X = stuQ2P.opCount;            
            xAxisLabel = "Depth (CC operators)";
        elseif(xAxis == "PSCount")
            axisDataMTQ3P.X = stuMTQ3P.PSCount;
            axisDataQ3P.X = stuQ3P.PSCount;
            axisDataMTQ2P.X = stuMTQ2P.PSCount;
            axisDataQ2P.X = stuQ2P.PSCount;            
            xAxisLabel = "Depth (Pauli rotations)";
        else
            error("Requested X axis unknown");
        end
        % Define Y axis according to specified option in string yAxis
        if(yAxis == "energy")
            axisDataMTQ3P.Y = stuMTQ3P.energy;
            axisDataQ3P.Y = stuQ3P.energy;
            axisDataMTQ2P.Y = stuMTQ2P.energy;
            axisDataQ2P.Y = stuQ2P.energy;            
            yAxisLabel = "Energy (Ha)";
        elseif(yAxis == "fidelity")
            axisDataMTQ3P.Y = log10(1 - stuMTQ3P.fidelity);
            axisDataQ3P.Y = log10(1 - stuQ3P.fidelity);
            axisDataMTQ2P.Y = log10(1 - stuMTQ2P.fidelity);
            axisDataQ2P.Y = log10(1 - stuQ2P.fidelity);            
            yAxisLabel = "1 - Fidelity";
        elseif(yAxis == "accuracy")
            axisDataMTQ3P.Y = log10(stuMTQ3P.accuracy);
            axisDataQ3P.Y = log10(stuQ3P.accuracy);
            axisDataMTQ2P.Y = log10(stuMTQ2P.accuracy);
            axisDataQ2P.Y = log10(stuQ2P.accuracy);            
            yAxisLabel = "Accuracy";
        else
            error("Requested Y axis unknown");
        end
    end % for s = 1 : noSamples


    figure
    hold on
    hQ3P = plot(axisDataQ3P(s).X, axisDataQ3P(s).Y, '-', 'Color', palette(2,:), 'LineWidth', 1.8, 'Marker','s');
    hMTQ3P = plot(axisDataMTQ3P(s).X, axisDataMTQ3P(s).Y, '-', 'Color', palette(4,:),'LineWidth', 1.8, 'Marker','s');
    hQ2P = plot(axisDataQ2P(s).X, axisDataQ2P(s).Y, '--', 'Color', palette(2,:), 'LineWidth', 1.8, 'Marker','s');
    hMTQ2P = plot(axisDataMTQ2P(s).X, axisDataMTQ2P(s).Y, '--', 'Color', palette(4,:),'LineWidth', 1.8, 'Marker','s');
    if(yAxis == "fidelity" || yAxis == "accuracy")
            yt = yticks; % Works
            yt = unique(floor(yt));
            yticks(yt); % Works for particular case
        yticklabels(num2str(10 .^ yt(:)));
    end
    hQ.Bar.LineWidth = 0.8;   % error bar stems
    hQ.Cap.LineWidth = 0.8;   % caps    
    hMTQ.Bar.LineWidth = 0.8;   % error bar stems
    hMTQ.Cap.LineWidth = 0.8;   % caps
    hold off

    xl = xlabel(xAxisLabel);
    yl = ylabel(yAxisLabel);
    legend([hQ2P, hMTQ2P, hQ3P, hMTQ3P], {"QITE 2P", "MT-QITE 2P", "QITE 3P", "MT-QITE 3P"});
    xl.FontSize = 22;
    yl.FontSize = 22;
    lgd.FontSize = 22;
    ax = gca;
    ax.FontSize = 22;
end