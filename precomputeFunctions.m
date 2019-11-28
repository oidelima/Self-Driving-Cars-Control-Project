function precomputeFunctions()
    load TestTrack.mat
    centerLine = TestTrack.cline;
    for i=1:8
    centerLine = subdivideTrack(centerLine);
    end
    global gradient_precision
    global J1data J1satdata J2data J2satdata Jdata dJdata F_totaldata 
    global J1vals J1satvals J2vals J2satvals Jvals dJvals F_totalvals
    
    xlim1 = 200;
    xlim2 = 1500;
    ylim1 = -200;
    ylim2 = 825;
    [xGrid, yGrid] = meshgrid(xlim1:2:xlim2,ylim1:2:ylim2);
    xStacked = reshape(xGrid, 1,[]);
    yStacked = reshape(yGrid, 1,[]);
    dataStacked = costFunctionTrack([xStacked; yStacked], centerLine, 1);
    dataGrid = reshape(dataStacked, size(xGrid));
    Jdata{1} = xGrid;
    Jdata{2} = yGrid;
    Jdata{3} = dataGrid;
    disp("Computed cost function table.");
    
    ddataStacked = torGradient(@costFunctionTrack, [xStacked; yStacked], gradient_precision(2), centerLine, 1);
    ddataGrid1 = reshape(ddataStacked(1,:), size(xGrid));
    ddataGrid2 = reshape(ddataStacked(2,:), size(xGrid));
    dJdata{1} = xGrid;
    dJdata{2} = yGrid;
    dJdata{3} = ddataGrid1;
    dJdata{4} = ddataGrid2;
    disp("Computed cost function gradient table.");
end