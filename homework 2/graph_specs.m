function graph_specs(xsize,ysize,fontsize)
    % size in points
    xsize= xsize*96;  ysize= ysize*96;

    % Graph font
    set(groot,'defaultTextInterpreter','Latex');
    set(groot,'defaultGraphplotInterpreter','Latex');
    % Legend font
    set(groot, 'defaultLegendInterpreter', 'Latex');
    % Legend location
    set(groot, 'defaultLegendLocation', 'SouthOutside');
    set(groot, 'defaultLegendOrientation', 'horizontal');
    % Legend box
    set(groot, 'defaultLegendBox', 'off')
    

    % Line width
    set(groot, 'defaultLineLineWidth',1.5);
    set(groot, 'defaultConstantlineLineWidth',1.25);
    % Marker size
    set(groot, 'defaultLineMarkerSize', 2);
    set(groot, 'defaultGraphplotMarkerSize', 1);
    
    % Line labels
    set( groot, 'defaultConstantlineInterpreter','latex');
    set( groot, 'defaultConstantlineLabelHorizontalAlignment','left');
    
    % Grids
    set(groot,'defaultAxesXGrid','on');
    set(groot,'defaultAxesYGrid','on');
    % Axes box
    set(groot,'defaultAxesBox','on');

    % Tiled Layouts
    set(groot,'defaultTiledLayoutTileSpacing','loose');
    set(groot,'defaultTiledLayoutPadding','loose');
    
    % Font sizes
    set( groot, 'defaultAxesFontSize', fontsize/1.5 );
    set( groot, 'defaultAxesLabelFontSizeMultiplier', 1.5 );
    set( groot, 'defaultAxesTitleFontSizeMultiplier', 1.5 );

    % Figure size
    size=   get(0,'screensize');
    pos=    [ (size(3)-xsize)/2 (size(4)-ysize)/2  xsize ysize ];
    set(groot,'defaultFigurePosition',pos);


end