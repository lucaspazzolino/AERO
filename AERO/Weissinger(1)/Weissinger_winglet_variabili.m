%% Weissinger winglet
close all
clear 
clc

%% Dati

U_Inf_Mag = 150; % 150 m/s
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag; % angolo di side-slip
rho = 1.225;

% Aereo 1
Rootchord1 = 1.625;
DihedralAngle1 = [0]; % [°]
SweepAngle1 = [0]; % [°]
TaperRatio1 = [0.672];  
Span1 = [11];
LEPosition_X1 = [0];
LEPosition_Y1 = [0];
LEPosition_Z1 = [0];
SemiSpanwiseDiscr1 = [10];
ChordwiseDiscr1 = [10];

% Aereo 2 - A10
Rootchord2 = 1.4;
DihedralAngle2 = [4]; % [°]
SweepAngle2 = [0]; % [°]
TaperRatio2 = [0];  
Span2 = [3.4];
LEPosition_X2 = [30];
LEPosition_Y2 = [0];
LEPosition_Z2 = [0];
SemiSpanwiseDiscr2 = [8];
ChordwiseDiscr2 = [8];

%% Selezione caso studio
fprintf('Selezionare il caso studio: \n')
fprintf('1 = Curva Cl-alpha primo corpo indisturbato \n')
fprintf('2 = Curva Cl-alpha secondo corpo indisturbato \n')
fprintf('3 = Curva Cl-alpha secondo corpo al variare dello sweep angle primo aereo \n')
fprintf('4 = Curva Cl-alpha secondo corpo al variare del dihedral angle primo aereo\n')
fprintf('5 = Curva Cl-alpha secondo corpo al variare del taper ratio primo aereo\n')
fprintf('6 = Curva Cl-alpha secondo corpo al variare del dihedral delle winglet \n \n')
fprintf('7 = Studio equilibrio \n \n')

value_input = input('Inserisci un numero: ');

if value_input == 1
    % Primo corpo indisturbato
    config.NCorpi = 1;
    
    config.RootChord = [Rootchord1];
    config.DihedralAngle = [DihedralAngle1]; % [°]
    config.SweepAngle = [SweepAngle1]; % [°]
    config.TaperRatio = [TaperRatio1]; 
    % config.AspectRatio = []; 
    config.Span = [Span1];
    config.LEPosition_X = [LEPosition_X1];
    config.LEPosition_Y = [LEPosition_Y1];
    config.LEPosition_Z = [LEPosition_Z1];
    
    config.RotationAngle_X = [0];
    % config.RotationAngle_Y = [5];
    config.RotationAngle_Z = [0];
    
    % Discretization options
    config.SemiSpanwiseDiscr = [SemiSpanwiseDiscr1];
    config.ChordwiseDiscr = [ChordwiseDiscr1];

    varing = config.SweepAngle; 
    alpha_vec = linspace(1,10,10);
elseif value_input == 2  
    % Secondo corpo indisturbato
    config.NCorpi = 1;
    
    config.RootChord = [Rootchord2];
    config.DihedralAngle = [DihedralAngle2]; % [°]
    config.SweepAngle = [SweepAngle2]; % [°]
    config.TaperRatio = [TaperRatio2]; 
    % config.AspectRatio = []; 
    config.Span = [Span2];
    config.LEPosition_X = [LEPosition_X2];
    config.LEPosition_Y = [LEPosition_Y2];
    config.LEPosition_Z = [LEPosition_Z2];
    
    config.RotationAngle_X = [0];
    % config.RotationAngle_Y = [5];
    config.RotationAngle_Z = [0];
    
    % Discretization options
    config.SemiSpanwiseDiscr = [SemiSpanwiseDiscr2];
    config.ChordwiseDiscr = [ChordwiseDiscr2];

    varing = config.SweepAngle;
    alpha_vec = linspace(1,10,10);
elseif value_input == 3  
    % variare dello sweep angle
    varing = linspace(15,45,10);
    alpha_vec = linspace(7,9,10);

    config.NCorpi = 2;
    
    config.RootChord = [Rootchord1,Rootchord2];
    config.DihedralAngle = [DihedralAngle1, DihedralAngle2]; % [°]
    % config.SweepAngle = [SweepAngle1,SweepAngle2]; % [°]
    config.TaperRatio = [TaperRatio1,TaperRatio2]; 
    % config.AspectRatio = []; 
    config.Span = [Span1,Span2];
    config.LEPosition_X = [LEPosition_X1,LEPosition_X2];
    config.LEPosition_Y = [LEPosition_Y1,LEPosition_Y2];
    config.LEPosition_Z = [LEPosition_Z1,LEPosition_Z2];
    
    config.RotationAngle_X = [0,0];
    % config.RotationAngle_Y = [5];
    config.RotationAngle_Z = [0,0];
    
    % Discretization options
    config.SemiSpanwiseDiscr = [SemiSpanwiseDiscr1,SemiSpanwiseDiscr2];
    config.ChordwiseDiscr = [ChordwiseDiscr1,ChordwiseDiscr2];
elseif value_input == 4
    % variare del dihedral angle
    varing = linspace(0,15,10);
    alpha_vec = linspace(7,9,10);

    config.NCorpi = 2;
    
    config.RootChord = [Rootchord1,Rootchord2];
    % config.DihedralAngle = [DihedralAngle1, DihedralAngle2]; % [°]
    config.SweepAngle = [SweepAngle1,SweepAngle2]; % [°]
    config.TaperRatio = [TaperRatio1,TaperRatio2]; 
    % config.AspectRatio = []; 
    config.Span = [Span1,Span2];
    config.LEPosition_X = [LEPosition_X1,LEPosition_X2];
    config.LEPosition_Y = [LEPosition_Y1,LEPosition_Y2];
    config.LEPosition_Z = [LEPosition_Z1,LEPosition_Z2];
    
    config.RotationAngle_X = [0,0];
    % config.RotationAngle_Y = [5];
    config.RotationAngle_Z = [0,0];
    
    % Discretization options
    config.SemiSpanwiseDiscr = [SemiSpanwiseDiscr1,SemiSpanwiseDiscr2];
    config.ChordwiseDiscr = [ChordwiseDiscr1,ChordwiseDiscr2];
elseif value_input == 5  
    % variare del taper ratio
    varing = linspace(0.1,1,10);
    alpha_vec = linspace(7,9,10);

    config.NCorpi = 2;
    
    config.RootChord = [Rootchord1,Rootchord2];
    config.DihedralAngle = [DihedralAngle1, DihedralAngle2]; % [°]
    config.SweepAngle = [SweepAngle1,SweepAngle2]; % [°]
    % config.TaperRatio = [TaperRatio1,TaperRatio2]; 
    % config.AspectRatio = []; 
    config.Span = [Span1,Span2];
    config.LEPosition_X = [LEPosition_X1,LEPosition_X2];
    config.LEPosition_Y = [LEPosition_Y1,LEPosition_Y2];
    config.LEPosition_Z = [LEPosition_Z1,LEPosition_Z2];
    
    config.RotationAngle_X = [0,0];
    % config.RotationAngle_Y = [5];
    config.RotationAngle_Z = [0,0];
    
    % Discretization options
    config.SemiSpanwiseDiscr = [SemiSpanwiseDiscr1,SemiSpanwiseDiscr2];
    config.ChordwiseDiscr = [ChordwiseDiscr1,ChordwiseDiscr2];
elseif value_input == 6  
    % variare del dihedral winglet
    varing = linspace(0,45,10);
    alpha_vec = linspace(0,25,15);

    config.NCorpi = 2;
    
    config.RootChord = [Rootchord1,Rootchord2];
    config.DihedralAngle = [DihedralAngle1, DihedralAngle2]; % [°]
    config.SweepAngle = [SweepAngle1,SweepAngle2]; % [°]
    config.TaperRatio = [TaperRatio1,TaperRatio2]; 
    % config.AspectRatio = []; 
    config.Span = [Span1,Span2];
    config.LEPosition_X = [LEPosition_X1,LEPosition_X2];
    config.LEPosition_Y = [LEPosition_Y1,LEPosition_Y2];
    config.LEPosition_Z = [LEPosition_Z1,LEPosition_Z2];
    
    config.RotationAngle_X = [0,0];
    % config.RotationAngle_Y = [5];
    config.RotationAngle_Z = [0,0];
    
    % Discretization options
    config.SemiSpanwiseDiscr = [SemiSpanwiseDiscr1,SemiSpanwiseDiscr2];
    config.ChordwiseDiscr = [ChordwiseDiscr1,ChordwiseDiscr2];
elseif value_input == 7  
    % variare dell'angolo d'attacco
    varing = 7.11; % angolo attacco primo corpo
    alpha_vec = 7.975; % angolo attacco secondo corpo

    config.NCorpi = 2;
    
    config.RootChord = [Rootchord1,Rootchord2];
    config.DihedralAngle = [DihedralAngle1, DihedralAngle2]; % [°]
    config.SweepAngle = [SweepAngle1,SweepAngle2]; % [°]
    config.TaperRatio = [TaperRatio1,TaperRatio2]; 
    % config.AspectRatio = []; 
    config.Span = [Span1,Span2];
    config.LEPosition_X = [LEPosition_X1,LEPosition_X2];
    config.LEPosition_Y = [LEPosition_Y1,LEPosition_Y2];
    config.LEPosition_Z = [LEPosition_Z1,LEPosition_Z2];
    
    config.RotationAngle_X = [0,0];
    % config.RotationAngle_Y = [5];
    config.RotationAngle_Z = [0,0];
    
    % Discretization options
    config.SemiSpanwiseDiscr = [SemiSpanwiseDiscr1,SemiSpanwiseDiscr2];
    config.ChordwiseDiscr = [ChordwiseDiscr1,ChordwiseDiscr2];
end

Cl_corpo1 = zeros(length(varing),length(alpha_vec));
Cl_corpo2 = zeros(length(varing),length(alpha_vec));
Cd_corpo1 = zeros(length(varing),length(alpha_vec));
Cd_corpo2 = zeros(length(varing),length(alpha_vec));
Efficienza1 = zeros(length(varing),length(alpha_vec));
Efficienza2 = zeros(length(varing),length(alpha_vec));
I1 = 0;
for variable = varing
    I1 = I1+1;
    if value_input == 1 || value_input == 2
        config.SweepAngle = variable;
    elseif value_input == 3
        config.SweepAngle = [variable, SweepAngle2];
    elseif value_input == 4
        config.DihedralAngle = [variable DihedralAngle2];
    elseif value_input == 5
        config.TaperRatio = [variable TaperRatio2];
    end
    
    I2 = 0;
    for alpha = alpha_vec
        I2 = I2+1;
        if value_input == 1 || value_input == 2
            config.RotationAngle_Y = alpha;
        elseif value_input == 3
            config.RotationAngle_Y = [7.2, alpha];
        elseif value_input == 4
            config.RotationAngle_Y = [7.2, alpha];
        elseif value_input == 5
            config.RotationAngle_Y = [7.2, alpha];
        elseif value_input == 6
            config.RotationAngle_Y = [7.2, alpha];
        elseif value_input == 7
            config.RotationAngle_Y = [variable, alpha];
        end
        if value_input ~=2
            if value_input == 6
                config2.DihedralAngle = variable;
            else
                config2.DihedralAngle = [75]; 
            end
        end
    
        %% Preliminary computations - BODY
        
        % Computing the span
        config.SemiSpan = config.Span./2;
        % Computing the surface
        config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
        config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
        % Computing the Tip chord
        config.TipChord = config.RootChord .* config.TaperRatio;
        
        % Compute MAC
        config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));
        
        % Create the geometry structure
        
        ControlPoints = cell(config.NCorpi, 1);
        InducedPoints = cell(config.NCorpi, 1);
        Normals = cell(config.NCorpi, 1);
        InfiniteVortices = cell(config.NCorpi, 1);
        Vortices = cell(config.NCorpi, 1);
        internalMesh = cell(config.NCorpi, 1);
        WingExtremes = cell(config.NCorpi, 1);
        
        % ControlPoints2 = cell(config.NCorpi, 1);
        % InducedPoints2 = cell(config.NCorpi, 1);
        % Normals2 = cell(config.NCorpi, 1);
        % InfiniteVortices2 = cell(config.NCorpi, 1);
        % Vortices2 = cell(config.NCorpi, 1);
        % internalMesh2 = cell(config.NCorpi, 1);
        % WingExtremes2 = cell(config.NCorpi, 1);
        
        Dihedral_vec = zeros(config.NCorpi,2); 
        
        if value_input ~=2
            for iCorpo = 1:1   
                config.SemiSpanwiseDiscr(1) = SemiSpanwiseDiscr1;
                [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
                %% Preliminary computations - WINGLET POSITIVA
                config2.RootChord = [config.TipChord(iCorpo)];
                % config2.DihedralAngle = [75]; % [°]
                config2.SweepAngle = [60]; % [°]
                % config.AspectRatio = []; 
                config2.Span = [2];
                config2.LEPosition_X = internalMesh{1, 1}{1, 1}.LEtip(1);
                config2.LEPosition_Y = internalMesh{1, 1}{1, 1}.LEtip(2);
                config2.LEPosition_Z = internalMesh{1, 1}{1, 1}.LEtip(3);
                
                config2.RotationAngle_X = [0];
                config2.RotationAngle_Y = [alpha];
                config2.RotationAngle_Z = [0];
                
                % Discretization options
                config2.SemiSpanwiseDiscr = [1];
                config2.ChordwiseDiscr = [config.ChordwiseDiscr(iCorpo)];
                
                config2.TaperRatio = [0.5];
                config2.SemiSpan = config2.Span./2;
                % Computing the surface
                config2.Surface = 2 * (config2.SemiSpan .* config2.RootChord .* ( 1 + config2.TaperRatio ) ./ 2);
                config2.SurfaceProjected = config2.Surface .* cosd(config2.DihedralAngle);
                % Computing the Tip chord
                config2.TipChord = config2.RootChord .* config2.TaperRatio;
                
                % Compute MAC
                config2.MAC = (2/3) .* config2.RootChord .* ( (1 + config2.TaperRatio + config2.TaperRatio.^2)./(1 + config2.TaperRatio));
                [ControlPoints2{iCorpo}, InducedPoints2{iCorpo}, Normals2{iCorpo}, InfiniteVortices2{iCorpo}, Vortices2{iCorpo}, internalMesh2{iCorpo}, WingExtremes2{iCorpo}] = createStructure(config2, iCorpo);
                ControlPoints2{iCorpo}(:, 2) = []; 
                Normals2{iCorpo}(:, 2) = []; 
                InfiniteVortices2{iCorpo}(:, 2) = []; 
                InducedPoints2{iCorpo}(:, 2) = []; 
                internalMesh2{iCorpo}(:, 2) = []; 
                % WingExtremes2{iCorpo}(:, 2) = []; 
                Vortices2{iCorpo}(:, 2) = [];
                
                ControlPoints{iCorpo} = [ControlPoints2{iCorpo} ControlPoints{iCorpo}]; 
                InducedPoints{iCorpo} = [InducedPoints2{iCorpo} InducedPoints{iCorpo}];
                Normals{iCorpo} = [Normals2{iCorpo} Normals{iCorpo}];
                Vortices{iCorpo} = [Vortices2{iCorpo} Vortices{iCorpo}];
                InfiniteVortices{iCorpo} = [InfiniteVortices2{iCorpo} InfiniteVortices{iCorpo}];
                internalMesh{iCorpo} = [internalMesh2{iCorpo} internalMesh{iCorpo}];
                WingExtremes{iCorpo} = [WingExtremes2{iCorpo} WingExtremes{iCorpo}];
                
                Dihedral_vec(1,1) = config2.DihedralAngle;
                
                %% Preliminary computations - WINGLET NEGATIVA
                config2.DihedralAngle = [-180-config2.DihedralAngle]; % [°] 
                config2.LEPosition_X = internalMesh{1, 1}{1, end}.LEtip(1);
                config2.LEPosition_Y = internalMesh{1, 1}{1, end}.LEtip(2);
                config2.LEPosition_Z = internalMesh{1, 1}{1, end}.LEtip(3);
                
                % Computing the surface
                config2.Surface = 2 * (config2.SemiSpan .* config2.RootChord .* ( 1 + config2.TaperRatio ) ./ 2);
                config2.SurfaceProjected = config2.Surface .* cosd(config2.DihedralAngle);
                % Computing the Tip chord
                config2.TipChord = config2.RootChord .* config2.TaperRatio;
                
                % Compute MAC
                config2.MAC = (2/3) .* config2.RootChord .* ( (1 + config2.TaperRatio + config2.TaperRatio.^2)./(1 + config2.TaperRatio));
                [ControlPoints2{iCorpo}, InducedPoints2{iCorpo}, Normals2{iCorpo}, InfiniteVortices2{iCorpo}, Vortices2{iCorpo}, internalMesh2{iCorpo}, WingExtremes2{iCorpo}] = createStructure(config2, iCorpo);
                ControlPoints2{iCorpo}(:, 2) = []; 
                Normals2{iCorpo}(:, 2) = []; 
                InfiniteVortices2{iCorpo}(:, 2) = []; 
                InducedPoints2{iCorpo}(:, 2) = []; 
                internalMesh2{iCorpo}(:, 2) = []; 
                % WingExtremes2{iCorpo}(:, 2) = []; 
                Vortices2{iCorpo}(:, 2) = [];
                
                ControlPoints{iCorpo} = [ControlPoints{iCorpo} ControlPoints2{iCorpo}]; 
                InducedPoints{iCorpo} = [InducedPoints{iCorpo} InducedPoints2{iCorpo}];
                Normals{iCorpo} = [Normals{iCorpo} Normals2{iCorpo}];
                Vortices{iCorpo} = [Vortices{iCorpo} Vortices2{iCorpo}];
                InfiniteVortices{iCorpo} = [InfiniteVortices{iCorpo} InfiniteVortices2{iCorpo}];
                internalMesh{iCorpo} = [internalMesh{iCorpo} internalMesh2{iCorpo}];
                WingExtremes{iCorpo} = [WingExtremes{iCorpo} WingExtremes2{iCorpo}];
                
                Dihedral_vec(1,2) = config2.DihedralAngle;
            end
            
            for iCorpo = 2:config.NCorpi
                [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
                Dihedral_vec(iCorpo,:) = config.DihedralAngle(iCorpo);
            end
            
            config.SemiSpanwiseDiscr(1) = SemiSpanwiseDiscr1+1;
%             if value_input == 2
%                 config.SemiSpanwiseDiscr(1) = SemiSpanwiseDiscr2+1;
%             else
%                 config.SemiSpanwiseDiscr(1) = SemiSpanwiseDiscr1+1;
%             end
        else
            iCorpo = 1;
            [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);
            Dihedral_vec(iCorpo,:) = config.DihedralAngle(iCorpo);
        end
            
        %% PLOT
        if (value_input == 1 || value_input == 2) && alpha == alpha_vec(1)
            for iCorpo = 1:config.NCorpi
            
                Mesh = internalMesh{iCorpo};
                [alt,lung] = size(Mesh);
            
                MPointsX = zeros(alt,lung);
                MPointsY = zeros(alt,lung);
                MPointsZ = zeros(alt,lung);  
                for iPan = 1:alt
                    A1 = 0;
                    for jPan = 1:lung/2
                        A1 = A1+1;
                        MPointsX(iPan,jPan) = Mesh{iPan,jPan}.LEtip(1);
                        MPointsY(iPan,jPan) = Mesh{iPan,jPan}.LEtip(2);
                        MPointsZ(iPan,jPan) = Mesh{iPan,jPan}.LEtip(3);
                    
                    end
                    A1 = A1+1;
                    MPointsX(iPan,A1) = Mesh{iPan,jPan}.LERoot(1);
                    MPointsY(iPan,A1) = Mesh{iPan,jPan}.LERoot(2);
                    MPointsZ(iPan,A1) = Mesh{iPan,jPan}.LERoot(3);
                    for jPan = lung/2+1:lung
                        A1 = A1+1;
                        MPointsX(iPan,A1) = Mesh{iPan,jPan}.LEtip(1);
                        MPointsY(iPan,A1) = Mesh{iPan,jPan}.LEtip(2);
                        MPointsZ(iPan,A1) = Mesh{iPan,jPan}.LEtip(3);
                    end
                end
                A1 = 0;
                for jPan = 1:lung/2
                    A1 = A1+1;
                    MPointsX(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(1);
                    MPointsY(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(2);
                    MPointsZ(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(3);
            
                end
                A1 = A1+1;
                MPointsX(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(1);
                MPointsY(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(2);
                MPointsZ(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(3);
                for jPan = lung/2+1:lung
                    A1 = A1+1;
                    MPointsX(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(1);
                    MPointsY(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(2);
                    MPointsZ(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(3);
                end
            
                surf(MPointsX,MPointsY,MPointsZ,'FaceColor','interp')
                hold on
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                colorbar;
                grid on;
                title('Profili alari');
                axis equal
            end
        end
        
            
        
        %% Matrices initialization
        
        NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
        matriceA = zeros(NPanelsTot, NPanelsTot);
        TermineNoto = zeros(NPanelsTot, 1);
        
        %% Construction of the matrix
        
        rowIndex = 0;
        for iCorpo = 1:config.NCorpi
            
            % Cycle on all of its chordwise panels
            for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
                % Cycle on all of its spanwise panels
                for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
                    
                    % Update row index
                    rowIndex = rowIndex + 1;
           
                    columnIndex = 0;
                    
                    ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
                    NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
                    
                    
                    for jCorpo = 1:config.NCorpi
                        
                        % Cycle on all of its chordwise panels
                        for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                            % Cycle on all of its spanwise panels
                            for SpanPanel_j = 1:2*(config.SemiSpanwiseDiscr(jCorpo))
                                
                                % Update column index
                                columnIndex = columnIndex + 1;
                                
                                % Compute the influence induced by first
                                % semi-infinite vortex
                                Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                                Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                                U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2,1);
        
                                
                                % Compute the influence induced by finite vortex
                                Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                                Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                                U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2,1);
        
                                
                                % Compute the influence induced by second
                                % semi-infinite vortex
                                Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                                Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                                U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2,1);
        
                                
                                matriceA(rowIndex, columnIndex) = dot(U, NormalHere);
                               
                                
                            end
                        end
        
                    end
                    
                
                    
                end
            end
        end
        
        %% Costruzione del termine noto
        rowIndex = 0;
        for iCorpo = 1:config.NCorpi
            
            % Cycle on all of its chordwise panels
            for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
                % Cycle on all of its spanwise panels
                for SpanPanel_i = 1:2*(config.SemiSpanwiseDiscr(iCorpo))
                    
                    % Update row index
                    rowIndex = rowIndex + 1;
          
                    NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
                    
                    TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
                    
                end
            end
        end
        
        %% Solve the linear system
        
        Solution = linsolve(matriceA, TermineNoto);
        
        Gamma = cell(config.NCorpi, 1);
        
        rowIndex = 0;
        for iCorpo = 1:config.NCorpi
            
            Gamma{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), (config.SemiSpanwiseDiscr(iCorpo))*2 );
            
             % Cycle on all of its chordwise panels
            for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
                % Cycle on all of its spanwise panels
                for SpanPanel_i = 1:2*(config.SemiSpanwiseDiscr(iCorpo))
                    
                    % Update row index
                    rowIndex = rowIndex + 1;
                    
                    Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
                end
                
            end
        %     Gamma{iCorpo}(:, end) = -Gamma{iCorpo}(:, end);
        end
        
        
        %% Compute the 2D and 3D Lift
        % Punti di collocazione
        quarterPoints = cell(config.NCorpi,1);
        
        for i = 1:config.NCorpi
            Mesh = internalMesh{i};
            Pos_LE = zeros(3,(config.SemiSpanwiseDiscr(i))*2);
            Pos_TE = zeros(3,(config.SemiSpanwiseDiscr(i))*2);
            for j = 1:(config.SemiSpanwiseDiscr(i))*2
                Pos_LE1 = Mesh{1,j}.LERoot(:);
                Pos_LE2 = Mesh{1,j}.LEtip(:);
                Pos_LE(1:3,j) = (Pos_LE1+Pos_LE2)./2;
                Pos_TE1 = Mesh{end,j}.TERoot(:);
                Pos_TE2 = Mesh{end,j}.TEtip(:);
                Pos_TE(1:3,j) = (Pos_TE1+Pos_TE2)./2;
            end
            quarterPoints{i} = Pos_LE+1/4*(Pos_TE-Pos_LE);
        end
        
        
        L2D = cell(config.NCorpi,1); 
        Cl2D = cell(config.NCorpi,1);
        corda = cell(config.NCorpi,1);
        
        for iCorpo = 1:config.NCorpi
            L2D{iCorpo} = zeros(1,(config.SemiSpanwiseDiscr(iCorpo))*2);
            Cl2D{iCorpo} = zeros(1,(config.SemiSpanwiseDiscr(iCorpo))*2);
        
            Mesh = internalMesh{iCorpo};
            for kcol = 1:1
                for kriga = 1:config.ChordwiseDiscr(iCorpo)
                    MatrixGamma = Gamma{iCorpo};
                    L2D{iCorpo}(1,kcol) = L2D{iCorpo}(1,kcol)+rho*U_Inf_Mag*MatrixGamma(kriga,kcol)*cosd(Dihedral_vec(iCorpo,1));
                end
                Pos_LE1 = Mesh{1,kcol}.LERoot(:);
                Pos_LE2 = Mesh{1,kcol}.LEtip(:);
                Pos_LE(1:3,kcol) = (Pos_LE1+Pos_LE2)./2;
                Pos_TE1 = Mesh{end,kcol}.TERoot(:);
                Pos_TE2 = Mesh{end,kcol}.TEtip(:);
                Pos_TE(1:3,kcol) = (Pos_TE1+Pos_TE2)./2;
                corda_var{iCorpo}(:,kcol) = Pos_TE(:,kcol)-Pos_LE(:,kcol);
                corda{iCorpo}(kcol) = norm(corda_var{iCorpo}(:,kcol));
            end
            for kcol = 2:(config.SemiSpanwiseDiscr(iCorpo))*2-1
                for kriga = 1:config.ChordwiseDiscr(iCorpo)
                    MatrixGamma = Gamma{iCorpo};
                    L2D{iCorpo}(1,kcol) = L2D{iCorpo}(1,kcol)+rho*U_Inf_Mag*MatrixGamma(kriga,kcol)*cosd(config.DihedralAngle(iCorpo));
                end
                Pos_LE1 = Mesh{1,kcol}.LERoot(:);
                Pos_LE2 = Mesh{1,kcol}.LEtip(:);
                Pos_LE(1:3,kcol) = (Pos_LE1+Pos_LE2)./2;
                Pos_TE1 = Mesh{end,kcol}.TERoot(:);
                Pos_TE2 = Mesh{end,kcol}.TEtip(:);
                Pos_TE(1:3,kcol) = (Pos_TE1+Pos_TE2)./2;
                corda_var{iCorpo}(:,kcol) = Pos_TE(:,kcol)-Pos_LE(:,kcol);
                corda{iCorpo}(kcol) = norm(corda_var{iCorpo}(:,kcol));
            end
            for kcol = (config.SemiSpanwiseDiscr(iCorpo))*2:(config.SemiSpanwiseDiscr(iCorpo))*2
                for kriga = 1:config.ChordwiseDiscr(iCorpo)
                    MatrixGamma = Gamma{iCorpo};
                    L2D{iCorpo}(1,kcol) = L2D{iCorpo}(1,kcol)+rho*U_Inf_Mag*MatrixGamma(kriga,kcol)*cosd(Dihedral_vec(iCorpo,2));
                end
                Pos_LE1 = Mesh{1,kcol}.LERoot(:);
                Pos_LE2 = Mesh{1,kcol}.LEtip(:);
                Pos_LE(1:3,kcol) = (Pos_LE1+Pos_LE2)./2;
                Pos_TE1 = Mesh{end,kcol}.TERoot(:);
                Pos_TE2 = Mesh{end,kcol}.TEtip(:);
                Pos_TE(1:3,kcol) = (Pos_TE1+Pos_TE2)./2;
                corda_var{iCorpo}(:,kcol) = Pos_TE(:,kcol)-Pos_LE(:,kcol);
                corda{iCorpo}(kcol) = norm(corda_var{iCorpo}(:,kcol));
            end
    
            coeff = (1/2*rho*U_Inf_Mag^2*corda{iCorpo});
            Cl2D{iCorpo} = L2D{iCorpo}./coeff;
    
            if (value_input == 1 || value_input == 2) && alpha == alpha_vec(1)
                Matrixcontpoints = ControlPoints{iCorpo};
                ControlPointsY = zeros(1,(config.SemiSpanwiseDiscr(iCorpo))*2);
                for iPan = 1:(config.SemiSpanwiseDiscr(iCorpo))*2
                    ControlPointsY(1,iPan) = Matrixcontpoints{1,iPan}.Coords(2);
                end
                estremosx = Matrixcontpoints{1,1}.Coords(2);
                estremodx = Matrixcontpoints{1,end}.Coords(2);
                xx = linspace(estremosx,estremodx,(config.SemiSpanwiseDiscr(iCorpo))*20);
                yy = spline(ControlPointsY,Cl2D{iCorpo},xx);
                figure()
                plot(ControlPointsY,Cl2D{iCorpo},'o',xx,yy);
                title('curva Cl2D del corpo ',iCorpo);
                xlabel('Apertura')
                ylabel('Cl')
                grid on
            end
        
            %Mesh = internalMesh{iCorpo};
            [alt,lung] = size(Mesh);
        
            MPointsX = zeros(alt,lung);
            MPointsY = zeros(alt,lung);
            MPointsZ = zeros(alt,lung);
        
            for iPan = 1:alt
                A1 = 0;
                for jPan = 1:lung/2
                    A1 = A1+1;
                    MPointsX(iPan,jPan) = Mesh{iPan,jPan}.LEtip(1);
                    MPointsY(iPan,jPan) = Mesh{iPan,jPan}.LEtip(2);
                    MPointsZ(iPan,jPan) = Mesh{iPan,jPan}.LEtip(3);
        
                end
                A1 = A1+1;
                MPointsX(iPan,A1) = Mesh{iPan,jPan}.LERoot(1);
                MPointsY(iPan,A1) = Mesh{iPan,jPan}.LERoot(2);
                MPointsZ(iPan,A1) = Mesh{iPan,jPan}.LERoot(3);
                for jPan = lung/2+1:lung
                    A1 = A1+1;
                    MPointsX(iPan,A1) = Mesh{iPan,jPan}.LEtip(1);
                    MPointsY(iPan,A1) = Mesh{iPan,jPan}.LEtip(2);
                    MPointsZ(iPan,A1) = Mesh{iPan,jPan}.LEtip(3);
                end
            end
            A1 = 0;
            for jPan = 1:lung/2
                A1 = A1+1;
                MPointsX(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(1);
                MPointsY(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(2);
                MPointsZ(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(3);
        
            end
            A1 = A1+1;
            MPointsX(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(1);
            MPointsY(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(2);
            MPointsZ(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(3);
            for jPan = lung/2+1:lung
                A1 = A1+1;
                MPointsX(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(1);
                MPointsY(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(2);
                MPointsZ(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(3);
            end
            if i==1 && value_input ~= 2
                MatrixGamma(:,end) = -MatrixGamma(:,end);
            else
                MatrixGamma(:,end) = MatrixGamma(:,end);
            end
            if (value_input == 1||value_input == 2) && alpha == alpha_vec(1)
                figure()
                surf(MPointsX,MPointsY,MPointsZ,MatrixGamma)
                hold on
                plot3(quarterPoints{iCorpo}(1,:),quarterPoints{iCorpo}(2,:),quarterPoints{iCorpo}(3,:),'o','MarkerSize',10,'Color','r')
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                colorbar;
                grid on;
                title('Distribuzione 3D di Gamma sul corpo',iCorpo);
                axis equal
            end
        
            % sistemo il segno
            L2D{iCorpo}(1,kcol) = -L2D{iCorpo}(1,kcol);
        end
        
        
        %% Compute 2D and 3D induced drag
        % Velocità indotta
        V_ind = cell(config.NCorpi,1);
        
        for iCorpo = 1:config.NCorpi
            for SpanPanel_i = 1:2*(config.SemiSpanwiseDiscr(iCorpo))
                U = 0;           
                ControlPointHere = (quarterPoints{iCorpo}(:,SpanPanel_i))';
                
                for jCorpo = 1:config.NCorpi
                    
                    % Cycle on all of its chordwise panels
                    for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                        % Cycle on all of its spanwise panels
                        for SpanPanel_j = 1:2*(config.SemiSpanwiseDiscr(jCorpo))
        
                            % Compute the influence induced by first
                            % semi-infinite vortex
                            Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                            Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                            U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2, Gamma{jCorpo}(ChordPanel_j,SpanPanel_j));
                            
                            % Compute the influence induced by second
                            % semi-infinite vortex
                            Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                            Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                            U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2, Gamma{jCorpo}(ChordPanel_j,SpanPanel_j));
                        end
                    end
                end
                V_ind{iCorpo}{SpanPanel_i} = U;
            end
        end
        
        %Alpha indotto e drag
        alpha_ind = cell(config.NCorpi,1);
        D2D = cell(config.NCorpi,1);
        Cd2D = cell(config.NCorpi,1);
        
        for i = 1:config.NCorpi  
            D2D{i} = zeros(1,(config.SemiSpanwiseDiscr(i))*2);
            for j = 1:1
                ez = [sind(config.RotationAngle_Y(i)); -sind(Dihedral_vec(i,1)); cosd(config.RotationAngle_Y(i))*cosd(Dihedral_vec(i,1))];
                ez_norm = norm(ez);
                ez = ez./ez_norm;
                alpha_ind{i}{j} = atan(dot(V_ind{i}{j},ez)/U_Inf_Mag);
                D2D{i}(j) = -L2D{i}(j)*sin(alpha_ind{i}{j});
            end
            for j = 2:(config.SemiSpanwiseDiscr(i))
                ez = [sind(config.RotationAngle_Y(i)); -sind(config.DihedralAngle(i)); cosd(config.RotationAngle_Y(i))*cosd(config.DihedralAngle(i))];
                ez_norm = norm(ez);
                ez = ez./ez_norm;
                alpha_ind{i}{j} = atan(dot(V_ind{i}{j},ez)/U_Inf_Mag);
                D2D{i}(j) = -L2D{i}(j)*sin(alpha_ind{i}{j});
            end
            for j = config.SemiSpanwiseDiscr(i)+1:(config.SemiSpanwiseDiscr(i))*2-1
                ez = [sind(config.RotationAngle_Y(i)); +sind(config.DihedralAngle(i)); cosd(config.RotationAngle_Y(i))*cosd(config.DihedralAngle(i))];      
                ez_norm = norm(ez);
                ez = ez./ez_norm;
                alpha_ind{i}{j} = atan(dot(V_ind{i}{j},ez)/U_Inf_Mag);
                D2D{i}(j) = -L2D{i}(j)*sin(alpha_ind{i}{j});
            end
            for j = (config.SemiSpanwiseDiscr(i))*2:(config.SemiSpanwiseDiscr(i))*2
                if i ~= 1 || (i == 1 && value_input == 2)
                    Dihedral_vec(i,2) = -config.DihedralAngle(i);
                end
                ez = [sind(config.RotationAngle_Y(i)); -sind(Dihedral_vec(i,2)); cosd(config.RotationAngle_Y(i))*cosd(Dihedral_vec(i,2))];
                ez_norm = norm(ez);
                ez = ez./ez_norm;
                alpha_ind{i}{j} = atan(dot(V_ind{i}{j},ez)/U_Inf_Mag);
                if i == 1  && value_input ~= 2
                    D2D{i}(j) = -L2D{i}(j)*sin(alpha_ind{i}{j});
                else
                    D2D{i}(j) = L2D{i}(j)*sin(alpha_ind{i}{j});
                end
            end
    
            coeff = (1/2*rho*U_Inf_Mag^2*corda{i});
            Cd2D{i} = D2D{i}./coeff;
    
            if (value_input == 1 || value_input == 2) && alpha == alpha_vec(1)
                Matrixcontpoints = ControlPoints{i};
                ControlPointsY = zeros(1,(config.SemiSpanwiseDiscr(i))*2);
                for iPan = 1:(config.SemiSpanwiseDiscr(i))*2
                    ControlPointsY(1,iPan) = Matrixcontpoints{1,iPan}.Coords(2);
                end
            %     estremi1 = WingExtremes{iCorpo}(:,1);
            %     estremi2 = WingExtremes{iCorpo}(:,3);
                estremosx = Matrixcontpoints{1,1}.Coords(2);
                estremodx = Matrixcontpoints{1,end}.Coords(2);
                xx = linspace(estremosx,estremodx,(config.SemiSpanwiseDiscr(i))*20);
                   
        
                yy = spline(ControlPointsY,Cd2D{i},xx);
                figure()
                plot(ControlPointsY,Cd2D{i},'o',xx,yy);
                title('curva Cd2D del corpo ',i);
                xlabel('Apertura')
                ylabel('Cd')
                grid on
            end
        end
        
        %% LTOT e DTOT
        Lift = zeros(1,config.NCorpi);
        Drag = zeros(1,config.NCorpi);
        for i = 1:config.NCorpi
            Mesh = internalMesh{i};
            deltaB = zeros(1,(config.SemiSpanwiseDiscr(i))*2);
            for j = 1:(config.SemiSpanwiseDiscr(i))*2
                Pos_LE1 = Mesh{1,j}.LERoot(:);
                Pos_LE2 = Mesh{1,j}.LEtip(:);
                deltaB(j) = norm(Pos_LE1-Pos_LE2);
            end
            Lift(i) = sum(deltaB.*L2D{i});
            Drag(i) = sum(deltaB.*D2D{i});
        
            Cl_tot = Lift./(1/2*rho*config.Surface*U_Inf_Mag^2);
            Cd_tot = Drag./(1/2*rho*config.Surface*U_Inf_Mag^2);

            if (value_input == 1 || value_input == 2)
                Lift_corpo1(I1,I2) = Lift(i);
                Drag_corpo1(I1,I2) = Drag(i);
                Cl_corpo1(I1,I2) = Cl_tot;
                Cd_corpo1(I1,I2) = Cd_tot;
            elseif (value_input ~= 1 && value_input ~= 2)
                Lift_corpo1(I1,I2) = Lift(1);
                Drag_corpo1(I1,I2) = Drag(1);

                Lift_corpo2(I1,I2) = Lift(2);
                Drag_corpo2(I1,I2) = Drag(2);

                Cl_corpo2(I1,I2) = Cl_tot(2);
                Cl_corpo1(I1,I2) = Cl_tot(1);

                Cd_corpo2(I1,I2) = Cd_tot(2);
                Cd_corpo1(I1,I2) = Cd_tot(1);

                Efficienza1(I1,I2) = Cl_corpo1(I1,I2)/Cd_corpo1(I1,I2);
                Efficienza2(I1,I2) = Cl_corpo2(I1,I2)/Cd_corpo2(I1,I2);
            end
        end
    end
end

%% PLOT Cl-alpha
figure()
for i = 1:length(varing)
    plot(alpha_vec, Cl_corpo1(i,:),'-o');
    hold on
    if value_input ~= 1 && value_input ~= 2
        legendTexts{i} = ['Variabile primo corpo = ', num2str(varing(i))];
        legend(legendTexts)
    end
end
grid on
xlabel('alpha')
ylabel('Cl')
title('Curva Cl-alpha corpo 1')

if value_input ~= 1 && value_input ~= 2
    figure()
    for i = 1:length(varing)
        plot(alpha_vec, Cl_corpo2(i,:),'-o');
        hold on
    end
    legend(legendTexts)
    grid on
    xlabel('alpha')
    ylabel('Cl')
    title('Curva Cl-alpha corpo 2')
end

%% PLOT polare
figure()
for i = 1:length(varing)
    plot(Cd_corpo1(i,:), Cl_corpo1(i,:),'-o');
    hold on
end
if value_input ~= 1 && value_input ~= 2
    legend(legendTexts)
end
grid on
xlabel('Cd')
ylabel('Cl')
title('Curva polare corpo 1')

if value_input ~= 1 && value_input ~= 2
    figure()
    for i = 1:length(varing)
        plot(Cd_corpo2(i,:), Cl_corpo2(i,:),'-o');
        hold on
    end
    if value_input ~= 1 && value_input ~= 2
        legend(legendTexts)
    end
    grid on
    xlabel('Cd')
    ylabel('Cl')
    title('Curva polare corpo 2')
end

    %% PLOT
if (value_input ~= 1 && value_input ~= 2)
    figure()
    for iCorpo = 1:config.NCorpi
    
        Mesh = internalMesh{iCorpo};
        [alt,lung] = size(Mesh);
    
        MPointsX = zeros(alt,lung);
        MPointsY = zeros(alt,lung);
        MPointsZ = zeros(alt,lung);
    
        for iPan = 1:alt
            A1 = 0;
            for jPan = 1:lung/2
                A1 = A1+1;
                MPointsX(iPan,jPan) = Mesh{iPan,jPan}.LEtip(1);
                MPointsY(iPan,jPan) = Mesh{iPan,jPan}.LEtip(2);
                MPointsZ(iPan,jPan) = Mesh{iPan,jPan}.LEtip(3);
            
            end
            A1 = A1+1;
            MPointsX(iPan,A1) = Mesh{iPan,jPan}.LERoot(1);
            MPointsY(iPan,A1) = Mesh{iPan,jPan}.LERoot(2);
            MPointsZ(iPan,A1) = Mesh{iPan,jPan}.LERoot(3);
            for jPan = lung/2+1:lung
                A1 = A1+1;
                MPointsX(iPan,A1) = Mesh{iPan,jPan}.LEtip(1);
                MPointsY(iPan,A1) = Mesh{iPan,jPan}.LEtip(2);
                MPointsZ(iPan,A1) = Mesh{iPan,jPan}.LEtip(3);
            end
        end
        A1 = 0;
        for jPan = 1:lung/2
            A1 = A1+1;
            MPointsX(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(1);
            MPointsY(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(2);
            MPointsZ(iPan+1,jPan) = Mesh{iPan,jPan}.TEtip(3);
    
        end
        A1 = A1+1;
        MPointsX(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(1);
        MPointsY(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(2);
        MPointsZ(iPan+1,A1) = Mesh{iPan,jPan}.TERoot(3);
        for jPan = lung/2+1:lung
            A1 = A1+1;
            MPointsX(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(1);
            MPointsY(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(2);
            MPointsZ(iPan+1,A1) = Mesh{iPan,jPan}.TEtip(3);
        end
    
        surf(MPointsX,MPointsY,MPointsZ,'FaceColor','interp')
        hold on
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        colorbar;
        grid on;
        title('Profili alari');
        axis equal
    end
end

% %% Plot lungo l'apertura
% if value_input == 1|| value_input == 2
%     Matrixcontpoints = ControlPoints{iCorpo};
%     ControlPointsY = zeros(1,(config.SemiSpanwiseDiscr(iCorpo))*2);
%     for iPan = 1:(config.SemiSpanwiseDiscr(iCorpo))*2
%         ControlPointsY(1,iPan) = Matrixcontpoints{1,iPan}.Coords(2);
%     end
% %     estremi1 = WingExtremes{iCorpo}(:,1);
% %     estremi2 = WingExtremes{iCorpo}(:,3);
% %     estremosx = estremi1{3}.LE(2);
% %     estremodx = estremi2{3}.LE(2); 
% 
%     estremosx = Matrixcontpoints{1,1}.Coords(2);
%     estremodx = Matrixcontpoints{1,end}.Coords(2);
%     xx = linspace(estremosx,estremodx,(config.SemiSpanwiseDiscr(iCorpo))*20);
% 
%     yy = spline(ControlPointsY,Cl2D{iCorpo},xx);
%     figure()
%     plot(ControlPointsY,Cl2D{iCorpo},'o',xx,yy);
%     title('curva Cl2D del corpo ',iCorpo);
%     xlabel('Profilo')
%     ylabel('Cl')
%     grid on
% 
%     yy = spline(ControlPointsY,Cd2D{iCorpo},xx);
%     figure()
%     plot(ControlPointsY,Cd2D{iCorpo},'o',xx,yy);
%     title('curva Cd2D del corpo ',iCorpo);
%     xlabel('Profilo')
%     ylabel('Cd')
%     grid on
% end