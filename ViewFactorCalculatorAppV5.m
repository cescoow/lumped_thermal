% Monte Carlo based View Factor calculator. Algorithm is based on ray
% tracing and ray intersection between emissed rays and target part (or
% world). Parts are loaded from STL files. Code is in experimental status!
% Validation shall be made with other tools (i.e. Fluent)
% Created by Francesco Lena. October 2024
% Any issues please contact me: francescorossilena@gmail.com

classdef ViewFactorCalculatorAppV5 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        GridLayout                 matlab.ui.container.GridLayout
        ControlPanel               matlab.ui.container.Panel
        LoadBodyAButton            matlab.ui.control.Button
        LoadBodyBButton            matlab.ui.control.Button
        LoadObstacleButton         matlab.ui.control.Button
        NumRaysLabel               matlab.ui.control.Label
        NumRaysEditField           matlab.ui.control.NumericEditField
        CalculateButton            matlab.ui.control.Button
        ViewFactorLabel            matlab.ui.control.Label
        ViewFactorObstacleLabel    matlab.ui.control.Label
        EstimatedTimeLabel         matlab.ui.control.Label
        BodyAFileLabel             matlab.ui.control.Label
        BodyBFileLabel             matlab.ui.control.Label
        ObstacleFileLabel          matlab.ui.control.Label
        ProgressLabel              matlab.ui.control.Label
        ShowRaysCheckBox           matlab.ui.control.CheckBox
        SaveDataButton             matlab.ui.control.Button
        LoadDataButton             matlab.ui.control.Button
        UIAxes                     matlab.ui.control.UIAxes
        CalculationModeDropDown    matlab.ui.control.DropDown  
        WatermarkLabel             matlab.ui.control.Label    
    end

    properties (Access = private)
        TR_A % Triangulation for Body A
        TR_B % Triangulation for Body B
        TR_Obstacle % Triangulation for Obstacle
        areasA
        cumulativeAreasA
        totalAreaA
        trianglesA
        numRays % Number of rays for Monte Carlo simulation
        isBodyALoaded = false
        isBodyBLoaded = false
        isObstacleLoaded = false
        rayOrigins   % Origins of rays for visualization
        rayEndpoints % Endpoints of rays for visualization
        rayStatuses  % Status of rays: 0 = not plotted, 1 = hit B, 2 = hit Obstacle, 3 = blocked by A
        viewFactor   % Calculated view factor between A and B
        viewFactorObstacle % Calculated view factor between A and Obstacle
        calculationMode = 'Between 2 Parts' % Default Property
    end

    methods (Access = private)

        % Button pushed function: LoadBodyAButton
        function LoadBodyAButtonPushed(app, ~)
            [file, path] = uigetfile('*.stl', 'Select STL file for Body A');
            if isequal(file, 0)
                return;
            end
            filename = fullfile(path, file);
            app.BodyAFileLabel.Text = sprintf('Body A File: %s', file);
            app.TR_A = stlread(filename);
            app.isBodyALoaded = true;
            app.prepareBodyA();
            app.updatePlot();
        end

        % Button pushed function: LoadBodyBButton
        function LoadBodyBButtonPushed(app, ~)
            [file, path] = uigetfile('*.stl', 'Select STL file for Body B');
            if isequal(file, 0)
                return;
            end
            filename = fullfile(path, file);
            app.BodyBFileLabel.Text = sprintf('Body B File: %s', file);
            app.TR_B = stlread(filename);
            app.isBodyBLoaded = true;
            app.updatePlot();
        end

        % Button pushed function: LoadObstacleButton
        function LoadObstacleButtonPushed(app, ~)
            [file, path] = uigetfile('*.stl', 'Select STL file for Obstacle');
            if isequal(file, 0)
                return;
            end
            filename = fullfile(path, file);
            app.ObstacleFileLabel.Text = sprintf('Obstacle File: %s', file);
            app.TR_Obstacle = stlread(filename);
            app.isObstacleLoaded = true;
            app.updatePlot();
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, ~)
            % Validate loaded bodies based on calculation mode
            if ~app.isBodyALoaded
                uialert(app.UIFigure, 'Please load Body A STL file before calculating.', 'Error');
                return;
            end
            switch app.calculationMode
                case 'Between 2 Parts'
                    if ~app.isBodyBLoaded
                        uialert(app.UIFigure, 'Please load Body B STL file before calculating.', 'Error');
                        return;
                    end
                case 'Between 2 Parts with Obstacle'
                    if ~app.isBodyBLoaded
                        uialert(app.UIFigure, 'Please load Body B STL file before calculating.', 'Error');
                        return;
                    end
                    if ~app.isObstacleLoaded
                        uialert(app.UIFigure, 'Please load Obstacle STL file before calculating.', 'Error');
                        return;
                    end
                otherwise
                    % Between Part and World requires only Body A
            end
            app.numRays = app.NumRaysEditField.Value;
            app.CalculateButton.Enable = 'off';
            app.ViewFactorLabel.Text = 'View Factor (A to B): Calculating...';
            app.ViewFactorObstacleLabel.Text = 'View Factor (A to Obstacle): Calculating...';
            drawnow;

            % Estimate time (rough estimation)
            estimatedTime = app.numRays * 0.001; % Adjust this multiplier based on performance
            app.EstimatedTimeLabel.Text = sprintf('Estimated Time: %.2f seconds', estimatedTime);
            drawnow;

            % Initialize rays data if visualization is enabled
            if app.ShowRaysCheckBox.Value
                maxRaysToStore = 1000; % Limit the number of rays to store for visualization
                raysStored = 0;
                app.rayOrigins = zeros(maxRaysToStore, 3);
                app.rayEndpoints = zeros(maxRaysToStore, 3);
                app.rayStatuses = zeros(maxRaysToStore, 1); % 1 = hit B, 2 = hit Obstacle, 3 = blocked by A
            else
                app.rayOrigins = [];
                app.rayEndpoints = [];
                app.rayStatuses = [];
            end

            % Run the Monte Carlo simulation
            hitCount = 0;
            hitObstacleCount = 0;
            app.ProgressLabel.Text = 'Progress: 0%';
            tic;
            for rayIndex = 1:app.numRays
                % Sample a random point and normal on Body A
                [pointA, normalA] = app.sampleRandomPointOnSurface();

                % Generate a random direction weighted by the Lambertian cosine law
                dir = app.cosineWeightedDirection(normalA);

                % Check for occlusions with Body A itself (Self-Occlusion)
                [occludedByA, ~] = app.rayMeshIntersect(pointA, dir, app.TR_A.Points, app.TR_A.ConnectivityList, 1e-6);

                if occludedByA
                    % Ray is blocked by Body A
                    if app.ShowRaysCheckBox.Value && raysStored < maxRaysToStore
                        raysStored = raysStored + 1;
                        app.rayOrigins(raysStored, :) = pointA;
                        % Extend the ray to the occlusion point on Body A
                        app.rayEndpoints(raysStored, :) = pointA + dir * 1e-3; % Small step to visualize blockage
                        app.rayStatuses(raysStored) = 3; % Blocked by A
                    end
                    continue; % Do not count this ray
                end

                if strcmp(app.calculationMode, 'Between 2 Parts')
                    % Calculate view factor between two parts
                    % Cast the ray and check for intersection with Body B
                    [intersectB, t_B] = app.rayMeshIntersect(pointA, dir, app.TR_B.Points, app.TR_B.ConnectivityList);

                    % If the ray hits Body B, count it
                    if intersectB
                        hitCount = hitCount + 1;
                        % Store the ray data if visualization is enabled
                        if app.ShowRaysCheckBox.Value && raysStored < maxRaysToStore
                            raysStored = raysStored + 1;
                            app.rayOrigins(raysStored, :) = pointA;
                            % Use the actual intersection point with Body B
                            app.rayEndpoints(raysStored, :) = pointA + dir * t_B;
                            app.rayStatuses(raysStored) = 1; % Hit B
                        end
                    end

                elseif strcmp(app.calculationMode, 'Between 2 Parts with Obstacle')
                    % Cast the ray and check for intersection with Obstacle first
                    [intersectObstacle, t_Obstacle] = app.rayMeshIntersect(pointA, dir, app.TR_Obstacle.Points, app.TR_Obstacle.ConnectivityList);

                    if intersectObstacle
                        % The ray is blocked by the Obstacle
                        hitObstacleCount = hitObstacleCount + 1;
                        if app.ShowRaysCheckBox.Value && raysStored < maxRaysToStore
                            raysStored = raysStored + 1;
                            app.rayOrigins(raysStored, :) = pointA;
                            app.rayEndpoints(raysStored, :) = pointA + dir * t_Obstacle;
                            app.rayStatuses(raysStored) = 2; % Hit Obstacle
                        end
                        continue; % Do not count this ray towards F_A_to_B
                    end

                    % If not blocked by Obstacle, check intersection with Body B
                    [intersectB, t_B] = app.rayMeshIntersect(pointA, dir, app.TR_B.Points, app.TR_B.ConnectivityList);

                    if intersectB
                        hitCount = hitCount + 1;
                        % Store the ray data if visualization is enabled
                        if app.ShowRaysCheckBox.Value && raysStored < maxRaysToStore
                            raysStored = raysStored + 1;
                            app.rayOrigins(raysStored, :) = pointA;
                            app.rayEndpoints(raysStored, :) = pointA + dir * t_B;
                            app.rayStatuses(raysStored) = 1; % Hit B
                        end
                    end

                else
                    % Calculate view factor between part and world
                    % If the ray is not occluded by Body A, it reaches the world
                    hitCount = hitCount + 1;
                    % Store the ray data if visualization is enabled
                    if app.ShowRaysCheckBox.Value && raysStored < maxRaysToStore
                        raysStored = raysStored + 1;
                        app.rayOrigins(raysStored, :) = pointA;
                        % Extend the ray to a large distance for visualization
                        app.rayEndpoints(raysStored, :) = pointA + dir * 100; % Arbitrary large distance
                        app.rayStatuses(raysStored) = 1; % Reaches World
                    end
                end

                % Update progress every 1% of completion
                if mod(rayIndex, max(floor(app.numRays / 100),1)) == 0
                    progressPercent = (rayIndex / app.numRays) * 100;
                    app.ProgressLabel.Text = sprintf('Progress: %.0f%%', progressPercent);
                    drawnow;
                end
            end
            
            progressPercent = (rayIndex / app.numRays) * 100;
            app.ProgressLabel.Text = sprintf('Progress: %.0f%%', progressPercent);
            drawnow;
                    
            elapsedTime = toc;
            viewFactor = hitCount / app.numRays;
            app.viewFactor = viewFactor; % Store the calculated view factor

            if strcmp(app.calculationMode, 'Between 2 Parts with Obstacle')
                viewFactorObstacle = hitObstacleCount / app.numRays;
                app.viewFactorObstacle = viewFactorObstacle; % Store the calculated view factor to Obstacle
                app.ViewFactorObstacleLabel.Text = sprintf('View Factor (A to Obstacle): %.6f', viewFactorObstacle);
            else
                app.viewFactorObstacle = NaN;
                app.ViewFactorObstacleLabel.Text = 'View Factor (A to Obstacle): N/A';
            end

            app.ViewFactorLabel.Text = sprintf('View Factor (A to B): %.6f', viewFactor);
            app.EstimatedTimeLabel.Text = sprintf('Elapsed Time: %.2f seconds', elapsedTime);
            app.CalculateButton.Enable = 'on';

            % Trim the rays arrays if needed
            if app.ShowRaysCheckBox.Value && raysStored > 0
                app.rayOrigins = app.rayOrigins(1:raysStored, :);
                app.rayEndpoints = app.rayEndpoints(1:raysStored, :);
                app.rayStatuses = app.rayStatuses(1:raysStored);
            else
                app.rayOrigins = [];
                app.rayEndpoints = [];
                app.rayStatuses = [];
            end

            % Update the plot to show the rays if visualization is enabled
            app.updatePlot();
        end

        function prepareBodyA(app)
            % Preprocess surface of Body A
            [app.areasA, app.cumulativeAreasA, app.totalAreaA, app.trianglesA] = app.preprocessSurface(app.TR_A.Points, app.TR_A.ConnectivityList);
        end

        function updatePlot(app)
            % Clear axes
            cla(app.UIAxes);
            hold(app.UIAxes, 'on');
            axis(app.UIAxes, 'equal');
            xlabel(app.UIAxes, 'X');
            ylabel(app.UIAxes, 'Y');
            zlabel(app.UIAxes, 'Z');
            view(app.UIAxes, 3);
            grid(app.UIAxes, 'on');

            % Plot Body A
            if app.isBodyALoaded
                patch(app.UIAxes, 'Faces', app.TR_A.ConnectivityList, 'Vertices', app.TR_A.Points, ...
                    'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Body A');
            end

            % Plot Body B
            if strcmp(app.calculationMode, 'Between 2 Parts') && app.isBodyBLoaded
                patch(app.UIAxes, 'Faces', app.TR_B.ConnectivityList, 'Vertices', app.TR_B.Points, ...
                    'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'DisplayName', 'Body B');
            elseif strcmp(app.calculationMode, 'Between 2 Parts with Obstacle') && app.isBodyBLoaded
                patch(app.UIAxes, 'Faces', app.TR_B.ConnectivityList, 'Vertices', app.TR_B.Points, ...
                    'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.25, 'DisplayName', 'Body B');
            end

            % Plot Obstacle
            if strcmp(app.calculationMode, 'Between 2 Parts with Obstacle') && app.isObstacleLoaded
                patch(app.UIAxes, 'Faces', app.TR_Obstacle.ConnectivityList, 'Vertices', app.TR_Obstacle.Points, ...
                    'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'DisplayName', 'Obstacle');
            end

            % Plot the rays if visualization is enabled
            if app.ShowRaysCheckBox.Value && ~isempty(app.rayOrigins)
                % Separate rays based on their status
                raysHitB = app.rayStatuses == 1;
                raysHitObstacle = app.rayStatuses == 2;
                raysBlockedByA = app.rayStatuses == 3;

                % Plot rays that hit B or reach the world
                if any(raysHitB)
                    X = [app.rayOrigins(raysHitB,1), app.rayEndpoints(raysHitB,1)]';
                    Y = [app.rayOrigins(raysHitB,2), app.rayEndpoints(raysHitB,2)]';
                    Z = [app.rayOrigins(raysHitB,3), app.rayEndpoints(raysHitB,3)]';
                    plot3(app.UIAxes, X, Y, Z, 'g-', 'LineWidth', 1, 'DisplayName', 'Hit B / World');
                end

                % Plot rays that hit the Obstacle
                if any(raysHitObstacle)
                    X = [app.rayOrigins(raysHitObstacle,1), app.rayEndpoints(raysHitObstacle,1)]';
                    Y = [app.rayOrigins(raysHitObstacle,2), app.rayEndpoints(raysHitObstacle,2)]';
                    Z = [app.rayOrigins(raysHitObstacle,3), app.rayEndpoints(raysHitObstacle,3)]';
                    plot3(app.UIAxes, X, Y, Z, 'r-', 'LineWidth', 1, 'DisplayName', 'Hit Obstacle');
                end

                % Plot rays that are blocked by Body A (Self-Occlusion)
                if any(raysBlockedByA)
                    X = [app.rayOrigins(raysBlockedByA,1), app.rayEndpoints(raysBlockedByA,1)]';
                    Y = [app.rayOrigins(raysBlockedByA,2), app.rayEndpoints(raysBlockedByA,2)]';
                    Z = [app.rayOrigins(raysBlockedByA,3), app.rayEndpoints(raysBlockedByA,3)]';
                    plot3(app.UIAxes, X, Y, Z, 'm-', 'LineWidth', 1, 'DisplayName', 'Blocked by A');
                end
            end
            camlight(app.UIAxes);
            lighting(app.UIAxes, 'gouraud');
            hold(app.UIAxes, 'off');
        end

        % Value changed function: ShowRaysCheckBox
        function ShowRaysCheckBoxValueChanged(app, ~)
            % Update the plot when the checkbox value changes
            app.updatePlot();
        end

        % Button pushed function: SaveDataButton
        function SaveDataButtonPushed(app, ~)
            [file, path] = uiputfile('*.mat', 'Save Simulation Data');
            if isequal(file, 0)
                return;
            end
            filename = fullfile(path, file);
            % Save relevant data
            TR_A = app.TR_A;
            TR_B = app.TR_B;
            TR_Obstacle = app.TR_Obstacle;
            numRays = app.numRays;
            viewFactor = app.viewFactor;
            viewFactorObstacle = app.viewFactorObstacle;
            rayOrigins = app.rayOrigins;
            rayEndpoints = app.rayEndpoints;
            rayStatuses = app.rayStatuses;
            calculationMode = app.calculationMode;
            save(filename, 'TR_A', 'TR_B', 'TR_Obstacle', 'numRays', 'viewFactor', 'viewFactorObstacle', 'rayOrigins', 'rayEndpoints', 'rayStatuses', 'calculationMode');
            uialert(app.UIFigure, 'Simulation data saved successfully.', 'Save Data');
        end

        % Button pushed function: LoadDataButton
        function LoadDataButtonPushed(app, ~)
            [file, path] = uigetfile('*.mat', 'Load Simulation Data');
            if isequal(file, 0)
                return;
            end
            filename = fullfile(path, file);
            % Load data
            data = load(filename);
            app.TR_A = data.TR_A;
            app.isBodyALoaded = true;
            app.BodyAFileLabel.Text = 'Body A File: Loaded from data';

            if isfield(data, 'TR_B') && ~isempty(data.TR_B)
                app.TR_B = data.TR_B;
                app.isBodyBLoaded = true;
                app.BodyBFileLabel.Text = 'Body B File: Loaded from data';
            else
                app.TR_B = [];
                app.isBodyBLoaded = false;
                app.BodyBFileLabel.Text = 'Body B File: None';
            end

            if isfield(data, 'TR_Obstacle') && ~isempty(data.TR_Obstacle)
                app.TR_Obstacle = data.TR_Obstacle;
                app.isObstacleLoaded = true;
                app.ObstacleFileLabel.Text = 'Obstacle File: Loaded from data';
            else
                app.TR_Obstacle = [];
                app.isObstacleLoaded = false;
                app.ObstacleFileLabel.Text = 'Obstacle File: None';
            end

            if isfield(data, 'numRays')
                app.numRays = data.numRays;
                app.NumRaysEditField.Value = app.numRays;
            else
                app.numRays = 1000;
                app.NumRaysEditField.Value = 1000;
            end

            if isfield(data, 'viewFactor')
                app.viewFactor = data.viewFactor;
                app.ViewFactorLabel.Text = sprintf('View Factor (A to B): %.6f', app.viewFactor);
            else
                app.viewFactor = NaN;
                app.ViewFactorLabel.Text = 'View Factor (A to B): ';
            end

            if isfield(data, 'viewFactorObstacle')
                app.viewFactorObstacle = data.viewFactorObstacle;
                if strcmp(app.calculationMode, 'Between 2 Parts with Obstacle')
                    app.ViewFactorObstacleLabel.Text = sprintf('View Factor (A to Obstacle): %.6f', app.viewFactorObstacle);
                else
                    if strcmp(app.calculationMode, 'Between Part and World')
                        app.ViewFactorObstacleLabel.Text = sprintf('View Factor (A to World): %.6f', app.viewFactorObstacle);
                    else
                    app.ViewFactorObstacleLabel.Text = 'View Factor (A to Obstacle): N/A';
                    end
                end
            else
                app.viewFactorObstacle = NaN;
                app.ViewFactorObstacleLabel.Text = 'View Factor (A to Obstacle): N/A';
            end

            if isfield(data, 'rayOrigins') && isfield(data, 'rayEndpoints') && isfield(data, 'rayStatuses')
                app.rayOrigins = data.rayOrigins;
                app.rayEndpoints = data.rayEndpoints;
                app.rayStatuses = data.rayStatuses;
            else
                app.rayOrigins = [];
                app.rayEndpoints = [];
                app.rayStatuses = [];
            end

            if isfield(data, 'calculationMode')
                app.calculationMode = data.calculationMode;
                app.CalculationModeDropDown.Value = app.calculationMode;
            else
                app.calculationMode = 'Between 2 Parts';
                app.CalculationModeDropDown.Value = 'Between 2 Parts';
            end

            % Prepare Body A
            app.prepareBodyA();

            % Update UI based on calculation mode
            app.updateCalculationModeUI();

            % Update plot
            app.updatePlot();
            uialert(app.UIFigure, 'Simulation data loaded successfully.', 'Load Data');
        end

        % Preprocess surface function
        function [areas, cumulativeAreas, totalArea, triangles] = preprocessSurface(app, vertices, faces)
            numFaces = size(faces, 1);
            areas = zeros(numFaces, 1);
            triangles(numFaces).v1 = [];
            triangles(numFaces).v2 = [];
            triangles(numFaces).v3 = [];
            triangles(numFaces).normal = [];

            for i = 1:numFaces
                vert_idx = faces(i, :);
                v1 = vertices(vert_idx(1), :);
                v2 = vertices(vert_idx(2), :);
                v3 = vertices(vert_idx(3), :);

                % Compute area
                area = 0.5 * norm(cross(v2 - v1, v3 - v1));
                areas(i) = area;

                % Store triangle data
                triangles(i).v1 = v1;
                triangles(i).v2 = v2;
                triangles(i).v3 = v3;
                normal = cross(v2 - v1, v3 - v1);
                if norm(normal) ~= 0
                    normal = normal / norm(normal);
                end
                triangles(i).normal = normal;
            end
            cumulativeAreas = cumsum(areas);
            totalArea = cumulativeAreas(end);
        end

        % Sample random point on surface function
        function [point, normal] = sampleRandomPointOnSurface(app)
            % Randomly select a triangle based on area weighting
            r = rand() * app.totalAreaA;
            triangleIndex = find(app.cumulativeAreasA >= r, 1);
            triangle = app.trianglesA(triangleIndex);

            % Generate random barycentric coordinates
            sqrt_r1 = sqrt(rand());
            r2 = rand();
            u = 1 - sqrt_r1;
            v = sqrt_r1 * (1 - r2);
            w = sqrt_r1 * r2;

            % Compute the random point on the triangle
            point = u * triangle.v1 + v * triangle.v2 + w * triangle.v3;
            normal = triangle.normal;
        end

        % Cosine weighted direction function
        function dir = cosineWeightedDirection(app, normal)
            % Generate a random direction in the hemisphere around the normal
            % using Lambertian cosine weighting
            % Generate random angles
            r1 = rand();
            r2 = rand();
            theta = acos(sqrt(1 - r1));
            phi = 2 * pi * r2;

            % Spherical to Cartesian coordinates
            x = sin(theta) * cos(phi);
            y = sin(theta) * sin(phi);
            z = cos(theta);

            % Create a coordinate system
            [T, B] = app.buildOrthonormalBasis(normal);
            dir = x * T + y * B + z * normal;
            dir = dir / norm(dir); % Ensure the direction is normalized
        end

        % Build orthonormal basis function
        function [T, B] = buildOrthonormalBasis(app, N)
            % Build an orthonormal basis (Tangent, Bitangent, Normal)
            if abs(N(1)) > abs(N(2))
                invLen = 1 / sqrt(N(1)^2 + N(3)^2);
                T = [-N(3) * invLen, 0, N(1) * invLen];
            else
                invLen = 1 / sqrt(N(2)^2 + N(3)^2);
                T = [0, N(3) * invLen, -N(2) * invLen];
            end
            B = cross(N, T);
        end

        % Ray-mesh intersection function
        function [intersect, t_min] = rayMeshIntersect(app, rayOrigin, rayDirection, vertices, faces, epsilon)
            if nargin < 6
                epsilon = 1e-6;
            end
            intersect = false;
            t_min = inf;

            % Use custom ray-triangle intersection function
            [intersects, t] = app.intersectRayTriangleMesh(rayOrigin, rayDirection, vertices, faces);
            if any(intersects)
                t_candidates = t(intersects);
                t_positive = t_candidates(t_candidates > epsilon);
                if ~isempty(t_positive)
                    t_min = min(t_positive);
                    intersect = true;
                end
            end
        end

        % Custom ray-triangle intersection function using Möller–Trumbore algorithm
        function [intersects, t] = intersectRayTriangleMesh(app, rayOrigin, rayDirection, vertices, faces)
            numFaces = size(faces, 1);
            intersects = false(numFaces, 1);
            t = inf(numFaces, 1);

            % Triangle vertices
            V0 = vertices(faces(:, 1), :);
            V1 = vertices(faces(:, 2), :);
            V2 = vertices(faces(:, 3), :);

            % Edges
            edge1 = V1 - V0;
            edge2 = V2 - V0;

            % Begin calculations
            h = cross(repmat(rayDirection, numFaces, 1), edge2, 2);
            a = dot(edge1, h, 2);
            valid = abs(a) > 1e-8;
            f = zeros(numFaces,1);
            f(valid) = 1 ./ a(valid);
            s = repmat(rayOrigin, numFaces, 1) - V0;
            u = zeros(numFaces,1);
            u(valid) = f(valid) .* dot(s(valid, :), h(valid, :), 2);
            q = cross(s(valid, :), edge1(valid, :), 2);
            v = zeros(numFaces,1);
            v(valid) = f(valid) .* dot(repmat(rayDirection, sum(valid), 1), q, 2);
            temp_t = zeros(numFaces,1);
            temp_t(valid) = f(valid) .* dot(edge2(valid, :), q, 2);

            % Check for valid intersections
            valid_intersections = valid & (u >= 0) & (u <= 1) & (v >= 0) & ((u + v) <= 1) & (temp_t > 0);
            intersects(valid_intersections) = true;
            t(valid_intersections) = temp_t(valid_intersections);
        end

        % Value changed function: CalculationModeDropDown
        function CalculationModeDropDownValueChanged(app, event)
            app.calculationMode = app.CalculationModeDropDown.Value;
            app.updateCalculationModeUI();
        end

        % Update UI elements based on calculation mode
        function updateCalculationModeUI(app)
            switch app.calculationMode
                case 'Between 2 Parts'
                    app.LoadBodyBButton.Enable = 'on';
                    app.BodyBFileLabel.Enable = 'on';
                    app.LoadObstacleButton.Enable = 'off';
                    app.ObstacleFileLabel.Enable = 'off';
                    app.ObstacleFileLabel.Text = 'Obstacle File: N/A';
                    app.ViewFactorObstacleLabel.Visible = 'off';
                case 'Between 2 Parts with Obstacle'
                    app.LoadBodyBButton.Enable = 'on';
                    app.BodyBFileLabel.Enable = 'on';
                    app.LoadObstacleButton.Enable = 'on';
                    app.ObstacleFileLabel.Enable = 'on';
                    if app.isObstacleLoaded
                        % Extract filename from label text
                        parts = split(app.ObstacleFileLabel.Text, ': ');
                        if numel(parts) == 2
                            filename = parts{2};
                            app.ObstacleFileLabel.Text = sprintf('Obstacle File: %s', filename);
                        end
                    else
                        app.ObstacleFileLabel.Text = 'Obstacle File: None';
                    end
                    app.ViewFactorObstacleLabel.Visible = 'on';
                otherwise
                    % Between Part and World
                    app.LoadBodyBButton.Enable = 'off';
                    app.BodyBFileLabel.Enable = 'off';
                    app.BodyBFileLabel.Text = 'Body B File: N/A';
                    app.LoadObstacleButton.Enable = 'off';
                    app.ObstacleFileLabel.Enable = 'off';
                    app.ObstacleFileLabel.Text = 'Obstacle File: N/A';
                    app.ViewFactorObstacleLabel.Visible = 'off';
            end
            app.updatePlot();
        end

    end

    % Component initialization
    methods (Access = private)

        % Create components
     function createComponents(app)

        % Create UIFigure and components
        app.UIFigure = uifigure('Visible', 'off');
        app.UIFigure.Position = [100 100 1000 700];
        app.UIFigure.Name = 'Monte Carlo View Factor Calculator V5';

        % Create GridLayout
        app.GridLayout = uigridlayout(app.UIFigure, [1, 2]);
        app.GridLayout.ColumnWidth = {'1.2x', '2x'};
        app.GridLayout.RowHeight = {'1x'};

        % Create ControlPanel
        app.ControlPanel = uipanel(app.GridLayout);
        app.ControlPanel.Title = 'Controls';
        app.ControlPanel.FontSize = 14;
        app.ControlPanel.Layout.Row = 1;
        app.ControlPanel.Layout.Column = 1;

        % Create CalculationModeDropDown
        app.CalculationModeDropDown = uidropdown(app.ControlPanel);
        app.CalculationModeDropDown.Items = {'Between 2 Parts', 'Between Part and World', 'Between 2 Parts with Obstacle'};
        app.CalculationModeDropDown.ValueChangedFcn = createCallbackFcn(app, @CalculationModeDropDownValueChanged, true);
        app.CalculationModeDropDown.Position = [20 600 200 22]; % Y position lowered by 150
        app.CalculationModeDropDown.Value = 'Between 2 Parts';

        % Create LoadBodyAButton
        app.LoadBodyAButton = uibutton(app.ControlPanel, 'push');
        app.LoadBodyAButton.ButtonPushedFcn = createCallbackFcn(app, @LoadBodyAButtonPushed, true);
        app.LoadBodyAButton.Position = [20 550 150 30]; % Y position lowered by 150
        app.LoadBodyAButton.Text = 'Load Body A';

        % Create LoadBodyBButton
        app.LoadBodyBButton = uibutton(app.ControlPanel, 'push');
        app.LoadBodyBButton.ButtonPushedFcn = createCallbackFcn(app, @LoadBodyBButtonPushed, true);
        app.LoadBodyBButton.Position = [20 510 150 30]; % Y position lowered by 150
        app.LoadBodyBButton.Text = 'Load Body B';

        % Create LoadObstacleButton
        app.LoadObstacleButton = uibutton(app.ControlPanel, 'push');
        app.LoadObstacleButton.ButtonPushedFcn = createCallbackFcn(app, @LoadObstacleButtonPushed, true);
        app.LoadObstacleButton.Position = [20 470 150 30]; % Y position lowered by 150
        app.LoadObstacleButton.Text = 'Load Obstacle';
        app.LoadObstacleButton.Enable = 'off'; % Initially disabled

        % Create BodyAFileLabel
        app.BodyAFileLabel = uilabel(app.ControlPanel);
        app.BodyAFileLabel.Position = [180 550 300 22]; % Y position lowered by 150
        app.BodyAFileLabel.Text = 'Body A File: None';

        % Create BodyBFileLabel
        app.BodyBFileLabel = uilabel(app.ControlPanel);
        app.BodyBFileLabel.Position = [180 510 300 22]; % Y position lowered by 150
        app.BodyBFileLabel.Text = 'Body B File: None';

        % Create ObstacleFileLabel
        app.ObstacleFileLabel = uilabel(app.ControlPanel);
        app.ObstacleFileLabel.Position = [180 470 300 22]; % Y position lowered by 150
        app.ObstacleFileLabel.Text = 'Obstacle File: N/A';
        app.ObstacleFileLabel.Enable = 'off'; % Initially disabled

        % Create NumRaysLabel
        app.NumRaysLabel = uilabel(app.ControlPanel);
        app.NumRaysLabel.HorizontalAlignment = 'right';
        app.NumRaysLabel.Position = [20 420 60 22]; % Y position lowered by 150
        app.NumRaysLabel.Text = 'Num Rays';

        % Create NumRaysEditField
        app.NumRaysEditField = uieditfield(app.ControlPanel, 'numeric');
        app.NumRaysEditField.Limits = [1 Inf];
        app.NumRaysEditField.RoundFractionalValues = 'on';
        app.NumRaysEditField.ValueDisplayFormat = '%.0f';
        app.NumRaysEditField.Position = [95 420 100 22]; % Y position lowered by 150
        app.NumRaysEditField.Value = 1000;

        % Create ShowRaysCheckBox
        app.ShowRaysCheckBox = uicheckbox(app.ControlPanel);
        app.ShowRaysCheckBox.Text = 'Show Rays';
        app.ShowRaysCheckBox.Position = [20 380 100 22]; % Y position lowered by 150
        app.ShowRaysCheckBox.Value = false;
        app.ShowRaysCheckBox.ValueChangedFcn = createCallbackFcn(app, @ShowRaysCheckBoxValueChanged, true);

        % Create CalculateButton
        app.CalculateButton = uibutton(app.ControlPanel, 'push');
        app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
        app.CalculateButton.Position = [20 340 150 30]; % Y position lowered by 150
        app.CalculateButton.Text = 'Calculate';

        % Create ViewFactorLabel
        app.ViewFactorLabel = uilabel(app.ControlPanel);
        app.ViewFactorLabel.Position = [20 300 300 22]; % Y position lowered by 150
        app.ViewFactorLabel.Text = 'View Factor (A to B): ';

        % Create ViewFactorObstacleLabel
        app.ViewFactorObstacleLabel = uilabel(app.ControlPanel);
        app.ViewFactorObstacleLabel.Position = [20 270 300 22]; % Y position lowered by 150
        app.ViewFactorObstacleLabel.Text = 'View Factor (A to Obstacle): N/A';
        app.ViewFactorObstacleLabel.Visible = 'off'; % Initially hidden

        % Create EstimatedTimeLabel
        app.EstimatedTimeLabel = uilabel(app.ControlPanel);
        app.EstimatedTimeLabel.Position = [20 240 300 22]; % Y position lowered by 150
        app.EstimatedTimeLabel.Text = 'Estimated Time: ';

        % Create ProgressLabel
        app.ProgressLabel = uilabel(app.ControlPanel);
        app.ProgressLabel.Position = [20 210 300 22]; % Y position lowered by 150
        app.ProgressLabel.Text = 'Progress: 0%';

        % Create SaveDataButton
        app.SaveDataButton = uibutton(app.ControlPanel, 'push');
        app.SaveDataButton.ButtonPushedFcn = createCallbackFcn(app, @SaveDataButtonPushed, true);
        app.SaveDataButton.Position = [20 160 100 30]; % Y position lowered by 150
        app.SaveDataButton.Text = 'Save Data';

        % Create LoadDataButton
        app.LoadDataButton = uibutton(app.ControlPanel, 'push');
        app.LoadDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDataButtonPushed, true);
        app.LoadDataButton.Position = [130 160 100 30]; % Y position lowered by 150
        app.LoadDataButton.Text = 'Load Data';

        % Create UIAxes
        app.UIAxes = uiaxes(app.GridLayout);
        title(app.UIAxes, 'Bodies Visualization')
        xlabel(app.UIAxes, 'X')
        ylabel(app.UIAxes, 'Y')
        zlabel(app.UIAxes, 'Z')
        app.UIAxes.Layout.Row = 1;
        app.UIAxes.Layout.Column = 2;
        app.UIAxes.Box = 'on';
        app.UIAxes.View = [45 30];
        grid(app.UIAxes, 'on');

        % Show the figure after all components are created
        app.UIFigure.Visible = 'on';

        %% Add Watermark Label
        % Create WatermarkLabel
        app.WatermarkLabel = uilabel(app.UIFigure);
        app.WatermarkLabel.Text = 'Created By: Francesco Lena';
        app.WatermarkLabel.FontSize = 10;
        app.WatermarkLabel.FontColor = [0.6 0.6 0.6]; % Light gray color
        app.WatermarkLabel.HorizontalAlignment = 'right';
        app.WatermarkLabel.VerticalAlignment = 'bottom';
        app.WatermarkLabel.Position = [app.UIFigure.Position(3)-200, 5, 190, 20];
        % Callback for figure resizing to adjust watermark
        app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @UIFigureSizeChanged, true);
    end


        % Callback function for UIFigure size changes
        function UIFigureSizeChanged(app, ~)
            % Adjust the position of the watermark label when the figure is resized
            figWidth = app.UIFigure.Position(3);
            figHeight = app.UIFigure.Position(4);
            % Set the position relative to the bottom right corner
            app.WatermarkLabel.Position = [figWidth - 200, 5, 190, 20];
        end
    end

    % App initialization and construction
    methods (Access = public)

        % Construct app
        function app = ViewFactorCalculatorAppV5

            % Create components and register app
            createComponents(app);

            % Initialize properties
            app.isBodyALoaded = false;
            app.isBodyBLoaded = false;
            app.isObstacleLoaded = false;
            app.rayOrigins = [];
            app.rayEndpoints = [];
            app.rayStatuses = [];
            app.viewFactor = NaN;
            app.viewFactorObstacle = NaN;

            % Update UI based on initial calculation mode
            app.updateCalculationModeUI();

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
