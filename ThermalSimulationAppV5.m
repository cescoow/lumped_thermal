% Matlab GUI wrapper for the Lumped Elements thermal simulation
% Created by Francesco Lena. October 2024
% Any issues please contact: francescorossilena@gmail.com

classdef ThermalSimulationAppV5 < matlab.apps.AppBase
    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        ElementsTab                     matlab.ui.container.Tab
        ElementsUITable                 matlab.ui.control.Table
        AddElementButton                matlab.ui.control.Button
        RemoveElementButton             matlab.ui.control.Button
        ConnectionsTab                  matlab.ui.container.Tab
        ConnectionsUITable              matlab.ui.control.Table
        AddConnectionButton             matlab.ui.control.Button
        RemoveConnectionButton          matlab.ui.control.Button
        RadiationConnectionsTab         matlab.ui.container.Tab
        RadiationConnectionsUITable     matlab.ui.control.Table
        AddRadiationConnectionButton    matlab.ui.control.Button
        RemoveRadiationConnectionButton matlab.ui.control.Button
        HeatLoadsTab                    matlab.ui.container.Tab
        HeatLoadsUITable                matlab.ui.control.Table
        AddHeatLoadButton               matlab.ui.control.Button
        RemoveHeatLoadButton            matlab.ui.control.Button
        BoundaryConditionsTab           matlab.ui.container.Tab
        BoundaryConditionsUITable       matlab.ui.control.Table
        AddBoundaryConditionButton      matlab.ui.control.Button
        RemoveBoundaryConditionButton   matlab.ui.control.Button
        SimulationTab                   matlab.ui.container.Tab
        RunSimulationButton             matlab.ui.control.Button
        SimulationTimeEditFieldLabel    matlab.ui.control.Label
        SimulationTimeEditField         matlab.ui.control.NumericEditField
        AmbientTemperatureEditFieldLabel matlab.ui.control.Label
        AmbientTemperatureEditField     matlab.ui.control.NumericEditField
        NumberOfTimeStepsEditFieldLabel matlab.ui.control.Label
        NumberOfTimeStepsEditField      matlab.ui.control.NumericEditField
        SaveDataButton                  matlab.ui.control.Button
        LoadDataButton                  matlab.ui.control.Button
        ResultsAxes                     matlab.ui.control.UIAxes
        SystemGraphAxes                 matlab.ui.control.UIAxes
    end

    properties (Access = private)
        elements                   % Struct array of elements
        connections                % Array of connections
        radiationConnections       % Array of radiation connections
        heatLoads                  % Array of heat loads
        boundaryConditions         % Array of boundary conditions
    end

    methods (Access = private)

        function startupFcn(app)
            % Initialize data with all required fields, including 'initial_temp'
            app.elements = struct('mass', {}, 'cp', {}, 'name', {}, 'g', {}, 'area', {}, 'epsilon', {}, 'initial_temp', {});
            app.connections = [];
            app.radiationConnections = [];
            app.heatLoads = [];
            app.boundaryConditions = [];

            % Initialize tables
            app.ElementsUITable.Data = {};
            app.ConnectionsUITable.Data = {};
            app.RadiationConnectionsUITable.Data = {};
            app.HeatLoadsUITable.Data = {};
            app.BoundaryConditionsUITable.Data = {};
        end

        % --- Element Functions ---

        function AddElementButtonPushed(app, ~)
            % Prompt user for element properties
            prompt = {'Name:', 'Mass:', 'Specific Heat (cp):', 'Conductance (g):', 'Area:', 'Emissivity (epsilon):', 'Initial Temp (K):'};
            dlgtitle = 'Add Element';
            dims = [1 35];
            definput = {'ElementName', '0.1', '400', '0.3', '0.1', '0.05', '293.15'};
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if ~isempty(answer)
                % Create new element
                newElement = struct('mass', str2double(answer{2}), 'cp', str2double(answer{3}), ...
                    'name', answer{1}, 'g', str2double(answer{4}), ...
                    'area', str2double(answer{5}), 'epsilon', str2double(answer{6}), ...
                    'initial_temp', str2double(answer{7}));
                
                % Append the new element ensuring consistent structure
                if isempty(app.elements)
                    app.elements = newElement;
                else
                    app.elements(end+1) = newElement;
                end
                
                % Update table
                app.updateElementsTable();
            end
        end

        function RemoveElementButtonPushed(app, ~)
            % Remove selected element
            idx = app.ElementsUITable.Selection;
            if ~isempty(idx)
                app.elements(idx(1)) = [];
                app.updateElementsTable();
            end
        end

        function updateElementsTable(app)
            % Update the elements table data
            numElements = length(app.elements);
            if numElements == 0
                app.ElementsUITable.Data = {};
                return;
            end
            data = cell(numElements, 7);
            for i = 1:numElements
                data{i, 1} = app.elements(i).name;
                data{i, 2} = app.elements(i).mass;
                data{i, 3} = app.elements(i).cp;
                data{i, 4} = app.elements(i).g;
                data{i, 5} = app.elements(i).area;
                data{i, 6} = app.elements(i).epsilon;
                data{i, 7} = app.elements(i).initial_temp;
            end
            app.ElementsUITable.Data = data;
        end

        function ElementsUITableCellEdit(app, event)
            idx = event.Indices(1);
            col = event.Indices(2);
            newData = event.NewData;

            switch col
                case 1
                    app.elements(idx).name = newData;
                case 2
                    app.elements(idx).mass = newData;
                case 3
                    app.elements(idx).cp = newData;
                case 4
                    app.elements(idx).g = newData;
                case 5
                    app.elements(idx).area = newData;
                case 6
                    app.elements(idx).epsilon = newData;
                case 7
                    app.elements(idx).initial_temp = newData;
            end
        end

        % --- Connection Functions ---

        function AddConnectionButtonPushed(app, ~)
            % Prompt user for connection properties
            if length(app.elements) < 2
                uialert(app.UIFigure, 'At least two elements are required to add a connection.', 'Error');
                return;
            end
            prompt = {'Element 1 Index:', 'Element 2 Index:', 'Interface Conductance:'};
            dlgtitle = 'Add Connection';
            dims = [1 35];
            definput = {'1', '2', '0.5'};
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if ~isempty(answer)
                % Validate indices
                elem1 = str2double(answer{1});
                elem2 = str2double(answer{2});
                if elem1 < 1 || elem1 > length(app.elements) || elem2 < 1 || elem2 > length(app.elements) || elem1 == elem2
                    uialert(app.UIFigure, 'Invalid element indices.', 'Error');
                    return;
                end
                % Add connection
                conn = [elem1, elem2, str2double(answer{3})];
                app.connections(end+1, :) = conn;
                app.updateConnectionsTable();
            end
        end

        function RemoveConnectionButtonPushed(app, ~)
            % Remove selected connection
            idx = app.ConnectionsUITable.Selection;
            if ~isempty(idx)
                app.connections(idx(1), :) = [];
                app.updateConnectionsTable();
            end
        end

        function updateConnectionsTable(app)
            % Update the connections table data
            numConnections = size(app.connections, 1);
            if numConnections == 0
                app.ConnectionsUITable.Data = {};
                return;
            end
            data = cell(numConnections, 3);
            for i = 1:numConnections
                data{i, 1} = app.connections(i, 1);
                data{i, 2} = app.connections(i, 2);
                data{i, 3} = app.connections(i, 3);
            end
            app.ConnectionsUITable.Data = data;
        end

        function ConnectionsUITableCellEdit(app, event)
            idx = event.Indices(1);
            col = event.Indices(2);
            newData = event.NewData;

            app.connections(idx, col) = newData;
        end

        % --- Radiation Connection Functions ---

        function AddRadiationConnectionButtonPushed(app, ~)
            % Prompt user for radiation connection properties
            if length(app.elements) < 2
                uialert(app.UIFigure, 'At least two elements are required to add a radiation connection.', 'Error');
                return;
            end
            prompt = {'Element 1 Index:', 'Element 2 Index:', 'View Factor (F12):'};
            dlgtitle = 'Add Radiation Connection';
            dims = [1 35];
            definput = {'1', '2', '1'};
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if ~isempty(answer)
                % Validate indices
                elem1 = str2double(answer{1});
                elem2 = str2double(answer{2});
                if elem1 < 1 || elem1 > length(app.elements) || elem2 < 1 || elem2 > length(app.elements) || elem1 == elem2
                    uialert(app.UIFigure, 'Invalid element indices.', 'Error');
                    return;
                end
                % Add radiation connection
                radConn = [elem1, elem2, str2double(answer{3})];
                app.radiationConnections(end+1, :) = radConn;
                app.updateRadiationConnectionsTable();
            end
        end

        function RemoveRadiationConnectionButtonPushed(app, ~)
            % Remove selected radiation connection
            idx = app.RadiationConnectionsUITable.Selection;
            if ~isempty(idx)
                app.radiationConnections(idx(1), :) = [];
                app.updateRadiationConnectionsTable();
            end
        end

        function updateRadiationConnectionsTable(app)
            % Update the radiation connections table data
            numRadConns = size(app.radiationConnections, 1);
            if numRadConns == 0
                app.RadiationConnectionsUITable.Data = {};
                return;
            end
            data = cell(numRadConns, 3);
            for i = 1:numRadConns
                data{i, 1} = app.radiationConnections(i, 1);
                data{i, 2} = app.radiationConnections(i, 2);
                data{i, 3} = app.radiationConnections(i, 3);
            end
            app.RadiationConnectionsUITable.Data = data;
        end

        function RadiationConnectionsUITableCellEdit(app, event)
            idx = event.Indices(1);
            col = event.Indices(2);
            newData = event.NewData;

            app.radiationConnections(idx, col) = newData;
        end

        % --- Heat Load Functions ---

        function AddHeatLoadButtonPushed(app, ~)
            % Prompt user for heat load properties
            if isempty(app.elements)
                uialert(app.UIFigure, 'No elements available to apply heat load.', 'Error');
                return;
            end
            prompt = {'Element Index:', 'Heat Load (W):'};
            dlgtitle = 'Add Heat Load';
            dims = [1 35];
            definput = {'1', '100'};
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if ~isempty(answer)
                % Validate element index
                elemIdx = str2double(answer{1});
                if elemIdx < 1 || elemIdx > length(app.elements)
                    uialert(app.UIFigure, 'Invalid element index.', 'Error');
                    return;
                end
                % Add heat load
                heatLoad = [elemIdx, str2double(answer{2})];
                app.heatLoads(end+1, :) = heatLoad;
                app.updateHeatLoadsTable();
            end
        end

        function RemoveHeatLoadButtonPushed(app, ~)
            % Remove selected heat load
            idx = app.HeatLoadsUITable.Selection;
            if ~isempty(idx)
                app.heatLoads(idx(1), :) = [];
                app.updateHeatLoadsTable();
            end
        end

        function updateHeatLoadsTable(app)
            % Update the heat loads table data
            numHeatLoads = size(app.heatLoads, 1);
            if numHeatLoads == 0
                app.HeatLoadsUITable.Data = {};
                return;
            end
            data = cell(numHeatLoads, 2);
            for i = 1:numHeatLoads
                data{i, 1} = app.heatLoads(i, 1);
                data{i, 2} = app.heatLoads(i, 2);
            end
            app.HeatLoadsUITable.Data = data;
        end

        function HeatLoadsUITableCellEdit(app, event)
            idx = event.Indices(1);
            col = event.Indices(2);
            newData = event.NewData;

            app.heatLoads(idx, col) = newData;
        end

        % --- Boundary Condition Functions ---

        function AddBoundaryConditionButtonPushed(app, ~)
            % Prompt user for boundary condition properties
            if isempty(app.elements)
                uialert(app.UIFigure, 'No elements available to apply boundary condition.', 'Error');
                return;
            end
            prompt = {'Element Index:', 'Conductance to Temperature:', 'Temperature (K):'};
            dlgtitle = 'Add Boundary Condition';
            dims = [1 35];
            definput = {'1', '1', '333.15'};
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            if ~isempty(answer)
                % Validate element index
                elemIdx = str2double(answer{1});
                if elemIdx < 1 || elemIdx > length(app.elements)
                    uialert(app.UIFigure, 'Invalid element index.', 'Error');
                    return;
                end
                % Add boundary condition
                bc = [elemIdx, str2double(answer{2}), str2double(answer{3})];
                app.boundaryConditions(end+1, :) = bc;
                app.updateBoundaryConditionsTable();
            end
        end

        function RemoveBoundaryConditionButtonPushed(app, ~)
            % Remove selected boundary condition
            idx = app.BoundaryConditionsUITable.Selection;
            if ~isempty(idx)
                app.boundaryConditions(idx(1), :) = [];
                app.updateBoundaryConditionsTable();
            end
        end

        function updateBoundaryConditionsTable(app)
            % Update the boundary conditions table data
            numBC = size(app.boundaryConditions, 1);
            if numBC == 0
                app.BoundaryConditionsUITable.Data = {};
                return;
            end
            data = cell(numBC, 3);
            for i = 1:numBC
                data{i, 1} = app.boundaryConditions(i, 1);
                data{i, 2} = app.boundaryConditions(i, 2);
                data{i, 3} = app.boundaryConditions(i, 3);
            end
            app.BoundaryConditionsUITable.Data = data;
        end

        function BoundaryConditionsUITableCellEdit(app, event)
            idx = event.Indices(1);
            col = event.Indices(2);
            newData = event.NewData;

            app.boundaryConditions(idx, col) = newData;
        end

        % --- Simulation Functions ---

        function RunSimulationButtonPushed(app, ~)
            if isempty(app.elements)
                uialert(app.UIFigure, 'No elements defined for simulation.', 'Error');
                return;
            end

            % Total simulation time
            t_total = app.SimulationTimeEditField.Value; % Total simulation time in seconds

            % Get the number of time steps from the user input
            num_steps = max(1, round(app.NumberOfTimeStepsEditField.Value));
            dt = t_total / num_steps;

            % Time vector for the full simulation
            t_full_sim = linspace(0, t_total, num_steps + 1);

            % Number of bodies
            num_bodies = length(app.elements);

            % Build the heat capacity matrix
            masses = [app.elements.mass];  % Extract mass from each element
            specific_heats = [app.elements.cp];  % Extract specific heat from each element
            E = app.buildHeatCapacity(masses, specific_heats);

            % Initialize the load matrix L and input vector u
            L_full = zeros(num_bodies, 0);  % Start with an empty L matrix
            u_full = zeros(0, length(t_full_sim));  % Empty u matrix

            % Apply heat loads
            for i = 1:size(app.heatLoads, 1)
                [L_full, u_full] = app.apply_heat(L_full, u_full, app.heatLoads(i, 1), app.heatLoads(i, 2), t_full_sim);
            end

            % Apply boundary conditions
            for i = 1:size(app.boundaryConditions, 1)
                [L_full, u_full] = app.apply_temperature(L_full, u_full, zeros(num_bodies), app.boundaryConditions(i, 1), ...
                    app.boundaryConditions(i, 2), app.boundaryConditions(i, 3), t_full_sim);
            end

            % Initialize state variables with individual initial temperatures
            X0 = [app.elements.initial_temp]; % Initial temperatures
            t_all = [];    % Time vector to store simulation results
            X_all = [];    % Temperature matrix to store simulation results

            % Start the time-stepping loop
            for step = 1:num_steps
                % Current time
                t_current = (step - 1) * dt;

                % Time indices for this time step
                t_indices = step : step + 1;

                % Time vector for this step, shifted to start from zero
                t_step = t_full_sim(t_indices) - t_current;

                % Extract input vector u for this time step
                u_step = u_full(:, t_indices);

                % Rebuild K_rad based on current temperatures
                K_rad = zeros(num_bodies);
                K = app.buildHeatTransfer(num_bodies);  % Start with a fresh conductance matrix

                % Add conduction connections
                for i = 1:size(app.connections, 1)
                    K = app.connect(K, app.elements, app.connections(i, 1), app.connections(i, 2), app.connections(i, 3));
                end

                % Add radiation connections
                for i = 1:size(app.radiationConnections, 1)
                    body_idx = app.radiationConnections(i, 1);
                    shield_idx = app.radiationConnections(i, 2);
                    F12 = app.radiationConnections(i, 3);

                    % Use current temperatures
                    T_body = X0(body_idx);
                    T_shield = X0(shield_idx);

                    % Calculate radiative conductance
                    [K, K_rad] = app.connect_rad_to_body(K, K_rad, app.elements, body_idx, shield_idx, F12, T_body, T_shield);
                end

                % Apply boundary conditions (conductances to temperature)
                for i = 1:size(app.boundaryConditions, 1)
                    body_idx = app.boundaryConditions(i, 1);
                    conductance_to_temperature = app.boundaryConditions(i, 2);
                    K(body_idx, body_idx) = K(body_idx, body_idx) + conductance_to_temperature;
                end

                % Build the state-space system for this time step
                A = -inv(E) * K;
                B = inv(E) * L_full;  % Input matrix B is based on the load matrix L_full
                C = eye(num_bodies);  % Output all temperatures
                D = zeros(num_bodies, size(L_full, 2));  % Zero matrix since there's no direct input-output coupling

                % Create state-space system
                sys = ss(A, B, C, D);

                % Simulate the system over this time step
                [~, ~, X_step] = lsim(sys, u_step', t_step, X0);

                % Update initial condition for the next time step
                X0 = X_step(end, :);

                % Store results
                t_all = [t_all; t_current + t_step'];
                X_all = [X_all; X_step];
            end

            % Get the final temperatures for each body
            final_temperatures = X0;

            % Plot the transient temperature response over time
            cla(app.ResultsAxes);
            plot(app.ResultsAxes, t_all, X_all, 'LineWidth', 2);
            title(app.ResultsAxes, 'Temperature Response Over Time', 'FontSize', 14);
            xlabel(app.ResultsAxes, 'Time [s]', 'FontSize', 12);
            ylabel(app.ResultsAxes, 'Temperature [K]', 'FontSize', 12);
            legend(app.ResultsAxes, {app.elements.name}, 'Interpreter', 'none', 'Location', 'best');
            grid(app.ResultsAxes, 'on');
            app.ResultsAxes.FontSize = 12;

            % Plot the system graph with final temperatures
            app.plot_system_graph(K, K_rad, final_temperatures, app.elements);
        end

        function plot_system_graph(app, K, K_rad, final_temperatures, elements)
            num_bodies = length(final_temperatures);

            % Create a graph object based on the connections in K
            G = graph(K~=0 & ~eye(size(K)), 'upper');

            % Prepare node labels combining element names and temperatures
            node_labels = arrayfun(@(i) sprintf('%s\nT=%.2fK', elements(i).name, final_temperatures(i)), ...
                1:num_bodies, 'UniformOutput', false);

            % Extract the actual edges from the graph object G
            edges = G.Edges.EndNodes;

            % Initialize edge colors (default all black for normal connections)
            edge_colors = repmat({'k'}, size(edges, 1), 1);

            % Update edge colors to red ('r') for radiation connections (those in K_rad)
            for i = 1:size(edges, 1)
                if K_rad(edges(i, 1), edges(i, 2)) ~= 0
                    edge_colors{i} = 'r';  % Set color to red for radiation connections
                end
            end

            % Extract the corresponding conductance values from K for the edges
            edge_conductances = arrayfun(@(r, c) K(r, c), edges(:,1), edges(:,2));
            edge_labels = arrayfun(@(x) sprintf('%.2f', abs(x)), edge_conductances, 'UniformOutput', false);

            % Plot the graph with node labels and edge labels
            cla(app.SystemGraphAxes);
            h = plot(app.SystemGraphAxes, G, 'Layout', 'force', 'NodeLabel', node_labels, 'EdgeLabel', edge_labels, 'LineWidth', 2);

            % Apply the edge colors to the plot
            for i = 1:size(edges, 1)
                highlight(h, edges(i, 1), edges(i, 2), 'EdgeColor', edge_colors{i}, 'LineWidth', 2);
            end

            % Enhance node appearance
            h.MarkerSize = 8;
            h.NodeColor = [0 0.4470 0.7410];
            h.NodeFontSize = 12;
            h.EdgeFontSize = 12;

            title(app.SystemGraphAxes, 'System Connections', 'FontSize', 14);
        end

        % --- Utility Functions ---

        % Function to build the heat capacity matrix
        function E = buildHeatCapacity(~, masses, specific_heats)
            heat_capacities = masses .* specific_heats;
            E = diag(heat_capacities);
        end

        % Function to initialize the heat transfer matrix
        function K = buildHeatTransfer(~, num_bodies)
            K = zeros(num_bodies);
        end

        % Modified connect function to include 1/2 of body conductances
        function K = connect(~, K, elements, body_1, body_2, interface_conductance)
            total_conductance = 0.5 * elements(body_1).g + 0.5 * elements(body_2).g + interface_conductance;
            K(body_1, body_1) = K(body_1, body_1) + total_conductance;
            K(body_2, body_2) = K(body_2, body_2) + total_conductance;
            K(body_1, body_2) = K(body_1, body_2) - total_conductance;
            K(body_2, body_1) = K(body_2, body_1) - total_conductance;
        end

        % Function to apply a heat load to a body
        function [L, u] = apply_heat(~, L, u, body_idx, heat_load, t_sim)
            num_bodies = size(L, 1);
            new_L_column = zeros(num_bodies, 1);
            new_L_column(body_idx) = 1;
            L = [L, new_L_column];
            new_u_row = heat_load * ones(1, length(t_sim));
            u = [u; new_u_row];
        end

        % Function to apply a temperature boundary condition to a body
        function [L, u] = apply_temperature(~, L, u, ~, body_idx, conductance_to_temperature, temperature, t_sim)
            num_bodies = size(L, 1);
            new_L_column = zeros(num_bodies, 1);
            new_L_column(body_idx) = conductance_to_temperature;
            L = [L, new_L_column];
            new_u_row = temperature * ones(1, length(t_sim));
            u = [u; new_u_row];
            % Do not modify K here; handle it in the loop
        end

        % Function to calculate radiative conductance between two bodies
        function [K, K_rad] = connect_rad_to_body(~, K, K_rad, elements, body_idx, shield_idx, F12, T_body, T_shield)
            epsilon_body = elements(body_idx).epsilon;
            epsilon_shield = elements(shield_idx).epsilon;
            area_body = elements(body_idx).area;
            area_shield = elements(shield_idx).area;
            T0 = (T_body + T_shield) / 2; % Reference temperature

            G_rad = ThermalSimulationAppV5.radiativeConductanceTwoBodies(epsilon_body, epsilon_shield, area_body, area_shield, F12, T0);

            % Update K matrix and K_rad for separate plotting
            K(body_idx, body_idx) = K(body_idx, body_idx) + G_rad;
            K(shield_idx, shield_idx) = K(shield_idx, shield_idx) + G_rad;
            K(body_idx, shield_idx) = K(body_idx, shield_idx) - G_rad;
            K(shield_idx, body_idx) = K(shield_idx, body_idx) - G_rad;

            % Store the radiation connection in K_rad for plotting
            K_rad(body_idx, shield_idx) = G_rad;
            K_rad(shield_idx, body_idx) = G_rad;
        end

        % --- Save and Load Functions ---

        function SaveDataButtonPushed(app, ~)
            % Prompt user to select file
            [file, path] = uiputfile('*.mat', 'Save Data As');
            if isequal(file, 0)
                return; % User canceled
            end
            % Save data to file
            elements = app.elements;
            connections = app.connections;
            radiationConnections = app.radiationConnections;
            heatLoads = app.heatLoads;
            boundaryConditions = app.boundaryConditions;
            simulationTime = app.SimulationTimeEditField.Value;
            ambientTemperature = app.AmbientTemperatureEditField.Value;
            numberOfTimeSteps = app.NumberOfTimeStepsEditField.Value;
            save(fullfile(path, file), 'elements', 'connections', 'radiationConnections', ...
                'heatLoads', 'boundaryConditions', 'simulationTime', 'ambientTemperature', 'numberOfTimeSteps');
        end

        function LoadDataButtonPushed(app, ~)
            % Prompt user to select file
            [file, path] = uigetfile('*.mat', 'Load Data');
            if isequal(file, 0)
                return; % User canceled
            end
            % Load data from file
            data = load(fullfile(path, file));
            if isfield(data, 'elements')
                app.elements = data.elements;
                % Ensure all elements have 'initial_temp' field
                for i = 1:length(app.elements)
                    if ~isfield(app.elements(i), 'initial_temp')
                        app.elements(i).initial_temp = 293.15; % Default initial temp
                    end
                end
            else
                app.elements = struct('mass', {}, 'cp', {}, 'name', {}, 'g', {}, 'area', {}, 'epsilon', {}, 'initial_temp', {});
            end
            if isfield(data, 'connections')
                app.connections = data.connections;
            else
                app.connections = [];
            end
            if isfield(data, 'radiationConnections')
                app.radiationConnections = data.radiationConnections;
            else
                app.radiationConnections = [];
            end
            if isfield(data, 'heatLoads')
                app.heatLoads = data.heatLoads;
            else
                app.heatLoads = [];
            end
            if isfield(data, 'boundaryConditions')
                app.boundaryConditions = data.boundaryConditions;
            else
                app.boundaryConditions = [];
            end
            if isfield(data, 'simulationTime')
                app.SimulationTimeEditField.Value = data.simulationTime;
            else
                app.SimulationTimeEditField.Value = 1000; % Default value if not saved
            end
            if isfield(data, 'ambientTemperature')
                app.AmbientTemperatureEditField.Value = data.ambientTemperature;
            else
                app.AmbientTemperatureEditField.Value = 293.15; % Default ambient temp
            end
            if isfield(data, 'numberOfTimeSteps')
                app.NumberOfTimeStepsEditField.Value = data.numberOfTimeSteps;
            else
                app.NumberOfTimeStepsEditField.Value = 100; % Default value if not saved
            end

            % Update tables
            app.updateElementsTable();
            app.updateConnectionsTable();
            app.updateRadiationConnectionsTable();
            app.updateHeatLoadsTable();
            app.updateBoundaryConditionsTable();
        end

    end

    % --- Static Methods ---
    methods (Static)
        % Static function to calculate radiative conductance between two bodies
        function G_rad = radiativeConductanceTwoBodies(epsilon1, epsilon2, A1, A2, F12, T0)
            % Stefan-Boltzmann constant (W/m^2*K^4)
            sigma = 5.670374419e-8;

            % Compute the radiative conductance G_rad between two surfaces
            G_rad = (sigma * A1 * F12 * 4 * T0^3) / ((1/epsilon1) + (A1/(A2 * epsilon2)) - 1);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1200 720];
            app.UIFigure.Name = 'Lumped Thermal Simulation App';

            % Set UIFigure properties for better aesthetics
            app.UIFigure.Color = [0.94, 0.94, 0.94];

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [10 10 350 680];

            % Create ElementsTab
            app.ElementsTab = uitab(app.TabGroup);
            app.ElementsTab.Title = 'Elements';

            % Create ElementsUITable
            app.ElementsUITable = uitable(app.ElementsTab);
            app.ElementsUITable.ColumnName = {'Name'; 'Mass'; 'cp'; 'g'; 'Area'; 'Epsilon'; 'Initial Temp (K)'};
            app.ElementsUITable.ColumnEditable = true(1,7);
            app.ElementsUITable.CellEditCallback = createCallbackFcn(app, @ElementsUITableCellEdit, true);
            app.ElementsUITable.Position = [10 90 330 550];
            app.ElementsUITable.FontSize = 12;

            % Create AddElementButton
            app.AddElementButton = uibutton(app.ElementsTab, 'push');
            app.AddElementButton.ButtonPushedFcn = createCallbackFcn(app, @AddElementButtonPushed, true);
            app.AddElementButton.Position = [10 20 150 50];
            app.AddElementButton.Text = 'Add Element';
            app.AddElementButton.FontSize = 12;

            % Create RemoveElementButton
            app.RemoveElementButton = uibutton(app.ElementsTab, 'push');
            app.RemoveElementButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveElementButtonPushed, true);
            app.RemoveElementButton.Position = [190 20 150 50];
            app.RemoveElementButton.Text = 'Remove Element';
            app.RemoveElementButton.FontSize = 12;

            % Create ConnectionsTab
            app.ConnectionsTab = uitab(app.TabGroup);
            app.ConnectionsTab.Title = 'Connections';

            % Create ConnectionsUITable
            app.ConnectionsUITable = uitable(app.ConnectionsTab);
            app.ConnectionsUITable.ColumnName = {'Element 1'; 'Element 2'; 'Interface Conductance'};
            app.ConnectionsUITable.ColumnEditable = true(1,3);
            app.ConnectionsUITable.CellEditCallback = createCallbackFcn(app, @ConnectionsUITableCellEdit, true);
            app.ConnectionsUITable.Position = [10 90 330 550];
            app.ConnectionsUITable.FontSize = 12;

            % Create AddConnectionButton
            app.AddConnectionButton = uibutton(app.ConnectionsTab, 'push');
            app.AddConnectionButton.ButtonPushedFcn = createCallbackFcn(app, @AddConnectionButtonPushed, true);
            app.AddConnectionButton.Position = [10 20 150 50];
            app.AddConnectionButton.Text = 'Add Connection';
            app.AddConnectionButton.FontSize = 12;

            % Create RemoveConnectionButton
            app.RemoveConnectionButton = uibutton(app.ConnectionsTab, 'push');
            app.RemoveConnectionButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveConnectionButtonPushed, true);
            app.RemoveConnectionButton.Position = [190 20 150 50];
            app.RemoveConnectionButton.Text = 'Remove Connection';
            app.RemoveConnectionButton.FontSize = 12;

            % Create RadiationConnectionsTab
            app.RadiationConnectionsTab = uitab(app.TabGroup);
            app.RadiationConnectionsTab.Title = 'Radiation Connections';

            % Create RadiationConnectionsUITable
            app.RadiationConnectionsUITable = uitable(app.RadiationConnectionsTab);
            app.RadiationConnectionsUITable.ColumnName = {'Element 1'; 'Element 2'; 'View Factor (F12)'};
            app.RadiationConnectionsUITable.ColumnEditable = true(1,3);
            app.RadiationConnectionsUITable.CellEditCallback = createCallbackFcn(app, @RadiationConnectionsUITableCellEdit, true);
            app.RadiationConnectionsUITable.Position = [10 90 330 550];
            app.RadiationConnectionsUITable.FontSize = 12;

            % Create AddRadiationConnectionButton
            app.AddRadiationConnectionButton = uibutton(app.RadiationConnectionsTab, 'push');
            app.AddRadiationConnectionButton.ButtonPushedFcn = createCallbackFcn(app, @AddRadiationConnectionButtonPushed, true);
            app.AddRadiationConnectionButton.Position = [10 20 150 50];
            app.AddRadiationConnectionButton.Text = 'Add Radiation Connection';
            app.AddRadiationConnectionButton.FontSize = 12;

            % Create RemoveRadiationConnectionButton
            app.RemoveRadiationConnectionButton = uibutton(app.RadiationConnectionsTab, 'push');
            app.RemoveRadiationConnectionButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveRadiationConnectionButtonPushed, true);
            app.RemoveRadiationConnectionButton.Position = [190 20 150 50];
            app.RemoveRadiationConnectionButton.Text = 'Remove Radiation Connection';
            app.RemoveRadiationConnectionButton.FontSize = 12;

            % Create HeatLoadsTab
            app.HeatLoadsTab = uitab(app.TabGroup);
            app.HeatLoadsTab.Title = 'Heat Loads';

            % Create HeatLoadsUITable
            app.HeatLoadsUITable = uitable(app.HeatLoadsTab);
            app.HeatLoadsUITable.ColumnName = {'Element Index'; 'Heat Load (W)'};
            app.HeatLoadsUITable.ColumnEditable = true(1,2);
            app.HeatLoadsUITable.CellEditCallback = createCallbackFcn(app, @HeatLoadsUITableCellEdit, true);
            app.HeatLoadsUITable.Position = [10 90 330 550];
            app.HeatLoadsUITable.FontSize = 12;

            % Create AddHeatLoadButton
            app.AddHeatLoadButton = uibutton(app.HeatLoadsTab, 'push');
            app.AddHeatLoadButton.ButtonPushedFcn = createCallbackFcn(app, @AddHeatLoadButtonPushed, true);
            app.AddHeatLoadButton.Position = [10 20 150 50];
            app.AddHeatLoadButton.Text = 'Add Heat Load';
            app.AddHeatLoadButton.FontSize = 12;

            % Create RemoveHeatLoadButton
            app.RemoveHeatLoadButton = uibutton(app.HeatLoadsTab, 'push');
            app.RemoveHeatLoadButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveHeatLoadButtonPushed, true);
            app.RemoveHeatLoadButton.Position = [190 20 150 50];
            app.RemoveHeatLoadButton.Text = 'Remove Heat Load';
            app.RemoveHeatLoadButton.FontSize = 12;

            % Create BoundaryConditionsTab
            app.BoundaryConditionsTab = uitab(app.TabGroup);
            app.BoundaryConditionsTab.Title = 'Boundary Conditions';

            % Create BoundaryConditionsUITable
            app.BoundaryConditionsUITable = uitable(app.BoundaryConditionsTab);
            app.BoundaryConditionsUITable.ColumnName = {'Element Index'; 'Conductance'; 'Temperature (K)'};
            app.BoundaryConditionsUITable.ColumnEditable = true(1,3);
            app.BoundaryConditionsUITable.CellEditCallback = createCallbackFcn(app, @BoundaryConditionsUITableCellEdit, true);
            app.BoundaryConditionsUITable.Position = [10 90 330 550];
            app.BoundaryConditionsUITable.FontSize = 12;

            % Create AddBoundaryConditionButton
            app.AddBoundaryConditionButton = uibutton(app.BoundaryConditionsTab, 'push');
            app.AddBoundaryConditionButton.ButtonPushedFcn = createCallbackFcn(app, @AddBoundaryConditionButtonPushed, true);
            app.AddBoundaryConditionButton.Position = [10 20 150 50];
            app.AddBoundaryConditionButton.Text = 'Add Boundary Condition';
            app.AddBoundaryConditionButton.FontSize = 12;

            % Create RemoveBoundaryConditionButton
            app.RemoveBoundaryConditionButton = uibutton(app.BoundaryConditionsTab, 'push');
            app.RemoveBoundaryConditionButton.ButtonPushedFcn = createCallbackFcn(app, @RemoveBoundaryConditionButtonPushed, true);
            app.RemoveBoundaryConditionButton.Position = [190 20 150 50];
            app.RemoveBoundaryConditionButton.Text = 'Remove Boundary Condition';
            app.RemoveBoundaryConditionButton.FontSize = 12;

            % Create SimulationTab
            app.SimulationTab = uitab(app.TabGroup);
            app.SimulationTab.Title = 'Simulation';

            % Create SimulationTimeEditFieldLabel
            app.SimulationTimeEditFieldLabel = uilabel(app.SimulationTab);
            app.SimulationTimeEditFieldLabel.HorizontalAlignment = 'right';
            app.SimulationTimeEditFieldLabel.Position = [10 260 120 22];
            app.SimulationTimeEditFieldLabel.Text = 'Simulation Time (s)';
            app.SimulationTimeEditFieldLabel.FontSize = 12;

            % Create SimulationTimeEditField
            app.SimulationTimeEditField = uieditfield(app.SimulationTab, 'numeric');
            app.SimulationTimeEditField.Position = [140 260 180 22];
            app.SimulationTimeEditField.Value = 1000;
            app.SimulationTimeEditField.FontSize = 12;

            % Create AmbientTemperatureEditFieldLabel
            app.AmbientTemperatureEditFieldLabel = uilabel(app.SimulationTab);
            app.AmbientTemperatureEditFieldLabel.HorizontalAlignment = 'right';
            app.AmbientTemperatureEditFieldLabel.Position = [10 220 120 22];
            app.AmbientTemperatureEditFieldLabel.Text = 'Ambient Temp (K)';
            app.AmbientTemperatureEditFieldLabel.FontSize = 12;

            % Create AmbientTemperatureEditField
            app.AmbientTemperatureEditField = uieditfield(app.SimulationTab, 'numeric');
            app.AmbientTemperatureEditField.Position = [140 220 180 22];
            app.AmbientTemperatureEditField.Value = 293.15;
            app.AmbientTemperatureEditField.FontSize = 12;

            % Create NumberOfTimeStepsEditFieldLabel
            app.NumberOfTimeStepsEditFieldLabel = uilabel(app.SimulationTab);
            app.NumberOfTimeStepsEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberOfTimeStepsEditFieldLabel.Position = [10 180 120 22];
            app.NumberOfTimeStepsEditFieldLabel.Text = 'Time Steps';
            app.NumberOfTimeStepsEditFieldLabel.FontSize = 12;

            % Create NumberOfTimeStepsEditField
            app.NumberOfTimeStepsEditField = uieditfield(app.SimulationTab, 'numeric');
            app.NumberOfTimeStepsEditField.Position = [140 180 180 22];
            app.NumberOfTimeStepsEditField.Value = 100;  % Default number of time steps
            app.NumberOfTimeStepsEditField.FontSize = 12;
            app.NumberOfTimeStepsEditField.Limits = [1 Inf];
            app.NumberOfTimeStepsEditField.RoundFractionalValues = true;  % Ensure integer values

            % Create RunSimulationButton
            app.RunSimulationButton = uibutton(app.SimulationTab, 'push');
            app.RunSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @RunSimulationButtonPushed, true);
            app.RunSimulationButton.Position = [10 120 310 50];
            app.RunSimulationButton.Text = 'Run Simulation';
            app.RunSimulationButton.FontSize = 14;
            app.RunSimulationButton.BackgroundColor = [0.47 0.67 0.19];
            app.RunSimulationButton.FontColor = [1 1 1];

            % Create SaveDataButton
            app.SaveDataButton = uibutton(app.SimulationTab, 'push');
            app.SaveDataButton.ButtonPushedFcn = createCallbackFcn(app, @SaveDataButtonPushed, true);
            app.SaveDataButton.Position = [10 60 150 30];
            app.SaveDataButton.Text = 'Save Data';
            app.SaveDataButton.FontSize = 12;

            % Create LoadDataButton
            app.LoadDataButton = uibutton(app.SimulationTab, 'push');
            app.LoadDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDataButtonPushed, true);
            app.LoadDataButton.Position = [170 60 150 30];
            app.LoadDataButton.Text = 'Load Data';
            app.LoadDataButton.FontSize = 12;

            % Create ResultsAxes
            app.ResultsAxes = uiaxes(app.UIFigure);
            title(app.ResultsAxes, 'Temperature Response Over Time', 'FontSize', 14)
            xlabel(app.ResultsAxes, 'Time [s]', 'FontSize', 12)
            ylabel(app.ResultsAxes, 'Temperature [K]', 'FontSize', 12)
            app.ResultsAxes.Position = [380 360 800 350];
            app.ResultsAxes.FontSize = 12;
            app.ResultsAxes.GridColor = [0.15 0.15 0.15];
            app.ResultsAxes.GridAlpha = 0.3;
            grid(app.ResultsAxes, 'on');

            % Create SystemGraphAxes
            app.SystemGraphAxes = uiaxes(app.UIFigure);
            title(app.SystemGraphAxes, 'System Connections', 'FontSize', 14)
            app.SystemGraphAxes.Position = [380 10 800 350];
            app.SystemGraphAxes.FontSize = 12;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    methods (Access = public)

        % Construct app
        function app = ThermalSimulationAppV5

            % Create UIFigure and components
            createComponents(app)

            % Execute the startup function
            startupFcn(app)
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
