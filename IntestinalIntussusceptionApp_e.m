%% IntestinalIntussusceptionApp_e.m
% MATLAB App for Intussusception Simulation with Enhanced Elasticity Control
% This script creates an App for simulating intestinal motility and intussusception formation
% Modified version with spectrum distribution of peristaltic speeds and elasticity control
% Creator: Ming-Yu Hsieh
% Chung Shan Medical University Hospital

classdef IntestinalIntussusceptionApp_e < matlab.apps.AppBase

    % Properties for UI components
    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        SimulationPanel            matlab.ui.container.Panel
        RunButton                  matlab.ui.control.Button
        ResetButton                matlab.ui.control.Button
        StatusLabel                matlab.ui.control.Label
        TableResultsButton         matlab.ui.control.Button
        TabGroup                   matlab.ui.container.TabGroup
        BasicParamsTab             matlab.ui.container.Tab
        GastroParamsTab            matlab.ui.container.Tab
        MotilityParamsTab          matlab.ui.container.Tab  % Motility spectrum tab
        ElasticityParamsTab        matlab.ui.container.Tab  % New elasticity tab
        AnatomyParamsTab           matlab.ui.container.Tab
        ResultsTab                 matlab.ui.container.Tab
        
        % Basic parameters
        SimTimeEditField           matlab.ui.control.NumericEditField
        SimTimeEditFieldLabel      matlab.ui.control.Label
        TimeStepEditField          matlab.ui.control.NumericEditField
        TimeStepEditFieldLabel     matlab.ui.control.Label
        GastroOnsetEditField       matlab.ui.control.NumericEditField
        GastroOnsetEditFieldLabel  matlab.ui.control.Label
        
        % Gastroenteritis parameters
        ProximalSpeedEditField           matlab.ui.control.NumericEditField
        ProximalSpeedEditFieldLabel      matlab.ui.control.Label
        TerminalSpeedEditField           matlab.ui.control.NumericEditField
        TerminalSpeedEditFieldLabel      matlab.ui.control.Label
        BacterialFactorEditField         matlab.ui.control.NumericEditField
        BacterialFactorEditFieldLabel    matlab.ui.control.Label
        ValveResistanceEditField         matlab.ui.control.NumericEditField
        ValveResistanceEditFieldLabel    matlab.ui.control.Label
        WaveAmplitudeEditField           matlab.ui.control.NumericEditField
        WaveAmplitudeEditFieldLabel      matlab.ui.control.Label
        WaveFrequencyEditField           matlab.ui.control.NumericEditField
        WaveFrequencyEditFieldLabel      matlab.ui.control.Label
        
        % Motility spectrum parameters
        SpeedTransitionTypeDropDown      matlab.ui.control.DropDown
        SpeedTransitionTypeLabel         matlab.ui.control.Label
        SpeedGradientEditField           matlab.ui.control.NumericEditField
        SpeedGradientEditFieldLabel      matlab.ui.control.Label
        SpeedVariabilityEditField        matlab.ui.control.NumericEditField
        SpeedVariabilityEditFieldLabel   matlab.ui.control.Label
        NoiseLevelEditField              matlab.ui.control.NumericEditField
        NoiseLevelEditFieldLabel         matlab.ui.control.Label
        SpeedPatternDropDown             matlab.ui.control.DropDown
        SpeedPatternLabel                matlab.ui.control.Label
        
         % Normal state speed control element
        NormalProximalSpeedEditField        matlab.ui.control.NumericEditField
        NormalProximalSpeedEditFieldLabel   matlab.ui.control.Label
        NormalTerminalSpeedEditField        matlab.ui.control.NumericEditField
        NormalTerminalSpeedEditFieldLabel   matlab.ui.control.Label

        % New elasticity parameters
        NormalProximalElasticityEditField      matlab.ui.control.NumericEditField
        NormalProximalElasticityEditFieldLabel matlab.ui.control.Label
        NormalDistalElasticityEditField        matlab.ui.control.NumericEditField
        NormalDistalElasticityEditFieldLabel   matlab.ui.control.Label
        GastroProximalElasticityEditField      matlab.ui.control.NumericEditField
        GastroProximalElasticityEditFieldLabel matlab.ui.control.Label
        GastroDistalElasticityEditField        matlab.ui.control.NumericEditField
        GastroDistalElasticityEditFieldLabel   matlab.ui.control.Label
        ElasticityPatternDropDown              matlab.ui.control.DropDown
        ElasticityPatternLabel                 matlab.ui.control.Label
        ElasticityGradientEditField            matlab.ui.control.NumericEditField
        ElasticityGradientEditFieldLabel       matlab.ui.control.Label
        ElasticityVariabilityEditField         matlab.ui.control.NumericEditField
        ElasticityVariabilityEditFieldLabel    matlab.ui.control.Label
        
        % Anatomy parameters
        IntestineLengthEditField         matlab.ui.control.NumericEditField
        IntestineLengthEditFieldLabel    matlab.ui.control.Label
        TerminalIleumStartEditField      matlab.ui.control.NumericEditField
        TerminalIleumStartEditFieldLabel matlab.ui.control.Label
        IntestineRadiusEditField         matlab.ui.control.NumericEditField
        IntestineRadiusEditFieldLabel    matlab.ui.control.Label
        CriticalRatioEditField           matlab.ui.control.NumericEditField
        CriticalRatioEditFieldLabel      matlab.ui.control.Label
        
        % Results display
        ResultsTextArea                  matlab.ui.control.TextArea
        ResultPlotArea                   matlab.ui.control.UIAxes
        SaveResultsButton                matlab.ui.control.Button
        SpeedProfilePlotArea             matlab.ui.control.UIAxes  % Speed profile plot
        ElasticityProfilePlotArea        matlab.ui.control.UIAxes  % New elasticity profile plot
    end
    
    % Properties for simulation results
    properties (Access = private)
        results               % Simulation results
        params                % Simulation parameters
    end
    
    methods (Access = public) %Change from Private to Public so that this program can be connected with Sweep.
        
        function out = run_simulation_withParams(app, proxSpeed, distSpeed, proxElastic, distElastic)
                app.ProxSpeedField.Value   = proxSpeed;
                app.DistSpeedField.Value   = distSpeed;
                app.ProxElasticField.Value = proxElastic;
                app.DistElasticField.Value = distElastic;
                out = app.run_simulation();
            end

        function params = getParametersFromUI(app)
            % Get all parameters from UI fields
            
            % Basic parameters
            params.simulation_time = app.SimTimeEditField.Value;
            params.dt = app.TimeStepEditField.Value;
            params.num_steps = round(params.simulation_time / params.dt);
            params.gastroenteritis_onset = app.GastroOnsetEditField.Value;
            
            % Intestinal geometry parameters
            params.intestine_length = app.IntestineLengthEditField.Value;
            params.terminal_ileum_start = app.TerminalIleumStartEditField.Value;
            params.intestine_radius = app.IntestineRadiusEditField.Value;
            params.num_segments = 200;  % Fixed value
            params.segment_length = params.intestine_length / params.num_segments;
            
            % Normal state parameters (fixed values)
            params.normal.proximal_speed = app.NormalProximalSpeedEditField.Value;
            params.normal.terminal_speed = app.NormalTerminalSpeedEditField.Value;
            params.normal.wave_amplitude = 0.4;
            params.normal.wave_frequency = 0.1;
            
            % Gastroenteritis state parameters
            params.gastroenteritis.proximal_speed = app.ProximalSpeedEditField.Value;
            params.gastroenteritis.terminal_speed = app.TerminalSpeedEditField.Value;
            params.gastroenteritis.wave_amplitude = app.WaveAmplitudeEditField.Value;
            params.gastroenteritis.wave_frequency = app.WaveFrequencyEditField.Value;
            
            % Bacterial and ileocecal valve factors
            params.bacterial_factor = app.BacterialFactorEditField.Value;
            params.valve_resistance = app.ValveResistanceEditField.Value;
            
            % New elasticity parameters
            params.normal.proximal_elasticity = app.NormalProximalElasticityEditField.Value;
            params.normal.distal_elasticity = app.NormalDistalElasticityEditField.Value;
            params.gastroenteritis.proximal_elasticity = app.GastroProximalElasticityEditField.Value;
            params.gastroenteritis.distal_elasticity = app.GastroDistalElasticityEditField.Value;
            params.elasticity_pattern = app.ElasticityPatternDropDown.Value;
            params.elasticity_gradient = app.ElasticityGradientEditField.Value;
            params.elasticity_variability = app.ElasticityVariabilityEditField.Value;
            
            % Intestinal wall parameters
            params.wall_elasticity = 0.2;  % Base elasticity (will be modified by profile)
            params.wall_tension = 0.5;
            params.critical_diameter_ratio = app.CriticalRatioEditField.Value;
            
            % Intestinal content parameters (fixed values)
            params.content_density = 1.0;
            params.secretion_rate = 0.01;
            params.absorption_rate = 0.008;
            params.inflamed_absorption_rate = 0.003;
            
            % Motility spectrum parameters
            params.speed_transition_type = app.SpeedTransitionTypeDropDown.Value;
            params.speed_gradient = app.SpeedGradientEditField.Value;
            params.speed_variability = app.SpeedVariabilityEditField.Value;
            params.noise_level = app.NoiseLevelEditField.Value;
            params.speed_pattern = app.SpeedPatternDropDown.Value;
            
            % Visualization parameters
            params.visualization_interval = 20;
        end
      
        function intestine = initialize_model(app, params)
            try
                % Initialize the small intestine model
                intestine = struct();
                
                % Configure the small intestine geometry
                intestine.x = linspace(0, params.intestine_length, params.num_segments);
                intestine.radius = ones(1, params.num_segments) * params.intestine_radius;
                intestine.initial_radius = intestine.radius;
                
                % Initialize the content distribution
                intestine.content = zeros(1, params.num_segments);
                
                % Add initial contents
                intestine.content(1:10) = 0.8;  
                
                % Initialize the wall position
                intestine.wall_position = zeros(2, params.num_segments);  % [inner, outter]
                for i = 1:params.num_segments
                    intestine.wall_position(1, i) = intestine.x(i);  % Initialize inner wall positiion
                    intestine.wall_position(2, i) = intestine.x(i);  % Initialize outter wall positiion
                end
                
                % Set up small intestine regional marker
                intestine.is_terminal_ileum = intestine.x >= params.terminal_ileum_start;
                
                % Initialize the waveform count
                intestine.wave_count = 0;
                
                % Initialize the peristaltic spectrum
                params_copy = params;
                                
                % Generate the speed and elasticity distributions
                try
                    intestine.local_speed = app.generate_speed_profile(params_copy);
                catch e
                    fprintf('Error generating the speed file: %s\n', e.message);
                    % Use a simple linear distribution as an alternative
                    intestine.local_speed = linspace(params.normal.proximal_speed, params.normal.terminal_speed, params.num_segments);
                end
                
                intestine.initial_speed_profile = intestine.local_speed; 
                
                try
                    intestine.local_elasticity = app.generate_elasticity_profile(params, false);
                catch e
                    fprintf('Error generating the elasticity file: %s\n', e.message);
                    % Use a simple linear distribution as an alternative
                    intestine.local_elasticity = linspace(params.normal.proximal_elasticity, params.normal.distal_elasticity, params.num_segments);
                end
                
                intestine.initial_elasticity_profile = intestine.local_elasticity;  
                
                % Set the wave amplitude
                intestine.local_amplitude = ones(1, params.num_segments) * params.normal.wave_amplitude;
                
                % Initial state: Non-gastroenteritis
                intestine.is_gastroenteritis = false;
            catch e
                fprintf('Error initializing the model: %s\n', e.message);
                fprintf('Error details: %s\n', getReport(e));
                intestine = [];  % Returning null indicates initialization failure
            end
        end

        function speed_profile = generate_speed_profile(app, params)
            % Normalized position along intestine (0 to 1)
            norm_position = linspace(0, 1, params.num_segments);
            
            % Start with uniform speed
            speed_profile = ones(1, params.num_segments);
            
            % Generate pattern based on selected profile type
            switch params.speed_pattern
                case 'Gradual Decrease'
                    % Gradual decrease from proximal to distal
                    speed_profile = params.normal.proximal_speed - ...
                        (params.normal.proximal_speed - params.normal.terminal_speed) * norm_position;
                    
                case 'Segmented Regions'
                    % Multiple segments with different speeds
                    num_regions = 5;
                    region_length = 1/num_regions;
                    for i = 1:num_regions
                        segment_idx = norm_position >= (i-1)*region_length & norm_position < i*region_length;
                        segment_speed = params.normal.proximal_speed - ...
                            (i-1)*(params.normal.proximal_speed - params.normal.terminal_speed)/(num_regions-1);
                        speed_profile(segment_idx) = segment_speed;
                    end
                    
                case 'Physiological'
                    % More realistic physiological pattern with gradual transition
                    % Faster in duodenum, slowing in jejunum, even slower in ileum
                    duodenum_idx = norm_position <= 0.1;
                    jejunum_idx = norm_position > 0.1 & norm_position <= 0.6;
                    ileum_idx = norm_position > 0.6;
                    
                    speed_profile(duodenum_idx) = params.normal.proximal_speed * 1.1;  % Fastest in duodenum
                    
                    % Gradual transition in jejunum
                    jejunum_pos = (norm_position(jejunum_idx) - 0.1) / 0.5;  % Normalized within jejunum
                    speed_profile(jejunum_idx) = params.normal.proximal_speed - ...
                        jejunum_pos * (params.normal.proximal_speed - params.normal.terminal_speed) * 0.5;
                    
                    % Gradual transition in ileum
                    ileum_pos = (norm_position(ileum_idx) - 0.6) / 0.4;  % Normalized within ileum
                    start_speed = params.normal.proximal_speed - ...
                        0.5 * (params.normal.proximal_speed - params.normal.terminal_speed);
                    speed_profile(ileum_idx) = start_speed - ...
                        ileum_pos * (start_speed - params.normal.terminal_speed);
                    
                case 'Custom'
                    % Apply speed transition type according to selection
                    switch params.speed_transition_type
                        case 'Linear'
                            speed_profile = params.normal.proximal_speed - ...
                                (params.normal.proximal_speed - params.normal.terminal_speed) * norm_position;
                            
                        case 'Sigmoid'
                            % Sigmoid transition (smooth S-curve)
                            sigmoid = 1./(1 + exp(-params.speed_gradient * (norm_position*2 - 1) * 10));
                            speed_profile = params.normal.proximal_speed - ...
                                (params.normal.proximal_speed - params.normal.terminal_speed) * sigmoid;
                            
                        case 'Exponential'
                            % Exponential decay
                            decay_factor = 5 * params.speed_gradient;  % Adjusts decay rate
                            speed_profile = params.normal.terminal_speed + ...
                                (params.normal.proximal_speed - params.normal.terminal_speed) * exp(-decay_factor * norm_position);
        
                        case 'Step with Transition'
                            % For normal state, use linear gradient (no sharp drop)
                            speed_profile = params.normal.proximal_speed - ...
                                (params.normal.proximal_speed - params.normal.terminal_speed) * norm_position;
                    end
            end
            
            % Add random variability (spatial noise)
            if params.noise_level > 0
                noise = randn(1, params.num_segments) * params.noise_level;
                speed_profile = speed_profile .* (1 + noise);
                
                % Ensure speeds remain within reasonable bounds
                speed_profile = max(speed_profile, 0.5);  % Minimum speed
                speed_profile = min(speed_profile, params.normal.proximal_speed * 1.5);  % Maximum speed
            end
            
            % Add local variability (spatially correlated)
            if params.speed_variability > 0
                % Create smooth variations using filtered noise
                window_size = 15;  % Size of spatial correlation
                raw_noise = randn(1, params.num_segments + window_size*2) * params.speed_variability;
                smooth_noise = filter(ones(1, window_size) / window_size, 1, raw_noise);
                smooth_noise = smooth_noise(window_size+1:end-window_size);
                
                % Apply the smooth variations to the speed profile
                speed_profile = speed_profile .* (1 + smooth_noise);
                
                % Ensure speeds remain within reasonable bounds
                speed_profile = max(speed_profile, 0.5);  % Minimum speed
                speed_profile = min(speed_profile, params.normal.proximal_speed * 1.5);  % Maximum speed
            end
        end


        % NEW: Generate elasticity profile function
        function elasticity_profile = generate_elasticity_profile(app, params, is_gastroenteritis)
            % Generate spectrum-based elasticity profile along the intestine
            
            % Normalized position along intestine (0 to 1)
            norm_position = linspace(0, 1, params.num_segments);
            
            % Start with uniform elasticity
            elasticity_profile = ones(1, params.num_segments);
            
            % Set base elasticity values based on normal or gastroenteritis state
            if is_gastroenteritis
                proximal_elasticity = params.gastroenteritis.proximal_elasticity;
                distal_elasticity = params.gastroenteritis.distal_elasticity;
            else
                proximal_elasticity = params.normal.proximal_elasticity;
                distal_elasticity = params.normal.distal_elasticity;
            end
            
            % Generate pattern based on selected profile type
            switch params.elasticity_pattern
                case 'Uniform'
                    % Use proximal elasticity for all segments (baseline)
                    elasticity_profile = elasticity_profile * proximal_elasticity;
                    
                case 'Linear Gradient'
                    % Linear transition from proximal to distal
                    elasticity_profile = proximal_elasticity - ...
                        (proximal_elasticity - distal_elasticity) * norm_position;
                    
                case 'Sigmoid Transition'
                    % Sigmoid transition for more natural change
                    midpoint = params.terminal_ileum_start / params.intestine_length;
                    sigmoid = 1./(1 + exp(-params.elasticity_gradient * (norm_position - midpoint) * 10));
                    elasticity_profile = proximal_elasticity - ...
                        (proximal_elasticity - distal_elasticity) * sigmoid;
                    
                case 'Anatomical'
                    % Different elasticity for different anatomical regions
                    duodenum_idx = norm_position <= 0.1;
                    jejunum_idx = norm_position > 0.1 & norm_position <= 0.6;
                    ileum_idx = norm_position > 0.6 & norm_position < 0.9;
                    terminal_ileum_idx = norm_position >= 0.9;
                    
                    % Assign elasticity values based on region
                    elasticity_profile(duodenum_idx) = proximal_elasticity * 0.9;  % Slightly less elastic in duodenum
                    elasticity_profile(jejunum_idx) = proximal_elasticity;
                    elasticity_profile(ileum_idx) = proximal_elasticity * 0.7 + distal_elasticity * 0.3;  % Mix
                    elasticity_profile(terminal_ileum_idx) = distal_elasticity;
                    
                case 'Terminal Ileum Focus'
                    % Enhanced focus on terminal ileum with distinct properties
                    elasticity_profile = proximal_elasticity * ones(1, params.num_segments);
                    
                    % Calculate distance to terminal ileum start
                    terminal_ileum_norm = params.terminal_ileum_start / params.intestine_length;
                    transition_width = 0.1;  % Width of transition zone
                    
                    % Create smooth transition to terminal ileum
                    for i = 1:params.num_segments
                        dist_to_terminal = abs(norm_position(i) - terminal_ileum_norm);
                        if norm_position(i) >= terminal_ileum_norm
                            % We're in terminal ileum
                            elasticity_profile(i) = distal_elasticity;
                        elseif dist_to_terminal < transition_width
                            % We're in the transition zone
                            transition_factor = (transition_width - dist_to_terminal) / transition_width;
                            elasticity_profile(i) = proximal_elasticity * (1-transition_factor) + ...
                                                   distal_elasticity * transition_factor;
                        end
                    end
            end
            
            % Add spatial variability if requested
            if params.elasticity_variability > 0
                % Create smooth variations using filtered noise
                window_size = 10;  % Size of spatial correlation
                raw_noise = randn(1, params.num_segments + window_size*2) * params.elasticity_variability;
                smooth_noise = filter(ones(1, window_size) / window_size, 1, raw_noise);
                smooth_noise = smooth_noise(window_size+1:end-window_size);
                
                % Apply the smooth variations to the elasticity profile
                elasticity_profile = elasticity_profile .* (1 + smooth_noise);
                
                % Ensure elasticity remains within reasonable bounds
                min_elasticity = 0.05;  % Minimum elasticity
                max_elasticity = max(proximal_elasticity, distal_elasticity) * 1.5;  % Maximum elasticity
                elasticity_profile = max(elasticity_profile, min_elasticity);
                elasticity_profile = min(elasticity_profile, max_elasticity);
            end
            
            % Add extra changes at specific regions if in gastroenteritis
            if is_gastroenteritis
                % Identify regions of interest based on anatomy
                terminal_ileum_norm = params.terminal_ileum_start / params.intestine_length;
                
                % Add inflammation effect in terminal ileum (potentially reduces elasticity further)
                terminal_ileum_idx = norm_position >= terminal_ileum_norm;
                elasticity_profile(terminal_ileum_idx) = elasticity_profile(terminal_ileum_idx) * 0.8;  % 20% reduction due to inflammation
                
                % Create localized inflammation patches (random locations)
                num_patches = 3;
                for i = 1:num_patches
                    % Random patch center
                    patch_center = rand() * 0.8 + 0.1;  % Between 0.1 and 0.9
                    patch_width = 0.05 + rand() * 0.05;  % Between 0.05 and 0.1
                    
                    % Apply inflammation effect
                    for j = 1:params.num_segments
                        dist_to_center = abs(norm_position(j) - patch_center);
                        if dist_to_center < patch_width
                            effect = (patch_width - dist_to_center) / patch_width;  % 0 to 1 based on distance
                            elasticity_profile(j) = elasticity_profile(j) * (1 - 0.3 * effect);  % Up to 30% reduction
                        end
                    end
                end
            end
        end
        
        function intestine = validateIntestineStructure(app, intestine, params)
              % Check and ensure the intestine struct contains all required fields
              required_fields = {'x', 'radius', 'initial_radius', 'content', 'wall_position', ...
                                'wave_count', 'is_terminal_ileum', 'local_speed', ...
                                'initial_speed_profile', 'local_elasticity', ...
                                'initial_elasticity_profile', 'local_amplitude', 'is_gastroenteritis'};
            
            for i = 1:length(required_fields)
                field = required_fields{i};
                if ~isfield(intestine, field)
                    app.updateStatus(['Missing field in intestine structure: ' field '. Initializing...']);
                    
                    % Initialize the missing fields appropriately
                    switch field
                        case 'radius'
                            intestine.radius = ones(1, params.num_segments) * params.intestine_radius;
                        case 'initial_radius'
                            intestine.initial_radius = intestine.radius;
                        case 'content'
                            intestine.content = zeros(1, params.num_segments);
                            % Add initial content to the proximal segment
                            intestine.content(1:50) = 0.8;
                            % Introduce some content into the middle segment
                            intestine.content(51:65) = 0.5;
                        case 'wall_position'
                            intestine.wall_position = zeros(2, params.num_segments); % [inner wall, outer wall]
                            for j = 1:params.num_segments
                                intestine.wall_position(1, j) = intestine.x(j); % Inner wall position
                                intestine.wall_position(2, j) = intestine.x(j); % Outter wall position
                            end
                        case 'wave_count'
                            intestine.wave_count = 0;
                        case 'is_terminal_ileum'
                            intestine.is_terminal_ileum = intestine.x >= params.terminal_ileum_start;
                        case 'local_speed'
                            intestine.local_speed = app.generate_speed_profile(params);
                        case 'initial_speed_profile'
                            intestine.initial_speed_profile = intestine.local_speed;
                        case 'local_elasticity'
                            intestine.local_elasticity = app.generate_elasticity_profile(params, false);
                        case 'initial_elasticity_profile'
                            intestine.initial_elasticity_profile = intestine.local_elasticity;
                        case 'local_amplitude'
                            intestine.local_amplitude = ones(1, params.num_segments) * params.normal.wave_amplitude;
                        case 'is_gastroenteritis'
                            intestine.is_gastroenteritis = false;
                        case 'x'
                            intestine.x = linspace(0, params.intestine_length, params.num_segments);
                    end
                end
            end
        end

        function results = run_simulation(app, intestine, params)
            % Add progress reporting
            fprintf('Start simulation, total steps: %d\n', params.num_steps);
            
            % Preallocate results storage
            results.time = linspace(0, params.simulation_time, params.num_steps);
            results.radius = zeros(params.num_steps, params.num_segments);
            results.content = zeros(params.num_steps, params.num_segments);
            results.wall_position = zeros(params.num_steps, 2, params.num_segments);
            results.speed_profile = zeros(params.num_steps, params.num_segments);  % Store speed profiles over time
            results.elasticity_profile = zeros(params.num_steps, params.num_segments);  % NEW: Store elasticity profiles
            results.intussusception_detected = false;
            results.intussusception_time = NaN;
            results.intussusception_location = NaN;
            results.intussusception_details = [];  % Store detailed info about intussusception events
            
            % Store initial state
            results.radius(1, :) = intestine.radius;
            results.content(1, :) = intestine.content;
            results.wall_position(1, :, :) = intestine.wall_position;
            results.speed_profile(1, :) = intestine.local_speed;
            results.elasticity_profile(1, :) = intestine.local_elasticity;  % Store initial elasticity profile
            
            % Time iteration
            for step = 2:params.num_steps
                % Current time
                t = (step-1) * params.dt;
                
                % Check if gastroenteritis state has begun
                if t >= params.gastroenteritis_onset && ~intestine.is_gastroenteritis
                    intestine.is_gastroenteritis = true;
                    
                    % Update speed profile for gastroenteritis state
                    intestine.local_speed = app.update_gastroenteritis_speed_profile(intestine, params);
                    
                    % Update elasticity profile for gastroenteritis state - NEW
                    intestine.local_elasticity = app.generate_elasticity_profile(params, true);
                    
                    % Update wave amplitude (can follow a spectrum as well)
                    intestine.local_amplitude = app.update_wave_amplitude_profile(intestine, params);
                    
                    app.updateStatus(sprintf('Time %.1f seconds: Gastroenteritis state activated', t));
                end
                
                % Update intestinal state with elasticity consideration
                intestine = app.update_peristalsis(intestine, t, params);
                
                % Update content movement with elasticity effects
                intestine = app.update_contents(intestine, params);
                
                % Check for intussusception with elasticity factors
                [intestine, intussusception_detected, intussusception_loc, intussusception_details] = ...
                    app.check_intussusception(intestine, t, params);
                
                if intussusception_detected && isnan(results.intussusception_time)
                    results.intussusception_detected = true;
                    results.intussusception_time = t;
                    results.intussusception_location = intussusception_loc;
                    results.intussusception_details = intussusception_details;
                    app.updateStatus(sprintf('Time %.1f seconds: Intussusception detected at position %.1f cm!', t, intussusception_loc));
                end
                
                % Store current state
                results.radius(step, :) = intestine.radius;
                results.content(step, :) = intestine.content;
                results.wall_position(step, :, :) = intestine.wall_position;
                results.speed_profile(step, :) = intestine.local_speed;
                results.elasticity_profile(step, :) = intestine.local_elasticity;  % Store current elasticity profile
                
                % Show progress
                if mod(step, floor(params.num_steps/10)) == 0
                    app.updateStatus(sprintf('Simulation progress: %.1f%%', step/params.num_steps*100));
                    drawnow;
                end
            end
            fprintf('Simulation completed。intussusception Test: %d\n', results.intussusception_detected);
        end
        
        % This is the corrected version with proper syntax for the update_gastroenteritis_speed_profile function

        function updated_speed = update_gastroenteritis_speed_profile(app, intestine, params)
    % Update the speed profile during gastroenteritis state based on selected pattern
    
    % Start with the initial speed profile
    updated_speed = intestine.initial_speed_profile;
    
    % Calculate normalized position
    normalize_pos = intestine.x / params.intestine_length;
    
    % Get gastroenteritis state parameters
    proximal_speed = params.gastroenteritis.proximal_speed;
    terminal_speed = params.gastroenteritis.terminal_speed;
    
    % Apply the selected speed pattern
    switch params.speed_pattern
        case 'Gradual Decrease'
            % Simple linear decrease from proximal to distal
            updated_speed = proximal_speed - (proximal_speed - terminal_speed) * normalize_pos;
            
        case 'Segmented Regions'
            % Multiple segments with different speeds
            num_regions = 5;
            region_length = 1/num_regions;
            for i = 1:length(normalize_pos)
                pos = normalize_pos(i);
                region_idx = floor(pos / region_length) + 1;
                if region_idx > num_regions
                    region_idx = num_regions;
                end
                segment_speed = proximal_speed - (region_idx-1)*(proximal_speed - terminal_speed)/(num_regions-1);
                updated_speed(i) = segment_speed;
            end
            
        case 'Physiological'
            % More realistic physiological pattern
            for i = 1:length(normalize_pos)
                pos = normalize_pos(i);
                if pos <= 0.1
                    % Duodenum - highest speed
                    updated_speed(i) = proximal_speed * 1.1;
                elseif pos <= 0.6
                    % Jejunum - gradual decrease
                    jejunum_pos = (pos - 0.1) / 0.5;
                    updated_speed(i) = proximal_speed - jejunum_pos * (proximal_speed - terminal_speed) * 0.5;
                else
                    % Ileum - further decrease
                    ileum_pos = (pos - 0.6) / 0.4;
                    start_speed = proximal_speed - 0.5 * (proximal_speed - terminal_speed);
                    updated_speed(i) = start_speed - ileum_pos * (start_speed - terminal_speed);
                end
            end
            
        case 'Custom'
            % Apply speed transition type according to selection
            switch params.speed_transition_type
                case 'Linear'
                    updated_speed = proximal_speed - (proximal_speed - terminal_speed) * normalize_pos;
                    
                case 'Sigmoid'
                    % Sigmoid transition (smooth S-curve)
                    for i = 1:length(normalize_pos)
                        sigmoid = 1./(1 + exp(-params.speed_gradient * (normalize_pos(i)*2 - 1) * 10));
                        updated_speed(i) = proximal_speed - (proximal_speed - terminal_speed) * sigmoid;
                    end
                    
                case 'Exponential'
                    % Exponential decay
                    for i = 1:length(normalize_pos)
                        decay_factor = 5 * params.speed_gradient;
                        updated_speed(i) = terminal_speed + (proximal_speed - terminal_speed) * exp(-decay_factor * normalize_pos(i));
                    end
                    
                case 'Step with Transition'
                    % Define regions
                    proximal_idx = normalize_pos < 0.3; 
                    middle_idx = normalize_pos >= 0.3 & normalize_pos < 0.85; 
                    transition_idx = normalize_pos >= 0.85 & normalize_pos <= 0.95; 
                    distal_idx = normalize_pos > 0.95; 
                    
                    % Set speeds for each region 
                    updated_speed(proximal_idx) = proximal_speed * 1.05;
                    
                    % Middle region - create a gentle decline
                    mid_pos = (normalize_pos(middle_idx) - 0.3) / (0.85 - 0.3);
                    updated_speed(middle_idx) = proximal_speed * 1.05 - mid_pos * (proximal_speed * 0.1);
                    
                    % Create S-curve transition
                    norm_trans_pos = (normalize_pos(transition_idx) - 0.85) / (0.95 - 0.85);
                    sigmoid_trans = 1./(1 + exp(-15 * (norm_trans_pos*2 - 1)));
                    
                    % Calculate connection points
                    start_speed = proximal_speed * 1.05 - 0.1 * proximal_speed;
                    end_speed = terminal_speed * 0.6;
                    
                    % Apply transition
                    updated_speed(transition_idx) = start_speed - sigmoid_trans * (start_speed - end_speed);
                    
                    % Terminal ileum - use parameter value directly with base speed
                    updated_speed(distal_idx) = terminal_speed * 0.6;
            end
    end
    
    % Apply bacterial factor (if needed)
    for i = 1:params.num_segments
        if intestine.is_terminal_ileum(i)
            reduction = (1 - params.bacterial_factor) * 0.5;
            updated_speed(i) = updated_speed(i) * (1 - reduction);
        end
    end
    
    % Add spatial variability across ALL segments including terminal ileum
    if params.speed_variability > 0
        window_size = 10;
        raw_noise = randn(1, params.num_segments + window_size*2) * params.speed_variability;
        smooth_noise = filter(ones(1, window_size) / window_size, 1, raw_noise);
        smooth_noise = smooth_noise(window_size+1:end-window_size);
        
        % Apply noise to the entire profile INCLUDING terminal ileum
        updated_speed = updated_speed .* (1 + smooth_noise);
        
        % Special handling for terminal ileum to ensure it gets enough variation
        terminal_idx = normalize_pos > 0.95;
        if any(terminal_idx)
            % Get the terminal ileum base speed (mean)
            terminal_base = mean(updated_speed(terminal_idx));
            
            % Create additional variation specific to terminal ileum
            % Using a higher variability factor for terminal ileum
            terminal_variation = 1.5 * params.speed_variability;
            terminal_noise = randn(1, sum(terminal_idx)) * terminal_variation;
            
            % Apply additional smoothing to the terminal noise
            terminal_smooth_window = min(5, floor(sum(terminal_idx)/2));
            if terminal_smooth_window > 1
                terminal_smooth_kernel = ones(1, terminal_smooth_window) / terminal_smooth_window;
                terminal_noise_padded = [terminal_noise, zeros(1, terminal_smooth_window)];
                smooth_terminal_noise = filter(terminal_smooth_kernel, 1, terminal_noise_padded);
                terminal_noise = smooth_terminal_noise(1:sum(terminal_idx));
            end
            
            % Apply the noise to terminal ileum section
            updated_speed(terminal_idx) = terminal_base * (1 + terminal_noise);
        end
        
        % Limit to reasonable range for the entire profile
        min_speed = terminal_speed * 0.5;
        max_speed = proximal_speed * 1.15;
        updated_speed = max(min_speed, min(max_speed, updated_speed));
    end
end

        function amplitude_profile = update_wave_amplitude_profile(app, intestine, params)
            % Generate spectrum-based amplitude profile during gastroenteritis
            
            % Base amplitude from gastroenteritis parameter
            amplitude_profile = ones(1, params.num_segments) * params.gastroenteritis.wave_amplitude;
            
            % Add position-dependent variation
            normalize_pos = intestine.x / params.intestine_length;
            
            % Apply position-dependent amplitude modulation
            for i = 1:params.num_segments
                % Higher amplitude in proximal, lower in distal parts
                if intestine.is_terminal_ileum(i)
                    amplitude_profile(i) = amplitude_profile(i) * 0.9;
                end
                
                % Add extra modulation at transition zones (areas of higher risk)
                distance_to_terminal = abs(intestine.x(i) - params.terminal_ileum_start);
                if distance_to_terminal < 10
                    % Enhance amplitude near terminal ileum transition
                    boundary_factor = exp(-0.05 * distance_to_terminal^2);
                    amplitude_profile(i) = amplitude_profile(i) * (1 + 0.2 * boundary_factor);
                end
            end
            
            % Add spatial variation if desired
            if params.speed_variability > 0
                window_size = 8;
                raw_noise = randn(1, params.num_segments + window_size*2) * params.speed_variability * 0.3;
                smooth_noise = filter(ones(1, window_size) / window_size, 1, raw_noise);
                smooth_noise = smooth_noise(window_size+1:end-window_size);
                
                amplitude_profile = amplitude_profile .* (1 + smooth_noise);
                
                % Bound the amplitude to reasonable values
                min_amp = 0.3;
                max_amp = params.gastroenteritis.wave_amplitude * 1.3;
                amplitude_profile = max(min_amp, min(max_amp, amplitude_profile));
            end
        end
        
        function intestine = update_peristalsis(app, intestine, t, params)
            % Update intestinal peristalsis state with elasticity consideration
            
            % Determine if a new peristaltic wave needs to be generated
            if intestine.is_gastroenteritis
                wave_freq = params.gastroenteritis.wave_frequency;
            else
                wave_freq = params.normal.wave_frequency;
            end
            
            if mod(t, 1/wave_freq) < params.dt
                intestine.wave_count = intestine.wave_count + 1;
            end
            
            % Reset intestinal radius to initial state (as baseline)
            intestine.radius = intestine.initial_radius;
            
            % Calculate peristaltic waves with spectrum-based speeds and elasticity
            for i = 1:params.num_segments
                x = intestine.x(i);
                
                % Get local properties
                local_speed = intestine.local_speed(i);
                local_amplitude = intestine.local_amplitude(i);
                local_elasticity = intestine.local_elasticity(i); % NEW: Get local elasticity
                
                % Calculate superposition effect of all waves at current position
                contraction = 0;
                for wave_num = 1:intestine.wave_count
                    % Wave start time
                    wave_start_time = (wave_num - 1) / wave_freq;
                    
                    % Skip if wave hasn't started yet
                    if t < wave_start_time
                        continue;
                    end
                    
                    % Calculate wave travel based on local speed profile
                    % Use piecewise integration of speed profile for more accurate wave propagation
                    % This accounts for varying speeds along the intestine
                    wave_distance = app.calculate_wave_position(intestine, wave_start_time, t, i, params);
                    
                    % Skip if wave hasn't reached current position or has passed far beyond
                    if wave_distance < x - 10 || wave_distance > x + 20
                        continue;
                    end
                    
                    % Calculate wave intensity at current position (Gaussian distribution)
                    % MODIFIED: Adjust wave intensity based on elasticity
                    elasticity_factor = 1.0 / (0.5 + local_elasticity); % Less elastic tissue has stronger contractions
                    wave_intensity = local_amplitude * elasticity_factor * exp(-0.1 * (x - wave_distance)^2);
                end
                
                % Limit maximum contraction
                contraction = min(contraction, 0.8);
                
                % Apply contraction effect to intestinal radius
                intestine.radius(i) = intestine.initial_radius(i) * (1 - contraction);
                
                % Update wall position (simulate physical displacement during peristalsis)
                % MODIFIED: Make displacement dependent on elasticity
                base_displacement = contraction * 0.5;  % Base wall displacement proportional to contraction
                
                % Modify displacement based on elasticity: less elastic = less displacement
                elasticity_modifier = (0.5 + local_elasticity);  % Scale factor based on elasticity
                displacement = base_displacement * elasticity_modifier;
                
                % Add enhanced distortion at areas with high speed or elasticity gradient
                if intestine.is_gastroenteritis && i > 1 && i < params.num_segments
                    % Calculate local speed gradient
                    speed_gradient = abs(intestine.local_speed(i) - intestine.local_speed(i-1));
                    
                    % NEW: Calculate local elasticity gradient
                    elasticity_gradient = abs(intestine.local_elasticity(i) - intestine.local_elasticity(i-1));
                    
                    % Add distortion proportional to speed gradient
                    if speed_gradient > 0.05  % Threshold to detect significant changes
                        gradient_factor = min(speed_gradient * 5, 3);  % Limit the effect
                        displacement = displacement * (1 + gradient_factor);
                    end
                    
                    % NEW: Add distortion proportional to elasticity gradient
                    if elasticity_gradient > 0.03  % Threshold for elasticity change
                        e_gradient_factor = min(elasticity_gradient * 8, 4);  % Elasticity changes can have stronger effects
                        displacement = displacement * (1 + e_gradient_factor);
                    end
                    
                    % Still check terminal ileum boundary as it remains an anatomical transition
                    distance_to_terminal = abs(intestine.x(i) - params.terminal_ileum_start);
                    if distance_to_terminal < 10
                        boundary_factor = exp(-0.05 * distance_to_terminal^2);
                        displacement = displacement * (1.5 * boundary_factor + 1);
                    end
                end
                
                % Update wall position
                intestine.wall_position(1, i) = intestine.x(i) - displacement;  % Inner wall
                intestine.wall_position(2, i) = intestine.x(i) + displacement;  % Outer wall
            end
        end
        
        function wave_pos = calculate_wave_position(app, intestine, wave_start_time, current_time, segment_idx, params)
            % Calculate wave position using speed profile integration
            % This accounts for different speeds along the intestine path
            
            elapsed_time = current_time - wave_start_time;
            
            % If elapsed time is very short, use simplified calculation
            if elapsed_time < 0.1
                wave_pos = intestine.local_speed(1) * elapsed_time;
                return;
            end
            
            % For more accurate tracking, integrate speed over position
            % This simulates a wave traveling through segments with different speeds
            
            % Start with an approximate position based on average speed
            avg_speed = mean(intestine.local_speed);
            approx_distance = avg_speed * elapsed_time;
            
            % Bound to intestine length
            approx_distance = min(approx_distance, params.intestine_length);
            
            % Find closest segment to this approximate position
            [~, approx_segment] = min(abs(intestine.x - approx_distance));
            
            % Refine the position by time integration from proximal end
            cumulative_time = 0;
            cumulative_distance = 0;
            
            % Forward integration from start
            for i = 1:params.num_segments-1
                segment_length = intestine.x(i+1) - intestine.x(i);
                segment_speed = (intestine.local_speed(i) + intestine.local_speed(i+1)) / 2;
                segment_time = segment_length / segment_speed;
                
                if cumulative_time + segment_time > elapsed_time
                    % Wave is within this segment - calculate precise position
                    remaining_time = elapsed_time - cumulative_time;
                    wave_pos = intestine.x(i) + remaining_time * segment_speed;
                    return;
                end
                
                cumulative_time = cumulative_time + segment_time;
                cumulative_distance = intestine.x(i+1);
                
                % If we've gone far enough, stop
                if cumulative_distance > params.intestine_length || i+1 >= params.num_segments
                    break;
                end
            end
            
            % If we get here, the wave has reached the end of the intestine
            wave_pos = params.intestine_length;
        end
        
        function intestine = update_contents(app, intestine, params)
            % Update intestinal content distribution with elasticity effects
            
            % Calculate content flow
            new_content = intestine.content;
            
            for i = 1:params.num_segments-1
                % Calculate pressure difference between adjacent segments (caused by radius difference)
                pressure_diff = intestine.radius(i) - intestine.radius(i+1);
                
                % Consider speed difference as additional factor in flow
                % This models how a speed transition zone can lead to content accumulation
                speed_ratio = max(1, intestine.local_speed(i) / max(0.1, intestine.local_speed(i+1)));
                
                % NEW: Consider elasticity difference as factor in flow
                elasticity_ratio = max(0.5, intestine.local_elasticity(i) / max(0.05, intestine.local_elasticity(i+1)));
                
                % Enhanced flow model that accounts for pressure, speed, and elasticity differences
                flow_factor = pressure_diff * speed_ratio * elasticity_ratio;
                
                % Consider ileocecal valve resistance (near end of intestine)
                resistance_factor = 1.0;
                if i >= 0.9 * params.num_segments
                    resistance_factor = params.valve_resistance;
                end
                
                % Calculate flow rate
                flow_rate = flow_factor * intestine.content(i) * params.dt / resistance_factor;
                
                % Ensure not to flow out too much content
                flow_rate = min(flow_rate, intestine.content(i));
                flow_rate = max(flow_rate, 0);  % Prevent backward flow
                
                % Update content
                new_content(i) = new_content(i) - flow_rate;
                new_content(i+1) = new_content(i+1) + flow_rate;
            end
            
            % Simulate secretion and absorption
            for i = 1:params.num_segments
                % Base secretion
                secretion = params.secretion_rate * params.dt;
                
                % Determine absorption rate
                if intestine.is_gastroenteritis
                    absorption = params.inflamed_absorption_rate * params.dt * intestine.content(i);
                else
                    absorption = params.absorption_rate * params.dt * intestine.content(i);
                end
                
                % Update content
                new_content(i) = new_content(i) + secretion - absorption;
            end
            
            % Add proximal content input
            if rand < 0.1  % Randomly add content
                new_content(1) = new_content(1) + 0.2;
            end
            
            % Ensure content amount is not negative
            new_content = max(new_content, 0);
            
            % Update content
            intestine.content = new_content;
        end
        
        function [intestine, detected, location, details] = check_intussusception(app, intestine, current_time, params)
            % Check for intussusception with enhanced detection for spectrum-based model
            % and with elasticity consideration
            detected = false;
            location = NaN;
            details = struct();
            
            % Initialize details structure
            details.time = current_time;
            details.speed_gradient = [];
            details.elasticity_gradient = [];  % NEW: Track elasticity gradient
            details.diameter_ratio = [];
            details.wall_overlap = [];
            details.critical_region = false;
            
            % Check the entire intestine for potential intussusception regions
            % This is more appropriate for spectrum model where transitions can occur anywhere
            for i = 1:params.num_segments-5
                % Calculate key factors that contribute to intussusception risk
                
                % 1. Diameter ratio between adjacent segments
                diameter_ratio = intestine.radius(i) / max(0.01, intestine.radius(i+1));
                
                % 2. Wall position overlap (physical invagination)
                wall_overlap = intestine.wall_position(2, i) - intestine.wall_position(1, i+1);
                
                % 3. Speed gradient (major contributor in spectrum model)
                if i < params.num_segments-1
                    speed_gradient = abs(intestine.local_speed(i) - intestine.local_speed(i+1));
                else
                    speed_gradient = 0;
                end
                
                % 4. NEW: Elasticity gradient (new contributor to risk)
                if i < params.num_segments-1
                    elasticity_gradient = abs(intestine.local_elasticity(i) - intestine.local_elasticity(i+1));
                else
                    elasticity_gradient = 0;
                end
                
                % Store values for highest risk segment
                if diameter_ratio > 1.5 || wall_overlap > 0 || speed_gradient > 0.2 || elasticity_gradient > 0.05
                    details.speed_gradient = speed_gradient;
                    details.elasticity_gradient = elasticity_gradient;  % Store elasticity gradient
                    details.diameter_ratio = diameter_ratio;
                    details.wall_overlap = wall_overlap;
                    
                    % 5. Check if this is in a critical anatomical region (e.g., near terminal ileum)
                    distance_to_terminal = abs(intestine.x(i) - params.terminal_ileum_start);
                    if distance_to_terminal < 15
                        details.critical_region = true;
                    end
                    
                    % Multi-factor detection: consider all factors together
                    % This better reflects real intussusception which requires multiple conditions
                    % NEW: Include elasticity gradient in risk calculation
                    risk_score = diameter_ratio * 0.4 + ...
                        (wall_overlap > 0) * 0.3 + ...
                        speed_gradient * 3.0 + ...
                        elasticity_gradient * 4.0; % Elasticity gradient has strong effect
                    
                    % Add additional risk if in critical region
                    if details.critical_region
                        risk_score = risk_score * 1.2;
                    end
                    
                    % If combined factors exceed threshold, detect intussusception
                    if (diameter_ratio > params.critical_diameter_ratio && wall_overlap > 0.1) || risk_score > 2.5
                        detected = true;
                        location = intestine.x(i);
                        
                        % Store detection details
                        details.position = location;
                        details.segment_index = i;
                        details.risk_score = risk_score;
                        
                        % Simulate wall position changes after intussusception (proximal enters distal)
                        overlap_length = min(5, params.num_segments - i);
                        
                        for j = 1:overlap_length
                            if i+j <= params.num_segments
                                intestine.wall_position(1, i+j) = intestine.wall_position(1, i);
                            end
                        end
                        
                        break;
                    end
                end
            end
        end
        
        function updateStatus(app, message)
            % Update status message in the app
            app.StatusLabel.Text = message;
            currentValue = app.ResultsTextArea.Value;
            if isempty(currentValue)
                app.ResultsTextArea.Value = {message};
            else
                app.ResultsTextArea.Value = [currentValue; {message}];
            end
            drawnow;
        end
        
        function plotResults(app)
            % Plot simulation results in the app
            
            % Check if params is null
            if isempty(app.params) || ~isfield(app.params, 'intestine_length')
                app.updateStatus('Error: Simulation parameters not initialized correctly');
                return;
            end
            
            % Get results and params
            results = app.results;
            params = app.params;

            % Get results and params
            app.results = app.results;
            app.params = app.params;
            
            % Create new figure for advanced visualization
            fig = figure('Position', [100, 100, 1200, 800], 'Visible', 'on');
            
            % Display intestinal radius changes over time and position
            subplot(5, 1, 1);  % Changed from 4,1,1 to 5,1,1 to add elasticity profile plot
            imagesc(linspace(0, params.intestine_length, params.num_segments), ...
                    results.time, results.radius);
            hold on;
            
            % Mark terminal ileum position
            line([params.terminal_ileum_start, params.terminal_ileum_start], [0, params.simulation_time], ...
                 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
            
            % Mark gastroenteritis onset time
            line([0, params.intestine_length], [params.gastroenteritis_onset, params.gastroenteritis_onset], ...
                 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
            
            % If intussusception detected, mark position
            if results.intussusception_detected
                plot(results.intussusception_location, results.intussusception_time, 'ro', ...
                     'MarkerSize', 10, 'MarkerFaceColor', 'r');
                text(results.intussusception_location + 2, results.intussusception_time, ...
                     'Intussusception', 'Color', 'r', 'FontWeight', 'bold');
            end
            
            colormap('jet');
            colorbar;
            title('Intestinal Radius Changes');
            xlabel('Intestinal Position (cm)');
            ylabel('Time (seconds)');
            
            % Add labels
            text(params.terminal_ileum_start + 5, 10, 'Terminal Ileum', 'Color', 'r', 'FontWeight', 'bold');
            text(10, params.gastroenteritis_onset + 5, 'Gastroenteritis Onset', 'Color', 'g', 'FontWeight', 'bold');
            
            % Display content distribution
            subplot(5, 1, 2);
            imagesc(linspace(0, params.intestine_length, params.num_segments), ...
                    results.time, results.content);
            hold on;
            
            % Mark terminal ileum position
            line([params.terminal_ileum_start, params.terminal_ileum_start], [0, params.simulation_time], ...
                 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
            
            % Mark gastroenteritis onset time
            line([0, params.intestine_length], [params.gastroenteritis_onset, params.gastroenteritis_onset], ...
                 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
            
            % If intussusception detected, mark position
            if results.intussusception_detected
                plot(results.intussusception_location, results.intussusception_time, 'ro', ...
                     'MarkerSize', 10, 'MarkerFaceColor', 'r');
            end
            
            colormap('jet');
            colorbar;
            title('Intestinal Content Distribution');
            xlabel('Intestinal Position (cm)');
            ylabel('Time (seconds)');
            
            % Display speed profile evolution
            subplot(5, 1, 3);
            imagesc(linspace(0, params.intestine_length, params.num_segments), ...
                    results.time, results.speed_profile);
            hold on;
            
            % Mark terminal ileum position
            line([params.terminal_ileum_start, params.terminal_ileum_start], [0, params.simulation_time], ...
                 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
            
            % Mark gastroenteritis onset time
            line([0, params.intestine_length], [params.gastroenteritis_onset, params.gastroenteritis_onset], ...
                 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
            
            % If intussusception detected, mark position
            if results.intussusception_detected
                plot(results.intussusception_location, results.intussusception_time, 'ro', ...
                     'MarkerSize', 10, 'MarkerFaceColor', 'r');
            end
            
            colormap('jet');
            colorbar;
            title('Intestinal Speed Profile (Spectrum Distribution)');
            xlabel('Intestinal Position (cm)');
            ylabel('Time (seconds)');
            
            % NEW: Display elasticity profile evolution
            subplot(5, 1, 4);
            imagesc(linspace(0, params.intestine_length, params.num_segments), ...
                    results.time, results.elasticity_profile);
            hold on;
            
            % Mark terminal ileum position
            line([params.terminal_ileum_start, params.terminal_ileum_start], [0, params.simulation_time], ...
                 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
            
            % Mark gastroenteritis onset time
            line([0, params.intestine_length], [params.gastroenteritis_onset, params.gastroenteritis_onset], ...
                 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
            
            % If intussusception detected, mark position
            if results.intussusception_detected
                plot(results.intussusception_location, results.intussusception_time, 'ro', ...
                     'MarkerSize', 10, 'MarkerFaceColor', 'r');
            end
            
            colormap('jet');
            colorbar;
            title('Intestinal Elasticity Profile');
            xlabel('Intestinal Position (cm)');
            ylabel('Time (seconds)');
            
            % Display intestinal shape at specific time
            subplot(5, 1, 5);
            
            % If intussusception occurred, show intestinal shape at that time
            if results.intussusception_detected
                time_idx = find(results.time >= results.intussusception_time, 1);
            else
                time_idx = round(size(results.radius, 1) * 0.8);  % Show shape at later time point
            end
            
            % Draw intestinal shape
            x_vals = linspace(0, params.intestine_length, params.num_segments);
            y_upper = squeeze(results.wall_position(time_idx, 2, :))';  % Transpose to ensure correct dimension
            y_lower = squeeze(results.wall_position(time_idx, 1, :))' * -1;  % Flip lower boundary for visualization
            
            % Ensure vector dimensions are consistent
            if size(y_upper, 1) > 1
                y_upper = y_upper';
            end
            if size(y_lower, 1) > 1
                y_lower = y_lower';
            end
            
            % Fill intestinal area
            fill([x_vals, fliplr(x_vals)], [y_upper, fliplr(y_lower)], [0.8, 0.8, 0.9]);
            hold on;
            
            % Draw intestinal walls
            plot(x_vals, y_upper, 'b-', 'LineWidth', 2);
            plot(x_vals, y_lower, 'b-', 'LineWidth', 2);
            
            % Mark terminal ileum
            line([params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [min(y_lower)-1, max(y_upper)+1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
            text(params.terminal_ileum_start + 2, max(y_upper) + 0.5, 'Terminal Ileum', 'Color', 'r');
            
            % Mark intussusception location
            if results.intussusception_detected
                plot(results.intussusception_location, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
                text(results.intussusception_location + 2, 0.5, 'Intussusception', 'Color', 'r', 'FontWeight', 'bold');
            end
            
            title(sprintf('Intestinal Shape at Time %.1f seconds', results.time(time_idx)));
            xlabel('Intestinal Position (cm)');
            ylabel('Wall Position');
            axis([0, params.intestine_length, min(y_lower)-1, max(y_upper)+1]);
            grid on;
            
            % Add summary
            if results.intussusception_detected
                sgtitle(sprintf('Intussusception Simulation Results: Occurred at %.1f seconds at %.1f cm', ...
                               results.intussusception_time, results.intussusception_location), 'FontSize', 14);
            else
                sgtitle('Intussusception Simulation Results: No intussusception occurred during simulation time', 'FontSize', 14);
            end
            
            % Also display a simplified version in the app's plot area
            cla(app.ResultPlotArea);
            axes(app.ResultPlotArea);
            
            % If intussusception occurred, show intestinal shape at that time
            if results.intussusception_detected
                time_idx = find(results.time >= results.intussusception_time, 1);
            else
                time_idx = round(size(results.radius, 1) * 0.8);  % Show shape at later time point
            end
            
            % Draw intestinal shape
            x_vals = linspace(0, params.intestine_length, params.num_segments);
            y_upper = squeeze(results.wall_position(time_idx, 2, :))';
            y_lower = squeeze(results.wall_position(time_idx, 1, :))' * -1;
            
            % Ensure vector dimensions are consistent
            if size(y_upper, 1) > 1
                y_upper = y_upper';
            end
            if size(y_lower, 1) > 1
                y_lower = y_lower';
            end
            
            % Fill intestinal area
            fill(app.ResultPlotArea, [x_vals, fliplr(x_vals)], [y_upper, fliplr(y_lower)], [0.8, 0.8, 0.9]);
            hold(app.ResultPlotArea, 'on');
            
            % Draw intestinal walls
            plot(app.ResultPlotArea, x_vals, y_upper, 'b-', 'LineWidth', 2);
            plot(app.ResultPlotArea, x_vals, y_lower, 'b-', 'LineWidth', 2);
            
            % Mark terminal ileum
            line(app.ResultPlotArea, [params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [min(y_lower)-1, max(y_upper)+1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
            
            % Mark intussusception location
            if results.intussusception_detected
                plot(app.ResultPlotArea, results.intussusception_location, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
                text(app.ResultPlotArea, results.intussusception_location, 0.5, 'Intussusception', 'Color', 'r', 'FontWeight', 'bold');
            end
            
            title(app.ResultPlotArea, 'Intestinal Shape');
            xlabel(app.ResultPlotArea, 'Intestinal Position (cm)');
            ylabel(app.ResultPlotArea, 'Wall Position');
            grid(app.ResultPlotArea, 'on');
            
            % Plot speed profile in the app
            cla(app.SpeedProfilePlotArea);
            axes(app.SpeedProfilePlotArea);

            % Show pre and post gastroenteritis speed profiles
            pre_gastro_idx = find(results.time < params.gastroenteritis_onset, 1, 'last');
            if isempty(pre_gastro_idx)
                pre_gastro_idx = 1;
            end
            
            post_gastro_idx = find(results.time >= params.gastroenteritis_onset, 1);
            if isempty(post_gastro_idx)
                post_gastro_idx = size(results.speed_profile, 1);
            end
            
            plot(app.SpeedProfilePlotArea, x_vals, results.speed_profile(pre_gastro_idx, :), 'b-', 'LineWidth', 2);
            hold(app.SpeedProfilePlotArea, 'on');
            plot(app.SpeedProfilePlotArea, x_vals, results.speed_profile(post_gastro_idx, :), 'r-', 'LineWidth', 2);


            % Mark terminal ileum
            line(app.SpeedProfilePlotArea, [params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [0, max(max(results.speed_profile))*1.1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
                
            % Mark intussusception location if detected
            if results.intussusception_detected
                plot(app.SpeedProfilePlotArea, results.intussusception_location, ...
                     results.speed_profile(find(results.time >= results.intussusception_time, 1), ...
                     round(results.intussusception_location/params.intestine_length*params.num_segments)), ...
                     'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            end
            
            title(app.SpeedProfilePlotArea, 'Intestinal Speed Profile');
            xlabel(app.SpeedProfilePlotArea, 'Intestinal Position (cm)');
            ylabel(app.SpeedProfilePlotArea, 'Peristaltic Speed (cm/s)');
            legend(app.SpeedProfilePlotArea, {'Normal', 'Gastroenteritis'}, 'Location', 'best');
            grid(app.SpeedProfilePlotArea, 'on');
            
            % NEW: Plot elasticity profile in the app
            cla(app.ElasticityProfilePlotArea);
            axes(app.ElasticityProfilePlotArea);
            
            % Show pre and post gastroenteritis elasticity profiles
            plot(app.ElasticityProfilePlotArea, x_vals, results.elasticity_profile(pre_gastro_idx, :), 'b-', 'LineWidth', 2);
            hold(app.ElasticityProfilePlotArea, 'on');
            
            if post_gastro_idx <= size(results.elasticity_profile, 1)
                plot(app.ElasticityProfilePlotArea, x_vals, results.elasticity_profile(post_gastro_idx, :), 'r-', 'LineWidth', 2);
            end
            
            % Mark terminal ileum
            line(app.ElasticityProfilePlotArea, [params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [0, max(max(results.elasticity_profile))*1.1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
                
            % Mark intussusception location if detected
            if results.intussusception_detected
                intus_idx = find(results.time >= results.intussusception_time, 1);
                if ~isempty(intus_idx) && intus_idx <= size(results.elasticity_profile, 1)
                    loc_idx = round(results.intussusception_location/params.intestine_length*params.num_segments);
                    if loc_idx > 0 && loc_idx <= params.num_segments
                        plot(app.ElasticityProfilePlotArea, results.intussusception_location, ...
                             results.elasticity_profile(intus_idx, loc_idx), ...
                             'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                    end
                end
            end
            
            title(app.ElasticityProfilePlotArea, 'Intestinal Elasticity Profile');
            xlabel(app.ElasticityProfilePlotArea, 'Intestinal Position (cm)');
            ylabel(app.ElasticityProfilePlotArea, 'Elasticity');
            legend(app.ElasticityProfilePlotArea, {'Normal', 'Gastroenteritis'}, 'Location', 'best');
            grid(app.ElasticityProfilePlotArea, 'on');
        end
        
        function saveResults(app)
            % Save simulation results to a file
            [filename, pathname] = uiputfile('*.mat', 'Save Simulation Results');
            if ~isequal(filename, 0) && ~isequal(pathname, 0)
                results = app.results;
                params = app.params;
                save(fullfile(pathname, filename), 'results', 'params');
                app.updateStatus(['Results saved to ' fullfile(pathname, filename)]);
            end
        end
        
        function exportResultsTable(app)
            % Create a summary figure with 5 subplots focusing on spectrum-based analysis and elasticity
            figSummary = figure('Name', 'Intussusception Simulation Summary (Spectrum Model with Elasticity)', 'Position', [100, 100, 1200, 900]);
            
            % Get data needed for plots
            results = app.results;
            params = app.params;
            
            % Define colors for plots
            color_normal = [0.0, 0.6, 1.0];     % Blue for normal state
            color_gastro = [1.0, 0.4, 0.4];      % Red for gastroenteritis
            color_intussusception = [0.9, 0.0, 0.0]; % Darker red for intussusception
                        
            % Compute data for the plots
            x_intestine = linspace(0, params.intestine_length, params.num_segments);
            
            % Determine location type if intussusception occurred
            if results.intussusception_detected
                distance_to_terminal = abs(results.intussusception_location - params.terminal_ileum_start);
                if distance_to_terminal < 5
                    locationType = 'Terminal Ileum Junction';
                elseif results.intussusception_location < params.terminal_ileum_start
                    locationType = 'Proximal Small Intestine';
                else
                    locationType = 'Terminal Ileum';
                end
                
                % Create title
                titleStr = sprintf('Intussusception Simulation Results: Occurred at %.1f seconds at %.1f cm (%s)', ...
                    results.intussusception_time, results.intussusception_location, locationType);
            else
                titleStr = 'Intussusception Simulation Results: No intussusception occurred';
            end
            
            % Set figure title
            sgtitle(titleStr, 'FontSize', 14);
            
            % Plot 1: Speed Profile Evolution
            subplot(3, 2, 1);
            
            % Get speed profiles at key time points
            initial_idx = 1;  % Initial state
            
            % Find index just before gastroenteritis onset
            pre_gastro_idx = find(results.time < params.gastroenteritis_onset, 1, 'last');
            if isempty(pre_gastro_idx)
                pre_gastro_idx = 1;
            end
            
            % Find index right after gastroenteritis onset
            post_gastro_idx = find(results.time >= params.gastroenteritis_onset, 1);
            if isempty(post_gastro_idx)
                post_gastro_idx = size(results.speed_profile, 1);
            end
            
            % If intussusception occurred, find that index
            if results.intussusception_detected
                intus_idx = find(results.time >= results.intussusception_time, 1);
                if isempty(intus_idx)
                    intus_idx = size(results.speed_profile, 1);
                end
            else
                intus_idx = size(results.speed_profile, 1);  % Use final state
            end
            
            % Plot speed profiles
            plot(x_intestine, results.speed_profile(initial_idx, :), '-', 'Color', color_normal, 'LineWidth', 2);
            hold on;
            plot(x_intestine, results.speed_profile(post_gastro_idx, :), '-', 'Color', color_gastro, 'LineWidth', 2);
            
            if results.intussusception_detected && intus_idx ~= post_gastro_idx
                plot(x_intestine, results.speed_profile(intus_idx, :), '-', 'Color', color_intussusception, 'LineWidth', 2);
            end

            % Add terminal ileum line
            line([params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [0, max(max(results.speed_profile))*1.1], ...
                 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
            
            % Mark intussusception point if detected
            if results.intussusception_detected
                plot(results.intussusception_location, ...
                     results.speed_profile(intus_idx, round(results.intussusception_location/params.intestine_length*params.num_segments)), ...
                     'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
            end
            
            title('Peristaltic Speed Profile Evolution');
            xlabel('Intestinal Position (cm)');
            ylabel('Speed (cm/s)');
            if results.intussusception_detected && intus_idx ~= post_gastro_idx
                legend({'Initial', 'After Gastroenteritis', 'At Intussusception'}, 'Location', 'best');
            else
                legend({'Initial', 'After Gastroenteritis'}, 'Location', 'best');
            end
            grid on;
            
            % Plot 2: Elasticity Profile Evolution (NEW)
            subplot(3, 2, 2);
            
            % Plot elasticity profiles
            plot(x_intestine, results.elasticity_profile(initial_idx, :), '-', 'Color', color_normal, 'LineWidth', 2);
            hold on;
            plot(x_intestine, results.elasticity_profile(post_gastro_idx, :), '-', 'Color', color_gastro, 'LineWidth', 2);
            
            if results.intussusception_detected && intus_idx ~= post_gastro_idx
                plot(x_intestine, results.elasticity_profile(intus_idx, :), '-', 'Color', color_intussusception, 'LineWidth', 2);
            end
            
            % Add terminal ileum line
            line([params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [0, max(max(results.elasticity_profile))*1.1], ...
                 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
            
            % Mark intussusception point if detected
            if results.intussusception_detected
                plot(results.intussusception_location, ...
                     results.elasticity_profile(intus_idx, round(results.intussusception_location/params.intestine_length*params.num_segments)), ...
                     'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
            end
            
            title('Intestinal Elasticity Profile Evolution');
            xlabel('Intestinal Position (cm)');
            ylabel('Elasticity');
            if results.intussusception_detected && intus_idx ~= post_gastro_idx
                legend({'Initial', 'After Gastroenteritis', 'At Intussusception'}, 'Location', 'best');
            else
                legend({'Initial', 'After Gastroenteritis'}, 'Location', 'best');
            end
            grid on;
            
            % Plot 3: Speed and Elasticity Gradients (key factors in spectrum model)
            subplot(3, 2, 3);
            
            % Calculate speed and elasticity gradients at each position
            speed_gradients = zeros(params.num_segments-1, 1);
            elasticity_gradients = zeros(params.num_segments-1, 1);
            
            for i = 1:params.num_segments-1
                speed_gradients(i) = abs(results.speed_profile(post_gastro_idx, i+1) - results.speed_profile(post_gastro_idx, i));
                elasticity_gradients(i) = abs(results.elasticity_profile(post_gastro_idx, i+1) - results.elasticity_profile(post_gastro_idx, i));
            end
            
            % Plot speed gradient
            yyaxis left;
            bar(x_intestine(1:end-1), speed_gradients, 'FaceColor', [0.8, 0.4, 0.4], 'FaceAlpha', 0.6);
            ylabel('Speed Gradient (cm/s per segment)');
            
            % Plot elasticity gradient on secondary axis
            yyaxis right;
            bar(x_intestine(1:end-1), elasticity_gradients, 'FaceColor', [0.4, 0.4, 0.8], 'FaceAlpha', 0.6);
            ylabel('Elasticity Gradient');
            
            % Add terminal ileum line
            line([params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [0, max(elasticity_gradients)*1.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
            
            % Mark intussusception point if detected
            if results.intussusception_detected
                iloc_idx = min(length(speed_gradients), round(results.intussusception_location/params.intestine_length*params.num_segments));
                if iloc_idx > 0
                    yyaxis right;
                    plot(results.intussusception_location, elasticity_gradients(iloc_idx), ...
                         'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
                end
            end
            
            title('Peristaltic Speed and Elasticity Gradients (Risk Factors)');
            xlabel('Intestinal Position (cm)');
            legend({'Speed Gradient', 'Elasticity Gradient'}, 'Location', 'best');
            grid on;
            
            % Plot 4: Content Accumulation vs Elasticity Profile
            subplot(3, 2, 4);
            
            % Calculate content accumulation
            content_accumulation = zeros(1, params.num_segments);
            
            if results.intussusception_detected
                pre_intus_indices = results.time < results.intussusception_time;
                post_gastro_indices = results.time >= params.gastroenteritis_onset & results.time < results.intussusception_time;
                
                if sum(post_gastro_indices) > 0
                    for i = 1:params.num_segments
                        content_accumulation(i) = mean(results.content(post_gastro_indices, i)) - mean(results.content(pre_intus_indices, i));
                    end
                end
            else
                pre_gastro_indices = results.time < params.gastroenteritis_onset;
                post_gastro_indices = results.time >= params.gastroenteritis_onset;
                
                if sum(post_gastro_indices) > 0 && sum(pre_gastro_indices) > 0
                    for i = 1:params.num_segments
                        content_accumulation(i) = mean(results.content(post_gastro_indices, i)) - mean(results.content(pre_gastro_indices, i));
                    end
                end
            end
            
            % Create bar chart for content accumulation
            bar(x_intestine, content_accumulation, 'FaceColor', [0.6, 0.6, 0.8]);
            hold on;
            
            % Plot elasticity profile on secondary axis
            yyaxis right;
            plot(x_intestine, results.elasticity_profile(post_gastro_idx, :), '-', 'Color', [0.8, 0.4, 0.4], 'LineWidth', 2);
            ylabel('Elasticity');
            
            % Switch back to left axis for additional elements
            yyaxis left;
            ylabel('Content Accumulation');
            
            % Add terminal ileum line
            line([params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [min(content_accumulation)-0.1, max(content_accumulation)+0.1], ...
                 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
            
            % Mark intussusception point if detected
            if results.intussusception_detected
                plot(results.intussusception_location, content_accumulation(round(results.intussusception_location/params.intestine_length*params.num_segments)), ...
                     'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
            end
            
            title('Content Accumulation vs. Elasticity Profile');
            xlabel('Intestinal Position (cm)');
            legend({'Content Accumulation', 'Elasticity Profile'}, 'Location', 'best');
            grid on;
            
            % Plot 5: Enhanced Intussusception Risk Index (with elasticity)
            subplot(3, 2, [5, 6]);
            
            % Calculate risk index based on multiple factors including elasticity
            risk_index = zeros(1, params.num_segments);
            
            for i = 1:params.num_segments-1
                % Speed gradient factor
                speed_grad = abs(results.speed_profile(post_gastro_idx, min(i+1, params.num_segments)) - results.speed_profile(post_gastro_idx, i));
                
                % Elasticity gradient factor (new)
                elasticity_grad = abs(results.elasticity_profile(post_gastro_idx, min(i+1, params.num_segments)) - results.elasticity_profile(post_gastro_idx, i));
                
                % Content accumulation factor
                content_factor = max(0, content_accumulation(i));
                
                % Anatomical position factor (higher risk near terminal ileum)
                distance_to_terminal = abs(x_intestine(i) - params.terminal_ileum_start);
                position_factor = exp(-0.01 * distance_to_terminal^2);
                
                % Elasticity value factor (low elasticity increases risk)
                elasticity_value_factor = max(0, (0.5 - results.elasticity_profile(post_gastro_idx, i)) * 2);
                
                % Combine factors with appropriate weighting including elasticity
                risk_index(i) = speed_grad * 3.0 + ...% Speed gradient
                elasticity_grad * 4.0 + ...% Elasticity gradient (stronger effect)
                content_factor * 1.0 + ...% Content accumulation
                position_factor * 0.5 + ...% Anatomical position
                elasticity_value_factor * 2.0; % Low elasticity factor
            end
            
            % Normalize risk index to 0-1 range
            if max(risk_index) > 0
                risk_index = risk_index / max(risk_index);
            end
            
            % Plot risk index
            area(x_intestine, risk_index, 'FaceColor', [1, 0.8, 0.8], 'EdgeColor', 'r');
            hold on;
            
            % Overlay line plot
            plot(x_intestine, risk_index, 'r-', 'LineWidth', 1.5);
            
            % Add terminal ileum line
            line([params.terminal_ileum_start, params.terminal_ileum_start], ...
                 [0, max(risk_index)*1.1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
            
            % Mark intussusception point if detected
            if results.intussusception_detected
                plot(results.intussusception_location, risk_index(round(results.intussusception_location/params.intestine_length*params.num_segments)), ...
                     'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
                text(results.intussusception_location, risk_index(round(results.intussusception_location/params.intestine_length*params.num_segments)), ...
                     sprintf('Intussusception at %.1fs', results.intussusception_time), ...
                     'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
            end
            
            title('Enhanced Intussusception Risk Index (With Elasticity Factors)');
            xlabel('Intestinal Position (cm)');
            ylabel('Risk Index (normalized)');
            grid on;
            
            % Add export button
            exportBtn = uicontrol(figSummary, 'Style', 'pushbutton', 'String', 'Export as Image', ...
                                 'Position', [850, 10, 120, 25], 'Callback', @exportFigureToImage);

            % Add export data button
            exportDataBtn = uicontrol(figSummary, 'Style', 'pushbutton', 'String', 'Export Data to CSV', ...
                                     'Position', [710, 10, 120, 25], 'Callback', @exportDataToCSV);

function exportFigureToImage(src, event)
    [filename, pathname] = uiputfile({'*.png';'*.jpg';'*.pdf'}, 'Save Results Figure');
    if ~isequal(filename, 0) && ~isequal(pathname, 0)
        fullPath = fullfile(pathname, filename);
        
        try
            % Check if exportapp function is available (MATLAB R2020a or newer)
            if exist('exportapp', 'file') == 2
                % Use exportapp to properly capture UI components
                exportapp(figSummary, fullPath);
                app.updateStatus(['Results figure successfully exported to ' fullPath]);
            else
                % Alternative approach that doesn't try to copy axes
                % This avoids the "multiple coordinate systems" error
                
                % Use print function directly on the original figure
                print(figSummary, fullPath, '-dpng', '-r300');
                
                app.updateStatus(['Results figure exported to ' fullPath]);
                
                % Inform user about possible limitations
                uiwait(msgbox(['The figure has been exported, but UI components may not be included. ' ...
                      'For best results, consider updating to MATLAB R2020a or newer.'], ...
                      'Export Information', 'modal'));
            end
        catch exportErr
            % If even the fallback method fails, use the most direct approach
            try
                % Just capture the figure as-is with print
                print(figSummary, fullPath, '-dpng', '-r300');
                app.updateStatus(['Basic figure export completed to ' fullPath]);
            catch
                % Last resort error reporting
                errMsg = ['Could not export figure: ' exportErr.message];
                app.updateStatus(errMsg);
                errordlg(errMsg, 'Export Error');
            end
        end
    end
end

        function exportDataToCSV(src, event)
    [filename, pathname] = uiputfile('*.xlsx', 'Save Results Data as Excel');
    if ~isequal(filename, 0) && ~isequal(pathname, 0)
        fullPath = fullfile(pathname, filename);
        
        % Create table for simulation parameters and results
        paramNames = {'Simulation Time (s)', 'Gastroenteritis Onset (s)', ...
            'Proximal Speed (cm/s)', 'Terminal Ileum Speed (cm/s)', ...
            'Speed Transition Type', 'Speed Gradient', ...
            'Speed Variability', 'Noise Level', 'Speed Pattern', ...
            'Normal Proximal Elasticity', 'Normal Distal Elasticity', ...
            'Gastro Proximal Elasticity', 'Gastro Distal Elasticity', ...
            'Elasticity Pattern', 'Elasticity Gradient', 'Elasticity Variability', ...
            'Bacterial Factor', 'Valve Resistance', ...
            'Intestine Length (cm)', 'Terminal Ileum Start (cm)', ...
            'Critical Diameter Ratio'};
        
        paramValues = {params.simulation_time, params.gastroenteritis_onset, ...
            params.gastroenteritis.proximal_speed, params.gastroenteritis.terminal_speed, ...
            params.speed_transition_type, params.speed_gradient, ...
            params.speed_variability, params.noise_level, params.speed_pattern, ...
            params.normal.proximal_elasticity, params.normal.distal_elasticity, ...
            params.gastroenteritis.proximal_elasticity, params.gastroenteritis.distal_elasticity, ...
            params.elasticity_pattern, params.elasticity_gradient, params.elasticity_variability, ...
            params.bacterial_factor, params.valve_resistance, ...
            params.intestine_length, params.terminal_ileum_start, ...
            params.critical_diameter_ratio};
        
        if results.intussusception_detected
            % Get detailed information about the intussusception event
            details = results.intussusception_details;
            
            resultNames = {'Outcome', 'Intussusception Time (s)', 'Intussusception Location (cm)', ...
                'Location Type', 'Speed Gradient at Location', 'Elasticity Gradient at Location', ...
                'Diameter Ratio at Location', 'Wall Overlap at Location', 'Critical Region', 'Risk Score'};
            
            % Extract details if available
            if isfield(details, 'speed_gradient')
                speed_gradient_val = details.speed_gradient;
            else
                speed_gradient_val = NaN;
            end
            
            if isfield(details, 'elasticity_gradient')
                elasticity_gradient_val = details.elasticity_gradient;
            else
                elasticity_gradient_val = NaN;
            end
            
            if isfield(details, 'diameter_ratio')
                diameter_ratio_val = details.diameter_ratio;
            else
                diameter_ratio_val = NaN;
            end
            
            if isfield(details, 'wall_overlap')
                wall_overlap_val = details.wall_overlap;
            else
                wall_overlap_val = NaN;
            end
            
            if isfield(details, 'critical_region')
                critical_region_val = details.critical_region;
            else
                critical_region_val = false;
            end
            
            if isfield(details, 'risk_score')
                risk_score_val = details.risk_score;
            else
                risk_score_val = NaN;
            end
            
            resultValues = {'Intussusception Detected', results.intussusception_time, ...
                results.intussusception_location, locationType, ...
                speed_gradient_val, elasticity_gradient_val, ...
                diameter_ratio_val, wall_overlap_val, critical_region_val, risk_score_val};
        else
            resultNames = {'Outcome'};
            resultValues = {'No Intussusception Detected'};
        end
        
        allNames = [paramNames, resultNames];
        allValues = [paramValues, resultValues];
        
        Summary = table(allNames', allValues', 'VariableNames', {'Parameter_Result', 'Value'});
        
        % Calculate full risk index based on key factors (similar to plotResults but for Excel output)
        risk_index = zeros(1, params.num_segments);
        
        % Get time index at gastroenteritis state
        post_gastro_idx = find(results.time >= params.gastroenteritis_onset, 1);
        if isempty(post_gastro_idx)
            post_gastro_idx = size(results.speed_profile, 1);
        end
        
        % Calculate risk based on multiple factors including elasticity
        x_intestine = linspace(0, params.intestine_length, params.num_segments);
        
        for i = 1:params.num_segments-1
            % Speed gradient factor
            speed_grad = abs(results.speed_profile(post_gastro_idx, min(i+1, params.num_segments)) - results.speed_profile(post_gastro_idx, i));
            
            % Elasticity gradient factor 
            elasticity_grad = abs(results.elasticity_profile(post_gastro_idx, min(i+1, params.num_segments)) - results.elasticity_profile(post_gastro_idx, i));
            
            % Content accumulation factor (based on average content in gastroenteritis phase if available)
            if results.intussusception_detected
                pre_intus_indices = results.time < results.intussusception_time;
                post_gastro_indices = results.time >= params.gastroenteritis_onset & results.time < results.intussusception_time;
                
                if sum(post_gastro_indices) > 0 && sum(pre_intus_indices) > 0
                    content_factor = max(0, mean(results.content(post_gastro_indices, i)) - mean(results.content(pre_intus_indices, i)));
                else
                    content_factor = 0;
                end
            else
                pre_gastro_indices = results.time < params.gastroenteritis_onset;
                post_gastro_indices = results.time >= params.gastroenteritis_onset;
                
                if sum(post_gastro_indices) > 0 && sum(pre_gastro_indices) > 0
                    content_factor = max(0, mean(results.content(post_gastro_indices, i)) - mean(results.content(pre_gastro_indices, i)));
                else
                    content_factor = 0;
                end
            end
            
            % Anatomical position factor (higher risk near terminal ileum)
            distance_to_terminal = abs(x_intestine(i) - params.terminal_ileum_start);
            position_factor = exp(-0.01 * distance_to_terminal^2);
            
            % Elasticity value factor (low elasticity increases risk)
            elasticity_value_factor = max(0, (0.5 - results.elasticity_profile(post_gastro_idx, i)) * 2);
            
            % Combine factors with appropriate weighting including elasticity
            risk_index(i) = speed_grad * 3.0 + ...                 % Speed gradient
                            elasticity_grad * 4.0 + ...            % Elasticity gradient
                            content_factor * 1.0 + ...             % Content accumulation
                            position_factor * 0.5 + ...            % Anatomical position
                            elasticity_value_factor * 2.0;         % Low elasticity factor
        end
        
        % Normalize risk index to 0-10 range (easier to interpret)
        if max(risk_index) > 0
            risk_index = risk_index / max(risk_index) * 10;
        end
        
        % Create 5cm intervals for the risk table
        interval_size = 5; % 5 cm intervals
        num_intervals = ceil(params.intestine_length / interval_size);
        
        % Arrays to hold interval data
        interval_start = zeros(num_intervals, 1);
        interval_end = zeros(num_intervals, 1);
        interval_risk = zeros(num_intervals, 1);
        interval_risk_category = cell(num_intervals, 1);
        
        % Calculate average risk for each interval
        for i = 1:num_intervals
            interval_start(i) = (i-1) * interval_size;
            interval_end(i) = min(i * interval_size, params.intestine_length);
            
            % Find segments within this interval
            segment_indices = x_intestine >= interval_start(i) & x_intestine < interval_end(i);
            if any(segment_indices)
                interval_risk(i) = mean(risk_index(segment_indices));
            else
                interval_risk(i) = 0;
            end
            
            % Assign risk category based on risk value
            if interval_risk(i) < 2
                interval_risk_category{i} = 'Very Low';
            elseif interval_risk(i) < 4
                interval_risk_category{i} = 'Low';
            elseif interval_risk(i) < 6
                interval_risk_category{i} = 'Moderate';
            elseif interval_risk(i) < 8
                interval_risk_category{i} = 'High';
            else
                interval_risk_category{i} = 'Very High';
            end
        end
        
        % Add anatomical landmarks information
        anatomical_region = cell(num_intervals, 1);
        for i = 1:num_intervals
            mid_point = (interval_start(i) + interval_end(i)) / 2;
            
            if mid_point < params.intestine_length * 0.1
                anatomical_region{i} = 'Duodenum';
            elseif mid_point < params.intestine_length * 0.6
                anatomical_region{i} = 'Jejunum';
            elseif mid_point < params.terminal_ileum_start
                anatomical_region{i} = 'Ileum';
            else
                anatomical_region{i} = 'Terminal Ileum';
            end
        end
        
        % Create risk table
        RiskDistribution = table(interval_start, interval_end, anatomical_region, interval_risk, interval_risk_category, ...
            'VariableNames', {'Start_Position_cm', 'End_Position_cm', 'Anatomical_Region', 'Risk_Value_0_10', 'Risk_Category'});
        
        % Export to Excel with two sheets
        writetable(Summary, fullPath, 'Sheet', 'Simulation Summary');
        writetable(RiskDistribution, fullPath, 'Sheet', 'Risk Distribution');
        
        % Add a third sheet with detailed risk factors for each interval
        % Calculate detailed risk factors for each interval
        speed_grad_factor = zeros(num_intervals, 1);
        elasticity_grad_factor = zeros(num_intervals, 1);
        content_accum_factor = zeros(num_intervals, 1);
        position_factor = zeros(num_intervals, 1);
        elasticity_value_factor = zeros(num_intervals, 1);
        
        for i = 1:num_intervals
            segment_indices = x_intestine >= interval_start(i) & x_intestine < interval_end(i);
            
            if any(segment_indices)
                % Calculate average gradients in the interval
                speed_grads = zeros(sum(segment_indices), 1);
                elasticity_grads = zeros(sum(segment_indices), 1);
                elasticity_values = zeros(sum(segment_indices), 1);
                
                idx = find(segment_indices);
                for j = 1:length(idx)
                    if idx(j) < params.num_segments
                        speed_grads(j) = abs(results.speed_profile(post_gastro_idx, idx(j)+1) - results.speed_profile(post_gastro_idx, idx(j)));
                        elasticity_grads(j) = abs(results.elasticity_profile(post_gastro_idx, idx(j)+1) - results.elasticity_profile(post_gastro_idx, idx(j)));
                    end
                    elasticity_values(j) = results.elasticity_profile(post_gastro_idx, idx(j));
                end
                
                speed_grad_factor(i) = mean(speed_grads) * 3.0;
                elasticity_grad_factor(i) = mean(elasticity_grads) * 4.0;
                
                % Calculate content accumulation
                if results.intussusception_detected
                    pre_intus_indices = results.time < results.intussusception_time;
                    post_gastro_indices = results.time >= params.gastroenteritis_onset & results.time < results.intussusception_time;
                    
                    if sum(post_gastro_indices) > 0 && sum(pre_intus_indices) > 0
                        content_before = mean(mean(results.content(pre_intus_indices, segment_indices)));
                        content_after = mean(mean(results.content(post_gastro_indices, segment_indices)));
                        content_accum_factor(i) = max(0, content_after - content_before) * 1.0;
                    end
                else
                    pre_gastro_indices = results.time < params.gastroenteritis_onset;
                    post_gastro_indices = results.time >= params.gastroenteritis_onset;
                    
                    if sum(post_gastro_indices) > 0 && sum(pre_gastro_indices) > 0
                        content_before = mean(mean(results.content(pre_gastro_indices, segment_indices)));
                        content_after = mean(mean(results.content(post_gastro_indices, segment_indices)));
                        content_accum_factor(i) = max(0, content_after - content_before) * 1.0;
                    end
                end
                
                % Calculate anatomical position factor
                mid_point = (interval_start(i) + interval_end(i)) / 2;
                distance_to_terminal = abs(mid_point - params.terminal_ileum_start);
                position_factor(i) = exp(-0.01 * distance_to_terminal^2) * 0.5;
                
                % Calculate elasticity value factor
                mean_elasticity = mean(elasticity_values);
                elasticity_value_factor(i) = max(0, (0.5 - mean_elasticity) * 2) * 2.0;
            end
        end
        
        % Create detailed risk factors table
        DetailedRiskFactors = table(interval_start, interval_end, anatomical_region, ...
            speed_grad_factor, elasticity_grad_factor, content_accum_factor, ...
            position_factor, elasticity_value_factor, interval_risk, ...
            'VariableNames', {'Start_Position_cm', 'End_Position_cm', 'Anatomical_Region', ...
            'Speed_Gradient_Factor', 'Elasticity_Gradient_Factor', 'Content_Accumulation_Factor', ...
            'Anatomical_Position_Factor', 'Elasticity_Value_Factor', 'Total_Risk_Value'});
        
        writetable(DetailedRiskFactors, fullPath, 'Sheet', 'Detailed Risk Factors');
        
        app.updateStatus(['Results data exported to ' fullPath]);
    end
end

        end
    end
    
    % App initialization and construction
    methods (Access = private)
        
        % Create UI components
        function createComponents(app)
            % Create main figure
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100, 100, 900, 650];
            app.UIFigure.Name = 'Intestinal Intussusception Simulation App (Spectrum Model with Elasticity)';
            app.UIFigure.Resize = 'off';
            
            % Create simulation panel
            app.SimulationPanel = uipanel(app.UIFigure);
            app.SimulationPanel.Title = 'Simulation Control';
            app.SimulationPanel.Position = [10, 560, 880, 80];
            
            % Create run button
            app.RunButton = uibutton(app.SimulationPanel, 'push');
            app.RunButton.Position = [20, 40, 100, 22];
            app.RunButton.Text = 'Run Simulation';
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            
            % Create reset button
            app.ResetButton = uibutton(app.SimulationPanel, 'push');
            app.ResetButton.Position = [130, 40, 100, 22];
            app.ResetButton.Text = 'Reset Parameters';
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            
            % Create status label
            app.StatusLabel = uilabel(app.SimulationPanel);
            app.StatusLabel.Position = [20, 10, 840, 22];
            app.StatusLabel.Text = 'Ready to run simulation. Set parameters and click Run.';
            
            % Create save results button
            app.SaveResultsButton = uibutton(app.SimulationPanel, 'push');
            app.SaveResultsButton.Position = [240, 40, 100, 22];
            app.SaveResultsButton.Text = 'Save Results';
            app.SaveResultsButton.ButtonPushedFcn = createCallbackFcn(app, @SaveResultsButtonPushed, true);
            app.SaveResultsButton.Enable = 'off';
            
            % Create table results button
            app.TableResultsButton = uibutton(app.SimulationPanel, 'push');
            app.TableResultsButton.Position = [350, 40, 130, 22];
            app.TableResultsButton.Text = 'Show Results Table';
            app.TableResultsButton.ButtonPushedFcn = createCallbackFcn(app, @TableResultsButtonPushed, true);
            app.TableResultsButton.Enable = 'off';
            
            % Create tab group
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [10, 10, 880, 540];
            
            % Create basic parameters tab
            app.BasicParamsTab = uitab(app.TabGroup);
            app.BasicParamsTab.Title = 'Basic Parameters';
            
            % Create simulation time field
            app.SimTimeEditFieldLabel = uilabel(app.BasicParamsTab);
            app.SimTimeEditFieldLabel.Position = [30, 460, 200, 22];
            app.SimTimeEditFieldLabel.Text = 'Simulation Time (seconds):';
            
            app.SimTimeEditField = uieditfield(app.BasicParamsTab, 'numeric');
            app.SimTimeEditField.Position = [240, 460, 100, 22];
            app.SimTimeEditField.Value = 600;
            
            % Create time step field
            app.TimeStepEditFieldLabel = uilabel(app.BasicParamsTab);
            app.TimeStepEditFieldLabel.Position = [30, 420, 200, 22];
            app.TimeStepEditFieldLabel.Text = 'Time Step (seconds):';
            
            app.TimeStepEditField = uieditfield(app.BasicParamsTab, 'numeric');
            app.TimeStepEditField.Position = [240, 420, 100, 22];
            app.TimeStepEditField.Value = 0.05;
            
            % Create gastroenteritis onset field
            app.GastroOnsetEditFieldLabel = uilabel(app.BasicParamsTab);
            app.GastroOnsetEditFieldLabel.Position = [30, 380, 200, 22];
            app.GastroOnsetEditFieldLabel.Text = 'Gastroenteritis Onset (seconds):';
            
            app.GastroOnsetEditField = uieditfield(app.BasicParamsTab, 'numeric');
            app.GastroOnsetEditField.Position = [240, 380, 100, 22];
            app.GastroOnsetEditField.Value = 30;
            
            % Create gastroenteritis parameters tab
            app.GastroParamsTab = uitab(app.TabGroup);
            app.GastroParamsTab.Title = 'Gastroenteritis Parameters';
            
            % Create proximal speed field
            app.ProximalSpeedEditFieldLabel = uilabel(app.GastroParamsTab);
            app.ProximalSpeedEditFieldLabel.Position = [30, 460, 280, 22];
            app.ProximalSpeedEditFieldLabel.Text = 'Proximal Intestine Speed during Gastroenteritis (cm/s):';
            
            app.ProximalSpeedEditField = uieditfield(app.GastroParamsTab, 'numeric');
            app.ProximalSpeedEditField.Position = [320, 460, 100, 22];
            app.ProximalSpeedEditField.Value = 4.0;
            
            % Create terminal speed field
            app.TerminalSpeedEditFieldLabel = uilabel(app.GastroParamsTab);
            app.TerminalSpeedEditFieldLabel.Position = [30, 420, 280, 22];
            app.TerminalSpeedEditFieldLabel.Text = 'Terminal Ileum Speed during Gastroenteritis (cm/s):';
            
            app.TerminalSpeedEditField = uieditfield(app.GastroParamsTab, 'numeric');
            app.TerminalSpeedEditField.Position = [320, 420, 100, 22];
            app.TerminalSpeedEditField.Value = 1.2;
            
            % Create wave amplitude field
            app.WaveAmplitudeEditFieldLabel = uilabel(app.GastroParamsTab);
            app.WaveAmplitudeEditFieldLabel.Position = [30, 380, 280, 22];
            app.WaveAmplitudeEditFieldLabel.Text = 'Peristaltic Wave Amplitude during Gastroenteritis:';
            
            app.WaveAmplitudeEditField = uieditfield(app.GastroParamsTab, 'numeric');
            app.WaveAmplitudeEditField.Position = [320, 380, 100, 22];
            app.WaveAmplitudeEditField.Value = 0.6;
            
            % Create wave frequency field
            app.WaveFrequencyEditFieldLabel = uilabel(app.GastroParamsTab);
            app.WaveFrequencyEditFieldLabel.Position = [30, 340, 280, 22];
            app.WaveFrequencyEditFieldLabel.Text = 'Peristaltic Wave Frequency during Gastroenteritis (Hz):';
            
            app.WaveFrequencyEditField = uieditfield(app.GastroParamsTab, 'numeric');
            app.WaveFrequencyEditField.Position = [320, 340, 100, 22];
            app.WaveFrequencyEditField.Value = 0.15;
            
            % Create bacterial factor field
            app.BacterialFactorEditFieldLabel = uilabel(app.GastroParamsTab);
            app.BacterialFactorEditFieldLabel.Position = [30, 300, 280, 22];
            app.BacterialFactorEditFieldLabel.Text = 'Bacterial Inhibition Factor (0-1):';
            
            app.BacterialFactorEditField = uieditfield(app.GastroParamsTab, 'numeric');
            app.BacterialFactorEditField.Position = [320, 300, 100, 22];
            app.BacterialFactorEditField.Value = 0.8;
            
            % Create valve resistance field
            app.ValveResistanceEditFieldLabel = uilabel(app.GastroParamsTab);
            app.ValveResistanceEditFieldLabel.Position = [30, 260, 280, 22];
            app.ValveResistanceEditFieldLabel.Text = 'Ileocecal Valve Resistance Factor:';
            
            app.ValveResistanceEditField = uieditfield(app.GastroParamsTab, 'numeric');
            app.ValveResistanceEditField.Position = [320, 260, 100, 22];
            app.ValveResistanceEditField.Value = 1.5;
            
            % Create motility spectrum parameters tab
            app.MotilityParamsTab = uitab(app.TabGroup);
            app.MotilityParamsTab.Title = 'Motility Spectrum';
            
            % Create normal proximal speed field
            app.NormalProximalSpeedEditFieldLabel = uilabel(app.MotilityParamsTab);
            app.NormalProximalSpeedEditFieldLabel.Position = [30, 260, 200, 22];
            app.NormalProximalSpeedEditFieldLabel.Text = 'Normal Proximal Speed (cm/s):';

            app.NormalProximalSpeedEditField = uieditfield(app.MotilityParamsTab, 'numeric');
            app.NormalProximalSpeedEditField.Position = [240, 260, 100, 22];
            app.NormalProximalSpeedEditField.Value = 2.5;  % Predict value 

            % Create normal terminal speed field
            app.NormalTerminalSpeedEditFieldLabel = uilabel(app.MotilityParamsTab);
            app.NormalTerminalSpeedEditFieldLabel.Position = [30, 220, 200, 22];
            app.NormalTerminalSpeedEditFieldLabel.Text = 'Normal Terminal Speed (cm/s):';

            app.NormalTerminalSpeedEditField = uieditfield(app.MotilityParamsTab, 'numeric');
            app.NormalTerminalSpeedEditField.Position = [240, 220, 100, 22];
            app.NormalTerminalSpeedEditField.Value = 2.0;  % Predict value 

            % Create speed pattern dropdown
            app.SpeedPatternLabel = uilabel(app.MotilityParamsTab);
            app.SpeedPatternLabel.Position = [30, 460, 200, 22];
            app.SpeedPatternLabel.Text = 'Speed Distribution Pattern:';
            
            app.SpeedPatternDropDown = uidropdown(app.MotilityParamsTab);
            app.SpeedPatternDropDown.Items = {'Gradual Decrease', 'Segmented Regions', 'Physiological', 'Custom'};
            app.SpeedPatternDropDown.Position = [240, 460, 180, 22];
            app.SpeedPatternDropDown.Value = 'Physiological';
            
            % Create speed transition type dropdown (for Custom pattern)
            app.SpeedTransitionTypeLabel = uilabel(app.MotilityParamsTab);
            app.SpeedTransitionTypeLabel.Position = [30, 420, 200, 22];
            app.SpeedTransitionTypeLabel.Text = 'Speed Transition Type (Custom):';
            
            app.SpeedTransitionTypeDropDown = uidropdown(app.MotilityParamsTab);
            app.SpeedTransitionTypeDropDown.Items = {'Linear', 'Sigmoid', 'Exponential', 'Step with Transition'};
            app.SpeedTransitionTypeDropDown.Position = [240, 420, 180, 22];
            app.SpeedTransitionTypeDropDown.Value = 'Sigmoid';
            
            % Create speed gradient field
            app.SpeedGradientEditFieldLabel = uilabel(app.MotilityParamsTab);
            app.SpeedGradientEditFieldLabel.Position = [30, 380, 200, 22];
            app.SpeedGradientEditFieldLabel.Text = 'Speed Gradient Factor:';
            
            app.SpeedGradientEditField = uieditfield(app.MotilityParamsTab, 'numeric');
            app.SpeedGradientEditField.Position = [240, 380, 100, 22];
            app.SpeedGradientEditField.Value = 1.0;
            
            % Create speed variability field
            app.SpeedVariabilityEditFieldLabel = uilabel(app.MotilityParamsTab);
            app.SpeedVariabilityEditFieldLabel.Position = [30, 340, 200, 22];
            app.SpeedVariabilityEditFieldLabel.Text = 'Local Speed Variability (0-0.5):';
            
            app.SpeedVariabilityEditField = uieditfield(app.MotilityParamsTab, 'numeric');
            app.SpeedVariabilityEditField.Position = [240, 340, 100, 22];
            app.SpeedVariabilityEditField.Value = 0.1;
            
            % Create noise level field
            app.NoiseLevelEditFieldLabel = uilabel(app.MotilityParamsTab);
            app.NoiseLevelEditFieldLabel.Position = [30, 300, 200, 22];
            app.NoiseLevelEditFieldLabel.Text = 'Random Noise Level (0-0.2):';
            
            app.NoiseLevelEditField = uieditfield(app.MotilityParamsTab, 'numeric');
            app.NoiseLevelEditField.Position = [240, 300, 100, 22];
            app.NoiseLevelEditField.Value = 0.05;
            
            % Create elasticity parameters tab (NEW)
            app.ElasticityParamsTab = uitab(app.TabGroup);
            app.ElasticityParamsTab.Title = 'Elasticity Parameters';
            
            % Create normal proximal elasticity field
            app.NormalProximalElasticityEditFieldLabel = uilabel(app.ElasticityParamsTab);
            app.NormalProximalElasticityEditFieldLabel.Position = [30, 460, 280, 22];
            app.NormalProximalElasticityEditFieldLabel.Text = 'Normal Proximal Intestine Elasticity:';
            
            app.NormalProximalElasticityEditField = uieditfield(app.ElasticityParamsTab, 'numeric');
            app.NormalProximalElasticityEditField.Position = [320, 460, 100, 22];
            app.NormalProximalElasticityEditField.Value = 0.4;
            
            % Create normal distal elasticity field
            app.NormalDistalElasticityEditFieldLabel = uilabel(app.ElasticityParamsTab);
            app.NormalDistalElasticityEditFieldLabel.Position = [30, 420, 280, 22];
            app.NormalDistalElasticityEditFieldLabel.Text = 'Normal Distal (Terminal Ileum) Elasticity:';
            
            app.NormalDistalElasticityEditField = uieditfield(app.ElasticityParamsTab, 'numeric');
            app.NormalDistalElasticityEditField.Position = [320, 420, 100, 22];
            app.NormalDistalElasticityEditField.Value = 0.3;
            
            % Create gastroenteritis proximal elasticity field
            app.GastroProximalElasticityEditFieldLabel = uilabel(app.ElasticityParamsTab);
            app.GastroProximalElasticityEditFieldLabel.Position = [30, 380, 280, 22];
            app.GastroProximalElasticityEditFieldLabel.Text = 'Gastroenteritis Proximal Intestine Elasticity:';
            
            app.GastroProximalElasticityEditField = uieditfield(app.ElasticityParamsTab, 'numeric');
            app.GastroProximalElasticityEditField.Position = [320, 380, 100, 22];
            app.GastroProximalElasticityEditField.Value = 0.25;
            
            % Create gastroenteritis distal elasticity field
            app.GastroDistalElasticityEditFieldLabel = uilabel(app.ElasticityParamsTab);
            app.GastroDistalElasticityEditFieldLabel.Position = [30, 340, 280, 22];
            app.GastroDistalElasticityEditFieldLabel.Text = 'Gastroenteritis Distal (Terminal) Elasticity:';
            
            app.GastroDistalElasticityEditField = uieditfield(app.ElasticityParamsTab, 'numeric');
            app.GastroDistalElasticityEditField.Position = [320, 340, 100, 22];
            app.GastroDistalElasticityEditField.Value = 0.15;
            
            % Create elasticity pattern dropdown
            app.ElasticityPatternLabel = uilabel(app.ElasticityParamsTab);
            app.ElasticityPatternLabel.Position = [30, 300, 280, 22];
            app.ElasticityPatternLabel.Text = 'Elasticity Distribution Pattern:';
            
            app.ElasticityPatternDropDown = uidropdown(app.ElasticityParamsTab);
            app.ElasticityPatternDropDown.Items = {'Uniform', 'Linear Gradient', 'Sigmoid Transition', 'Anatomical', 'Terminal Ileum Focus'};
            app.ElasticityPatternDropDown.Position = [320, 300, 180, 22];
            app.ElasticityPatternDropDown.Value = 'Terminal Ileum Focus';
            
            % Create elasticity gradient field
            app.ElasticityGradientEditFieldLabel = uilabel(app.ElasticityParamsTab);
            app.ElasticityGradientEditFieldLabel.Position = [30, 260, 280, 22];
            app.ElasticityGradientEditFieldLabel.Text = 'Elasticity Gradient Factor:';
            
            app.ElasticityGradientEditField = uieditfield(app.ElasticityParamsTab, 'numeric');
            app.ElasticityGradientEditField.Position = [320, 260, 100, 22];
            app.ElasticityGradientEditField.Value = 1.2;
            
            % Create elasticity variability field
            app.ElasticityVariabilityEditFieldLabel = uilabel(app.ElasticityParamsTab);
            app.ElasticityVariabilityEditFieldLabel.Position = [30, 220, 280, 22];
            app.ElasticityVariabilityEditFieldLabel.Text = 'Elasticity Local Variability (0-0.3):';
            
            app.ElasticityVariabilityEditField = uieditfield(app.ElasticityParamsTab, 'numeric');
            app.ElasticityVariabilityEditField.Position = [320, 220, 100, 22];
            app.ElasticityVariabilityEditField.Value = 0.15;
            
            % Create anatomy parameters tab
            app.AnatomyParamsTab = uitab(app.TabGroup);
            app.AnatomyParamsTab.Title = 'Anatomy Parameters';
            
            % Create intestine length field
            app.IntestineLengthEditFieldLabel = uilabel(app.AnatomyParamsTab);
            app.IntestineLengthEditFieldLabel.Position = [30, 460, 200, 22];
            app.IntestineLengthEditFieldLabel.Text = 'Intestine Length (cm):';
            
            app.IntestineLengthEditField = uieditfield(app.AnatomyParamsTab, 'numeric');
            app.IntestineLengthEditField.Position = [240, 460, 100, 22];
            app.IntestineLengthEditField.Value = 200;  % Updated to 200 cm
            
            % Create terminal ileum start field
            app.TerminalIleumStartEditFieldLabel = uilabel(app.AnatomyParamsTab);
            app.TerminalIleumStartEditFieldLabel.Position = [30, 420, 200, 22];
            app.TerminalIleumStartEditFieldLabel.Text = 'Terminal Ileum Start Position (cm):';
            
            app.TerminalIleumStartEditField = uieditfield(app.AnatomyParamsTab, 'numeric');
            app.TerminalIleumStartEditField.Position = [240, 420, 100, 22];
            app.TerminalIleumStartEditField.Value = 180;  % Updated to 180 cm
            
            % Create intestine radius field
            app.IntestineRadiusEditFieldLabel = uilabel(app.AnatomyParamsTab);
            app.IntestineRadiusEditFieldLabel.Position = [30, 380, 200, 22];
            app.IntestineRadiusEditFieldLabel.Text = 'Intestine Radius (cm):';
            
            app.IntestineRadiusEditField = uieditfield(app.AnatomyParamsTab, 'numeric');
            app.IntestineRadiusEditField.Position = [240, 380, 100, 22];
            app.IntestineRadiusEditField.Value = 1.0;
            
            % Create critical ratio field
            app.CriticalRatioEditFieldLabel = uilabel(app.AnatomyParamsTab);
            app.CriticalRatioEditFieldLabel.Position = [30, 340, 200, 22];
            app.CriticalRatioEditFieldLabel.Text = 'Critical Diameter Ratio for Intussusception:';
            
            app.CriticalRatioEditField = uieditfield(app.AnatomyParamsTab, 'numeric');
            app.CriticalRatioEditField.Position = [240, 340, 100, 22];
            app.CriticalRatioEditField.Value = 1.8;
            
            % Create results tab
            app.ResultsTab = uitab(app.TabGroup);
            app.ResultsTab.Title = 'Results';
            
            % Create results text area
            app.ResultsTextArea = uitextarea(app.ResultsTab);
            app.ResultsTextArea.Position = [10, 10, 350, 180];
            app.ResultsTextArea.Value = {'Ready to run simulation.'};
            
            % Create result plot area
            app.ResultPlotArea = uiaxes(app.ResultsTab);
            app.ResultPlotArea.Position = [370, 10, 500, 180];
            app.ResultPlotArea.XLim = [0 200];  % Updated for 200 cm intestine
            app.ResultPlotArea.YLim = [-5 5];
            title(app.ResultPlotArea, 'Intestinal Shape');
            xlabel(app.ResultPlotArea, 'Intestinal Position (cm)');
            ylabel(app.ResultPlotArea, 'Wall Position');
            
            % Create speed profile plot area
            app.SpeedProfilePlotArea = uiaxes(app.ResultsTab);
            app.SpeedProfilePlotArea.Position = [10, 200, 430, 320];
            app.SpeedProfilePlotArea.XLim = [0 200];  % Updated for 200 cm intestine
            app.SpeedProfilePlotArea.YLim = [0 5];
            title(app.SpeedProfilePlotArea, 'Intestinal Speed Profile');
            xlabel(app.SpeedProfilePlotArea, 'Intestinal Position (cm)');
            ylabel(app.SpeedProfilePlotArea, 'Peristaltic Speed (cm/s)');
            
            % Create elasticity profile plot area
            app.ElasticityProfilePlotArea = uiaxes(app.ResultsTab);
            app.ElasticityProfilePlotArea.Position = [450, 200, 420, 320];
            app.ElasticityProfilePlotArea.XLim = [0 200];  % Updated for 200 cm intestine
            app.ElasticityProfilePlotArea.YLim = [0 0.5];
            title(app.ElasticityProfilePlotArea, 'Intestinal Elasticity Profile');
            xlabel(app.ElasticityProfilePlotArea, 'Intestinal Position (cm)');
            ylabel(app.ElasticityProfilePlotArea, 'Elasticity');
            
            % Make the UI visible
            app.UIFigure.Visible = 'on';
        end
    end
    
    % App creation and deletion
    methods (Access = public)
        
        % Construct app
        function app = IntestinalIntussusceptionApp_e
            % Create components
            createComponents(app)
            
            % Initialize app
            app.ResetButtonPushed();
            
            % Register the app with App Designer
            registerApp(app, app.UIFigure)
        end
        
        % Code that executes before app deletion
        function delete(app)
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
    
    % Callbacks that handle component events
    methods (Access = private)
        
        % Button pushed function: TableResultsButton
        function TableResultsButtonPushed(app, event)
            app.exportResultsTable();
        end
        
        % Button pushed function: SaveResultsButton
        function SaveResultsButtonPushed(app, event)
            app.saveResults();
        end
        
        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
        % Disable UI during simulation
         app.RunButton.Enable = 'off';
         app.SaveResultsButton.Enable = 'off';
         app.TableResultsButton.Enable = 'off';
         app.StatusLabel.Text = 'Starting simulation...';
         drawnow;
    
        % Get parameters from UI
        app.params = app.getParametersFromUI();
        app.updateStatus('Parameters loaded. Initializing model...');
        app.updateStatus(sprintf('With the speed pattern: %s', app.params.speed_pattern));
        if strcmp(app.params.speed_pattern, 'Custom')
            app.updateStatus(sprintf('Custom mode type: %s', app.params.speed_transition_type));
        end
    
        % Initialize model (execute only once)
        intestine = app.initialize_model(app.params);
    
        % Validate structure and continue
        intestine = app.validateIntestineStructure(intestine, app.params);
        app.updateStatus('Model initialized. Running simulation...');
    
        % Run simulation (execute only once)
        app.results = app.run_simulation(intestine, app.params);
    
        % Update UI with results
       if app.results.intussusception_detected
           app.updateStatus(sprintf('Simulation complete. Intussusception detected at %.1f seconds at %.1f cm!', ...
                                app.results.intussusception_time, app.results.intussusception_location));
      else
          app.updateStatus('Simulation complete. No intussusception detected.');
     end
    
        % Plot results
        app.updateStatus('Generating visualizations...');
        app.plotResults();
    
        % Generate results table
        app.updateStatus('Creating results summary table...');
        app.exportResultsTable();
    
        % Enable UI
        app.RunButton.Enable = 'on';
        app.SaveResultsButton.Enable = 'on';
        app.TableResultsButton.Enable = 'on';
        app.StatusLabel.Text = 'Simulation complete.';
    end
     
        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            % Reset parameters to default values
            app.SimTimeEditField.Value = 600;
            app.TimeStepEditField.Value = 0.1;
            app.GastroOnsetEditField.Value = 30;
            
            app.ProximalSpeedEditField.Value = 3.5;
            app.TerminalSpeedEditField.Value = 1.5;
            app.WaveAmplitudeEditField.Value = 0.6;
            app.WaveFrequencyEditField.Value = 0.15;
            app.BacterialFactorEditField.Value = 0.8;
            app.ValveResistanceEditField.Value = 1.5;
            
            app.IntestineLengthEditField.Value = 300; 
            app.TerminalIleumStartEditField.Value = 270;  % 90% of length
            app.IntestineRadiusEditField.Value = 1.0;
            app.CriticalRatioEditField.Value = 1.8;
            
            % Reset normal parameters
            app.NormalProximalSpeedEditField.Value = 2.5;
            app.NormalTerminalSpeedEditField.Value = 2.0;
            
            % Reset motility spectrum parameters
            app.SpeedPatternDropDown.Value = 'Custom';
            app.SpeedTransitionTypeDropDown.Value = 'Step with Transition';
            app.SpeedGradientEditField.Value = 1.0;
            app.SpeedVariabilityEditField.Value = 0.1;
            app.NoiseLevelEditField.Value = 0.05;
            
            % Reset elasticity parameters (NEW)
            app.NormalProximalElasticityEditField.Value = 0.4;
            app.NormalDistalElasticityEditField.Value = 0.35;
            app.GastroProximalElasticityEditField.Value = 0.35;
            app.GastroDistalElasticityEditField.Value = 0.15;
            app.ElasticityPatternDropDown.Value = 'Terminal Ileum Focus';
            app.ElasticityGradientEditField.Value = 1.5;
            app.ElasticityVariabilityEditField.Value = 0.15;
            
            % Clear results
            app.ResultsTextArea.Value = {'Parameters reset to default values.'};
            cla(app.ResultPlotArea);
            cla(app.SpeedProfilePlotArea);
            cla(app.ElasticityProfilePlotArea);  % Clear elasticity plot
            title(app.ResultPlotArea, 'Intestinal Shape');
            title(app.SpeedProfilePlotArea, 'Intestinal Speed Profile');
            title(app.ElasticityProfilePlotArea, 'Intestinal Elasticity Profile');
            
            app.SaveResultsButton.Enable = 'off';
            app.TableResultsButton.Enable = 'off';
            app.StatusLabel.Text = 'Parameters reset to default values.';
        end
    end
end

%% Create and run the app
% app = IntestinalIntussusceptionApp_e;
