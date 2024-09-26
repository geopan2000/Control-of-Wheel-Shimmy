% Exercise 4.6

    mdl = 'coursework_model_template';
% Initialize parameters:
    sigma_rel = 0.2;        
    k_values = 5000:5000:100000;% values for loop (counter)
    z_initial=0; %e linearisation should always be performed at δ=0; 
    eigenvaluesEx6 = [];
    k_for_eigenvalues = [];
    speed=10; %constant speed

% Range the stiffness value k:
    for k = k_values
        [A, B, C, D] = linmod(mdl);% Linearize the model
        eigs_current = eig(A);% Compute eigenvalues of A
        eigenvaluesEx6 = [eigenvaluesEx6; eigs_current];% Store the eigenvalues in diff column 
        % of the matrix eigenvaluesEX6 for every interation.
        % Store k for each eigenvalue
        k_for_eigenvalues = [k_for_eigenvalues; k * ones(length(eigs_current), 1)];
    end

% Plot eigenvalues on complex plane
    figure;
    scatter(real(eigenvaluesEx6), imag(eigenvaluesEx6), 36, k_for_eigenvalues, 'filled');
    cb = colorbar; % colorbar to identify the value of k for every eig.
    colormap('hot');
    ylabel(cb, 'Stiffness (N/m)');
    ylim([0 inf]);
    xlabel('Real Part');
    ylabel('Imaginary Part');
    title('Eigenvalues for various stiffness');
    grid on;
    clim([min(k_values) max(k_values)]); % Color axis range for k values
    colormap jet;

%% 
%Exercise 4.7
% Initialize parameters
    sim_time = 4; % 4 seconds of simulation
    set_param(mdl, 'StopTime', num2str(sim_time));
    z_initial = 0.01;
    
% Simulation k = 30000
    k=30000;
    [A, B, C, D] = linmod(mdl);
    sim_out1 = sim(mdl); 
% Extract output data for the 1st simulation
    y_pin1_data = y_pin.signals.values; 
    time1 = y_pin.time;
    z_forks1_data = z_forks.signals.values;

% Simulation with k = 100000
    k=100000;
    [A, B, C, D] = linmod(mdl);
    sim_out2 = sim(mdl); % Run

% Extract output data for the 2nd simulation
    y_pin2_data = y_pin.signals.values; 
    time2 = y_pin.time;
    z_forks2_data = z_forks.signals.values;

% Lateral displacement k = 30000
    figure;
    subplot(2, 1, 1);
    plot(time1, y_pin1_data, 'r-');
    xlabel('Time (s)');
    ylabel('Lateral Displacement y (m)');
    legend('k = 30000 N/m');
    title('Lateral Displacement y vs. Time');

% Rotation k = 30000
    subplot(2, 1, 2);
    plot(time1, z_forks1_data, 'r-');
    xlabel('Time (s)');
    ylabel('z Rotation (rad)');
    legend('k = 30000 N/m');
    title('z Rotation vs. Time');

% Lateral displacement k=100000
    figure;
    subplot(2, 1, 1);
    plot(time2, y_pin2_data, 'b-');
    xlabel('Time (s)');
    ylabel('Lateral Displacement y (m)');
    legend('k = 100000 N/m');
    title('Lateral Displacement y vs. Time');

% Rotation k=100000
    subplot(2, 1, 2);
    plot(time2,z_forks2_data, 'b-');
    xlabel('Time (s)');
    ylabel('z Rotation (rad)');
    legend('k = 100000 N/m');
    title('z Rotation vs. Time');
%%
% Exercise 4.8

% Initialize parameters
    k = 100000;
    z_initial=0;
    speed_values = 1:1:50; %Range speed (loop counter)
    eigenvalues_all = [];
    speeds_for_eigenvalues = []; % To keep track of which speed corresponds to which eigenvalue

% Loop for range the speed values
    for speed = speed_values
        [A, B, C, D] = linmod(mdl);
        eigs_current = eig(A);
        eigenvalues_all = [eigenvalues_all; eigs_current];
        speeds_for_eigenvalues = [speeds_for_eigenvalues; repmat(speed, length(eigs_current), 1)];
    end

% Plot eigenvalues
    figure;
    scatter(real(eigenvalues_all), imag(eigenvalues_all), 30, speeds_for_eigenvalues, 'filled');
    colorbar;
    cb = colorbar; % Create a handle for the colorbar
    ylabel(cb, 'Forward Speed (m/s)'); % Add label to the colorbar
    xlabel('Real Part');
    ylabel('Imaginary Part');
    title('Eigenvalues for Various Speeds');
    xlim([-20 inf]); % Limit for the real part
    ylim([0 inf]);   % Limit for the imaginary part
    grid on;
    clim([min(speed_values) max(speed_values)]); % Color axis range for speed values
    colormap jet; % Set the colormap

 %% Exercise 4.9
    speed = 10;
    k = 100000;
    sigma_values = 0.02:0.02:1;
    eigenvalues_all = [];
    sigma_for_eigenvalues = [];

    % Loop to range sigma values
    for sigma_rel = sigma_values
        [A, B, C, D] = linmod(mdl);
        eigs_current = eig(A);
        eigenvalues_all = [eigenvalues_all; eigs_current];
        sigma_for_eigenvalues = [sigma_for_eigenvalues; sigma_rel * ones(length(eigs_current), 1)];
    end

% Plot eigenvalues on the complex plane
    figure;
    scatter(real(eigenvalues_all), imag(eigenvalues_all), 30, sigma_for_eigenvalues, 'filled');
    xlabel('Real Part');
    ylabel('Imaginary Part');
    title('Eigenvalues for Various relaxation lengths');
    xlim([-20, 2]);
    ylim([0 inf])
    grid on;
    cb = colorbar;
    ylabel(cb, 'Relaxation length (m)');
    clim([min(sigma_values) max(sigma_values)]);
    colormap jet;

% Initialize parameters
    sigma_rel = 1;
    z_initial = 0.01;
    [A, B, C, D] = linmod(mdl);
    sim(mdl); %run

% Extract output
    y_pin3_data = y_pin.signals.values; 
    time3 = y_pin.time;
    z_forks3_data = z_forks.signals.values;

% Plot Lat dis and rotation
    figure;
    subplot(2, 1, 1);
    plot(time3, y_pin3_data);
    title('y King-Pin Displacement vs. Time');
    xlabel('Time (s)');
    ylabel('Displacement (m)');
   
    subplot(2, 1, 2);
    plot(time3, z_forks3_data);
    title('z Fork Rotation vs. Time');
    xlabel('Time (s)');
    ylabel('Rotation (rad)');

%% 
% 4.10
% Intialize parameters
    mdl = 'coursework_model_template';
    sigma_rel = 0.2;
    k=60000;
    speed=10;
    z_initial=0;
    [A, B, C, D] = linmod(mdl);
    eigs_current = eig(A) %print eigvalues

 % simulate with delta=0
    z_initial = 0.1;
    [A, B, C, D] = linmod(mdl);
    sim_time = 4; % 4 seconds of simulation
    set_param(mdl, 'StopTime', num2str(sim_time));
    sim_out3 = sim(mdl); 

% Extract output
    y_pin1_data = y_pin.signals.values; 
    time1 = y_pin.time;
    z_forks1_data = z_forks.signals.values;

% Lateral displacement y
    figure;
    subplot(2, 1, 1);
    plot(time1, y_pin1_data, 'r-');
    xlabel('Time (s)');
    ylabel('Lateral Displacement y (m)');
    legend('k = 60000 N/m');
    title('Lateral Displacement y vs. Time');

% z rotation of forks for both simulations
    subplot(2, 1, 2);
    plot(time1, z_forks1_data, 'r-');
    xlabel('Time (s)');
    ylabel('z Rotation (rad)');
    legend('k = 60000 N/m');
    title('z Rotation vs. Time');
    hold off;

% open the new changed template which has the control feedback loop
    mdl = 'coursework_model_template_changed';
% Initialize parameters 
    sigma_rel = 0.2;
    k=60000;
    speed=10;
    z_initial=0;
    lamda=0;

 % Specify the input and output points for linearization
    in = linio('coursework_model_template_changed/Gain2', 1, 'input'); % path of input
    out = linio('coursework_model_template_changed/PS-Simulink Converter1', 1, 'output'); % path of out
    io = [in, out];
    sys1 = linearize(mdl, io);%linirization with respect to a specific in and out
    Open_loop_P=pole(sys1) %Validation of marginal stability
    Open_loop_Z=zero(sys1) %Validation of marginal stability
    [GM,PM,~,~] = margin(sys1); 

    f1=figure;
    nyquist(sys1);
    xlim([-0.01 0.001]);%focus 
    hold off;

    lamda=-1*GM %Gm for marginal stability and 5*GM quick stabilization
    in2 = linio('coursework_model_template_changed/Gain2', 1, 'input');
    out2 = linio('coursework_model_template_changed/PS-Simulink Converter1', 1, 'output');
    io2 = [in2, out2];
    sys2= linearize(mdl, io2);
    Close_loop_P=pole(sys2)%print poles to validate stability
    Close_loop_Z=zero(sys2)%print zeroes to validate stability
    [GM,PM,~,~] = margin(sys2) 
    f2=figure;
    nyquist(sys2);
    hold off;

% Plots
    z_initial = 0.1;
    sim_time = 4; % 4 seconds of simulation
    set_param(mdl, 'StopTime', num2str(sim_time));
    sim_out4 = sim(mdl); 

% Extract output data
    y_pin4_data = y_pin.signals.values; 
    time4 = y_pin.time;
    z_forks4_data = z_forks.signals.values;

% Lateral displacement
    figure;
    subplot(2, 1, 1);
    plot(time4, y_pin4_data, 'r-');
    xlabel('Time (s)');
    ylabel('Lateral Displacement y (m)');
    legend('k = 60000 N/m');
    title('Lateral Displacement y vs. Time');

% Rotation
    subplot(2, 1, 2);
    plot(time4, z_forks4_data, 'r-');
    xlabel('Time (s)');
    ylabel('z Rotation (rad)');
    legend('k = 60000 N/m');
    title('z Rotation vs. Time');
    hold off;
