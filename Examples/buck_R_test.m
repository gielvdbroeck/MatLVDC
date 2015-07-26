% BUCK_R_TEST.M
%
% Simulation example for MatLVDC v2.0 with a unipolar network that consists of one line that connects a
% voltage source and a CO component (Buck_R, i.e. a buck converter that regulates the output voltage
% across a resistor at a constant value)
%
% Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last modification: July 23, 2015

clear all;

%% Define the network topology
% The first step comprises setting up the network topology and the branch parameters
% Consult the user guide for more information about the different column headers

% Bus data
                %BUSDC_I        %BUSDC_CONFIG
bus_data = [    1               1;
                2               1];
                
% Branch data
                %F_BUS      %T_BUS      %BR_CONFIG      %BR_R1      %BR_R2      %BR_Rn      %BR_L1      %BR_L2      %BR_Ln      %BR_C1          %BR_C2
branch_data = [ 1           2           1               0.1         0           0.1         0.5e-5      0           0.5e-5      0               0; % LINE 1-2
              ];

%% Create a DC System
% A DC System object is created, which enables the analysis features of MatLVDC. Upon construction, the
% network topology should be included in the function call.

% Define the nominal network voltage between the phase conductor and the neutral conductor
nominalVoltage = 380;

% Create the DC System object
dc_system = DCSystem(bus_data, branch_data, nominalVoltage);
disp('DC system initialized');

%% Add components to the network
% First create the components and add them subsequently to the DC system at the specified bus number. The
% component connection (po, on, pon) determines how the component will be connected to the network.

C_bus = 1e-5;
component{1} = Voltage_source('V source BUS 1',Connection.po);
%component{2} = Current_source('I source BUS 1',Connection.po);
component{2} = Buck_R('BUS 2 CPL', Connection.po, [30 470e-6 10e-3 1.4 7.9 1.9 10 15]);
component{3} = Capacitor('CAP BUS 2',Connection.po,C_bus);

dc_system.connectComponent(component,[1 2 2]);

disp('Components added');

%% Show the network's external outputs, inputs and state-variables
% The external outputs need to be set, as well as the initial state for simulation. The index column
% corresponds to the index that has been added to the quantity in MatLVDC
labels = dc_system.getLabels();

%% Run a steady-state simulation
% The steadyState function of the DCSystem is used to determine the steady state

% Define the external inputs to the DC system
U_ext = [       380     300      ];

% Run the analysis function
X_bfs = dc_system.steadyState_bfs(U_ext,1e-6);
X = dc_system.steadyState(U_ext,X_bfs,[]);

%% Linearized system
[A, sysPoles] = dc_system.linearSys(X,U_ext);

%% Run dynamic simulation

% Initial state
x0 = X_bfs;

% External inputs
u_ext = SignalBuilder();
u_ext.newSignal(component{1}.u_ext{1}.index,@constant,380); % voltage source value
%u_ext.newSignal(component{2}.u_ext{1}.index,@constant,300);
u_ext.newSignal(component{2}.u_ext{1}.index,@ramp_signal,[0.05 0.055 300 360]);

disp('Starting simulation...');
options = odeset('AbsTol',1e-3,'RelTol',1e-4,'Stats','on');
resultSet = dc_system.simulate(@ode15s, options, [0.04 0.1], u_ext, x0);
disp('Simulation finished');

%% Plot results
figure;
title('Linearized system poles');
scatter(real(sysPoles),imag(sysPoles),'x');
grid on;

figure;
subplot(2,1,1);
p = plot(resultSet.t,resultSet.x);
set(p,{'DisplayName'},labels.x);
xlabel('Time [s]');
legend('-DynamicLegend');
grid on;

subplot(2,1,2);
p = plot(resultSet.t,resultSet.y_ext);
set(p,{'DisplayName'},labels.y_ext);
xlabel('Time [s]');
legend('-DynamicLegend');
grid on;