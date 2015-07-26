%% Constant power DC network
% The example demonstrates how to run a steady-state simulation of a constant-power DC network. The
% example network consists of 6 unipolar DC busses (plus and neutral conductors) which are
% interconnected by 5 branches. The network is depicted below (busses are indicated by numbers between
% brackets; branches are indicated by -branch no.-). At BUS 1 voltage sources are connected, at all other
% busses constant power loads (or generation if the power is set to negative values).
%
%                         (4)
%                          |
%                          4
%                          |
% (1) -1- (2) -2- (3) -3- (5)
%                          |
%                          5
%                          |
%                         (6)
%
% Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be)
% - Last modification: July 19, 2015

%% Define the network topology
% The first step comprises setting up the network topology and the branch parameters
% Consult the user guide for more information about the different column headers

% Bus data
                %BUSDC_I        %BUSDC_CONFIG
bus_data = [    1               1;
                2               1;
                3               1;
                4               1;
                5               1;
                6               1];

% Branch data
                %F_BUS      %T_BUS      %BR_CONFIG      %BR_R1      %BR_R2      %BR_Rn      %BR_L1      %BR_L2      %BR_Ln      %BR_C1          %BR_C2
branch_data = [ 1           2           1               0.01        0.01        0.01        0.5e-3      0.5e-3      0.5e-3      0               0; % branch 1
                2           3           1               0.01        0.01        0.01        0.5e-3      0.5e-3      0.5e-3      0               0; % branch 2
                3           5           1               0.01        0.01        0.01        0.5e-3      0.5e-3      0.5e-3      0               0; % branch 3
                5           4           1               0.01        0.01        0.01        0.5e-3      0.5e-3      0.5e-3      0               0; % branch 4
                5           6           1               0.01        0.01        0.01        0.5e-3      0.5e-3      0.5e-3      0               0; % branch 5
              ];

%% Create a DC System
% A DC System object is created, which enables the analysis features of MatLVDC. Upon construction, the
% network topology should be included in the function call.

% Define the nominal network voltage between the phase conductor and the neutral conductor
nominalVoltage = 380;

% Create the DC System object
dc_system_cpl = DCSystem(bus_data, branch_data, nominalVoltage);

%% Add components to the network
% First create the components and add them subsequently to the DC system at the specified bus number. The
% component configuration (1=unipolar positive, 2=unipolar negative, 3=bipolar) determines how the
% component will be connected to the network.

% Initialize cell array that will contain the Component objects
comp = cell(1,6);
busses = zeros(6,1);

% Constant voltage source at BUS 1
% Note: BUS 1 is always by definition the slack bus and needs a voltage controlling component connected to it
% Note: the Name class allows you to use human readable strings instead of the numbering convention used in the
% software
comp{1} = Voltage_source('BUS 1 VOLTAGE SOURCE', Connection.po);
% Indicate that component 1 should be connected to BUS 1
busses(1) = 1;

% Constant power load capacitance
C = 1e-3;

% Constant power loads/generators at BUS 2..6
for k=2:6
    comp{k} = CPL_C(['BUS ' num2str(k) ' CPL_C'], Connection.po, C);
    % Indicate that component k should be connected to BUS k
    busses(k) = k;
end

% Add components to the DC system
dc_system_cpl.connectComponent(comp, busses);

%% Display quantities
labels = dc_system_cpl.getLabels();

%% Run a steady-state simulation
% The steadyState function of the DCSystem is used to determine the steady state

% Define the external inputs to the dc_system_cpl (i.e. BUS 1 voltage and BUS 2..6 the power output of the components)
            % BUS 1 [V]     % BUS 2 [W] % BUS 3 [W] % BUS 4 [W] % BUS 5 [W] % BUS 6 [W]
u_ext = [   380             2e3         2e3         2e3         2e3         2e3     ];

% Run the analysis function
X = dc_system_cpl.steadyState_bfs(u_ext,1e-6);
X = dc_system_cpl.steadyState(u_ext,X,[]);

%% Linearized system
[A, sysPoles] = dc_system_cpl.linearSys(X,u_ext);

%% Run a time-domain simulation
%x0 = [0; 0; 0; 0; 0; 380; 380; 380; 380; 380];
x0 = X;
t_span = [0 0.3];

% External inputs: the user can custom define input signals to the external inputs
u_ext = SignalBuilder();
u_ext.newSignal(comp{1}.u_ext{1}.index,@constant,380);
u_ext.newSignal(comp{2}.u_ext{1}.index,@ramp_signal,[0.1 0.15 2e3 -2e3]);
%u_ext.newSignal(comp{2}.u_ext{1}.index,@constant,2e3);
u_ext.newSignal(comp{3}.u_ext{1}.index,@constant,2e3);
u_ext.newSignal(comp{4}.u_ext{1}.index,@constant,2e3);
u_ext.newSignal(comp{5}.u_ext{1}.index,@constant,2e3);
u_ext.newSignal(comp{6}.u_ext{1}.index,@constant,2e3);

resultSet = dc_system_cpl.simulate(@ode15s, [], t_span, u_ext, x0);

%% Plot results
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
p = plot(resultSet.t,resultSet.x);
set(p,{'DisplayName'},labels.x);
xlabel('Time [s]');
%legend(labels.x,'Location','eastoutside');
grid on;

subplot(2,1,2);
p = plot(resultSet.t,resultSet.y_ext);
set(p,{'DisplayName'},labels.y_ext);
xlabel('Time [s]');
%legend('-DynamicLegend');
grid on;