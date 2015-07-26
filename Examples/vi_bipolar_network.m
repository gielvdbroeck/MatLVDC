% VI_BIPOLAR_NETWORK.M
%
% Simulation example for MatLVDC v2.0 with a bipolar network that consists of one line that connects a
% voltage and current source.
%
% Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last modification: July 16, 2015

clear all

%% Define the network topology
% The first step comprises setting up the network topology and the branch parameters
% Consult the user guide for more information about the different column headers

% Bus data
                %BUSDC_I        %BUSDC_CONFIG
bus_data = [    1               3;
                2               3];
                
% Branch data
                %F_BUS      %T_BUS      %BR_CONFIG      %BR_R1      %BR_R2      %BR_Rn      %BR_L1      %BR_L2      %BR_Ln      %BR_C1          %BR_C2
branch_data = [ 1           2           3               0.1         0.1         0.1         0.5e-3      0.5e-3      0.5e-3      0               0; % LINE 1-2
              ];

%% Create a DC System
% A DC System object is created, which enables the analysis features of MatLVDC. Upon construction, the
% network topology should be included in the function call.

% Define the nominal network voltage between the phase conductor and the neutral conductor
nominalVoltage = 380;

% Create the DC System object
dc_system_vi = DCSystem(bus_data, branch_data, nominalVoltage);

%% Add components to the network
% First create the components and add them subsequently to the DC system at the specified bus number. The
% component connection (po, on, pon) determines how the component will be connected to the network.

component{1} = Voltage_source('V source BUS 1+',Connection.po);
component{2} = Voltage_source('V source BUS 1-',Connection.on);
component{3} = Current_source('I source BUS 2+',Connection.po);
component{4} = Voltage_source('V source BUS 2+',Connection.po);
component{5} = Voltage_source('V source BUS 2-',Connection.on);

dc_system_vi.connectComponent(component,[1 1 2 2 2]);

%% Run dynamic simulation

% Initial state
x0 = [0 0];

% External inputs
u_ext = SignalBuilder();
u_ext.newSignal(component{1}.u_ext{1}.index,@constant,380); % voltage source value
u_ext.newSignal(component{2}.u_ext{1}.index,@constant,380); % voltage source value
u_ext.newSignal(component{3}.u_ext{1}.index,@constant,10); % current source value
u_ext.newSignal(component{4}.u_ext{1}.index,@step_signal,[0.1 380 360]); % voltage source value
u_ext.newSignal(component{5}.u_ext{1}.index,@constant,380); % voltage source value

resultSet = dc_system_vi.simulate(@ode23tb, [], [0 1], u_ext, x0);

%% Plot results
p = plot(resultSet.t,resultSet.x);
ylabel('Current [A]');
xlabel('Time [s]');
legend('LINE 1-2 +','LINE 1-2 -');
grid on;
