% MatLVDC is an open source Matlab toolbox to analyze DC distribution network dynamics.
% It features time-domain analysis, steady-state analysis and small-signal linearization.
%
% Copyright (C) 2015 Giel Van den Broeck, Tuan Dat Mai, Johan Driesen
% Contact: giel.vandenbroeck@energyville.be
%
% Scientific publications derived from this toolbox should refer to publications mentioned
% in README
%
% This program is free software: you can redistribute it and/or modify it under the terms of the
% GNU General Public License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version. This program is distributed in the hope
% that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program.
% If not, see http://www.gnu.org/licenses/.

classdef DCSystem < handle
    % DCSYSTEM The class DCSYSTEM represents a DC system to be analyzed
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 16, 2015
    
    properties
        % Array containing the busses in the DC system
        busses;
        % Object representing the DC electrical network
        network;
    end
    
    properties(Access=private)
        % Structure containing the labels of the states, internal/external inputs, internal/external outputs
        labels;
        % Counter object that keeps track of the state-variables and I/O
        % indexes
        counter = Counter();
    end
    
    methods(Access=private)
        function [] = registerElement(obj, element)
            % REGISTERELEMENT Register an element to the DC system
            
            %% Obtain indexes
            x_index = []; u_int_index = []; y_int_index = []; u_ext_index = []; y_ext_index = [];
            
            if isa(element,'Component')
                    [x_index, u_int_index, y_int_index, u_ext_index, y_ext_index] = obj.counter.getComponentIndexes(element);
            elseif isa(element,'Network')
                    [x_index, u_int_index, y_int_index, y_ext_index] = obj.counter.getNetworkIndexes(element);
            end
            
            %% Add indexes to the element
            for k=1:length(x_index), element.x{k}.index = x_index(k); end
            for k=1:length(u_int_index), element.u_int{k}.index = u_int_index(k); end
            for k=1:length(y_int_index), element.y_int{k}.index = y_int_index(k); end
            for k=1:length(u_ext_index), element.u_ext{k}.index = u_ext_index(k); end
            for k=1:length(y_ext_index), element.y_ext{k}.index = y_ext_index(k); end
        end
        
        function [] = addBusses(obj, bus_data)
            % ADDBUSSES Add busses to the DC system
            
            for k=1:size(bus_data,1)
                obj.busses{k} = Bus(bus_data(k,:));
                %disp(obj.busses{k}.label);
            end
        end
        
        function interconnectionM = interconnectionMatrices(obj)
            % GETINTERCONNECTIONMATRICES sets up the interconnection matrices H0, H1, H2 and H3
            %
            % u_N = H0 y_v
            % v = H1 y_v
            % i = H2 y_N + H3 y_c
            %
            % (Consult the paper for the details)
            
            for k=1:length(obj.busses)
                %% NETWORK components
                % u_N = H0 y_v
                [u_N_index_p, u_N_index_n] = obj.network.getIndex(obj.busses{k});
                [y_v_index_p, y_v_index_n] = obj.busses{k}.getIndex('y_v');
                % Positive terminal
                interconnectionM.H0(u_N_index_p, y_v_index_p) = 1;
                % Negative terminal
                interconnectionM.H0(u_N_index_n, y_v_index_n) = 1;
                                
                %% CO components
                % v = H1 y_v
                [v_index_p, v_index_n] = obj.busses{k}.getIndex('v');
                interconnectionM.H1(v_index_p,y_v_index_p) = 1;
                interconnectionM.H1(v_index_n,y_v_index_n) = 1;
                
                %% VO components
                % i = H2 y_N + H3 y_c
                [i_index_p, i_index_n] = obj.busses{k}.getIndex('i');
                [y_c_index_p, y_c_index_n] = obj.busses{k}.getIndex('y_c');
                [y_N_index_p, y_N_index_n] = obj.network.getIndex(obj.busses{k});
                % Positive terminal
                interconnectionM.H2(i_index_p,y_N_index_p) = 1;
                interconnectionM.H3(i_index_p,y_c_index_p) = -1;
                % Negative terminal
                interconnectionM.H2(i_index_n,y_N_index_n) = 1;
                interconnectionM.H3(i_index_n,y_c_index_n) = -1;
            end
        end
        
        function [dx] = diffEqns(obj, t, x, u_ext, interconnectionM)
            % DIFFEQNS returns the DC system differential equations which are a function of the current
            % state, time and the external inputs. The interconnection matrices should be provided by the
            % structure interconnectionM (these can be set-up with the interconnectionMatrices method).
            %
            % INPUTS
            % t                 current time value
            % x                 state vector
            % u_ext             external input vector
            % interconnectionM  structure containing the interconnection matrices (can be set-up using
            %                   the interconnectionMatrices method in the DC System class)
            
            %% Segment state vectors
            x_ = x; clear x;
            xv = obj.counter.getNumberOfStates('vo');
            xc = obj.counter.getNumberOfStates('co');
            xn = obj.counter.getNumberOfStates('network');
            x.vo = x_(1:xv);
            x.co = x_(xv+1:xv+xc);
            x.network = x_(xv+xc+1:xv+xc+xn);
            
            %% Calculate subsystem internal outputs
            % Structure y_int (y_int.network, y_int.vo, y_int.co) is calculated by method in DC System
            y_int = obj.calcIntOutputs(t, x, u_ext);
            
            %% Calculate subsystem internal inputs
            % Structure u_int (u_int.network, u_int.vo, u_int.co) is calculated by method in DC System
            % Based upon interconnection matrices
            u_int = obj.calcIntInputs(y_int, interconnectionM);
            
            %% Obtain bus differential equations
            dx = struct('vo',zeros(xv,1),'co',zeros(xc,1),'network',zeros(xn,1));
            for b=1:length(obj.busses)
                bus = obj.busses{b};
                [index, dx_] = bus.diffEqns(t,x,u_int,u_ext);
                dx.vo(index.vo) = dx_.vo;
                dx.co(index.co) = dx_.co;
            end
            
            %% Obtain network differential equations
            [index, dx_] = obj.network.diffEqns(t,x.network,u_int.network,u_ext);
            dx.network(index) = dx_;
            
            %% Join state vectors
            dx = [dx.vo; dx.co; dx.network];
        end
        
        function u_int = calcIntInputs(obj, y_int, interconnectionM)
            % CALCINTINPUTS calculates the internal input values of the DC system
            %
            % INPUTS
            % y_int             structure containing the internal outputs (y_int.network, y_int.vo,
            % y_int.co)
            % interconnectionM  structure containing the interconnection matrices (can be set-up using
            %                   the interconnectionMatrices method in the DC System class)
            %
            % OUTPUTS
            % u_int             structure containing the internal input vectors (u_int.vo, u_int.co,
            %                   u_int.network)
            
            u_int.network = interconnectionM.H0 * transp(y_int.vo);
            u_int.co = interconnectionM.H1 * transp(y_int.vo);
            if ~isempty(y_int.co)
                u_int.vo = interconnectionM.H2 * transp(y_int.network) + interconnectionM.H3 * transp(y_int.co);
            else
                % If no CO components are connected
                u_int.vo = interconnectionM.H2 * transp(y_int.network);
            end
        end
        
        function y_int = calcIntOutputs(obj, t, x, u_ext)
            % CALCINTOUTPUTS calculates the internal output values of the DC system
            %
            % INPUTS
            % t                 current time value
            % x                 state vector
            % u_ext             external input vector
            %
            % OUTPUTS
            % y_int             structure containing the internal outputs (y_int.vo, y_int.co,
            %                   y_int.network)
            
            %% Calculate component internal outputs
            for b=1:length(obj.busses)
                [index, y_int_] = obj.busses{b}.calcIntOutputs(t, x, u_ext);
                y_int.vo(index.vo) = y_int_.vo;
                y_int.co(index.co) = y_int_.co;
            end
            
            %% Calculate network internal outputs
            [index, y_int_] = obj.network.calcIntOutputs(t, x.network);
            y_int.network(index) = y_int_.network;
        end
        
        function y_ext = calcExtOutputs(obj, t, x, u_ext, interconnectionM)
            % CALCEXTOUTPUTS calculates the external output values of the DC system
            %
            % INPUTS
            % t                 current time value
            % x                 state vector
            % u_ext             external input vector
            %
            % OUTPUTS
            % y_ext             structure containing the external outputs (y_ext.vo, y_ext.co,
            %                   y_ext.network)
            
            %% Segment state-vector
            xv = obj.counter.getNumberOfStates('vo');
            xc = obj.counter.getNumberOfStates('co');
            xn = obj.counter.getNumberOfStates('network');
            x_.vo = x(1:xv);
            x_.co = x(xv+1:xv+xc);
            x_.network = x(xv+xc+1:xv+xc+xn);
            
            %% Calculate subsystem internal outputs
            % Structure y_int (y_int.network, y_int.vo, y_int.co) is calculated by method in DC System
            y_int = obj.calcIntOutputs(t, x_, u_ext);
            
            %% Calculate subsystem internal inputs
            % Structure u_int (u_int.network, u_int.vo, u_int.co) is calculated by method in DC System
            % Based upon interconnection matrices
            u_int = obj.calcIntInputs(y_int, interconnectionM);
            
            %% Calculate component external outputs
            y_ext = zeros(obj.counter.getNumberOfExtOutputs,1);
            for b=1:length(obj.busses)
                [index, y_ext_] = obj.busses{b}.calcExtOutputs(t, x_, u_int, u_ext);
                y_ext(index) = y_ext_;
            end
            
            %% Calculate network external outputs
            [index, y_ext_] = obj.network.calcExtOutputs(t, x_.network, u_int);
            y_ext(index) = y_ext_;
        end
        
        function [] = validate(obj)
            % VALIDATE checks the DCSystem prior to solving for most common mistakes
            
            for k=1:length(obj.busses)
                obj.busses{k}.validate;
            end
            
        end
    end
    
    methods
        
        function obj = DCSystem(bus_data, branch_data, nominalVoltage)
            % DCSYSTEM The class DCSYSTEM represents a DC system to be analyzed
            
            disp('--------------------------------------------------');
            disp('MatLVDC  Copyright (C) 2015 KU Leuven, Research Group ELECTA');
            disp('This program comes with ABSOLUTELY NO WARRANTY; for details consult LICENCE file');
            disp('This is free software and you are welcome to redistribute it under GNU GPL licence.')
            disp('--------------------------------------------------');
            disp('Scientific publications derived from this toolbox should refer to publications mentioned in README');
            disp('--------------------------------------------------');
            disp('Developed by Giel Van den Broeck, Tuan Dat Mai, Johan Driesen');
            disp('--------------------------------------------------');
            
            %% Initialize counter to keep track of the quantity indexes
            obj.counter = Counter();
            
            %% Create the Electrical Network object
            obj.network = Network(bus_data, branch_data, nominalVoltage);
            % Register the DC network element to the DC system
            obj.registerElement(obj.network);
            
            %% Add the busses to the network
            obj.addBusses(bus_data);
        end
        
        function resultSet = simulate(obj, odesolver, solverOpts, t_span, u_ext, x0)
            % SIMULATE runs a dynamic simulation of the DC system
            %
            % INPUTS
            % odesolver         odesolver to use
            % solverOpts        odesolver options
            % t_span            time span of the simulation
            % u_ext             SignalBuilder object containing the external inputs as a function of time
            % x0                Initial conditions
            %
            % OUTPUTS
            % resultSet         structure containing the simulation results
            
            %% Validate DC system
            obj.validate;
            
            %% Check length of initial conditions vector x0
            if length(x0)~=obj.counter.getNumberOfStates()
                error('The number of initial conditions does not match with the number of state-variables of the DC system');
            end
            
            %% Check number of external inputs
            if length(u_ext.signalFcn)~=obj.counter.getNumberOfExtInputs
                error('The number of external inputs does not match the number of DC system external inputs.');
            end
            
            %% Set-up the interconnection matrices
            interconnM = obj.interconnectionMatrices();
            
            %% Run odesolver
            if ~isempty(solverOpts)
                % odesolver is a function callback according to ode23tb, ode45 or ode15s which is declared in the main routine.
                [resultSet.t, resultSet.x] = odesolver(@(t,x) obj.diffEqns(t, x, u_ext.getInputs(t), interconnM), t_span, x0, solverOpts);
            elseif length(t_span)>2
                % Use fixed time step solver
                resultSet.x = odesolver(@(t,x) obj.diffEqns(t, x, u_ext.getInputs(t), interconnM), t_span, x0);
                resultSet.t = t_span;
            else
                [resultSet.t, resultSet.x] = odesolver(@(t,x) obj.diffEqns(t, x, u_ext.getInputs(t), interconnM), t_span, x0);
            end
            
            %% Calculate external outputs
            resultSet.y_ext = zeros(size(resultSet.t,1),obj.counter.getNumberOfExtOutputs);
            for k=1:length(resultSet.t)
                resultSet.y_ext(k,:) = obj.calcExtOutputs(resultSet.t(k), resultSet.x(k,:), u_ext.getInputs(resultSet.t(k)), interconnM);
            end
        end
        
        function [] = connectComponent(obj, component, bus_no)
            % CONNECTCOMPONENT connect a component to the DC system at the corresponding bus
            %
            % INPUTS
            % component         the component(s) that should be added to the bus
            % bus_no            number referring to the bus to connect the component to
            %
            % Note: component can contain an array of components that will be connected to the bus number
            % specified by the array of bus_no.
            
            %% Verify function inputs
            if length(component)~=length(bus_no)
                error('The number of components does not match with the number of specified busses the components should be connected to.');
            end
            
            %% Adjust function inputs when only one component specified
            % convert to cell array
            if length(bus_no)==1
                component_{1} = component;
                clear component;
                component = component_;
            end
            
            %% Add components to the appropriate bus and register them to the DC system
            for k=1:length(component)
                if bus_no(k)>length(obj.busses)
                    error(['Unable to connect component ' component.label ' to BUS ' num2str(bus_no(k))]);
                else
                    obj.busses{bus_no(k)}.connectComponent(component{k});
                    obj.registerElement(component{k});
                end
            end
            
        end
        
        function [lbls,index] = getLabels(obj)
            % GETLABELS displays the labels of the network and components that make up the DC system
            
            %% Get counter values
            xc = obj.counter.getNumberOfStates('co');
            xv = obj.counter.getNumberOfStates('vo');
            
            %% Network
            [lbls_,index_] = obj.network.getLabels();
            lbls.x.network = lbls_.x; % Because of offset in state-vector x = [xv; xc; xn]
            index.x.network = index_.x + xv + xc;
            lbls.u_ext = lbls_.u_ext;
            index.u_ext = index_.u_ext;
            lbls.y_ext = lbls_.y_ext;
            index.y_ext = index_.y_ext;
            
            %% Busses
            
            lbls.x.vo = [];
            index.x.vo = [];
            lbls.x.co = [];
            index.x.co = [];
            
            for b=1:length(obj.busses)
                [lbls_,index_] = obj.busses{b}.getLabels;
                lbls.x.vo = [lbls.x.vo; lbls_.x.vo];
                index.x.vo = [index.x.vo; index_.x.vo];
                lbls.x.co = [lbls.x.co; lbls_.x.co];
                index.x.co = [index.x.co; index_.x.co];
                lbls.u_ext = [lbls.u_ext; lbls_.u_ext];
                index.u_ext = [index.u_ext; index_.u_ext];
                lbls.y_ext = [lbls.y_ext; lbls_.y_ext];
                index.y_ext = [index.y_ext; index_.y_ext];
            end
            
            index.x.co = index.x.co + xv; % Because of offset in state-vector x = [xv; xc; xn]
            
            %% External inputs
            t = table(index.u_ext,lbls.u_ext,'VariableNames',{'Index' 'External_inputs'});
            t = sortrows(t,{'Index'},{'ascend'});
            disp(t);
            
            %% External outputs
            t = table(index.y_ext,lbls.y_ext,'VariableNames',{'Index' 'External_outputs'});
            t = sortrows(t,{'Index'},{'ascend'});
            disp(t);
            
            %% State-variables
            
            % Join state-variables
            index.x = [index.x.vo; index.x.co; index.x.network];
            lbls.x = [lbls.x.vo; lbls.x.co; lbls.x.network];
            
            t = table(index.x,lbls.x,'VariableNames',{'Index' 'States'});
            t = sortrows(t,{'Index'},{'ascend'});
            disp(t);
        end
        
        function [X] = steadyState(obj, u_ext, x0, solverOpts)
            % Calculate the steady state value of the state variables using
            % fsolve, given the external inputs
            %
            % INPUTS
            % u_ext             array specifying the external inputs
            % x0                array specifying the steady-state initial search point (if empty, the function will calculate an appropriate start point)
            % solverOpts        the solver options for fsolve
            
            %% Validate DC system
            obj.validate;
            
            %% Validate external inputs
            if length(u_ext)~=obj.counter.getNumberOfExtInputs
                error('The number of external inputs does not match the number of DC system external inputs.');
            end
            
            %% Solve numerically
            interconnectionM = obj.interconnectionMatrices();
            
            if ~isempty(solverOpts)
                X = fsolve(@(x) obj.diffEqns(0, x, u_ext, interconnectionM), x0, solverOpts);
            else
                X = fsolve(@(x) obj.diffEqns(0, x, u_ext, interconnectionM), x0);
            end
        end
        
        function [X] = steadyState_bfs(obj,u_ext,max_tol)
            % STEADYSTATE_BFS determines the steady-state of the DC network components using the
            % backward-forward sweep method. Note: this method is only applicable to radial distribution
            % networks. The first bus in the network is the slack bus and remains constant at the nominal
            % network voltage with respect to its neutral terminal.
            %
            % INPUTS
            % u_ext         vector containing the external input steady-state values
            % max_tol       specified tolerance for the DC bus voltage
            
            %% Validate DC system
            obj.validate;
            
            %% Verify function input
            if length(u_ext)~=obj.counter.getNumberOfExtInputs
                error('The number of external inputs does not match the number of DC system external inputs.');
            end
            
            %% Determine the sequence in which the busses should be handled
            % Determine the parent busses and branches
            % i.e. the bus directly above the bus under consideration and the interconnecting branch
            [sequence, parent_busses, parent_branches, incM] = obj.network.breadthFirstSearch();
                
            %% Initialize bus voltage matrix
            % [ positive voltage        neutral voltage         negative voltage ] with respect to bus 1
            % neutral voltage terminal
            
            [noBusses, noBranches] = size(incM);
            
            bus_voltage = zeros(noBusses,3);
            bus_voltage(:,1) = obj.network.nominalVoltage;
            bus_voltage(:,2) = 0;
            bus_voltage(:,3) = -obj.network.nominalVoltage;
            
            %% BFS method
            tol = max_tol;
            iter = 0;
            max_iter = 1000;
            while tol>=max_tol && iter<=max_iter
                iter = iter + 1;
                if iter==max_iter
                    warning('Maximum number of iterations achieved. BFS algorithm stopped.');
                end
                %% Backward sweep
                % Initialize branch current matrix and voltage drop matrix to zero
                branch_current = zeros(noBranches,3);
                dv = branch_current;
                for b=1:noBusses-1
                    bus = obj.busses{sequence(end-b+1)};
                    branch = parent_branches(bus.number);
                    % Calculate the bus currents (ip=positive terminal bus current, in=negative terminal bus current, u_ext=vector containing steady-state external input values)
                    vp = bus_voltage(bus.number,1)-bus_voltage(bus.number,2);
                    vn = bus_voltage(bus.number,2)-bus_voltage(bus.number,3);
                    [ip, in] = bus.calcStStCurrent(vp, vn, u_ext);
                    % STEPS:
                    % 1. Calculate the branch current through the parent branch
                    % 2. Calculate the branch voltage drop
                    % 3. Update the parent bus voltage
                    switch obj.network.branch_data(branch,Name.getBranchField('BR_CONFIG'))
                        case Connection.po
                            % Branch segment containing a positive and neutral conductor
                            branch_current(branch,1) = incM(bus.number,branch) * (ip - incM(bus.number,[1:branch-1 branch+1:end]) * branch_current([1:branch-1 branch+1:end],1));
                            dv(branch,1) = obj.network.calcLineVoltageDrop(branch,'+',branch_current(branch,1));
                            bus_voltage(parent_busses(bus.number),1) = bus_voltage(bus.number,1) + dv(branch,1);
                        case Connection.on
                            % Branch segment containing a negative and neutral conductor
                            branch_current(branch,3) = incM(bus.number,branch) * (-in - incM(bus.number,[1:branch-1 branch+1:end]) * branch_current([1:branch-1 branch+1:end],3));
                            dv(branch,3) = obj.network.calcLineVoltageDrop(branch,'-',branch_current(branch,3));
                            bus_voltage(parent_busses(bus.number),3) = bus_voltage(bus.number,3) + dv(branch,3);
                        case Connection.pon
                            % Branch segment containing a positive, negative and neutral conductor
                            branch_current(branch,1) = incM(bus.number,branch) * (ip - incM(bus.number,[1:branch-1 branch+1:end]) * branch_current([1:branch-1 branch+1:end],1));
                            branch_current(branch,3) = incM(bus.number,branch) * (-in - incM(bus.number,[1:branch-1 branch+1:end]) * branch_current([1:branch-1 branch+1:end],3));
                            dv(branch,1) = obj.network.calcLineVoltageDrop(branch,'+',branch_current(branch,1));
                            dv(branch,3) = obj.network.calcLineVoltageDrop(branch,'-',branch_current(branch,3));
                            bus_voltage(parent_busses(bus.number),1) = bus_voltage(bus.number,1) + dv(branch,1);
                            bus_voltage(parent_busses(bus.number),3) = bus_voltage(bus.number,3) + dv(branch,3);
                    end
                    % Neutral branch voltage drop
                    branch_current(branch,2) = -branch_current(branch,1) - branch_current(branch,3);
                    dv(branch,2) = obj.network.calcLineVoltageDrop(branch,'o',branch_current(branch,2)); 
                    bus_voltage(parent_busses(bus.number),2) = bus_voltage(bus.number,2) + dv(branch,2);
                end
                
                %% Forward sweep
                % Update bus voltage 1
                tol = norm(bus_voltage(1,:) - obj.network.nominalVoltage*[1 0 -1]);
                bus_voltage(1,:) = obj.network.nominalVoltage*[1 0 -1];
                % Update other bus voltages
                for b=2:noBusses
                    bus = obj.busses{sequence(b)};
                    branch = parent_branches(bus.number);
                    switch obj.network.branch_data(branch,Name.getBranchField('BR_CONFIG'))
                        case Connection.po
                            bus_voltage(bus.number,1) = bus_voltage(parent_busses(bus.number),1) - dv(branch,1);
                        case Connection.on
                            bus_voltage(bus.number,3) = bus_voltage(parent_busses(bus.number),3) - dv(branch,3);
                        case Connection.pon
                            bus_voltage(bus.number,1) = bus_voltage(parent_busses(bus.number),1) - dv(branch,1);
                            bus_voltage(bus.number,3) = bus_voltage(parent_busses(bus.number),3) - dv(branch,3);
                    end
                    bus_voltage(bus.number,2) = bus_voltage(parent_busses(bus.number),2) - dv(branch,2);
                end
                disp(['Tolerance: ' num2str(tol)]);
            end
            
            
            
            %% Determine steady-state state-variables
            
            % Initialize state-vector
            xv = obj.counter.getNumberOfStates('vo');
            xc = obj.counter.getNumberOfStates('co');
            xn = obj.counter.getNumberOfStates('network');
            
            X = zeros(xv+xc+xn,1);
            
            % Determine steady-state network state-variables
            % i.e. the branch currents                
            for b=1:noBranches
                branch_index = obj.network.getBranchIndex(b,'+');
                if branch_index~=0, X(xv+xc+branch_index) = branch_current(b,1); end
                branch_index = obj.network.getBranchIndex(b,'-');
                if branch_index~=0, X(xv+xc+branch_index) = branch_current(b,3); end
            end
            
            % Determine steady-state component state-variables
            for b=1:noBusses
                vp = bus_voltage(b,1)-bus_voltage(b,2);
                vn = bus_voltage(b,2)-bus_voltage(b,3);
                [index, X_] = obj.busses{b}.calcStSt(vp,vn,u_ext);
                X(index.vo) = X_.vo;
                X(index.co + xv) = X_.co;
            end

        end
        
        function [A, sysPoles] = linearSys(obj, X, U_ext)
            % Return a linearized system matrix around the equilibrium
            % point X, that corresponds to the external input u_ext
            %
            % INPUTS
            % X         The equilibrium point around which the linearized system will be calculated
            % u_ext     The external inputs that are applied to the system
            %
            % OUTPUTS
            % A         The linearized system matrix (dx/dt = A x)
            % sysPoles  The poles of the system matrix (i.e. the eigenvalues of A)
            
            %% Convert u_ext to column vector if necessary
            if size(U_ext,1)<size(U_ext,2) 
                U_ext = transp(U_ext);
            end

            %% Verify whether X is an equilibrium point
            iM = obj.interconnectionMatrices();
            if norm(obj.diffEqns(0, X, U_ext, iM))>0.2
                warning(message('Error obtaining a linearized system: the specified point is not an equilibrium point.'));
            end
            
            %% Segment state vectors
            x_ = X; clear X;
            xv = obj.counter.getNumberOfStates('vo');
            xc = obj.counter.getNumberOfStates('co');
            xn = obj.counter.getNumberOfStates('network');
            X.vo = x_(1:xv);
            X.co = x_(xv+1:xv+xc);
            X.network = x_(xv+xc+1:xv+xc+xn);
            
            %% Calculate steady-state internal inputs
            
            % Calculate subsystem internal outputs
            % Structure y_int (y_int.network, y_int.vo, y_int.co) is calculated by method in DC System
            t = 0;
            Y_int = obj.calcIntOutputs(t, X, U_ext);
            
            % Calculate subsystem internal inputs
            % Structure u_int (u_int.network, u_int.vo, u_int.co) is calculated by method in DC System
            % Based upon interconnection matrices
            U_int = obj.calcIntInputs(Y_int, iM);
            
            %% Obtain component system function derivatives
            
            dxdx = struct('vo',[],'co',[],'network',[]);
            dxdu_int = struct('vo',[],'co',[],'network',[]);
            dy_intdx = struct('vo',[],'co',[],'network',[]);
            
            for b=1:length(obj.busses)
                [dxdx_, dxdu_int_, dy_intdx_] = obj.busses{b}.ddiffEqns(t,X,U_int,U_ext);
                dxdx.vo(dxdx_.vo.index.row,dxdx_.vo.index.col) = dxdx_.vo.value;
                dxdx.co(dxdx_.co.index.row,dxdx_.co.index.col) = dxdx_.co.value;
                dxdu_int.vo(dxdu_int_.vo.index.row,dxdu_int_.vo.index.col) = dxdu_int_.vo.value;
                dxdu_int.co(dxdu_int_.co.index.row,dxdu_int_.co.index.col) = dxdu_int_.co.value;
                dy_intdx.vo(dy_intdx_.vo.index.row,dy_intdx_.vo.index.col) = dy_intdx_.vo.value;
                dy_intdx.co(dy_intdx_.co.index.row,dy_intdx_.co.index.col) = dy_intdx_.co.value;
            end
            
            %% Set-up the system matrix
            
            An = obj.network.systemMatrix.A;
            Bn = obj.network.systemMatrix.B;
            Cn = obj.network.systemMatrix.C;
            A = [  dxdx.vo                              dxdu_int.vo * iM.H3 * dy_intdx.co       dxdu_int.vo * iM.H2 * Cn;
                   dxdu_int.co * iM.H1 * dy_intdx.vo    dxdx.co                                 zeros(xc,xn);
                   Bn * iM.H0 * dy_intdx.vo             zeros(xn,xc)                            An];
            
            %% Calculate the system poles
            sysPoles = eig(A);
        end
    end
end