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

classdef Network < Element
    % NETWORK The class NETWORK handles the network equations
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 16, 2015
    
    properties
        % Incidence matrix of the circuit	
        inc
        % Nominal DC voltage in the network between phase and neutral conductor
        nominalVoltage
        % Network data
        bus_data;
        branch_data;
        % Structure containing the system matrices (systemMatrix.A .B .C)
        systemMatrix
    end
    
    properties(Access=private)
        % Structure containing the network matrices (networkMatrix.L .R .V)
        networkMatrix
        % Matrix containing the index of the bus voltages in the input vector u_N
        bus_index;
        % Matrix containing the index of the branch currents in the state vector i_L
        branch_index;
    end
    
    methods(Static)
        function [inc] = incidenceMatrix(bus_data, branch_data)
            % Construct the incidence matrix based on the bus_data and the
            % branch_data
            [numBranches, ~] = size(branch_data);
            [numBusses, ~] = size(bus_data);
            inc = zeros(numBusses, numBranches);
            for k=1:numBranches
                F_bus = branch_data(k,1);
                T_bus = branch_data(k,2);
                inc(F_bus,k) = -1;
                inc(T_bus,k) = 1;
            end
        end
    end
    
    methods(Access=private)
        
        function [C] = busCurrentMatrix(obj)
            % BUSCURRENTMATRIX returns the matrix that relates the bus currents to the state-variable inductor currents
            
            %% Set-up incidence matrix
            % as if the DC distribution network is bipolar (and eliminate redundant rows and columns
            % subsequently)
            incM = Network.incidenceMatrix(obj.bus_data, obj.branch_data);
            
            Mp = [incM               zeros(size(incM))];
            Mn = [zeros(size(incM))  -incM];
            
            C_ = zeros(2*size(incM));
            C = C_;
            
            for k=1:size(Mp,1) % rows
                C_(2*k-1,:) = Mp(k,:);
                C_(2*k,:) = Mn(k,:);
            end
            
            no_branches = size(incM,2);
            for k=1:no_branches % columns
                C(:,2*k-1) = C_(:,k);
                C(:,2*k) = C_(:,k+no_branches);
            end

            %% Eliminate redundant bus currents
            remove = [];
            for b=1:size(obj.bus_data,1)
                switch obj.bus_data(b,Name.getBusField('BUSDC_CONFIG'))
                    case Connection.po
                        remove = [remove 2*b];
                    case Connection.on
                        remove = [remove 2*b-1];
                end
            end
            C(remove,:) = [];
            
            %% Eliminate redundant inductor currents
            remove = [];
            for b=1:size(obj.branch_data,1)
                switch obj.branch_data(b,Name.getBranchField('BR_CONFIG'))
                    case Connection.po
                        remove = [remove 2*b];
                    case Connection.on
                        remove = [remove 2*b-1];
                end
            end
            C(:,remove) = [];
        end
        
        function [L, R, V, bus_index, branch_index] = networkMatrices(obj)
            % NETWORKMATRICES sets up the network matrices of the DC distribution network
            %
            % L di/dt = R i + V v
            %
            % OUTPUTS
            % L             matrix containing the line inductances (conform equation above)
            % R             matrix containing the line resistances (conform equation above)
            % V             matrix containing the incidence matrix (conform equation above)
            % bus_index     matrix containing the bus number, positive terminal voltage index and
            %               negative terminal voltage index
            % branch_index  matrix containing the branch number, positive branch current index and 
            %               negative branch current index
            
            %% Initialize network matrices
            % Inductance matrix
            L = zeros(2*size(obj.branch_data,1));
            % Resistance matrix
            R = L;
            % Voltage matrix
            V = zeros(2*size(obj.branch_data,1),2*size(obj.bus_data,1));
            
            %% Set-up network matrices
            
            branch_index = zeros(size(obj.branch_data,1),3);
            
            for l=1:size(obj.branch_data,1)
                
                F_BUS = obj.branch_data(l, Name.getBranchField('F_BUS'));
                T_BUS = obj.branch_data(l, Name.getBranchField('T_BUS'));
                
                % Check feasibility
                if obj.branch_data(l,Name.getBranchField('BR_CONFIG'))==Connection.on && obj.bus_data(F_BUS,Name.getBusField('BUSDC_CONFIG'))==Connection.po
                    error(['Cannot connect branch ' num2str(l) ' to bus ' num2str(F_BUS)]);
                elseif obj.branch_data(l,Name.getBranchField('BR_CONFIG'))==Connection.on && obj.bus_data(T_BUS,Name.getBusField('BUSDC_CONFIG'))==Connection.po
                    error(['Cannot connect branch ' num2str(l) ' to bus ' num2str(T_BUS)]);
                elseif obj.branch_data(l,Name.getBranchField('BR_CONFIG'))==Connection.po && obj.bus_data(F_BUS,Name.getBusField('BUSDC_CONFIG'))==Connection.on
                    error(['Cannot connect branch ' num2str(l) ' to bus ' num2str(F_BUS)]);
                elseif obj.branch_data(l,Name.getBranchField('BR_CONFIG'))==Connection.po && obj.bus_data(T_BUS,Name.getBusField('BUSDC_CONFIG'))==Connection.on
                    error(['Cannot connect branch ' num2str(l) ' to bus ' num2str(T_BUS)]);
                end
                
                Lp = obj.branch_data(l, Name.getBranchField('BR_L1'));
                Ln = obj.branch_data(l, Name.getBranchField('BR_L2'));
                Lo = obj.branch_data(l, Name.getBranchField('BR_Ln'));
                Rp = obj.branch_data(l, Name.getBranchField('BR_R1'));
                Rn = obj.branch_data(l, Name.getBranchField('BR_R2'));
                Ro = obj.branch_data(l, Name.getBranchField('BR_Rn'));
                
                L(2*l-1,2*l-1) = -(Lp+Lo);
                L(2*l-1,2*l) = -Lo;
                L(2*l,2*l-1) = Lo;
                L(2*l,2*l) = Lo+Ln;
                
                R(2*l-1,2*l-1) = Rp+Ro;
                R(2*l-1,2*l) = Ro;
                R(2*l,2*l-1) = -Ro;
                R(2*l,2*l) = -(Ro+Rn);
                
                V(2*l-1,2*F_BUS-1) = -1;
                V(2*l-1,2*T_BUS-1) = 1;
                V(2*l,2*F_BUS) = -1;
                V(2*l,2*T_BUS) = 1;
            end
            
            %% Eliminate redundant branches and set-up branch index matrix
            branch_index_counter = 0;
            remove = [];
            for l=1:size(obj.branch_data,1)
                branch_index(l,1) = l;
                switch obj.branch_data(l,Name.getBranchField('BR_CONFIG'))
                    case Connection.po
                        remove = [remove 2*l];
                        branch_index_counter = branch_index_counter + 1;
                        branch_index(l,2) = branch_index_counter;
                    case Connection.on
                        remove = [remove 2*l-1];
                        branch_index_counter = branch_index_counter + 1;
                        branch_index(l,3) = branch_index_counter;
                    case Connection.pon
                        branch_index_counter = branch_index_counter + 1;
                        branch_index(l,2) = branch_index_counter;
                        branch_index_counter = branch_index_counter + 1;
                        branch_index(l,3) = branch_index_counter;
                end
            end
            
            L(remove,:) = [];
            L(:,remove) = [];
            R(remove,:) = [];
            R(:,remove) = [];
            V(remove,:) = [];
        
            %% Eliminate redundant busses and set-up bus index matrix
            remove_cols = [];
            bus_index = zeros(size(obj.bus_data,1),3);
            bus_index_counter = 0;
            for b=1:size(obj.bus_data,1)
                bus_index(b,1) = b;
                switch obj.bus_data(b,Name.getBusField('BUSDC_CONFIG'))
                    case Connection.po
                        remove_cols = [remove_cols 2*b];
                        bus_index_counter = bus_index_counter + 1;
                        bus_index(b,2) = bus_index_counter;
                    case Connection.on
                        remove_cols = [remove_cols 2*b-1];
                        bus_index_counter = bus_index_counter + 1;
                        bus_index(b,3) = bus_index_counter;
                    case Connection.pon
                        bus_index_counter = bus_index_counter + 1;
                        bus_index(b,2) = bus_index_counter;
                        bus_index_counter = bus_index_counter + 1;
                        bus_index(b,3) = bus_index_counter;
                end
            end
            V(:,remove_cols) = [];
        
        end
        
        function [A, B, C] = systemMatrices(obj)
            A = obj.networkMatrix.L\obj.networkMatrix.R;
            B = obj.networkMatrix.L\obj.networkMatrix.V;
            C = obj.busCurrentMatrix();
        end
        
    end
    
    methods
        
        function obj = Network(bus_data, branch_data, nominalVoltage)
            % NETWORK creates a new network object based upon the bus and branch data
            
            %% Create a new network object
            obj.bus_data = bus_data;
            obj.branch_data = branch_data;
            obj.nominalVoltage = nominalVoltage;
            obj.label = 'Network';
            
            %% Set-up equation matrices
            [obj.networkMatrix.L, obj.networkMatrix.R, obj.networkMatrix.V, obj.bus_index, obj.branch_index] = obj.networkMatrices();
            
            %% Set-up system matrices
            [obj.systemMatrix.A, obj.systemMatrix.B, obj.systemMatrix.C] = obj.systemMatrices();
            
            %% Create quantities
            % Based upon the bus_index and branch_index matrices, quantities are created
            for b=1:size(obj.branch_index,1)
                if obj.branch_index(b,2)~=0
                    % Create positive branch current quantity
                    obj.x{end+1} = Quantity(['i_{L' num2str(b) '+}'],'A');
                    obj.y_ext{end+1} = Quantity(['i_{L' num2str(b) '+}'],'A');
                    obj.y_ext{end+1} = Quantity(['P_{L' num2str(b) '+}'],'A');
                end
                if obj.branch_index(b,3)~=0
                    % Create negative branch current quantity
                    obj.x{end+1} = Quantity(['i_{L' num2str(b) '-}'],'A');
                    obj.y_ext{end+1} = Quantity(['i_{L' num2str(b) '-}'],'A');
                    obj.y_ext{end+1} = Quantity(['P_{L' num2str(b) '-}'],'A');
                end
                obj.y_ext{end+1} = Quantity(['i_{L' num2str(b) 'o}'],'A');
                obj.y_ext{end+1} = Quantity(['P_{L' num2str(b) 'o}'],'A');
            end
            
            % Bus current (y_N) and voltage (u_N) vectors
            for b=1:size(obj.bus_index,1)
                if obj.bus_index(b,2)~=0
                    % Create positive bus current/voltage quantity
                    obj.y_int{end+1} = Quantity(['y_{N' num2str(b) '+}'],'A');
                    obj.u_int{end+1} = Quantity(['u_{N' num2str(b) '+}'],'V');
                end
                if obj.bus_index(b,3)~=0
                    % Create negative bus current/voltage quantity
                    obj.y_int{end+1} = Quantity(['y_{N' num2str(b) '-}'],'A');
                    obj.u_int{end+1} = Quantity(['u_{N' num2str(b) '-}'],'V');
                end
            end
            
        end
        
        function [index, dx] = diffEqns(obj,t,x,u_int,u_ext)
            % DIFFEQNS calculates dx of the Network object, i.e. a vector that contains the
            % time-derivative of the branch current state-variables
            %
            % INPUTS
            % t         current time value
            % x         current state-variable vector
            % u_int     internal input vector, i.e. the bus terminal voltages
            % u_ext     external input vector (not used)
            
%             i_L_index = []; % REMOVE
%             for k=1:length(obj.x)
%                 i_L_index = [i_L_index; obj.x{k}.index]; %#ok<*AGROW>
%             end
            i_L_index = Quantity.getIndex(obj.x);
            i_L = x(i_L_index);
            
%             u_N_index = []; % REMOVE
%             for k=1:length(obj.u_int)
%                 u_N_index = [u_N_index; obj.u_int{k}.index];
%             end
            u_N_index = Quantity.getIndex(obj.u_int);
            u_N = u_int(u_N_index);
            
            dx = obj.systemMatrix.A * i_L + obj.systemMatrix.B * u_N;
            index = i_L_index;
        end
        
        function [index, y_int] = calcIntOutputs(obj,t,x)
            % CALCINTOUTPUTS calculates y_N of the Network object, i.e. a vector that contains the bus
            % currents as well as the index where they should be placed
            %
            % INPUTS
            % t         current time value
            % x         current state-variable vector
            % u_ext     external input vector (not used)
            
            i_L_index = [];
            for k=1:length(obj.x)
                i_L_index = [i_L_index; obj.x{k}.index]; %#ok<AGROW>
            end
            i_L = x(i_L_index);
            
            y_N_index = [];
            for k=1:length(obj.y_int)
                y_N_index = [y_N_index; obj.y_int{k}.index]; %#ok<AGROW>
            end
            
            if size(i_L,2)>size(i_L,1)
                % Convert i_L to column vector
                i_L = transp(i_L);
            end
            
            y_int.network = obj.systemMatrix.C * i_L;
            index = y_N_index;
        end
        
        function [index, y_ext] = calcExtOutputs(obj,t,x,u_int,u_ext)
            % CALCEXTOUTPUTS calculates the external outputs of the Network object (i.e. the branch
            % currents and losses)
            
            y_ext = zeros(0,1);
            for b=1:size(obj.branch_index,1)
                ip = 0;
                in = 0;
                if obj.branch_index(b,2)~=0
                    % Positive branch current
                    ip = x(obj.branch_index(b,2));
                    y_ext(end+1) = ip;
                    % Power loss: R*i^2
                    y_ext(end+1) = obj.branch_data(b,Name.getBranchField('BR_R1')) * x(obj.branch_index(b,2))^2;
                end
                if obj.branch_index(b,3)~=0
                    % Negative branch current
                    in = x(obj.branch_index(b,3));
                    y_ext(end+1) = in;
                    % Power loss: R*i^2
                    y_ext(end+1) = obj.branch_data(b,Name.getBranchField('BR_R2')) * x(obj.branch_index(b,3))^2;
                end
                % Neutral branch current
                y_ext(end+1) = -ip-in;
                % Power loss: R*i^2
                y_ext(end+1) = obj.branch_data(b,Name.getBranchField('BR_Rn')) * (ip+in)^2;
            end
            
            index = zeros(length(obj.y_ext),1);
            for k=1:length(index)
                index(k) = obj.y_ext{k}.index;
            end
        end
        
        function [index_p, index_n] = getIndex(obj,bus)
            % GETINDEX returns the index of the internal in and outputs for
            % the specified bus
            %
            % INPUTS
            % bus        the bus for which the index of the internal inputs
            %            should be returned
            %
            % OUTPUTS
            % index_p   indexes corresponding to the postive phase
            % index_n   indexes corresponding to the negative phase
            
            index_p = [];
            index_n = [];
            
            if obj.bus_index(bus.number,2)>0
                index_p = obj.bus_index(bus.number,2);
            end
            
            if obj.bus_index(bus.number,3)>0
                index_n = obj.bus_index(bus.number,3);
            end
        end
        
        function [index] = getBranchIndex(obj,branchNo,phase)
            % GETBRANCHINDEX returns the index of the branch current in the state vector
            %
            % INPUTS
            % branchNo      the number of the branch
            % phase         the positive(+) or negative(-) branch current
            %
            % OUTPUTS
            % index         the index of the branch current in the state-vector
            %               Note: if index=0, the branch segment is non-existing (e.g. when the branch is
            %               unipolar, it contains no negative or positive branch segment)
            
            switch phase
                case '+'
                    col = 2;
                case '-'
                    col = 3;
            end
            
            index = obj.branch_index(branchNo, col);
        end
        
        function [sequence, parent_busses, parent_branches, incM] = breadthFirstSearch(obj)
            % BREADTHFIRSTSEARCH returns the sequence in which the busses should be traversed in the
            % backward-forward sweep steady-state solution algorithm.
            
            %% Incidence matrix
            % Unoriented incidence matrix
            incM = obj.incidenceMatrix(obj.bus_data,obj.branch_data);
            [no_busses, no_branches] = size(incM);
            
            %% Adjacency matrix
            adjM = zeros(no_busses);
            for k=1:no_branches
                connected_busses = find(abs(incM(:,k))==1);
                adjM(connected_busses(1),connected_busses(2))=1;
                adjM(connected_busses(2),connected_busses(1))=1;
            end
            
            %% Traverse tree
            % Starting from bus 1, add the children of bus 1 to the sequence until all busses added to
            % the sequence
            sequence = zeros(no_busses,1);
            sequence(1) = 1;
            index = 2;
            for b=1:no_busses
                for h=(b+1):no_busses
                    if (adjM(b,h)==1)
                        % Bus h is a child of Bus c
                        if index>no_busses
                            % If the sequence vector contains more elements than the number of busses,
                            % some busses were added twice, indicating that the network is not radial
                            error('The BFS algorithm cannot only be applied to radial networks');
                        else
                            sequence(index) = h;
                            index = index + 1;
                        end
                    end
                end
            end
            
            parent_branches = zeros(no_busses,1);
            parent_busses = zeros(no_busses,1);
            incM_ = incM;
            sequence(1) = 1;
            index = 2;
            next_bus = 1;
            for j=1:no_busses
                b = next_bus(1);
                % List branches connected to the current bus
                branch_index = find(abs(incM_(b,:))==1);
                % For each branch: determine the child bus (cb) and set its parent bus and parent branch
                for k=1:length(branch_index)
                    cb = find(incM_(:,branch_index(k))==-(incM_(b,branch_index(k))));
                    next_bus(end+1) = cb;
                    sequence(index) = cb;
                    index = index + 1;
                    parent_busses(cb) = b;
                    parent_branches(cb) = branch_index(k);
                end
                incM_(:,branch_index) = 0;
                next_bus(1) = [];
            end
        end
        
        function [Dv] = calcLineVoltageDrop(obj,branchNo,phase,current)
            % CALCLINEVOLTAGEDROP calculates the steady-state voltage drop along the specified line
            % segment, provided the branch current
            %
            % INPUTS
            % branchNo      the number of the branch
            % phase         the positive(+), negative(-) or neutral(o) branch voltage drop
            % current       the branch current to calculate the voltage drop            
            
            switch phase
                case '+'
                    R = obj.branch_data(branchNo,Name.getBranchField('BR_R1'));
                case '-'
                    R = obj.branch_data(branchNo,Name.getBranchField('BR_R2'));
                case 'o'
                    R = obj.branch_data(branchNo,Name.getBranchField('BR_Rn'));
            end
            
            Dv = R * current;            
        end
        
    end
end

