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

classdef Bus < Element
    % BUS The class BUS constitutes the busses in the DC network
    % 
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 16, 2015
    
    properties
        % The bus number
        number
    end
    
    properties(SetAccess=private)
        % Bus connection that defines the bus terminals
        connection = Connection.getDefaultValue;
        % Array containing the VO components
        vo_components = cell(0,1);
        % Array containing the CO components
        co_components = cell(0,1);
    end
    
    methods
        
        function obj = Bus(bus_data)
            % Create a new bus object and collect the branch data
            %
            % INPUTS
            % bus_data          a matrix structure as defined in the manual that contains the data of the bus
            obj.connection = bus_data(Name.getBusField('BUSDC_CONFIG'));
            obj.label = ['BUS ' num2str(bus_data(Name.getBusField('BUSDC_I')))];
            obj.number = bus_data(Name.getBusField('BUSDC_I'));
            obj.vo_components = cell(0,1);
            obj.co_components = cell(0,1);
        end
        
        function [index, dx] = diffEqns(obj,t,x,u_int,u_ext)
            % DIFFEQNS returns the differential equations of the components connected to the Bus
            %
            % INPUTS
            % t                 current time value
            % x                 state vector
            % u_int             internal input vector
            % u_ext             external input vector
            %
            % OUTPUTS
            % index             structure (.vo .co) containing the indexes of the VO and CO components in
            %                   the state vector derivative dx
            % dx                structure (.vo .co) containing the state-vector derivative values
            
            index = struct('vo',[],'co',[]);
            dx = struct('vo',[],'co',[]);
            
            %% Differential equations of VO components
            for c=1:length(obj.vo_components)
                dx_= obj.vo_components{c}.diffEqns(t,x.vo(obj.vo_components{c}.getIndex('x')),u_int.vo(obj.vo_components{c}.getIndex('u_int')),u_ext(obj.vo_components{c}.getIndex('u_ext')));
                index_ = obj.vo_components{c}.getIndex('x');
                index.vo = [index.vo; index_]; %#ok<*AGROW>
                dx.vo = [dx.vo; dx_];
            end
            
            %% Differential equations of CO components
            for c=1:length(obj.co_components)
                dx_ = obj.co_components{c}.diffEqns(t,x.co(obj.co_components{c}.getIndex('x')),u_int.co(obj.co_components{c}.getIndex('u_int')),u_ext(obj.co_components{c}.getIndex('u_ext')));
                index_ = obj.co_components{c}.getIndex('x');
                index.co = [index.co; index_];
                dx.co = [dx.co; dx_];
            end
        end
        
        function [dxdx, dxdu_int, dy_intdx] = ddiffEqns(obj,t,x,u_int,u_ext)
            % DIFFEQNS returns the second-derivative equations of the components connected to the Bus as
            % a Quantity object
            %
            % INPUTS
            % t                 current time value
            % x                 structure containing the state vectors (x.vo, x.co)
            % u_int             internal input vector
            % u_ext             external input vector
            
            dxdx.vo = Quantity('dxdx',''); dxdx.co = Quantity('dxdx',''); dxdx_ = Quantity('dxdx','');
            dxdu_int.vo = Quantity('dxdu_int',''); dxdu_int.co = Quantity('dxdu_int',''); dxdu_int_ = Quantity('dxdu_int','');
            dy_intdx.vo = Quantity('dy_intdx',''); dy_intdx.co = Quantity('dy_intdx',''); dy_intdx_ = Quantity('dy_intdx','');
            
            % Note: the index of each quantity contains a structure holding the row, respectively column, indexes
            % of each linearized subsystem matrix
            
            %% Equations of VO components
            x_index = [];
            u_int_index = [];
            y_int_index = [];
            for c=1:length(obj.vo_components)
                [dxdx_.value, dxdu_int_.value, dy_intdx_.value] = obj.vo_components{c}.ddiffEqns(t,x.vo(obj.vo_components{c}.getIndex('x')),u_int.vo(obj.vo_components{c}.getIndex('u_int')),u_ext(obj.vo_components{c}.getIndex('u_ext')));
                x_index_ = obj.vo_components{c}.getIndex('x');
                u_int_index_ = obj.vo_components{c}.getIndex('u_int');
                y_int_index_ = obj.vo_components{c}.getIndex('y_int');
                x_index = max([x_index x_index_]);
                u_int_index = max([u_int_index u_int_index_]);
                y_int_index = max([y_int_index y_int_index_]);
                dxdx.vo.value(x_index_,x_index_) = dxdx_.value;
                dxdu_int.vo.value(x_index_,u_int_index_) = dxdu_int_.value;
                dy_intdx.vo.value(y_int_index_,x_index_) = dy_intdx_.value;
            end
            dxdx.vo.index.row = 1:x_index;
            dxdx.vo.index.col = 1:x_index;
            dxdu_int.vo.index.row = 1:x_index;
            dxdu_int.vo.index.col = 1:u_int_index;
            dy_intdx.vo.index.row = 1:y_int_index;
            dy_intdx.vo.index.col = 1:x_index;
            
            %% Equations of CO components
            x_index = [];
            u_int_index = [];
            y_int_index = [];
            for c=1:length(obj.co_components)
                [dxdx_.value, dxdu_int_.value, dy_intdx_.value] = obj.co_components{c}.ddiffEqns(t,x.co(obj.co_components{c}.getIndex('x')),u_int.co(obj.co_components{c}.getIndex('u_int')),u_ext(obj.co_components{c}.getIndex('u_ext')));
                x_index_ = obj.co_components{c}.getIndex('x');
                u_int_index_ = obj.co_components{c}.getIndex('u_int');
                y_int_index_ = obj.co_components{c}.getIndex('y_int');
                x_index = max([x_index x_index_]);
                u_int_index = max([u_int_index u_int_index_]);
                y_int_index = max([y_int_index y_int_index_]);
                dxdx.co.value(x_index_,x_index_) = dxdx_.value;
                dxdu_int.co.value(x_index_,u_int_index_) = dxdu_int_.value;
                dy_intdx.co.value(y_int_index_,x_index_) = dy_intdx_.value;
            end
            dxdx.co.index.row = 1:x_index;
            dxdx.co.index.col = 1:x_index;
            dxdu_int.co.index.row = 1:x_index;
            dxdu_int.co.index.col = 1:u_int_index;
            dy_intdx.co.index.row = 1:y_int_index;
            dy_intdx.co.index.col = 1:x_index;
        end 
        
        function [index, y_ext] = calcExtOutputs(obj,t,x,u_int,u_ext)
            % CALCEXTOUTPUTS calculates the external output values of the components connected to the Bus
            %
            % INPUTS
            % t                 current time value
            % x                 state vector structure (.vo and .co)
            % u_int             internal input vector
            % u_ext             external input vector
            
            index = [];
            y_ext = [];
            
            components = [obj.vo_components; obj.co_components];
            
            for c=1:length(components)
                % CO or VO component
                if components{c}.type == Type.VO, s = 'vo';
                elseif components{c}.type == Type.CO, s = 'co'; end
                y_ext_ = components{c}.calcExtOutputs(t,x.(s)(components{c}.getIndex('x')),u_int.(s)(components{c}.getIndex('u_int')),u_ext(components{c}.getIndex('u_ext')));
                index_ = components{c}.getIndex('y_ext');
                index = [index; index_];
                y_ext = [y_ext; y_ext_];
            end
        end
        
        function [index, y_int] = calcIntOutputs(obj,t,x,u_ext)
            % CALCEXTOUTPUTS calculates the external output values of the components connected to the Bus
            %
            % INPUTS
            % t                 current time value
            % x                 state vector
            % u_ext             external input vector
            
            y_int.vo = [];
            y_int.co = [];
            
            index.vo = [];
            index.co = [];
            y_int.vo = [];
            y_int.co = [];
            
            for c=1:length(obj.co_components)
                y_int_ = obj.co_components{c}.calcIntOutputs(t, x.co(obj.co_components{c}.getIndex('x')), u_ext(obj.co_components{c}.getIndex('u_ext')));
                index_ = obj.co_components{c}.getIndex('y_int');
                index.co = [index.co; index_];
                y_int.co = [y_int.co; y_int_];
            end
            
            for c=1:length(obj.vo_components)
                y_int_ = obj.vo_components{c}.calcIntOutputs(t, x.vo(obj.vo_components{c}.getIndex('x')), u_ext(obj.vo_components{c}.getIndex('u_ext')));
                index_ = obj.vo_components{c}.getIndex('y_int');
                index.vo = [index.vo; index_];
                y_int.vo = [y_int.vo; y_int_];
            end
            
        end
        
        function [] = connectComponent(obj, component)
            % CONNECTCOMPONENT Connect a component to the bus
            %
            % INPUTS
            % component         the component that should be added to the bus
            %
            % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be)
            % - Last modification: July 14, 2015
            
            %% Verify feasible connection
            % Unipolar negative and unipolar positive components cannot be connected to unipolar
            % positive, resp. unipolar negative busses
            
            if component.connection==Connection.po && obj.connection==Connection.on
                error(['Positive terminal components ' component.label ' cannot be connected to negative terminal busses.']);
            elseif component.connection==Connection.on && obj.connection==Connection.po
                error(['Negative terminal component ' component.label ' cannot be connected to positive terminal busses.']);
            end
                
            %% Connect component
            switch component.type
                case Type.VO
                    % Check that only one VO component is connected to each terminal
                    for c=1:length(obj.vo_components)
                        if component.connection == obj.vo_components{c}.connection
                            error(['Unable to connect ' component.label ' to BUS ' num2str(obj.number)]);
                        elseif (obj.vo_components{c}.connection == Connection.pon) || (~isempty(obj.vo_components) && component.connection==Connection.pon)
                            % Trying to connect a VO component to a bus where already a bipolar VO
                            % component is connected -OR- trying to connect a bipolar VO component to a
                            % bus where already another VO component is connected
                            error(['Unable to connect ' component.label ' to BUS ' num2str(obj.number)]);
                        end
                    end
                    obj.vo_components{end+1,1} = component;
                case Type.CO
                    obj.co_components{end+1,1} = component;
            end
        end
        
        function [index_p, index_n] = getIndex(obj,id)
            % GETINDEX returns the index of the internal in and outputs, as
            % specified by id.
            %
            % id can be either:
            % v         internal inputs of CO components
            % i         internal inputs of VO components
            % y_c       internal outputs of CO components
            % y_v       internal outputs of VO components
            %
            % INPUTS
            % id        (as specified above)
            %
            % OUTPUTS
            % index_p   indexes corresponding to the postive phase
            % index_n   indexes corresponding to the negative phase
            index_p = [];
            index_n = [];
            switch id
                case 'v'
                    % internal inputs of CO components
                    for c=1:length(obj.co_components)
                        if obj.co_components{c}.connection==Connection.po
                            index_p = [index_p; obj.co_components{c}.u_int{1}.index];
                        elseif obj.co_components{c}.connection==Connection.on
                            index_n = [index_n; obj.co_components{c}.u_int{1}.index];
                        elseif obj.co_components{c}.connection==Connection.pon
                            index_p = [index_p; obj.co_components{c}.u_int{1}.index];
                            index_n = [index_n; obj.co_components{c}.u_int{2}.index];
                        end
                    end
                case 'i'
                    % internal inputs of VO components
                    for c=1:length(obj.vo_components)
                        if obj.vo_components{c}.connection==Connection.po
                            index_p = [obj.vo_components{c}.u_int{1}.index];
                        elseif obj.vo_components{c}.connection==Connection.on
                            index_n = [obj.vo_components{c}.u_int{1}.index];
                        elseif obj.vo_components{c}.connection==Connection.pon
                            index_p = [obj.vo_components{c}.u_int{1}.index];
                            index_n = [obj.vo_components{c}.u_int{2}.index];
                        end
                    end
                case 'y_c'
                    % internal outputs of CO components
                    for c=1:length(obj.co_components)
                        if obj.co_components{c}.connection==Connection.po
                            index_p = [index_p; obj.co_components{c}.y_int{1}.index];
                        elseif obj.co_components{c}.connection==Connection.on
                            index_n = [index_n; obj.co_components{c}.y_int{1}.index];
                        elseif obj.co_components{c}.connection==Connection.pon
                            index_p = [index_p; obj.co_components{c}.y_int{1}.index];
                            index_n = [index_n; obj.co_components{c}.y_int{2}.index];
                        end
                    end
                case 'y_v'
                    % internal outputs of VO components
                    for c=1:length(obj.vo_components)
                        if obj.vo_components{c}.connection==Connection.po
                            index_p = [index_p; obj.vo_components{c}.y_int{1}.index];
                        elseif obj.vo_components{c}.connection==Connection.on
                            index_n = [index_n; obj.vo_components{c}.y_int{1}.index];
                        elseif obj.vo_components{c}.connection==Connection.pon
                            index_p = [index_p; obj.vo_components{c}.y_int{1}.index];
                            index_n = [index_n; obj.vo_components{c}.y_int{2}.index];
                        end
                    end
            end
        end
        
        function [labels, index] = getLabels(obj)
            % GETLABELS returns the labels of the quantities (i.e. parameters, states, internal
            % inputs, external inputs, internal outputs and external outputs) of the components connected
            % to the bus
            vars = {'x' 'u_ext' 'y_ext' 'parameters'};
            labels.x.co = [];
            labels.x.vo = [];
            index.x.co = [];
            index.x.vo = [];
            for j=2:length(vars)
                labels.(vars{j}) = [];
                index.(vars{j}) = [];
            end
            components = [obj.vo_components; obj.co_components];
            for k=1:length(components)
                [comp_labels, comp_index] = getLabels@Element(components{k});
                % State-labels
                switch(components{k}.type)
                    case Type.CO
                        labels.x.co = [labels.x.co; comp_labels.x];
                        index.x.co = [index.x.co; comp_index.x];
                    case Type.VO
                        labels.x.vo = [labels.x.vo; comp_labels.x];
                        index.x.vo = [index.x.vo; comp_index.x];
                end
                % Other labels
                for j=2:length(vars)
                    labels.(vars{j}) = [labels.(vars{j}); comp_labels.(vars{j})];
                    index.(vars{j}) = [index.(vars{j}); comp_index.(vars{j})];
                end
            end
        end
        
        function [ip, in] = calcStStCurrent(obj, vp, vn, u_ext)
            % CALCSTSTCURRENT calculates the steady-state bus current of the components connected to the bus, provided the bus positive and negative terminal voltages vp and vn and the steady-state
            % values of the external inputs in the vector u_ext
            
            %% Initialize bus current values
            ip = 0;
            in = 0;
            
            %% Calculate component steady-state output currents
            components = [obj.vo_components; obj.co_components];
            for c=1:length(components)
                switch components{c}.connection
                    case Connection.po
                        ip = ip + components{c}.calcStStCurrent(vp,u_ext(components{c}.getIndex('u_ext')));
                    case Connection.on
                        in = in + components{c}.calcStStCurrent(vn,u_ext(components{c}.getIndex('u_ext')));
                    case Connection.pon
                        i = components{c}.calcStStCurrent([vp vn],u_ext(components{c}.getIndex('u_ext')));
                        if length(i)==2
                            ip = ip + i(1);
                            in = in + i(2);
                        else
                            error(['The function calcStStCurrent in component ' components{c}.label ' should return the positive and negative current in vector format.']);
                        end
                end
            end
        end
        
        function [index, X] = calcStSt(obj, vp, vn, u_ext)
            % CALCSTST calculate the steady-state quantities of the components connected to the bus, provided
            % the terminal voltages vp (positive phase), vn (negative phase) and the external inputs.
            
            components = [obj.vo_components; obj.co_components];
            X = struct('vo',[],'co',[]);
            index = struct('vo',[],'co',[]);
            for c=1:length(components)
                switch components{c}.connection
                    case Connection.po
                        X_ = components{c}.calcStSt(vp,u_ext(components{c}.getIndex('u_ext')));
                    case Connection.on
                        X_ = components{c}.calcStSt(vn,u_ext(components{c}.getIndex('u_ext')));
                    case Connection.pon
                        X_ = components{c}.calcStSt([vp vn],u_ext(components{c}.getIndex('u_ext')));
                end
                if ~isempty(X_)
                    switch components{c}.type
                        case Type.CO
                            X.co = [X.co; X_];
                            index.co = [index.co; components{c}.getIndex('x')];
                        case Type.VO
                            X.vo = [X.vo; X_];
                            index.vo = [index.vo; components{c}.getIndex('x')];
                    end
                end
            end
            X.vo = transp(X.vo);
            X.co = transp(X.co);
        end
        
        function [] = validate(obj)
            % VALIDATE checks the Bus prior to solving for most common mistakes
            
            %% Check whether VO components are properly connected
            
            % Array containing the number of vo_components connected at the positive (first element) and
            % negative terminal
            no_vo_components = zeros(2,1);
            
            for c=1:length(obj.vo_components)
                switch obj.vo_components{c}.connection
                    case Connection.po
                        no_vo_components(1) = no_vo_components(1) + 1;
                    case Connection.on
                        no_vo_components(2) = no_vo_components(2) + 1;
                    case Connection.pon
                        no_vo_components(1) = no_vo_components(1) + 1;
                        no_vo_components(2) = no_vo_components(2) + 1;
                end
            end
            
            if no_vo_components(1)==0 && obj.connection~=Connection.on
                error(['No VO component is connected to the positive terminal of ' obj.label]);
            elseif no_vo_components(2)==0 && obj.connection~=Connection.po
                error(['No VO component is connected to the negative terminal of ' obj.label]);
            elseif no_vo_components(1)>1 || no_vo_components(2)>1
                error(['Too many VO components are connected to the positive and/or terminal of ' obj.label]);
            end
            
        end
    end
end

