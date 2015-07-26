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

classdef Component < Element
    % COMPONENT The class COMPONENT represents a component in the DC network
    % 
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be)
    % - Last modification: July 14, 2015
    
    properties
        % Component configuration
        configuration = Configuration.getDefaultValue;
        % Component connection
        connection = Connection.getDefaultValue;
        % Component type
        type;
    end
    
    methods (Abstract)
        % Return the time-derivative values of the state-variables
        diffEqns(obj);
        % Return the second time-derivative values of the state-variables and internal inputs
        ddiffEqns(obj);
        % Return the internal input values
        calcIntOutputs(obj);
        % Return the external output values
        calcExtOutputs(obj);
        % Return the steady-state input current when providing the steady-state bus voltage (when the component is bipolar,
        % the function should return a vector containing two elements: the first element equals the positive
        % bus terminal current and the second element the negative bus terminal current)
        calcStStCurrent(obj);
        % Calculates the steady-state vector of the component, provided the voltage and external inputs.
        % When the component is bipolar, voltage should be a vector of which the first element contains
        % the positive terminal voltage and the second element the negative terminal voltage.
        % (when the component is bipolar,
        % the voltage vector v contains two elements: the first element equals the positive
        % bus terminal voltage and the second element the negative bus terminal voltage)
        calcStSt(obj,v,u_ext);
    end
    
    methods
        function [index] = getIndex(obj,id)
            % GETINDEX returns the index of the internal in- and outputs, as
            % specified by id.
            %
            % id can be either:
            % u_int         internal inputs
            % y_int         internal outputs
            % x             states
            % u_ext         external inputs
            % y_ext         external outputs
            %
            % INPUTS
            % id        (as specified above)
            %
            % OUTPUTS
            % index   indexes
            
            var = obj.(id);
            
            if ~isempty(var)
                index = zeros(length(var),1);
                for k=1:length(var)
                    index(k) = var{k}.index;
                end
            else
                index = [];
            end
        end
    end
    
end

