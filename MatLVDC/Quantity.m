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

classdef Quantity
    % QUANTITY Type definition of a quantity
    %   States, inputs and outputs are quantities in the DC system.
    %   Each quantity has an index, a unit and a label.
    %
    %   Created by: Giel Van den Broeck - Last update: July 22, 2015
    
    properties
        % Index of the quantity
        index = 0;
        % Unit of the quantity
        unit;
        % Label of the quantity
        label;
        % Value of the quantity
        value;
    end
    
    methods(Static)
        function index = getIndex(quantities)
            % GETINDEX returns the indexes of the quantities
            %
            % INPUTS
            % quantities        cell array of quantities of which the indexes should be returned
            %
            % OUTPUT
            % index             array of indexes corresponding to the input quantities
            
            index = zeros(length(quantities),1);
            for k=1:length(quantities)
                index(k) = quantities{k}.index;
            end
        end
        
        function index = getValues(quantities)
            % GETVALUES returns the values of the quantities as a vector
            %
            % INPUTS
            % quantities        cell array of quantities of which the values should be returned
            %
            % OUTPUT
            % index             vector of values corresponding to the input quantities
            
            index = zeros(length(quantities),1);
            for k=1:length(quantities)
                index(k) = quantities{k}.value;
            end
        end
    end
    
    methods
        function obj = Quantity(label, unit)
            % QUANTITY create a new quantity with the specified label and unit
            obj.index = [];
            obj.unit = unit;
            obj.label = label;
        end
        
        function ind = get.index(obj)
            ind = obj.index;
        end
    end
    
end

