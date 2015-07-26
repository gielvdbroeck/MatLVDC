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

classdef Element < handle
    % ELEMENT The parent class of all elements constituting a DC system
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 13, 2015
    
    properties
        % Parameter indexes of the element
        parameters;
        % State indexes of the element
        x;        
        % Internal input indexes of the element
        u_int;
        % External input indexes of the element
        u_ext;
        % Internal output indexes of the element
        y_int;
        % External output indexes of the element
        y_ext;
        % Element label
        label;
    end
    
    methods (Abstract)
        % Return the time-derivative values of the state-variables and their index
        diffEqns(obj);
        % Return the internal input values and their index
        calcIntOutputs(obj);
        % Return the external output values and their index
        calcExtOutputs(obj);
    end
    
    methods
        function [labels, index] = getLabels(obj)
            % GETLABELS returns the labels of the Element's quantities (i.e. parameters, states, internal
            % inputs, external inputs, internal outputs and external outputs)
            
            labels = []; index = [];
            
            %% States
            vars{1} = 'x';
            
            %% Parameters
            vars{2} = 'parameters';
            
            %% External inputs
            vars{3} = 'u_ext';
            
            %% External outputs
            vars{4} = 'y_ext';
            
            %% Get the labels
            for k=1:length(vars)
                labels_ = cell(length(obj.(vars{k})),1);
                index_ = zeros(length(obj.(vars{k})),1);
                for j=1:length(obj.(vars{k}))
                    labels_{j} = [obj.label ' - ' obj.(vars{k}){j}.label];
                    if isempty(obj.(vars{k}){j}.index)
                        index_(j) = 0;
                    else
                        index_(j) = obj.(vars{k}){j}.index;
                    end
                end
                labels.(vars{k}) = labels_;
                index.(vars{k}) = index_;
            end
        end
    end
end

