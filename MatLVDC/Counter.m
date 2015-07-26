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

classdef Counter < handle
    %COUNTER Object that keeps track of the number of state-variables and
    %in- and outputs
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) -
    % Last update: July 13, 2015
    
    properties(Access=private)
        % Structure containing counters of VO, CO component and network
        % internal input's index
        u_int = struct('vo',0,'co',0,'network',0);
        % Structure containing counters of VO, CO component and network
        % internal output's index
        y_int = struct('vo',0,'co',0,'network',0);
        % Counter containing external input's index
        u_ext = 0;
        % Counter containing external output's index
        y_ext = 0;
        % Structure containing state-variables
        x = struct('vo',0,'co',0,'network',0);
    end
    
    methods(Access=private)
        function [index,counter] = addVariables(obj,counter,noVariables)
            % ADDVARIABLES adds the specified number of variables
            % noVariables to the counter and returns the index(es) of the
            % added variables
            
            index = counter+1:counter+noVariables;
            counter = counter + noVariables;            
        end
    end
    
    methods
        
        function obj = Counter()
            % COUNTER initialize a new counter object
            obj.u_ext = 0;
            obj.y_ext = 0;
            obj.x = struct('vo',0,'co',0,'network',0);
            obj.y_int = struct('vo',0,'co',0,'network',0);
            obj.u_int = struct('vo',0,'co',0,'network',0);
        end
        
        function [x_index, u_int_index, y_int_index, u_ext_index, y_ext_index] = getComponentIndexes(obj,component)
            % GETCOMPONENTINDEXES adds indexes to a component and increment
            % the corresponding I/O and state-variable counters
            
            % Add state variables
            if component.type==Type.VO
                [x_index, obj.x.vo] = obj.addVariables(obj.x.vo,length(component.x));
            elseif component.type==Type.CO
                [x_index, obj.x.co] = obj.addVariables(obj.x.co,length(component.x));
            end
            
            % Add external inputs
            [u_ext_index, obj.u_ext] = obj.addVariables(obj.u_ext,length(component.u_ext));
            % Add external outputs
            [y_ext_index, obj.y_ext] = obj.addVariables(obj.y_ext,length(component.y_ext));
            
            % Add internal inputs and outputs
            noIntInputs = length(component.u_int);
            noIntOutputs = length(component.y_int);
            if component.type==Type.VO
                [u_int_index, obj.u_int.vo] = obj.addVariables(obj.u_int.vo,noIntInputs);
                [y_int_index, obj.y_int.vo] = obj.addVariables(obj.y_int.vo,noIntOutputs);
            elseif component.type==Type.CO
                [u_int_index, obj.u_int.co] = obj.addVariables(obj.u_int.co,noIntInputs);
                [y_int_index, obj.y_int.co] = obj.addVariables(obj.y_int.co,noIntOutputs);
            end
        end
        
        function [x_index, u_int_index, y_int_index, y_ext_index] = getNetworkIndexes(obj,network)
            % GETNETWORKINDEXES adds indexes to a network and increment
            % the corresponding I/O and state-variable counters
            %
            % Note: a network does not contain external inputs
            
            % Add state variables
            [x_index, obj.x.network] = obj.addVariables(obj.x.network,length(network.x));
            % Add external outputs
            [y_ext_index, obj.y_ext] = obj.addVariables(obj.y_ext,length(network.y_ext));
            
            % Add internal inputs and outputs
            noIntInputs = length(network.u_int);
            noIntOutputs = length(network.y_int);

            [u_int_index, obj.u_int.network] = obj.addVariables(obj.u_int.network,noIntInputs);
            [y_int_index, obj.y_int.network] = obj.addVariables(obj.y_int.network,noIntOutputs);
        end
        
        function n = getNumberOfStates(obj,varargin)
            % GETNUMBEROFSTATES returns the current counter value of state-variables
            % corresponding to id (id=vo, co or network)
            if nargin==2
                n = obj.x.(varargin{1});
            else
                n = obj.x.network + obj.x.vo + obj.x.co;
            end
        end
        
        function n = getNumberOfExtOutputs(obj)
            % GETNUMBEROFEXTOUTPUTS returns the current counter value of external outputs
            n = obj.y_ext;
        end
        
        function n = getNumberOfExtInputs(obj)
            % GETNUMBEROFEXTINPUTS returns the current counter value of external inputs
            n = obj.u_ext;
        end
        
    end
    
end

