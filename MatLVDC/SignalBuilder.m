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

classdef SignalBuilder < handle
    %SIGNALBUILDER create input signals for the simulation of DC Systems
    
    properties
        signalFcn       % The signal time functions to be used
        parameters      % The parameters of the signals
    end
    
    methods
        function obj = SignalBuilder
            obj.signalFcn = cell(0);
            obj.parameters = cell(0);
        end
        
        function [] = newSignal(obj, extInputNo, signalFcn, params)
            % INPUTS
            % extInputNo    the index of the external input where the signal will be applied
            % signalFcn     a function of the format fcn(t, parameters)
            % that returns the signal value at time t.
            % params        the signal parameters
            obj.signalFcn{extInputNo} = signalFcn;
            obj.parameters{extInputNo} = params;
        end        
        function inputs = getInputs(obj, time)
            % Return the input vector that comprises all input signals that
            % are defined in this SignalBuilder
            %inputs = zeros(length(obj.signalFcn),1);
            inputs = [];
            for k=1:length(obj.signalFcn)
                if ~isempty(obj.signalFcn{k})
                    inputs = [inputs; obj.signalFcn{k}(time, obj.parameters{k})];
                end
            end
        end
    end
    
end