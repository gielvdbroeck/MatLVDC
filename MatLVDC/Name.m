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

classdef Name
    % NAME Convert between human interpretable names and the application's numbering conventions
    %
    % e.g. 'unipolar pos' refers to a bus or branch configuration where only a positive and neutral
    % phase conductor exist. The negative phase conductor is not present. Similarly 'unipolar neg' is
    % defined. In MatLVDC a numbering convention is used to indicate this bus or branch configuration.
    % The STATIC class Name includes functions to convert between the names which are easier to interpret
    % and the numbering convention used in the software.
    %
    % Different aliases can as well be programmed
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be)
    % Last modification: June 14, 2015
    
    methods(Static)
        function configuration_number = configuration(configuration_name)
            % CONFIGURATION Return the configuration number that corresponds to the configuration name
            % Three different configurations exist:
            %       unipolar positive:  bus or branch contains a positive and neutral terminal
            %       unipolar negative:  bus or branch contains a neutral and negative terminal
            %       bipolar:            bus or branch contains a positive, neutral and negative terminal
            %
            % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be)
            % Last modification: June 14, 2015
            switch configuration_name
                case {'unipolar plus', 'unipolar_plus', 'unipolar pos', 'unipolar_pos'}
                    configuration_number = 1;
                case {'unipolar min', 'unipolar_min', 'unipolar neg', 'unipolar_neg'}
                    configuration_number = 2;
                case 'bipolar'
                    configuration_number = 3;
                otherwise
                    error(ErrorWarning.getErrorWarning(1));
            end
        end
        
        function [ colNo ] = getBusField( name )
            % GETBUSFIELD Return the corresponding column number that corresponds to the
            % column header that is given as an input argument
            %
            % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be)
            % Last modification: June 14, 2015
            switch name
                case 'BUSDC_I'
                    colNo = 1;
                case 'BUSDC_CONFIG'
                    colNo = 2;
                otherwise
                    error(['The bus field ' name 'does not exist']);
            end

        end
        
        function [ colNo ] = getBranchField(colHeader)
            % GETBRANCHFIELD Return the column number that corresponds to the column header of the branch
            % matrix

            switch colHeader
                case 'F_BUS'
                    colNo = 1;
                case 'T_BUS'
                    colNo = 2;
                case 'BR_CONFIG'
                    colNo = 3;
                case 'BR_R1'
                    colNo = 4;
                case 'BR_R2'
                    colNo = 5;
                case 'BR_Rn'
                    colNo = 6;
                case 'BR_L1'
                    colNo = 7;
                case 'BR_L2'
                    colNo = 8;
                case 'BR_Ln'
                    colNo = 9;
                case 'BR_C1'
                    colNo = 10;
                case 'BR_C2'
                    colNo = 11;
                otherwise
                    error(['The branch field ' colHeader 'does not exist']);
            end

        end
    end
    
end

