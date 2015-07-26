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

classdef Configuration < uint32
    % CONFIGURATION Class containing enumeration of possible component and
    % bus configurations.
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 13, 2015
    
    enumeration
      Unipolar  (1)
      Bipolar   (3)
      default   (0)
    end
    
    methods(Static)
        function retVal = getDefaultValue()
            retVal = Configuration.default;    
        end
    end
end

