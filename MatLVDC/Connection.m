classdef Connection < uint32
    % CONNECTION Class containing enumeration of possible component and
    % bus connections.
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 13, 2015
    
    enumeration
        % Unipolar connection between positive and neutral phase
        po  (1)
        % Unipolar connection between neutral and negative phase
        on  (2)
        % Bipolar connection between positive, neutral and negative phase
        pon (3)
        % Default value
        default (0)
    end
    
    methods(Static)
        function retVal = getDefaultValue()
            retVal = Connection.default;    
        end
    end
    
end

