classdef Current_source < Component
    % CURRENT_SOURCE Model of a unipolar current source
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 20, 2015
    
    methods
        function obj = Current_source(label,connection)
            obj.type = Type.CO;
            obj.connection = connection;
            obj.configuration = Configuration.Unipolar;
            obj.x = cell(0);
            obj.u_int{1} = Quantity('Voltage','V');
            obj.y_int{1} = Quantity('Current','A');
            obj.u_ext{1} = Quantity('Current','A');
            obj.y_ext = cell(0);
            obj.parameters = cell(0);
            obj.label = label;
        end
        
        function dx = diffEqns(obj,t,x,u_int,u_ext)
            dx = [];
        end
        
        function y_int = calcIntOutputs(obj,t,x,u_ext)
            y_int(1) = u_ext(1);
        end
        
        function [dxdx, dxdu_int, dy_intdx, dy_intdu_int] = ddiffEqns(obj,t,x,u_int,u_ext)
            % Time-derivative values of the state-variable functions and the internal output
            % functions
            dxdx = [];
            dxdu_int = [];
            dy_intdx(1) = 0;
            dy_intdu_int(1) = 0;
        end
        
        function y_ext = calcExtOutputs(obj,t,x,u_int,u_ext)
            y_ext = [];
        end
        
        function i = calcStStCurrent(obj,v,u_ext)
            i = u_ext(1);
        end
        
        function X = calcStSt(obj,v,u_ext);
            X = [];
        end
    end
    
end

