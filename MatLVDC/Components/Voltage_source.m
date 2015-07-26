classdef Voltage_source < Component
    % VOLTAGE_SOURCE Model of a unipolar voltage source
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 14, 2015
    
    methods
        function obj = Voltage_source(label, connection)
            obj.type = Type.VO;
            obj.connection = connection;
            obj.configuration = Configuration.Unipolar;
            obj.x = cell(0);
            obj.u_int{1} = Quantity('Current','A');
            obj.y_int{1} = Quantity('Voltage','V');
            obj.u_ext{1} = Quantity('Voltage','V');
            obj.y_ext{1} = Quantity('Voltage','V');
            obj.parameters = cell(0);
            obj.label = label;
        end
        
        function dx = diffEqns(obj,t,x,u_int,u_ext)
            dx = [];
        end
        
        function y_int = calcIntOutputs(obj,t,x,u_ext)
            %y_int(1) = u_ext(obj.u_ext{1}.index);
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
            % v = v
            %y_ext(1) = u_ext(obj.u_ext{1}.index);
            y_ext(1) = u_ext(1);
            % P = v * i
            % y_ext(2) = u_ext(obj.u_ext{1}.index) * u_int(obj.u_int{1}.index);
        end
        
        function i = calcStStCurrent(obj,v,u_ext)
            error(['Cannot calculate the steady-state current provided the bus voltage for Voltage source ' obj.label]);
        end
        
        function X = calcStSt(obj,v,u_ext);
            X = [];
        end
    end
    
end

