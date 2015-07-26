classdef Capacitor < Component
    % CAPACITOR Model of a capacitor
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 16, 2015
    
    methods
        function obj = Capacitor(label, connection, capacitance)
            obj.type = Type.VO;
            if connection==Connection.pon
                error([label ' cannot be connected as bipolar component']);
            else
                obj.connection = connection;
            end
            obj.configuration = Configuration.Unipolar;
            
            obj.x{1} = Quantity('Voltage','V');
            obj.u_int{1} = Quantity('Current','A');
            obj.y_int{1} = Quantity('Voltage','V');
            obj.u_ext = cell(0);
            obj.y_ext = cell(0);
            
            obj.parameters{1} = Quantity('Capacitance','C');
            
            if capacitance>0
                obj.parameters{1}.value = capacitance;
            else
                error('Invalid capacitance value');
            end
            
            obj.label = label;
        end
        
        function dx = diffEqns(obj,t,x,u_int,u_ext)
            dx = u_int(1)/obj.parameters{1}.value;
        end
        
        function y_int = calcIntOutputs(obj,t,x,u_ext)
            y_int(1) = x(1);
        end
        
        function [dxdx, dxdu_int, dy_intdx, dy_intdu_int] = ddiffEqns(obj,t,x,u_int,u_ext)
            % Time-derivative values of the state-variable functions and the internal output
            % functions
            dxdx(1) = 0;
            dxdu_int(1) = 1;
            dy_intdx(1) = 1;
            dy_intdu_int(1) = 0;
        end
        
        function y_ext = calcExtOutputs(obj,t,x,u_int,u_ext)
            y_ext = [];
        end
        
        function i = calcStStCurrent(obj,v,u_ext)
            % INPUTS
            % v         the steady-state voltage
            % OUTPUTS
            % i         the steady-state current
            i = 0;
        end
        
        function X = calcStSt(obj,v,u_ext);
            X(1) = v;
        end
    end
    
end

