classdef CPL_C < Component
    % CPL_C Model of a constant power load in parallel with a capacitor
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 16, 2015
    
    methods
        function obj = CPL_C(label, connection, capacitance)
            obj.type = Type.VO;
            obj.connection = connection;
            obj.configuration = Configuration.Unipolar;
            obj.x{1} = Quantity('Voltage','A');
            obj.u_int{1} = Quantity('Current','A');
            obj.y_int{1} = Quantity('Voltage','V');
            obj.u_ext{1} = Quantity('Power','W');
            obj.y_ext = cell(0);
            obj.parameters{1} = Quantity('Capacitance','F');
            if capacitance>0
                obj.parameters{1}.value = capacitance;
            else
                error(['The capacitance value specified for ' label ' is <0']);
            end
            obj.label = label;
        end
        
        function dx = diffEqns(obj,t,x,u_int,u_ext)
            %P = u_ext(obj.u_ext{1}.index);
            P = u_ext(1);
            %v = x(obj.x{1}.index);
            v = x(1);
            %i = u_int(obj.u_int{1}.index);
            i = u_int(1);
            C = obj.parameters{1}.value;
            if v==0
                warning(['Constant power load voltage is zero, which leads to a singularity in the differential equations of component ' obj.label]);
            else
                dx = (i-P/v)/C;
            end
        end
        
        function y_int = calcIntOutputs(obj,t,x,u_ext)
            %y_int(1) = x(obj.x{1}.index);
            y_int(1) = x(1);
        end
        
        function [dxdx, dxdu_int, dy_intdx] = ddiffEqns(obj,t,x,u_int,u_ext)
            % Time-derivative values of the state-variable functions and the internal output
            % functions
            P = u_ext(1);
            v = x(1);
            i = u_int(1);
            C = obj.parameters{1}.value;
            
            dxdx(1) = P/C/v^2;
            dxdu_int(1) = 1/C;
            dy_intdx(1) = 1;
            dy_intdu_int(1) = 0;
        end
        
        function y_ext = calcExtOutputs(obj,t,x,u_int,u_ext)
            y_ext = [];
        end
        
        function i = calcStStCurrent(obj,v,u_ext)
            % INPUTS
            % v         the steady-state terminal voltage
            % OUTPUTS
            % i         the steady-state current
            %P = u_ext(obj.u_ext{1}.index);
            P = u_ext(1);
            i = P/v;
        end
        
         function X = calcStSt(obj,v,u_ext);
            X = v;
        end
    end
    
end

