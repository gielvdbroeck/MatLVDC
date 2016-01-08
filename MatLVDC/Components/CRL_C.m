classdef CRL_C < Component
    % CRL_C Model of a constant resistance load in parallel with a capacitor
    %
    % Created by Shutang You (syou3@vols.utk.edu) - Last update: Jan 3, 2016
    
    methods
        function obj = CRL_C(label, connection, capacitance)
            obj.type = Type.VO;  %
            obj.connection = connection;
            obj.configuration = Configuration.Unipolar;
            obj.x{1} = Quantity('Voltage','V');  % the voltage over the load
            obj.u_int{1} = Quantity('Current','A'); % the current over the resistive and capacitive load
            obj.y_int{1} = Quantity('Voltage','V'); % the voltage over the load
            obj.u_ext{1} = Quantity('Resistance','Ω'); % load resistance value
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
            R = u_ext(1);  % the resistance of the load
            v = x(1);   % the voltage over the load
            i = u_int(1);  % the current over the resistive and capacitve load
            C = obj.parameters{1}.value;
            dx = (i-v/R)/C; % voltage deviation over time
        end
        
        function y_int = calcIntOutputs(obj,t,x,u_ext)
            y_int(1) = x(1);  % the voltage over the load
        end
        
        function [dxdx, dxdu_int, dy_intdx] = ddiffEqns(obj,t,x,u_int,u_ext)
            % Time-derivative values of the state-variable functions and the internal output
            % functions
            R = u_ext(1);  % the resistance of the load
            v = x(1);   % the voltage over the load
            i = u_int(1); % the current over the resistive and capacitve load
            C = obj.parameters{1}.value; 
            
            dxdx(1) = -1/(R*C);   
            dxdu_int(1) = 1/C;
            dy_intdx(1) = 1;
        end
        
        function y_ext = calcExtOutputs(obj,t,x,u_int,u_ext)
            y_ext = [];
        end
        
        function i = calcStStCurrent(obj,v,u_ext)
            % INPUTS
            % v         the steady-state terminal voltage
            % OUTPUTS
            % i         the steady-state current
            %R = u_ext(obj.u_ext{1}.index);
            R = u_ext(1);  
            i = v/R;
        end
        
        function X = calcStSt(obj,v,u_ext);
            X = v;
        end
        end
end
