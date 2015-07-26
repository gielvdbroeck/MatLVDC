classdef Buck_R < Component
    % BUCK_R Model of a buck converter (averaged model) driving a resistive load, of which the output
    % voltage is regulated by a cascaded control structure
    %
    % Created by Giel Van den Broeck (giel.vandenbroeck@energyville.be) - Last update: July 23, 2015
    
    methods
        function obj = Buck_R(label,connection,parameters)
            %% Component properties
            obj.type = Type.CO;
            obj.connection = connection;
            obj.configuration = Configuration.Unipolar;
            
            %% States
            %   x(1) = inductor current [A]
            %   x(2) = capacitor voltage [V]
            %   x(3) = current controller state
            %   x(4) = voltage controller state
            obj.x{1} = Quantity('inductor current','A');
            obj.x{2} = Quantity('capacitor voltage','V');
            obj.x{3} = Quantity('current controller state','');
            obj.x{4} = Quantity('voltage controller state','');
            
            %% Internal inputs
            %   u_int(1) = input voltage [V]
            obj.u_int{1} = Quantity('Voltage','V');
            
            %% Internal outputs
            %   y_int(1) = input current [A]
            obj.y_int{1} = Quantity('Current','A');
            
            %% External inputs
            %   u_ext(1) = load voltage reference [V]
            obj.u_ext{1} = Quantity('Voltage','V');
            
            %% External outputs
            obj.y_ext{1} = Quantity('Duty cycle','');
            obj.y_ext{2} = Quantity('Current setpoint','A');
            
            %% Parameters
            %   p(1) = R load [Ohm]
            %   p(2) = C [F]
            %   p(3) = L [H]
            %   p(4) = Kpi proportional current controller parameter
            %   p(5) = Kii integral current controller parameter
            %   p(6) = Kpv proportional voltage controller parameter
            %   p(7) = Kiv integral voltage controller parameter
            obj.parameters{1} = Quantity('R load','\Omega');
            obj.parameters{2} = Quantity('C','F');
            obj.parameters{3} = Quantity('L','H');
            obj.parameters{4} = Quantity('Kpi proportional current controller parameter','');
            obj.parameters{5} = Quantity('Kii integral current controller parameter','');
            obj.parameters{6} = Quantity('Kpv proportional voltage controller parameter','');
            obj.parameters{7} = Quantity('Kiv integral voltage controller parameter','');
            obj.parameters{8} = Quantity('Inductor current limit','');
            
            if length(parameters)~=length(obj.parameters)
                error(['Invalid parameters provided for' label]);
            else
                for k=1:length(parameters)
                    obj.parameters{k}.value = parameters(k);
                end
            end
            
            %% Component label
            obj.label = label;
        end
        
        function dx = diffEqns(obj,t,x,u_int,u_ext)
            
            %% Parameter declarations
            R = obj.parameters{1}.value;
            C = obj.parameters{2}.value;
            L = obj.parameters{3}.value;
            Kpi = obj.parameters{4}.value;
            Kii = obj.parameters{5}.value;
            Kpv = obj.parameters{6}.value;
            Kiv = obj.parameters{7}.value;
            i_max = obj.parameters{8}.value;

            %% Intermediate signals
            % Reference current
            i_ref = min(Kpv*(u_ext(1)-x(2)) + x(4),i_max);
            % Duty cycle (limited between 0 and 1)
            d = saturate(Kpi*(i_ref-x(1)) + x(3), 0, 1); 

            %% State equations
            dx(1) = 1/L * (d*u_int(1) - x(2));
            dx(2) = 1/C * (x(1) - x(2)/R);
            dx(3) = Kii * (i_ref - x(1));
            dx(4) = Kiv * (u_ext(1) - x(2));
        end
        
        function y_int = calcIntOutputs(obj,t,x,u_ext)
            Kpi = obj.parameters{4}.value;
            Kpv = obj.parameters{6}.value;
            i_max = obj.parameters{8}.value;
            i_ref = min(Kpv*(u_ext(1)-x(2)) + x(4),i_max);
            d = saturate(Kpi*(i_ref-x(1)) + x(3), 0, 1);
            y_int(1) = x(1)*d;
        end
        
        function [dxdx, dxdu_int, dy_intdx] = ddiffEqns(obj,t,x,u_int,u_ext)
            % Time-derivative values of the state-variable functions and the internal output
            % functions
            
            R = obj.parameters{1}.value;
            C = obj.parameters{2}.value;
            L = obj.parameters{3}.value;
            Kpi = obj.parameters{4}.value;
            Kii = obj.parameters{5}.value;
            Kpv = obj.parameters{6}.value;
            Kiv = obj.parameters{7}.value;
            i_max = obj.parameters{8}.value;
            
            dxdx = zeros(4);
            dxdx(1,1) = -Kpi*u_int(1)/L;
            dxdx(1,2) = -1/L-u_int(1)*Kpi*Kpv/L;
            dxdx(1,3) = u_int(1)/L;
            dxdx(1,4) = Kpi*u_int(1)/L;
            
            dxdx(2,1) = 1/C;
            dxdx(2,2) = -R/C;
            
            dxdx(3,1) = -Kii;
            dxdx(3,2) = -Kpv * Kii;
            dxdx(3,4) = Kii;
            
            dxdx(4,2) = -Kiv;
            
            dxdu_int = zeros(4,1);
            i_ref = min(Kpv*(u_ext(1)-x(2)) + x(4),i_max);
            d = Kpi*(i_ref-x(1)) + x(3);
            dxdu_int(1) = d/L;
            
            dy_intdx = zeros(1,4);
            dy_intdx(1) = ( Kpi * ( Kpv * (u_ext(1)-x(2)) + x(4) - x(1) ) + x(3) ) * u_int(1) - x(1) * Kpi * u_int(1);
            dy_intdx(2) = -x(1) * u_int(1) * Kpv * Kpi;
            dy_intdx(3) = x(1) * u_int(1);
            dy_intdx(4) = x(1) * u_int(1) * Kpi;
        end
        
        function y_ext = calcExtOutputs(obj,t,x,u_int,u_ext)
            Kpi = obj.parameters{4}.value;
            Kpv = obj.parameters{6}.value;
            i_max = obj.parameters{8}.value;
            i_ref = min(Kpv*(u_ext(1)-x(2)) + x(4),i_max);
            d = saturate(Kpi*(i_ref-x(1)) + x(3), 0, 1);
            y_ext(1) = d;
            y_ext(2) = i_ref;
        end
        
        function i = calcStStCurrent(obj,v,u_ext)
            R = obj.parameters{1}.value;
            v_ref = u_ext(1);
            i = v_ref^2/R/v;
        end
        
        function X = calcStSt(obj,v,u_ext);
            R = obj.parameters{1}.value;
            X(2) = u_ext(1);
            X(1) = X(2)/R;
            X(3) = X(2)/v;
            X(4) = X(1);
        end
    end
    
end

