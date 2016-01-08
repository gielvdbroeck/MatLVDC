function [ val ] = rectangular_pulse(time, parameters)
%RECTANGULAR_PULSE represents a rectangular pulse disturbance signal to be applied to the DC System object
%via the SignalBuilder class.
%
% INPUTS
% time          the time when the step_signal should be evaluated
% parameters    array containing the parameters of the step signal
%       parameters(1) = time when the pulse starts
%       parameters(2) = time when the pulse stops
%       parameters(3) = initial output value
%       parameters(4) = pulse output value

if time<parameters(1)||time>parameters(2)
    val = parameters(3);
else
    val = parameters(4);
end
end
