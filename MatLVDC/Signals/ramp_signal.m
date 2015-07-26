function [ val ] = ramp_signal(time, parameters)
%RAMP_SIGNAL represents a ramp signal to be applied to the DC System object
%via the SignalBuilder class.
%
% INPUTS
% time          the time when the step_signal should be evaluated
% parameters    array containing the parameters of the step signal
%       parameters(1) = time when the ramp starts
%       parameters(2) = time when the ramp settles
%       parameters(3) = initial output value
%       parameters(4) = final output value

if time<parameters(1)
    val = parameters(3);
elseif time>parameters(2)
    val = parameters(4);
else
    val = (parameters(4)-parameters(3))/(parameters(2)-parameters(1))*(time-parameters(1))+parameters(3);
end

end

