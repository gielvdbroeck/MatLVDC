function [ val ] = step_signal(time, parameters)
%STEP_SIGNAL represents a step signal to be applied to the DC System object
%via the SignalBuilder class.
%
% INPUTS
% time          the time when the step_signal should be evaluated
% parameters    array containing the parameters of the step signal
%       parameters(1) = time when the step takes place
%       parameters(2) = initial output value
%       parameters(3) = final output value

if time<parameters(1)
    val = parameters(2);
else
    val = parameters(3);
end

end

