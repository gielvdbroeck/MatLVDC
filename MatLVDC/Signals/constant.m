function [ val ] = constant(~, parameters)
%CONSTANT represents a constant signal to be applied to the DC System object
%via the SignalBuilder class.
%
% INPUTS
% time          the time when the step_signal should be evaluated
% parameters    array containing the parameters of the step signal
%       paramters(1) = constant signal value

val = parameters(1);

end

