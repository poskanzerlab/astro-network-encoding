function [ t1, t2, t3 ] = PointsMaxCurvature_4CoE_V2( Coefficients )
%% PointsMaxCurvature Find the frames of maximum curvature for a trace fit
%with a sigmoid curve (see createFit_sigmoid20210816) as described in
%Fedchyshyn and Wang, 2007. JPhysiology 581.2(581-602)
%   Set the fourth derivative of the sigmoid equation equal to 0 and solve
%   for X. There should be three answers(t1,t2 and t3), t1 and t3 are the points of
%   maximum curvature for the sigmoid curve. The point/frame that comes
%   first (t1 or t3) is the frame at which the trace will be designated as
%   beginning to change
%   This version of PointsMaxCurvature (_4CoE_V2) calculates the 4th
%   derivative within matlab instead of using the equation worked out
%   outside of matlab (20210928)
%% 
syms a b c d f(x) %create the symbolic variables a b and c and assign them the coefficients solved for in createFit_sigmoid20200221

f(x) = (a/(1+exp((b-x)/c)))+d;
fourthDerivative = diff(f,x,4);

t = solve(fourthDerivative == 0, x);
a = Coefficients(1);
b = Coefficients(2);
c = Coefficients(3);
d = Coefficients(4);
t_eval = subs(t);

t1 = double(t_eval(1));
t2 = double(t_eval(2));
t3 = double(t_eval(3));
end