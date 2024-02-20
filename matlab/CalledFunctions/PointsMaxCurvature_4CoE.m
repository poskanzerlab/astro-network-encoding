function [ t1, t2, t3 ] = PointsMaxCurvature_4CoE( Coefficients )
%% PointsMaxCurvature Find the frames of maximum curvature for a trace fit
%with a sigmoid curve (see createFit_sigmoid20210816) as described in
%Fedchyshyn and Wang, 2007. JPhysiology 581.2(581-602)
%   Set the fourth derivative of the sigmoid equation equal to 0 and solve
%   for X. There should be three answers(t1,tmid and t3), t1 and t3 are the points of
%   maximum curvature for the sigmoid curve. The point/frame that comes
%   first (t1 or t3) is the frame at which the trace will be designated as
%   beginning to change

%% 
syms a b c d x t1 tmid t3 %create the symbolic variables a b and c and assign them the coefficients solved for in createFit_sigmoid20200221
a = Coefficients(1);
b = Coefficients(2);
c = Coefficients(3);
d = Coefficients(4);

fourthDerivative = (a*exp((b-x)/d)*(exp(3*((b-x)/d))-11*exp(2*((b-x)/d))+11*exp((b-x)/d)-1))/...
                             (exp((b-x)/d) + 1)^5 == 0;

t = solve(fourthDerivative, x);
t1 = double(t(1));
t2 = double(t(2));
t3 = double(t(3));
end

