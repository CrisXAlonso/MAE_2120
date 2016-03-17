function [sigX] = Homework_5(f1, f2, d1, d2, l, dO, radO, radI)
%This function solves for the stress in the axial direction of a hollow
%cylinder under 4 point bending. The force 'f1' is the first force from the
%left of the cylinder acting a distance 'd1; from the left edge. the force 
%'f2' is the second force from the left of the cylinder acting a distance 
%'d2' from the left edge. The total length of the cynlinder is 'l' and the
%distance from the left edge of the point of examination is 'dO'. radO is
%the outer radius of the cylinder in centimeters and radI is the inner
%radius of the cylinder in centimeters. f1 and f2 must be provided in
%newtons and d1, d2, l, and dO in meters. The output sigX is in MPa.
%!!!!!!!!!!!!!!The point O must be between d1 and d2 !!!!!!!!!!!!!!!!!!!!!

if ~exist('f1','var') || isempty(f1)
    f1 = -200; % in newtons
end
if ~exist('f2','var') || isempty(f2)
    f2 = -200; % in newtons
end
if ~exist('d1','var') || isempty(d1)
    d1 = 3; % in meters
end
if ~exist('d2','var') || isempty(d2)
    d2 = 9; % in meters
end
if ~exist('l','var') || isempty(l)
    l = 12; % in meters
end
if ~exist('dO','var') || isempty(dO)
    dO = 6; % in meters
end
if ~exist('radO','var') || isempty(radO)
    %Creates a span of evenly spaced outer radii values
    radO = linspace(1.5, 2.5, 21); % in centimeters
end
if ~exist('radI','var') || isempty(radI)
    radI = 1; % in centimeters
end

if (dO < d1) || (dO > d2)
    error('The point O must be between d1 and d2!');
end

%Creates matrices to solve for the reaction forces.
A = [1 1; 0 l];
B = [-(f1+f2); -(d1*f1+d2*f2)];
%Solves for the reaction forces.
R = A\B;

radO = radO/100; %Converting to meters
radI = radI/100; %Converting to meters

%Function that calculates the stress along the axial direction due to the
%moment of bending using the equation sigma = My/I
    function sigX = bending_stress2(M, ri, ro)
        IVal = pi/4*(ro.^4-ri.^4);
        sigX = 1e-6*M.*ro./IVal;
    end

%Calculating the moment at point O
M = R(1,1)*dO + f1*(dO-d1);
%Calculating the stress at point 0
sigX = bending_stress2(M, radI, radO);

%Plot sigma_x vs. ri
figure(1); %Create figure #1
clf; %Clear current figure
h=plot(radO, sigX, 'g');
set(gca,'Box','on','LineWidth',2,...
'FontName','Helvetica',...
'FontSize',14);
%Set axis properties
set(h,'LineWidth',2);
%Set linewidth for curve
xlabel('r_o (m)');
ylabel('\sigma_x at O (MPa)');
title('Bending stress at Point O vs Outer Radius');
axis square; %Make axis box square


end