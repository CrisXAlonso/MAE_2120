function [] = Alonso_cxa2_Project1()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = 9.81; %Settting gravitatinal constant in m/s^2
rhoAir = 1.23; %Constant density of air in kg/m^3
phiAll = 5; %Allowable angle of twist in degrees
fileNameVec = {'Alonso_cxa2_parameters.txt', 'Alonso_cxa2_geometry.txt', ...
    'Alonso_cxa2_material.txt'}; %File names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that takes in a file name string and outputs the colum vector
    %with the floating point data in those text files
    function dataVec = readData(fileName)
        fileID = fopen(fileName); dataVec = fscanf(fileID, '%f'); ...
            fclose(fileID);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that solves for reaction forces on a mean with 2 supports and
    %a single force. Fa is the magnitude of the force with positive being
    %upward, la is the distance from the left end of the rod at which the
    %force is applied, lb is the distance from the left end of the rod
    %where the first support is found, and likewise of lc but for the
    %second support
    function reac = reacSolve(Fa,la, lb, lc)
        reac = zeros(2,1); 
        reac(2,1) = Fa*(lb-la)/(lc-lb); %Derived equations from force and
        reac(1,1) = -Fa-reac(2,1);      %Moment balance
                                       
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solving Shear %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that solves shear force at a particular x location on a beam
    %with 3 forces, where Fa is the leftmost force at a distance la from
    %the left edge, Fb is the next force from the left at a distance lb
    %from the left edge, and lc is the distance from the left edge where
    %the rightmost force is found.
    function V = shearSolve(Fa, la, Fb, lb, lc, x)
        V = zeros(1, length(x)); %Preallocation
        for i = 1:length(x)
            if (x(i) >= 0) && (x(i)<la) %No forces contributing to shear
                V(i) = 0;
            elseif (x(i) >= la) && (x(i) < lb) %One force cont. to shear
                V(i) = Fa;
            elseif (x(i) >= lb) && (x(i) < lc) %Two forces cont. to shear
                V(i) = Fa + Fb; 
            else  %Other end of beam, no forces contributing to shear
                V(i) = 0; 
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solving Moments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that solves the bending momnet at a particular x location on
    %a beam with 3 forces, where Fa is the leftmost force at a distance la 
    %from the left edge, Fb is the next force from the left at a distance lb
    %from the left edge, and lc is the distance from the left edge where
    %the rightmost force is found.
    function M = momentSolve(Fa, la, Fb, lb, lc, x)
        M = zeros(1, length(x)); %Preallocation
        for i = 1:length(x) %Iterating over the length of the rod
            if (x(i) >= 0) && (x(i)<la)
                M(i) = 0; %No forces contributing to moment
            elseif (x(i) >= la) && (x(i) < lb)
                M(i) = Fa*(x(i)-la); %One force contributing to moment
            elseif (x(i) >= lb) && (x(i) < lc)
                M(i) = Fa*(x(i)-la) + Fb*(x(i)-lb);  %Both cont. to moment
            else 
                M(i) = 0; %Other end of the beam, no forces cont. to moment
            end; 
        end; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Max Effective Stress without Shear %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that solves for the max effective stress in a rod ignoring
    %transverse shear, where x is the vector of x values, MzVec is the
    %vector of Mz values, MyVec is the vector of My values, ra is the
    %radius of the leftmost portion of the rod, rb is the radius of the
    %rightmost portion of the rod, la is the length at which the radius
    %change occurs measured from the left end, and T is the torque on that
    %portion of the rod.
    function sigmaH = effStressSolve(x, MzVec, MyVec, ra, rb, la, T)
        sigmaH = zeros(1,length(x)); %Preallocating
        for i = 1:length(x) %Iterating over length of beam
            xVal = x(i); %Setting current x Value
            if xVal < la; rVal = ra; else rVal = rb; end %Picking radius
            %Finding I J and Mz and My values
            I = (pi/4)*rVal^4; J = 2*I; Mz = MzVec(i); My = MyVec(i);
            %Finding total moment and normal stress along beam and 
            %shear due to torsion
            Mtot = sqrt(Mz^2+My^2); sigmaX = Mtot*rVal/I; tau = T*rVal/J;
            %Finding radius of mohrs circle
            rMohrs = sqrt((sigmaX/2)^2+tau^2);
            %Finding principle stresses
            sig1 = sigmaX/2+rMohrs; sig3 = sigmaX/2-rMohrs;
            %Calculating Von Mises stress
            sigmaH(i) = (1/sqrt(2))*sqrt(sig1^2+sig3^2+(sig3-sig1)^2);
        end   
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Max Effective Stress with Shear %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that solves for the max effective stress at K in a rod NOT 
    %ignoring transverse shear: x is the vector of x values, MzVec is the
    %vector of Mz values, MyVec is the vector of My values, ra is the
    %radius of the leftmost portion of the rod, rb is the radius of the
    %rightmost portion of the rod, la is the length at which the radius
    %change occurs measured from the left end, and T is the torque on that
    %portion of the rod.
    function sigmaH = effStressShearSolve(x, MzVec, MyVec, ...
            VzVec, VyVec, ra, rb, la, T)
        sigmaH = zeros(1,length(x)); %Preallocation
        for i = 1:length(x) %Iterating over length of beam
            xVal = x(i); %Setting x value
            %Picking approprate radius
            if xVal < la; rVal = ra; else rVal = rb; end
            %Calculating I, J, Mz and My values
            I = (pi/4)*rVal^4; J = 2*I; Mz = MzVec(i); My = MyVec(i);
            %Calculating total moment, Vz and Vy
            Mtot = sqrt(Mz^2+My^2); Vz = VzVec(i); Vy = VyVec(i);
            %Calculating total shear
            Vtot = sqrt(Vz^2+Vy^2);
            %Checking for appropriate sign of total shear
            if (sign(Vy) + sign(Vz) == 0) || (sign(Vy) + sign(Vz) == -1)
                Vtot = -Vtot;
            end
            %Finding angle from y axies of moment and shear vectors and
            %calculating the angle of difference beta
            thM = atan2(Mz,My); thV = atan2(Vz,Vy); beta = thM-thV;
            %Calculating normal stress along x and shear due to torsion
            sigmaX = -Mtot*rVal*abs(cos(beta))/I; tau = T*rVal/J;
            %Calculating transverse shear
            shearV = (2/3)*Vtot*(rVal^2)/J;
            %Calculating total shear
            totShear = shearV+tau;
            %Calculating Von Mises stress
            sigmaH(i) = (1/sqrt(2))*sqrt(2*sigmaX^2+6*(totShear^2));
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Checking Shear Assumption%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that checks that shear assumption made for two effective
    %stress measurements, one with transverse shear and one without, is
    %valid to ignore the transverse shear effects: noShear is the effective
    %stress vector of the calculation made without shear, yesShear is the
    %effective stress calcultion made with shear, and msg is the error
    %message to display if transverse shear cannot be overlooked.
    function [] = checkShear(noShear, yesShear, msg)
        once = false; %Sets a bool to ensure the msg only gets displayed once
        for i = 1:length(noShear) %Iterating over the beam
           %Checks if the stress value with shear is greater than the one
           %without and if so diplays the warning
           if (yesShear(i) > noShear(i)) && ~once
               warning(msg);
               once = true;
           end; 
        end; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Tau Max at H and K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that computes the maximum shear stress at points k and h along
%the beam where x is the vector of x values, MzVec is the vector of Mz
%values, MyVec is the vector of My values, VzVec is the vector of shear
%force in z, VyVec is the vector of shear force in z, ra is the radius of
%the leftmost portion of the beam, rb is the radius of the rightmost
%portion of the beam, la is the distance from the left edge of the beam at
%which this change in radius occurs, and T is the torque on the beam
function [tauH, tauK] = tauSolve(x, MzVec, MyVec, ...
            VzVec, VyVec, ra, rb, la, T)
        tauH = zeros(1,length(x)); %Preallocation
        tauK = zeros(1,length(x)); %Preallocation
        for i = 1:length(x)
            xVal = x(i); %Sets x value
            if xVal < la %Sets correct radius value
                rVal = ra; 
            else
                rVal = rb;
            end
            %Calculates I, J, Mz, and My values
            I = (pi/4)*rVal^4; J = 2*I; Mz = MzVec(i); My = MyVec(i);
            %Calculates total moment Vz and Vy
            Mtot = sqrt(Mz^2+My^2); Vz = VzVec(i); Vy = VyVec(i);
            Vtot = sqrt(Vz^2+Vy^2); %Calculates total shear force
            %Sets the proper sign on the shear forces
            if (sign(Vy) + sign(Vz) == 0) || (sign(Vy) + sign(Vz) == -1)
                Vtot = -Vtot;
            end
            %Finding angle from y axies of moment and shear vectors and
            %calculating the angle of difference beta 
            thM = atan2(Mz,My); thV = atan2(Vz,Vy); beta = thM-thV;
            %Calculating shear at point H
            tauH(i) = rVal/J*sqrt(Mtot^2+T^2);
            %Calculating shear at point K
            tauK(i) = rVal/J*sqrt((Mtot*cos(beta))^2+((2/3)*rVal*Vtot+T)^2);
        end; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Angle of Twist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function that calculated the angle of twist in degrees of a rod
    %subject to torsional loading: where x is the vector of x values along
    %the rod, T is the torsional load, E is the Young's modulus, v is the
    %poisson's ratio, ra is the radius of the leftmost section of the rod,
    %rb is the radius of the rightmost section of the rod, la is the
    %distance from the left of the rod at which this change in radius
    %occurs
    function phi = twist(x, T, E, v, ra, rb, la)
        %Preallocation, setting E to Pa from GPa, and creating saveVal
        phi = zeros(1, length(x)); E = E*10^9; saveVal = 0;
        for i = 1:length(x) %Iterating over length of beam
            xVal = x(i); %Setting x value
            G = E/(2*(1+v)); %Finding shear modulus
            %If on first section of beam, calculates the coresponding I, J,
            %torsional shear, shear strain, and angle of twist values
            if xVal < la
                rVal = ra; I = (pi/4)*rVal^4; J = 2*I; tau = T*rVal/J;
                gamma = tau/G; phi(i) = gamma*xVal/rVal; saveVal = phi(i);
            else
            %If on first section of beam, calculates the coresponding I, J,
            %torsional shear, shear strain, and angle of twist values
            %referenced from where the previous one left off
                rVal = rb; I = (pi/4)*rVal^4; J = 2*I; tau = T*rVal/J;
                gamma = tau/G; phi(i) = gamma*(xVal-la)/rVal + saveVal;
            end; %Converts to degrees
        end; phi = phi*360/(2*pi); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read text files and assign values to variables
paramVec = readData(fileNameVec{1}); % r_blade, V_wind, gamma, C_p, M_blades 
geomVec = readData(fileNameVec{2}); % l1 l2 l3 l4 l5 l6 l7 l8 d1 d2 d3 d4 r1 r2
matVec = readData(fileNameVec{3}); % E, v,  sigma0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Derived Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lTotLow = geomVec(1,1) + geomVec(2,1) + geomVec(3,1) + geomVec(4,1);
lTotHigh = geomVec(5,1) + geomVec(6,1) + geomVec(7,1) + geomVec(8,1);
W_blades = -g*paramVec(5,1);
xVecLow = linspace(0, lTotLow, 10000);
xVecHigh = linspace(0, lTotHigh, 10000);
Pwind = paramVec(4,1)*0.5*rhoAir*pi*paramVec(1,1)^2*paramVec(2,1)^3;
U = paramVec(2,1)*paramVec(3,1); f = U/(2*pi*paramVec(1,1));
Tin = Pwind/(2*pi*f); TinVec = zeros(1, length(xVecLow));
TinVec(1, :) = Tin; FgearLow = Tin/geomVec(13,1); FgearHigh = -FgearLow;
Tout = -Tin*geomVec(14,1)/geomVec(13,1);
ToutVec = zeros(1, length(xVecLow)); ToutVec(1, :) = Tout;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deriving Reaction Forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rVecLowXY = reacSolve(W_blades, 0, geomVec(1,1), geomVec(1,1)+geomVec(2,1)+...
    geomVec(3));
rVecLowXZ = reacSolve(FgearLow,lTotLow,geomVec(1,1),geomVec(1,1)+geomVec(2,1)+...
    geomVec(3));
rVecHighXZ = reacSolve(FgearHigh,0,geomVec(5,1),geomVec(5,1)+geomVec(6,1)+...
    geomVec(7,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deriving Shears
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VyVecLow = shearSolve(W_blades, 0, rVecLowXY(1,1), geomVec(1,1), ...
    geomVec(1,1)+geomVec(2,1)+geomVec(3), xVecLow);
VzVecLow = shearSolve(rVecLowXZ(1,1), geomVec(1,1), rVecLowXZ(2,1), ...
    geomVec(1,1)+geomVec(2,1)+geomVec(3,1), lTotLow ,xVecLow);
VzVecHigh = shearSolve(FgearHigh, 0, rVecHighXZ(1,1), geomVec(5,1), ...
    geomVec(5,1)+geomVec(6,1)+geomVec(7,1), xVecHigh);
VyVecHigh = zeros(1, length(xVecHigh));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deriving Moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MzVecLow = momentSolve(W_blades, 0, rVecLowXY(1,1), geomVec(1,1), ...
    geomVec(1,1)+geomVec(2,1)+geomVec(3), xVecLow);
MyVecLow = momentSolve(rVecLowXZ(1,1), geomVec(1,1), rVecLowXZ(2,1), ...
    geomVec(1,1)+geomVec(2,1)+geomVec(3,1), lTotLow ,xVecLow);
MyVecHigh = momentSolve(FgearHigh, 0, rVecHighXZ(1,1), geomVec(5,1), ...
    geomVec(5,1)+geomVec(6,1)+geomVec(7,1), xVecHigh);
MzVecHigh = zeros(1, length(xVecHigh));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deriving Max Eff Stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lowEffStressA = effStressSolve(xVecLow,MzVecLow, MyVecLow, geomVec(9,1)/2,...
    geomVec(10,1)/2, geomVec(1,1)+geomVec(2,1), Tin);
highEffStressA = effStressSolve(xVecHigh,MzVecHigh, MyVecHigh, geomVec(11,1)/2,...
    geomVec(12,1)/2, geomVec(5,1)+geomVec(6,1), Tout);
lowEffStressB = effStressShearSolve(xVecLow,MzVecLow, MyVecLow,VzVecLow, ...
    VyVecLow, geomVec(9,1)/2,...
    geomVec(10,1)/2, geomVec(1,1)+geomVec(2,1), Tin);
highEffStressB = effStressShearSolve(xVecHigh,MzVecHigh, MyVecHigh, ...
    VzVecHigh, VyVecHigh, geomVec(11,1)/2,...
    geomVec(12,1)/2, geomVec(5,1)+geomVec(6,1), Tout);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deriving Max Tau at H and K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lowTauH, lowTauK] = tauSolve(xVecLow,MzVecLow, MyVecLow,VzVecLow, ...
    VyVecLow, geomVec(9,1)/2,...
    geomVec(10,1)/2, geomVec(1,1)+geomVec(2,1), Tin);
[highTauH, highTauK] = tauSolve(xVecHigh,MzVecHigh, MyVecHigh, ...
    VzVecHigh, VyVecHigh, geomVec(11,1)/2,...
    geomVec(12,1)/2, geomVec(5,1)+geomVec(6,1), Tout);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deriving Angles of Twist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lowPhi = twist(xVecLow, Tin, matVec(1,1), matVec(2,1), geomVec(9,1)/2, ...
    geomVec(10,1)/2, geomVec(1,1)+geomVec(2,1));
highPhi = twist(xVecHigh, Tout, matVec(1,1), matVec(2,1), geomVec(11,1)/2, ...
    geomVec(12,1)/2, geomVec(5,1)+geomVec(6,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding Max and Min Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[maxVyLow, maxIVyLow] = max(VyVecLow);
str1 = ['(', num2str(round(xVecLow(maxIVyLow)*100)/100), ', ', ...
    num2str(maxVyLow), ')'];
[minVyLow, minIVyLow] = min(VyVecLow);
str2 = ['(', num2str(round(xVecLow(minIVyLow)*100)/100), ', ', ...
    num2str(minVyLow), ')'];
[maxVzLow, maxIVzLow] = max(VzVecLow);
str3 = ['(', num2str(round(xVecLow(maxIVzLow)*100)/100), ', ', ...
    num2str(maxVzLow), ')'];
[minVzLow, minIVzLow] = min(VzVecLow);
str4 = ['(', num2str(round(xVecLow(minIVzLow)*100)/100), ', ', ...
    num2str(minVzLow), ')'];
[maxVzHigh, maxIVzHigh] = max(VzVecHigh);
str5 = ['(', num2str(round(xVecHigh(maxIVzHigh)*100)/100), ', ', ...
    num2str(maxVzHigh), ')'];
[minVzHigh, minIVzHigh] = min(VzVecHigh);
str6 = ['(', num2str(round(xVecHigh(minIVzHigh)*100)/100), ', ', ...
    num2str(minVzHigh), ')'];
[maxMzLow, maxIMzLow] = max(MzVecLow);
str7 = ['(', num2str(round(xVecLow(maxIMzLow)*100)/100), ', ', ...
    num2str(maxMzLow), ')'];
[minMzLow, minIMzLow] = min(MzVecLow);
str8 = ['(', num2str(round(xVecLow(minIMzLow)*100)/100), ', ', ...
    num2str(minMzLow), ')'];
[maxMyLow, maxIMyLow] = max(MyVecLow);
str9 = ['(', num2str(round(xVecLow(maxIMyLow)*100)/100), ', ', ...
    num2str(maxMyLow), ')'];
[minMyLow, minIMyLow] = min(MyVecLow);
str10 = ['(', num2str(round(xVecLow(minIMyLow)*100)/100), ', ', ...
    num2str(minMyLow), ')'];
[maxMyHigh, maxIMyHigh] = max(MyVecHigh);
str11 = ['(', num2str(round(xVecHigh(maxIMyHigh)*100)/100), ', ', ...
    num2str(maxMyHigh), ')'];
[minMyHigh, minIMyHigh] = min(MyVecHigh);
str12 = ['(', num2str(round(xVecHigh(minIMyHigh)*100)/100), ', ', ...
    num2str(minMyHigh), ')'];
str13 = ['(0, ', num2str(Tin), ')']; str14 = ['(0, ', num2str(Tout), ')'];
[maxlowEffStressA, maxIlowEffStressA] = max(lowEffStressA);
str15 = ['(', num2str(round(xVecLow(maxIlowEffStressA)*100)/100), ', ', ...
    num2str(maxlowEffStressA/10^6), ')'];
[minlowEffStressA, minIlowEffStressA] = min(lowEffStressA);
str16 = ['(', num2str(round(xVecLow(minIlowEffStressA)*100)/100), ', ', ...
    num2str(minlowEffStressA/10^6), ')'];
[maxhighEffStressA, maxIhighEffStressA] = max(highEffStressA);
str17 = ['(', num2str(round(xVecHigh(maxIhighEffStressA)*100)/100), ', ', ...
    num2str(maxhighEffStressA/10^6), ')'];
[minhighEffStressA, minIhighEffStressA] = min(highEffStressA);
str18 = ['(', num2str(round(xVecHigh(minIhighEffStressA)*100)/100), ', ', ...
    num2str(minhighEffStressA/10^6), ')'];
[maxlowEffStressB, maxIlowEffStressB] = max(lowEffStressB);
str19 = ['(', num2str(round(xVecLow(maxIlowEffStressB)*100)/100), ', ', ...
    num2str(maxlowEffStressB/10^6), ')'];
[minlowEffStressB, minIlowEffStressB] = min(lowEffStressB);
str20 = ['(', num2str(round(xVecLow(minIlowEffStressB)*100)/100), ', ', ...
    num2str(minlowEffStressB/10^6), ')'];
[maxhighEffStressB, maxIhighEffStressB] = max(highEffStressB);
str21 = ['(', num2str(round(xVecHigh(maxIhighEffStressB)*100)/100), ', ', ...
    num2str(maxhighEffStressB/10^6), ')'];
[minhighEffStressB, minIhighEffStressB] = min(highEffStressB);
str22 = ['(', num2str(round(xVecHigh(minIhighEffStressB)*100)/100), ', ', ...
    num2str(minhighEffStressB/10^6), ')'];
[maxlowTauH, maxIlowTauH] = max(lowTauH);
str23 = ['(', num2str(round(xVecLow(maxIlowTauH)*100)/100), ', ', ...
    num2str(maxlowTauH/10^6), ')'];
[minlowTauH, minIlowTauH] = min(lowTauH);
str24 = ['(', num2str(round(xVecLow(minIlowTauH)*100)/100), ', ', ...
    num2str(minlowTauH/10^6), ')'];
[maxlowTauK, maxIlowTauK] = max(lowTauK);
str25 = ['(', num2str(round(xVecLow(maxIlowTauK)*100)/100), ', ', ...
    num2str(maxlowTauK/10^6), ')'];
[minlowTauK, minIlowTauK] = min(lowTauK);
str26 = ['(', num2str(round(xVecLow(minIlowTauK)*100)/100), ', ', ...
    num2str(minlowTauK/10^6), ')'];
[maxhighTauH, maxIhighTauH] = max(highTauH);
str27 = ['(', num2str(round(xVecHigh(maxIhighTauH)*100)/100), ', ', ...
    num2str(maxhighTauH/10^6), ')'];
[minhighTauH, minIhighTauH] = min(highTauH);
str28 = ['(', num2str(round(xVecHigh(minIhighTauH)*100)/100), ', ', ...
    num2str(minhighTauH/10^6), ')'];
[maxhighTauK, maxIhighTauK] = max(highTauK);
str29 = ['(', num2str(round(xVecHigh(maxIhighTauK)*100)/100), ', ', ...
    num2str(maxhighTauK/10^6), ')'];
[minhighTauK, minIhighTauK] = min(highTauK);
str30 = ['(', num2str(round(xVecHigh(minIhighTauK)*100)/100), ', ', ...
    num2str(minhighTauK/10^6), ')'];
[maxlowPhi, maxIlowPhi] = max(lowPhi);
str31 = ['(', num2str(round(xVecLow(maxIlowPhi)*100)/100), ', ', ...
    num2str(maxlowPhi), ')'];
[minlowPhi, minIlowPhi] = min(lowPhi);
str32 = ['(', num2str(round(xVecLow(minIlowPhi)*100)/100), ', ', ...
    num2str(minlowPhi), ')'];
[maxhighPhi, maxIhighPhi] = max(highPhi);
str33 = ['(', num2str(round(xVecHigh(maxIhighPhi)*100)/100), ', ', ...
    num2str(maxhighPhi), ')'];
[minhighPhi, minIhighPhi] = min(highPhi);
str34 = ['(', num2str(round(xVecHigh(minIhighPhi)*100)/100), ', ', ...
    num2str(minhighPhi), ')'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Checking Shear Assumption %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg1 = 'Cannot ignore transverse shear in low-speed shaft';
msg2 = 'Cannot ignore transverse shear in high-speed shaft';
checkShear(lowEffStressA, lowEffStressB, msg1);
checkShear(highEffStressA, highEffStressB, msg2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1 = figure(1);
set(fig1, 'Position', [50,1,1200, 800]);
s1 = subplot(2,1,1); hold on;
plot(xVecLow, VyVecLow);
plot(xVecLow(maxIVyLow), maxVyLow, 'xk', 'Markersize', 14);
text(xVecLow(maxIVyLow), maxVyLow, str1);
plot(xVecLow(minIVyLow), minVyLow, 'xk', 'Markersize', 14);
text(xVecLow(minIVyLow), minVyLow, str2);
title(s1, 'V_{y} Shear Diagram for Low Speed Shaft in XY Plane', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Shear Force (N)', 'Fontsize', 14);
s2 = subplot(2,1,2); hold on;
plot(xVecLow, MzVecLow);
plot(xVecLow(maxIMzLow), maxMzLow, 'xk', 'Markersize', 14);
text(xVecLow(maxIMzLow), maxMzLow, str7);
plot(xVecLow(minIMzLow), minMzLow, 'xk', 'Markersize', 14);
text(xVecLow(minIMzLow), minMzLow, str8);
title(s2, 'M_{z} Moment Diagram for Low Speed Shaft in XY Plane', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Bending Moment (N m)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure(2); set(fig2, 'Position', [50,1,1200, 800]); 
s3 = subplot(2,1,1); hold on; 
plot(xVecLow, VzVecLow); plot(xVecLow(maxIVzLow), maxVzLow, 'xk',...
    'Markersize', 14);
text(xVecLow(maxIVzLow), maxVzLow, str3);
plot(xVecLow(minIVzLow), minVzLow, 'xk', 'Markersize', 14);
text(xVecLow(minIVzLow), minVzLow, str4);
title(s3, 'V_{z} Shear Diagram for Low Speed Shaft in XZ Plane', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Shear Force (N)', 'Fontsize', 14);
s4 = subplot(2,1,2); hold on;
plot(xVecLow, MyVecLow);
plot(xVecLow(maxIMyLow), maxMyLow, 'xk', 'Markersize', 14);
text(xVecLow(maxIMyLow), maxMyLow, str9);
plot(xVecLow(minIMyLow), minMyLow, 'xk', 'Markersize', 14);
text(xVecLow(minIMyLow), minMyLow, str10);
title(s4, 'M_{y} Moment Diagram for Low Speed Shaft in XZ Plane', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Bending Moment (N m)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig3 = figure(3);
set(fig3, 'Position', [50,1,1200, 800]);
s5 = subplot(2,1,1); hold on;
plot(xVecHigh, VzVecHigh);
plot(xVecHigh(maxIVzHigh), maxVzHigh, 'xk', 'Markersize', 14);
text(xVecHigh(maxIVzHigh), maxVzHigh, str5);
plot(xVecHigh(minIVzHigh), minVzHigh, 'xk', 'Markersize', 14);
text(xVecHigh(minIVzHigh), minVzHigh, str6);
title(s5, 'V_{z} Shear Diagram for High Speed Shaft in XZ Plane', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Shear Force (N)', 'Fontsize', 14);
s6 = subplot(2,1,2); hold on;
plot(xVecHigh, MyVecHigh);
plot(xVecHigh(maxIMyHigh), maxMyHigh, 'xk', 'Markersize', 14);
text(xVecHigh(maxIMyHigh), maxMyHigh, str11);
plot(xVecHigh(minIMyHigh), minMyHigh, 'xk', 'Markersize', 14);
text(xVecHigh(minIMyHigh), minMyHigh, str12);
title(s6, 'M_{y} Moment Diagram for High Speed Shaft in XZ Plane', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Bending Moment (N m)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig4 = figure(4);
set(fig4, 'Position', [50,1,1200, 800]);
s7 = subplot(2,1,1); hold on;
plot(xVecLow, TinVec);
plot(0, Tin, 'xk', 'Markersize', 14);
text(0, Tin, str13);
title(s7, 'Torque Along Shaft for Low Speed Shaft', 'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Torque (N m)', 'Fontsize', 14);
s8 = subplot(2,1,2); hold on;
plot(xVecHigh, ToutVec);
plot(0, Tout, 'xk', 'Markersize', 14);
text(0, Tout, str14);
title(s8, 'Torque Along Shaft for High Speed Shaft', 'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Torque (N m)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig5 = figure(5);
set(fig5, 'Position', [50,1,1200, 800]);
s9 = subplot(2,1,1); hold on;
plot(xVecLow, lowEffStressA/10^6);
plot(xVecLow(maxIlowEffStressA), maxlowEffStressA/10^6, 'xk', ...
    'Markersize', 14);
plot(xVecLow(minIlowEffStressA), minlowEffStressA/10^6, 'xk', ...
    'Markersize', 14);
text(xVecLow(maxIlowEffStressA), maxlowEffStressA/10^6, str15);
text(xVecLow(minIlowEffStressA), minlowEffStressA/10^6, str16);
title(s9, 'Max Effective Stress Along Shaft for Low Speed Shaft (No V)', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Max Effective Stress (MPa)', 'Fontsize', 14);
s10 = subplot(2,1,2); hold on;
plot(xVecHigh, highEffStressA/10^6);
plot(xVecHigh(maxIhighEffStressA), maxhighEffStressA/10^6, 'xk', ...
    'Markersize', 14);
plot(xVecHigh(minIhighEffStressA), minhighEffStressA/10^6, 'xk', ...
    'Markersize', 14);
text(xVecHigh(maxIhighEffStressA), maxhighEffStressA/10^6, str17);
text(xVecHigh(minIhighEffStressA), minhighEffStressA/10^6, str18);
title(s10, 'Max Effective Stress Along Shaft for High Speed Shaft (No V)', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Max Effective Stress (MPa)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig6 = figure(6);
set(fig6, 'Position', [50,1,1200, 800]);
s11 = subplot(2,1,1); hold on;
plot(xVecLow, lowEffStressB/10^6);
plot(xVecLow(maxIlowEffStressB), maxlowEffStressB/10^6, 'xk', 'Markersize', ...
    14);
plot(xVecLow(minIlowEffStressB), minlowEffStressB/10^6, 'xk', 'Markersize', ...
    14);
text(xVecLow(maxIlowEffStressB), maxlowEffStressB/10^6, str19);
text(xVecLow(minIlowEffStressB), minlowEffStressB/10^6, str20);
title(s11, 'Max Effective Stress Along Shaft for Low Speed Shaft (With V)', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Max Effective Stress (MPa)', 'Fontsize', 14);
s12 = subplot(2,1,2); hold on;
plot(xVecHigh, highEffStressB/10^6);
plot(xVecHigh(maxIhighEffStressB), maxhighEffStressB/10^6, 'xk', ...
    'Markersize', 14);
plot(xVecHigh(minIhighEffStressB), minhighEffStressB/10^6, 'xk', ...
    'Markersize', 14);
text(xVecHigh(maxIhighEffStressB), maxhighEffStressB/10^6-5, str21);
text(xVecHigh(minIhighEffStressB), minhighEffStressB/10^6-5, str22);
title(s12, 'Max Effective Stress Along Shaft for High Speed Shaft (With V)', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Max Effective Stress (MPa)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig7 = figure(7);
set(fig7, 'Position', [50,1,1200, 800]);
s13 = subplot(2,1,1); hold on;
plot(xVecLow, lowTauH/10^6);
plot(xVecLow(maxIlowTauH), maxlowTauH/10^6, 'xk', 'Markersize', 14);
plot(xVecLow(minIlowTauH), minlowTauH/10^6, 'xk', 'Markersize', 14);
text(xVecLow(maxIlowTauH), maxlowTauH/10^6, str23);
text(xVecLow(minIlowTauH), minlowTauH/10^6, str24);
title(s13, 'Max Shear Stress At H Along Shaft for Low Speed Shaft', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Max Shear Stress (MPa)', 'Fontsize', 14);
s14 = subplot(2,1,2); hold on;
plot(xVecLow, lowTauK/10^6);
plot(xVecLow(maxIlowTauK), maxlowTauK/10^6, 'xk', 'Markersize', 14);
plot(xVecLow(minIlowTauK), minlowTauK/10^6, 'xk', 'Markersize', 14);
text(xVecLow(maxIlowTauK), maxlowTauK/10^6, str25);
text(xVecLow(minIlowTauK), minlowTauK/10^6, str26);
title(s14, 'Max Shear Stress At K Along Shaft for Low Speed Shaft', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Max Shear Stress (MPa)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig8 = figure(8);
set(fig8, 'Position', [50,1,1200, 800]);
s15 = subplot(2,1,1); hold on;
plot(xVecHigh, highTauH/10^6);
plot(xVecHigh(maxIhighTauH), maxhighTauH/10^6, 'xk', 'Markersize', 14);
plot(xVecHigh(minIhighTauH), minhighTauH/10^6, 'xk', 'Markersize', 14);
text(xVecHigh(maxIhighTauH), maxhighTauH/10^6, str27);
text(xVecHigh(minIhighTauH), minhighTauH/10^6, str28);
title(s15, 'Max Shear Stress At H Along Shaft for High Speed Shaft', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Max Shear Stress (MPa)', 'Fontsize', 14);
s16 = subplot(2,1,2); hold on;
plot(xVecHigh, highTauK/10^6);
plot(xVecHigh(maxIhighTauK), maxhighTauK/10^6, 'xk', 'Markersize', 14);
plot(xVecHigh(minIhighTauK), minhighTauK/10^6, 'xk', 'Markersize', 14);
text(xVecHigh(maxIhighTauK), maxhighTauK/10^6, str29);
text(xVecHigh(minIhighTauK), minhighTauK/10^6, str30);
title(s16, 'Max Shear Stress At K Along Shaft for High Speed Shaft', ...
    'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Max Shear Stress (MPa)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig9 = figure(9);
set(fig9, 'Position', [50,1,1200, 800]);
s17 = subplot(2,1,1); hold on;
plot(xVecLow, lowPhi);
plot(xVecLow(maxIlowPhi), maxlowPhi, 'xk', 'Markersize', 14);
plot(xVecLow(minIlowPhi), minlowPhi, 'xk', 'Markersize', 14);
text(xVecLow(maxIlowPhi), maxlowPhi, str31);
text(xVecLow(minIlowPhi), minlowPhi, str32);
title(s17, 'Angle of Twist on Low Speed Shaft', 'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Angle of Twist (deg)', 'Fontsize', 14);
s18 = subplot(2,1,2); hold on;
plot(xVecHigh, highPhi);
plot(xVecHigh(maxIhighPhi), maxhighPhi, 'xk', 'Markersize', 14);
plot(xVecHigh(minIhighPhi), minhighPhi, 'xk', 'Markersize', 14);
text(xVecHigh(maxIhighPhi)+0.03, maxhighPhi, str33);
text(xVecHigh(minIhighPhi)+0.03, minhighPhi, str34);
title(s18, 'Angle of Twist on High Speed Shaft', 'Fontsize', 20);
xlabel('Distance Along Shaft (m)', 'Fontsize', 14);
ylabel('Angle of Twist (deg)', 'Fontsize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finding max values for effective stresses in each beam to find S.F.
max1 = max(maxhighEffStressA/10^6, maxhighEffStressB/10^6);
max2 = max(maxlowEffStressA/10^6, maxlowEffStressB/10^6);
%Printing safety factors:
fprintf('The safety factor for stress in the low speed shaft is %6.4f\n',...
    matVec(3,1)/max2);
fprintf('The safety factor for stress in the high speed shaft is %6.4f\n',...
    matVec(3,1)/max1);
fprintf('The safety factor for twist in the low speed shaft is %6.4f\n', ...
    abs(phiAll/maxlowPhi));
fprintf('The safety factor for twist in the high speed shaft is %6.4f\n', ...
    abs(phiAll/minhighPhi));
end

