
forceVec = [0 19.04 38.53 58.81 65.63 67.60 67.16 67.75 68.63 69.43 70.02];
delLVec = [0 0.0794 0.16 0.2505 0.2815 0.2952 0.3944 0.6573 1.0951 1.5156 1.9534];
initL = 50.8;
initD = 8.56;

forceVec = forceVec*1000;
delLVec = delLVec/1000;
initL = initL/1000;
initD = initD/1000;

a = pi*(initD/2)^2;

stressVec = forceVec/a;
strainVec = delLVec/initL;

figure(1);
plot(strainVec, stressVec/10^6);

clear all;

stressVec = [0 202 403 587 785 822 836 832 829 828 864 897 912 918 915 899 ...
    871 831 772 689 574];
strainVec = [0 0.099 0.195 0.283 0.382 0.405 0.423 0.451 0.887 1.988 2.94 ...
    4.51 5.96 8.07 9.94 12.04 13.53 15.03 16.70 18.52 20.35];

strainVec = strainVec/100;

figure(2);
plot(strainVec, stressVec);

E = (stressVec(4) - stressVec(1))/(strainVec(4)-strainVec(1));
fprintf('The Youngs Modulus is %e MPa\n', E);
fprintf('The max tensile strength is %e MPa\n', max(stressVec));
fprintf('The percent elongation is %e \n', strainVec(end)*100);

ri = 9.09/2;
ri = ri/1000;
rf = 5.56/2;
rf = rf/1000;
ai = pi*ri^2;
af = pi*rf^2;

perRA = 100*(ai-af)/ai;
fprintf('The percent reduction in area is %e \n', perRA);

clear all;

engStressVec = [0 125 257 359 317 333 357 397 458 507 541 576 558 531 476 379];
engStrainVec = [0 0.0006 0.0012 0.0017 0.0035 0.007 0.01 0.017 0.03 0.05 0.079 0 0 0 0 0];
radVec = [6.32 0 0 0 0 0 0 0 0 0 5.99 5.72 5.33 5.08 4.45 3.5];
radVec = radVec/2000;
stressVec = zeros(1, length(engStressVec));
strainVec = zeros(1, length(engStrainVec));
strainVec(1:11) = exp(engStrainVec(1:11))-1;
stressVec(1:11) = engStressVec(1:11)./(1+strainVec(1:11));
areaVec = radVec.^2*pi;
Ainit = 0.00316^2*pi;
stressVec(12:16) = (engStressVec(12:16).*areaVec(12:16))/Ainit;
engStrainVec(12:16) = log(Ainit./areaVec(12:16));
strainVec(12:16) = exp(engStrainVec(12:16))-1;