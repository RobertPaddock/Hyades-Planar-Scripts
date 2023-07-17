%This script reads Hyades Planar simulations. It is a more basic script
%with the shock tracking etc. removed.

file = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Desktop Hyades Files\Dated\230601\File2.cdf'

figures = 1;
setlimits = 0;
Measured = 0;


%Import the key variables, and get in correct format
Radius  = ncread(file,'R');
Time  = ncread(file,'DumpTimes');
Radiuscm  = ncread(file,'Rcm');
Density  = ncread(file,'Rho'); %g/cm^3
ElecTemp  = ncread(file,'Te')*1000; %Converted to eV from keV
IonTemp  = ncread(file,'Ti')*1000; %Converted to eV from keV
IonTempKelvin = 1.16E4 .*IonTemp;
Pressure = ncread(file,'Pres').*0.0000001; %Converted to J/cm^3 from dyn/cm^2
Velocity = ncread(file, 'U')./100000; %Converted to km/s from cm/s
IonThermalEnergy  = ncread(file,'Eion').*10^-7; %Ion Thermal Energy, converted from erg to J.
ElectronThermalEnergy  = ncread(file,'Eelc').*10^-7; %Electron Thermal Energy, converted from erg to J.
%KineticEnergy  = ncread(file,'Ekint').*10^-7;
DepositedLaserPowerZones = ncread(file,'Deplas').*10^-7; %Deposited Laser Energy per zone, converted from erg to J.
DepositedLaserPower = sum(DepositedLaserPowerZones).'; %Was divided by 0.8 previously
DepositedEnergy = trapz(Time, DepositedLaserPower);
Volume =  ncread(file,'Vol'); %cm^3
Mass = Density.*Volume(2:end,:);
InternalEnergy = IonThermalEnergy+ElectronThermalEnergy;
ElectronDensity = ncread(file,'Dene'); % Electron density per cm3
try
Zions = ncread(file, 'Zions');

LaserDep =ncread(file,'Deplas');
catch 
end
Zbar = ncread(file, 'Zbar');
LaserEnergy  = ncread(file,'Elasin')*(10^-7); %Was divided by 0.8 previously
LaserPower=[0; diff(LaserEnergy)./diff(Time)];
%SimulationPower = LaserPower*0.8;

%Laser statistics calculated. Hyades has a 1cm^2 area, so multiply by area
%in cm^2 to get actual power. This also means that laser power is actually
%Laser Intensity.
AdjustedEnergy = (LaserEnergy) * 0.02 * 0.1;
AdjustedPower = (LaserPower) * 0.02 * 0.1;
TotalLaserPower = AdjustedPower;
SimulationPower = TotalLaserPower*0.8;
DepositedLaserPower = (DepositedLaserPower) * 0.02 * 0.1;
LaserIntensity = max(LaserPower);

% FoamDensity = 0.253;
FoamDensity = Density(end,1);

%Define shock front as last zone in increasing radius in the sim where
%there is a substantial difference between density now and original density
%- since density doesn't really change until shock reaches it. Doesn't work
%so well for gold/ start of al where there is preheat.
Shock = (Density - repmat(Density(:,1), 1, length(Time))) > 0.1;
%Shock = (Density ./ repmat(Density(:,1), 1, length(Time))) > 1.5;

%To avoid it picking up changing density in aluminised plastic files, I
%have also specified that these layers shouldn't be included.
FoamMaxRadius = max(Radius(:,1));
BulkRegions = Radius(1:end-1,1)<FoamMaxRadius;
Shock = Shock.*BulkRegions;
ShockIndex = arrayfun(@(x)find(Shock(:,x),1,'last'),1:size(Shock,2), 'UniformOutput',false);
tf = cellfun('isempty',ShockIndex); % true for empty cells
ShockIndex(tf) = {1};
ShockIndex= [ShockIndex{:}];
ShockIndexReshaped=sub2ind(size(Radius),ShockIndex,1:length(Time));
ShockRadius=Radius(ShockIndexReshaped);



figure;
DensityPlot = surf(Time.*10^9,Radius(1:end-1,:).*10000, Density);
%xlim([(MaxDensityTime-.5e-9) (MaxDensityTime+.5e-9)]);
%xlim([xlowerlim xupperlim]);
%ylim([0 0.01]);

%zlim([0 MaxDensity]);
colorbar;
colormap(jet);
title('Density Plot')
xlabel('Time (ns)');
ylabel('Radius (\mu m)')
zlabel('Mass Density (g/cm^3)');
shading interp
view(2)
set(gca,'colorscale','log')
caxis([0.01 20]);
ylim([0 110])

% hold on
% plot3(Time(1:length(ShockFrontRadius))*10^9,ShockFrontRadius.*10000, max(Density(:, 1:length(ShockFrontRadius))), 'b');
% plot3(Time(1:length(ShockRadius))*10^9,ShockRadius.*10000, max(Density(:, 1:length(ShockRadius))), 'w');
% hold off
% xline(Time(minShockinFoam)*10^9);
% xline(Time(minShockinQuartz)*10^9);
% xline(Time(maxShockinFoam)*10^9);

figure;
DensityPlot = surf(Time.*10^9,Radius(1:end-1,:).*10000, IonTempKelvin);

colorbar;
colormap(jet);
title('Ion Temp Plot');
xlabel('Time (ns)');
ylabel('Radius (cm)');
zlabel('Mass Density (g/cm^3)');
shading interp
set(gca,'colorscale','log')
view(2)

% figure;
% DensityPlot = surf(Time.*10^9,Radius(1:end-1,:), ElecTemp);
% %xlim([(MaxDensityTime-.5e-9) (MaxDensityTime+.5e-9)]);
% %xlim([xlowerlim xupperlim]);
% %ylim([0 0.01]);
% ylim([0 0.01]);
% %zlim([0 MaxDensity]);
% colorbar;
% colormap(jet);
% title('Elec Temp Plot')
% xlabel('Time (ns)');
% ylabel('Radius (cm)');
% zlabel('Mass Density (g/cm^3)');
% shading interp
% view(2)
% set(gca,'colorscale','log')

%Plot zoning
figure
MassDifference = 100*diff(Mass(:,1))./Mass(2:end,1);
plot(MassDifference)
ylabel('Percentage Mass Difference')
yyaxis right;
plot(Mass(:,1), ':');
ylabel('Zone Mass')
title('Zoning Plot')
xlabel('Zone')

figure
title('Laser Power and input energy')
yyaxis left
plot(Time.*10^9, TotalLaserPower./10^12)
ylabel('Power (TW)')
yyaxis right
plot(Time.*10^9, AdjustedEnergy)
ylabel('Energy (J)')

dim = [.2 .5 .3 .3];
str = "Laser Intensity = " + string(LaserIntensity./10^12) + " TW/cm^2";
annotation('textbox',dim,'String',str,'FitBoxToText','on');


% Zone log plot - good for tracking shocks
figure;
DensityPlot = surf(Time.*10^9,repmat([1:size(Radius)-1].',1,length(Time)), Density);
%xlim([(MaxDensityTime-.5e-9) (MaxDensityTime+.5e-9)]);
%xlim([xlowerlim xupperlim]);
%ylim([0 0.01]);
xlim([0 Time(end).*10^9]);
%zlim([0 MaxDensity]);
colorbar;
colormap(jet);
title('Density Plot')
xlabel('Time (ns)');
ylabel('Zones')
zlabel('Mass Density (g/cm^3)');
shading interp
view(2)
set(gca,'colorscale','log')
caxis([0.01 20]);
hold on
%plot3(Time*10^9,ShockRadius.*10000, max(Density), 'w');
hold off

% InverseDensityScaleLength = abs( diff(log(Density))./diff(Radius(1:end-1:end, :)));
InverseDensityScaleLength = abs( diff(log(Density)));

InverseDensityScaleLength(isnan(InverseDensityScaleLength)) = 0.001;
InverseDensityScaleLength((InverseDensityScaleLength==0)) = 0.001;
InverseDensityScaleLength(isinf(InverseDensityScaleLength)) = 0.001; 
Zones = repmat([1:size(Radius,1)-2].', 1, length(Time));

 figure
        surf(Time*10^9,Zones, InverseDensityScaleLength);
        xlim([0 max(Time)*10^9]);
        set(gca,'ColorScale','log');
        colormap(jet);
        title('Shock trajectories');
        xlabel('Time (ns)');
        ylabel('Cell number');
        zlabel('Mass Density (g/cm^3)');
        shading interp
        colorbar
                view(2)
                caxis([0.001 1]);

                figure
yyaxis left
plot(Time.*10^9, IonTemp(end, :).*11600, '-', 'LineWidth',2);
ylim([280 300]);
ylabel('Ion temperature (Kelvin)')
yyaxis right
plot(Time.*10^9, ElecTemp(end, :).*11600, '--', 'LineWidth',2);
ylim([280 300]);
ylabel('Electron temperature (Kelvin)')
xlabel('Time (ns)')
   

function txt = myupdatefcn(~,event_obj,Time, minShockinFoam)
% Customizes text of data tips
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');
txt = {['X: ',num2str(pos(1))],...
       ['Y: ',num2str(pos(2))],...
       ['Time: ',num2str(Time(I+minShockinFoam))]};
end

function [fitresult, gof] = createFit(x, y)
warning('off','all');
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

end




function [FoamApproximateUp, FoamApproximatePressure, Ps, RayleighLine, R, HugoniotPressure, Up, Up1, P1] = ExactImpedanceMatch(AverageQuartzShockFrontVelocity, AverageFoamShockFrontVelocity)

% AverageQuartzShockFrontVelocity = QuartzShockVelocity_MC;
% AverageFoamShockFrontVelocity = FoamShockVelocity_MC;

% AverageQuartzShockFrontVelocity = AllOverlapVelocities(1,:).';
% AverageFoamShockFrontVelocity = AllOverlapVelocities(2,:).';
% AverageQuartzShockFrontVelocity = [14.4545; 15; 18];
% AverageFoamShockFrontVelocity = [14.4183; 15; 18];

%Define initial quartz parameters
Rho0 = 2.65;
V0 = 1/Rho0; %V is specific volume, 1/density

HugoniotUp = 0:0.001:30;

%Linear Hugoniot used in this one
a = 1.754;
b = 1.862;
c = -3.364E-2;
d = 5.666E-4;
HugoniotUs = a + (b.*HugoniotUp) + (c.*HugoniotUp.^2) + (d.*HugoniotUp.^3);
R = [a,b,c,d];

% Cov = [2.097E-2 -6.159E-3 5.566E-4 -1.572E-5;...
%     -6.159E-3 1.877E-3 -1.742E-4 5.017E-6; ...
%     5.566E-4 -1.742E-4 1.650E-5 -4.834E-7; ...
%     -1.572E-5 5.017E-6 -4.834E-7 1.438E-8];
% 
% I = [1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 1];
% 
% Coeffs = [1.754  1.862 -3.364E-2 5.666E-4];
% 
% %Covariance matrix is not positive semi definite, but adding this tiny
% %number is. Don't think this has a huge impact, and adding model error has
% %no effect anyway./
%  R = mvnrnd(Coeffs,Cov+0.0000000001*I,1);
%  a = R(1); b=R(2); c=R(3); d=R(4);
%  
% HugoniotUs = a + (b.*HugoniotUp) + (c.*HugoniotUp.^2) + (d.*HugoniotUp.^3);



%R-H expressions for pressure and density
HugoniotPressure = Rho0 .* HugoniotUs .* HugoniotUp;
HugoniotRho = Rho0 ./ (1 - HugoniotUp ./ HugoniotUs);
HugoniotV = 1./HugoniotRho;

% figure
% plot(HugoniotUp, HugoniotUs);
% xlabel('Particle Velocity Up (km/s)')
% ylabel('Shock Velocity Us (km/s)')
% xlim([4 20])
% title('Quartz Hugoniot')

%Select point for release isentrope

% [~, HugoniotIndex] = min(abs(repmat(HugoniotUs, length(AverageQuartzShockFrontVelocity), 1) - AverageQuartzShockFrontVelocity).');
% % [~, HugoniotIndex2] = min(abs(HugoniotPressure - 600))
% Us1 = [HugoniotUs(HugoniotIndex).'];
% P1 = [HugoniotPressure(HugoniotIndex).'];
% Up1 = [HugoniotUp(HugoniotIndex).'];
% Rho1 = [HugoniotRho(HugoniotIndex).'];
% V1 = 1 ./ Rho1;
% eta = Rho1/Rho0;

Us1 = AverageQuartzShockFrontVelocity;
Up1 = interp1(HugoniotUs, HugoniotUp, Us1);
P1 = interp1(HugoniotUs, HugoniotPressure, Us1);
Rho1 = interp1(HugoniotUs, HugoniotRho, Us1);
V1 = 1./Rho1;
eta = Rho1/Rho0;

% figure
% plot(HugoniotUp, HugoniotUs);
% xlabel('Particle Velocity Up (km/s)')
% ylabel('Shock Velocity Us (km/s)')
% % ylim([0 10000])
% xlim([4 20])
% title('Quartz Hugoniot')
% hold on
% scatter(Up1, Us1)
% hold off

%Constants/parameters determined from Us_quartz
Gruneison = 0.64;

for i = 1:size(V1);
V(i,:)= [V0:-(V0-V1(i))/1000:V1(i)];
end
    

%Interpolate V and Up to get Up for the V points we use in our analysis.
%Produce the Hugoniot corresponding to those V points
Up_H = interp1(HugoniotV, HugoniotUp, V);
Us_H = a + (b.*Up_H) + (c.*Up_H.^2) + (d.*Up_H.^3);
PH = Rho0 .* Us_H .* Up_H;

%In the paper, E1-E0 has an integral over dummy vairable VPrime with V as a lowwer limit. These limits can be
%rearranged to give an integral between fixed limits, and V as an upper
%limit. This can then be done as an array function, usinng trapz-cumtrapz.
VPrime = V; %Used as dummy, but actually just equals V.
EnergyDiffIntegrand = ((VPrime ./ V1).^Gruneison) .* PH .* ...
    (1 - (Gruneison./2) .* ((V0./VPrime) - 1));
% EnergyDiffIntegral = -(trapz(VPrime, EnergyDiffIntegrand,2) - cumtrapz(VPrime, EnergyDiffIntegrand,2));
for i= 1:size(VPrime,1)
    EnergyDiffIntegral(i,:) = -(trapz(VPrime(i,:), EnergyDiffIntegrand(i,:)) - cumtrapz(VPrime(i,:), EnergyDiffIntegrand(i,:)));
end


EnergyDiff = (P1.*V0/2) .* ( (eta-1)./eta ) .* (V1./V).^Gruneison - ((V1./V).^Gruneison).*EnergyDiffIntegral;

Ps = PH .* (1 - (Gruneison./2) .* ((V0./V) - 1)) + (Gruneison./V).* EnergyDiff;

% figure
% plot(V.', PH.')
% hold on
% scatter(V1.', P1.')
% plot(V.', Ps.')
% ylabel('Pressure (GPa)')
% xlabel('Specific Volume')
% legend('Hugoniot', 'Release Isentrope')

    

%As above, integrates over dummy variable with P1 as lower limit. We can
%reasrrange to have a fixed integral - an integral with P as an upper
%limit.
for i= 1:size(V,1)
[dPdV(i,:)]= gradient(Ps(i,:), V(i,:));
end
% ReducedV = V(2:end);
% ReducedPs = Ps(2:end);
Cs = sqrt(-(V.^2) .* dPdV);
UpIntegrand = V./Cs;
% UpIntegral = -(trapz(Ps, UpIntegrand,2) - cumtrapz(Ps, UpIntegrand,2));
for i= 1:size(Ps,1)
UpIntegral(i,:) = -(trapz(Ps(i,:), UpIntegrand(i,:)) - cumtrapz(Ps(i,:), UpIntegrand(i,:)));
end
%The expression they give has a plus, but this doesn't seem to work - I
%need a - to recreate their plots.
Up = Up1 - UpIntegral;

% figure
% plot(HugoniotUp, HugoniotPressure);
% xlabel('Particle Velocity Up (km/s)')
% ylabel('Pressure (GPa)')
% hold on
% plot(Up.', Ps.');
% hold off
% xlim([6 30])
% ylim([0 1100])
% title('Knudsen and Desjarlis plot recreated')

% RayleighLine = 0.253 * FoamShockVelocity * HugoniotUp;
RayleighLine = 0.253 .* AverageFoamShockFrontVelocity .* Up;

% figure
% plot(HugoniotUp, HugoniotPressure);
% xlabel('Particle Velocity Up (km/s)')
% ylabel('Pressure (GPa)')
% hold on
% plot(Up.', Ps.');
% plot(HugoniotUp, 0.253 .* AverageFoamShockFrontVelocity .* HugoniotUp)
% hold off
% xlim([6 30])
% ylim([0 1100])
% legend('Hugoniot', 'Release Isentrope', 'Rayleigh Line')

[~, FoamHugoniotIndex] = min(abs(Ps - RayleighLine).');
% FoamHugoniotIndex = MinHugoniotWindowIndex+FoamHugoniotIndex;

for i = 1:length(FoamHugoniotIndex)
FoamApproximatePressure(i) = RayleighLine(i, FoamHugoniotIndex(i));
FoamApproximateUp(i) = Up(i, FoamHugoniotIndex(i));
end

FoamApproximateRho = 0.253 ./ (1 - FoamApproximateUp ./ AverageFoamShockFrontVelocity);
end


