%This script reads the results of a planar Hyades simulation. This is a
%more complicated version of the script, setup for my 4 layer target. It
%will attempt to do shock tracking in the final two layers, produce plots
%of shock velocity and other shock variables vs time, and estimate Hugoniot
%quantities based on the measured shock velocity. It will not work
%completely if the shock has not left the rear of the target.

file = '\\aldaq1.physics.ox.ac.uk\Archer\Robert\Desktop Hyades Files\Dated\230705\SiO2_RoomTemp\File3.cdf'

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
%in cm^2 to get actual power. This also means that the 'LaserPower'
%variable is actually laser intensity, while 'AdjustedPower' is the actual
%power.
AdjustedEnergy = (LaserEnergy) * pi() * (0.03)^2;
AdjustedPower = (LaserPower) * pi() * (0.03)^2;
TotalLaserPower = AdjustedPower;
SimulationPower = TotalLaserPower*0.8;
DepositedLaserPower = (DepositedLaserPower) * pi() * (0.03)^2;
LaserIntensity = max(LaserPower);

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

if max(abs(Density(:,1) - 2.2)<0.01)
    QuartzDensity = 2.2;
    disp('Amorphous quartz, 2.2 density')
elseif max(abs(Density(:,1) - 2.6)<0.01)
    QuartzDensity = 2.6;
    disp('Silica, 2.6 density')
    else
            QuartzDensity = 2.6;
    disp('Quartz type not identified')
    end

%Identify the different regions based on density. In aluminised plastic
%runs we have additional al and ch sections, so I'm specifying that only
%those in the bulk regions count
Quartzzones = (abs(Density(:,1) - QuartzDensity)<0.01);
Quartzzones = Quartzzones .* BulkRegions;
CHzones = abs(Density(:,1) - FoamDensity)<0.01;
CHzones = CHzones .* BulkRegions;
if setlimits == 1
    CHzones = (Radius(:,1) > FoamMinRadius).*(Radius(:,1) < FoamMaxRadius);
end
Goldzones = abs(Density(:,1) - 19.3)<0.01;
[~,minQuartzIndex] = find(Quartzzones.', 1, 'first');
[~,maxQuartzIndex] = find(Quartzzones.', 1, 'last');
[~,minFoamIndex] = find(CHzones.', 1, 'first');
[~,maxFoamIndex] = find(CHzones.', 1, 'last');
% maxFoamIndex = maxFoamIndex-1; %Newly added - suddenly it started erroring that max foam index was the last radius (and so didn't exist in foam etc.)
[~,minGoldIndex] = find(Goldzones.', 1, 'first');

QuartzMaxRadius = Radius(maxQuartzIndex, 1);
QuartzMinRadius = Radius(minQuartzIndex, 1);
FoamMinRadius = Radius(minFoamIndex, 1);

%Identify when the shock is in Al or Foam by comparing shock radius to
%above radii. Final condition ensures the shock has not yet reached the
%rear surface (as I had problems with the shock radius at times past this
%in some files).
ShockInQuartz=(QuartzMaxRadius>ShockRadius).*(ShockRadius>QuartzMinRadius).*(Velocity(end, :)<1);
ShockInFoam=(FoamMaxRadius>ShockRadius).*(ShockRadius>FoamMinRadius).*(Velocity(end, :)<1);

%Identify the indexes of when the shock enters and leaves the different
%regions.
[~,minShockinQuartz] = find(ShockInQuartz, 1, 'first');
[~,maxShockinQuartz] = find(ShockInQuartz, 1, 'last');
[~,minShockinFoam] = find(ShockInFoam, 1, 'first');
[~,maxShockinFoam] = find(ShockInFoam, 1, 'last');


%Take average of ion temp over the region when the shock leaves the region.
%The Foam calculation does not use all zones in the foam, as the shock
%front can jump from within the region to outside in one time step - thus
%some material is not shocked at the maxShockinFoam time, and erraneous
%results were gerenated.
QuartzBreakoutTemp = mean(IonTemp(minQuartzIndex:maxQuartzIndex, maxShockinQuartz));
FoamBreakoutTemp = mean(IonTemp(minFoamIndex:ShockIndex(maxShockinFoam), maxShockinFoam));

%Find the gradient around the shock radius. We see that the gradient is
%positive at the shock front, but is negative once it reaches the front and
%decays through material. Find when gradient first goes negative behind
%shock to define the shock front. We use shock as foot of the shock, and
%shock front as the peak - what we actually want to measure.
ShockPressureGradient = diff(Pressure,1,2);

ShockPressureSecondGradient = diff(ShockPressureGradient,1,2);
for i = minShockinQuartz:maxShockinFoam
l = find(flipud(ShockPressureGradient(1:ShockIndex(i),i)) < 0, 1);

ShockFrontIndex(i) = ShockIndex(i)-l;
ShockFrontRadius(i) = Radius(ShockFrontIndex(i), i);
end
ShockFrontVelocity = (gradient(ShockFrontRadius)./gradient(Time(1:length(ShockFrontRadius)).')/100000);

%We find that preheating of the gold layer causes expansion, and this
%throws off the shock tracking as to when the shock enters the quartz.
%However, we find that we can identify this more accurately using inverse
%density scale length. We can identify for each zone the time at which this
%parameter is maximum, which correlates well with the primary shock.
InverseDensityScaleLength = abs( diff(log(Density))./diff(Radius(1:end-1:end, :)));
InverseDensityScaleLength(isnan(InverseDensityScaleLength)) = 0.001;
InverseDensityScaleLength((InverseDensityScaleLength==0)) = 0.001;
InverseDensityScaleLength(isinf(InverseDensityScaleLength)) = 0.001; 
Zones = repmat([1:size(Radius,1)-2].', 1, length(Time));
[~, ScaleLengthMax] = max(InverseDensityScaleLength.');
figure
plot(Zones(:,1), ScaleLengthMax )
ShockInQuartzFromScaleLength = ScaleLengthMax(minQuartzIndex);
minShockinQuartz = ShockInQuartzFromScaleLength;

%Calculate time in foam from breakout times, and use this to calculate the
%measured velocity from breakout times. Then, figure out the intercept for
%y=mx+c so we can plot this later
TimeInFoam = Time(maxShockinFoam) - Time(minShockinFoam);
MeasuredFoamShockVelocity = (FoamMaxRadius - FoamMinRadius)*10^-5/TimeInFoam ;
MeasuredFoamTrajectoryIntercept = FoamMinRadius - ((MeasuredFoamShockVelocity/10^-5)*Time(minShockinFoam));

TimeInQuartz = Time(maxShockinQuartz) - Time(minShockinQuartz);
MeasuredQuartzShockVelocity = (QuartzMaxRadius - QuartzMinRadius)*10^-5/TimeInQuartz ;
MeasuredQuartzTrajectoryIntercept = QuartzMinRadius - ((MeasuredQuartzShockVelocity/10^-5)*Time(minShockinQuartz));

%Calculate particle velocity at each time by averaging velocity in all
%zones in region behind the shock front. Same for pressure (not quite
%constant behind shock, so approximating here). The NaN data point is
%because the slight boundary motion means that tracking the regions by
%radius doesn't quite work, and so when the shock first enters the foam
%it's at a lower radius than the minFoamIndex radius.
for j=1:minShockinQuartz-1
     IonTempTracked(j) = 0; 
    ShockFrontParticleVelocity(j)=0;
    ShockFrontPressure(j) = 0;
    ShockFrontDensity(j) = 0;
   
end
for j=minShockinQuartz:maxShockinQuartz
    IonTempTracked(j) = mean(IonTemp(minQuartzIndex:ShockIndex(j), j));
    MassAveragedTempTracked(j) = sum(Mass(minQuartzIndex:ShockIndex(j), j).*IonTemp(minQuartzIndex:ShockIndex(j), j))./sum(Mass(minQuartzIndex:ShockIndex(j), j));
    
    ShockFrontParticleVelocity(j)=Velocity(ShockFrontIndex(j),j);
    ShockFrontPressure(j) = Pressure(ShockFrontIndex(j),j);
    ShockFrontDensity(j) = Density(ShockFrontIndex(j),j);
end
for j=maxShockinQuartz+1:maxShockinFoam
    IonTempTracked(j) = mean(IonTemp(minFoamIndex:ShockIndex(j), j));
    MassAveragedTempTracked(j) = sum(Mass(minFoamIndex:ShockIndex(j), j).*IonTemp(minFoamIndex:ShockIndex(j), j))./sum(Mass(minFoamIndex:ShockIndex(j), j));
    ShockFrontParticleVelocity(j)=Velocity(ShockFrontIndex(j),j);
    ShockFrontPressure(j) = Pressure(ShockFrontIndex(j),j);
    ShockFrontDensity(j) = Density(ShockFrontIndex(j),j);
end

for j=maxShockinFoam+1:length(Time)
    IonTempTracked(j) = mean(IonTemp(minFoamIndex:size(IonTemp,1), j));
    MassAveragedTempTracked(j) = sum(Mass(minFoamIndex:size(IonTemp,1), j).*IonTemp(minFoamIndex:size(IonTemp,1), j))./sum(Mass(minFoamIndex:size(IonTemp,1), j));
end

%Ensure that the shock front is in the right material. If not, this will
%screw up the density calculation (and the other calculations too). Starts
%from above intervals -1, so that the first value in the quartz
%(minShockinQuartz) = 1.
ShockFrontValidQuartz = find(ShockFrontIndex(1, minShockinQuartz:maxShockinQuartz)>minQuartzIndex)+minShockinQuartz-1;
ShockFrontValidFoam = find(ShockFrontIndex(1, minShockinFoam:maxShockinFoam)>minFoamIndex)+maxShockinQuartz;

%Calculate length of each time interval (so that we can take time averages)
TimeSpacing = diff(Time);

%Take time averages to find mean of key variables
AverageQuartzShockFrontPressure = sum(TimeSpacing(ShockFrontValidQuartz).*ShockFrontPressure(ShockFrontValidQuartz).')./sum(TimeSpacing(ShockFrontValidQuartz)).*10^-3;;
AverageQuartzShockFrontParticleVelocity = sum(TimeSpacing(ShockFrontValidQuartz).*ShockFrontParticleVelocity(ShockFrontValidQuartz).')./sum(TimeSpacing(ShockFrontValidQuartz));
AverageQuartzShockFrontDensity = sum(TimeSpacing(ShockFrontValidQuartz).*ShockFrontDensity(ShockFrontValidQuartz).')./sum(TimeSpacing(ShockFrontValidQuartz));
AverageFoamShockFrontPressure = sum(TimeSpacing(ShockFrontValidFoam).*ShockFrontPressure(ShockFrontValidFoam).')./sum(TimeSpacing(ShockFrontValidFoam)).*10^-3;;
AverageFoamShockFrontParticleVelocity = sum(TimeSpacing(ShockFrontValidFoam).*ShockFrontParticleVelocity(ShockFrontValidFoam).')./sum(TimeSpacing(ShockFrontValidFoam));
AverageFoamShockFrontDensity = sum(TimeSpacing(ShockFrontValidFoam).*ShockFrontDensity(ShockFrontValidFoam).')./sum(TimeSpacing(ShockFrontValidFoam));



%Fit a linear curve to the shock front radius vs time for both quartz and
%foam, and use this to determine shock velocities.
[Fit1Variables,~] = createFit(Time(ShockFrontValidQuartz), ShockFrontRadius(ShockFrontValidQuartz));
Fit1Coeffs = coeffvalues(Fit1Variables);
AverageQuartzShockFrontVelocity = Fit1Coeffs(1)/100000;
[Fit2Variables,~] = createFit(Time(ShockFrontValidFoam), ShockFrontRadius(ShockFrontValidFoam));
Fit2Coeffs = coeffvalues(Fit2Variables);
AverageFoamShockFrontVelocity = Fit2Coeffs(1)/100000;

%Split the foam and quartz time periods into intervals, so we can generate
%profiles at each interval
NumofIntervals = 8;
FoamTimeIntervals = round((maxShockinFoam - minShockinFoam) / NumofIntervals);
QuartzTimeIntervals = round((maxShockinQuartz - minShockinQuartz) / NumofIntervals);
for i = 1:NumofIntervals
    FoamTimeIndex(i) = minShockinFoam + round(FoamTimeIntervals/2) + (i-1)*FoamTimeIntervals;
    QuartzTimeIndex(i) = minShockinQuartz + round(QuartzTimeIntervals/2) + (i-1)*QuartzTimeIntervals;
end



%Define Hugoniot for Quartz using data from Knudson and Desjarlis
%Fitting parameters
a = 6.26;
b = 1.20;
c = 2.56;
d = 0.37;
%Particle velocities for curve
HugoniotUp = -25:0.001:25;
%Principle Hugoniot
HugoniotUs = a + (b.*HugoniotUp) - (c.*HugoniotUp.*exp(-d.*HugoniotUp));

%Define a window for which the data is valid - this must cover the particle
%velocity of the quartz!!! Using quartz shock velocity, identify quartz
%particle velocity from Hugoniot
[~, MinHugoniotWindowIndex] = min(abs(HugoniotUp-4));
[~, MaxHugoniotWindowIndex] = min(abs(HugoniotUp-20));
[ QuartzApproximateUs, QuartzApproximateUsIndex ] = min( abs( HugoniotUs(MinHugoniotWindowIndex:MaxHugoniotWindowIndex)-AverageQuartzShockFrontVelocity ) );
QuartzApproximateUsIndex = QuartzApproximateUsIndex+ MinHugoniotWindowIndex;
QuartzApproximateUs = HugoniotUs(QuartzApproximateUsIndex);
QuartzApproximateUp = HugoniotUp(QuartzApproximateUsIndex);


%Reflect Hugoniot around the particle velocity of the Quartz for release
%isentrope
%Particle velocities
ReflectedHugoniotUp = -(HugoniotUp-(2*QuartzApproximateUp));
%Reflected Hugoniot
ReflectedHugoniotUs = a + (b.*ReflectedHugoniotUp) - (c.*ReflectedHugoniotUp.*exp(-d.*ReflectedHugoniotUp));
%Express in terms of pressure
HugoniotPressure = 2.6 * HugoniotUs .* HugoniotUp;
ReflectedHugoniotPressure = 2.6 * ReflectedHugoniotUs .* -(HugoniotUp-(2*QuartzApproximateUp));

if Measured==0
%Plot Rayleigh line for Foam. Find intercept, and thus foam parameters
RayleighLine = FoamDensity * AverageFoamShockFrontVelocity * HugoniotUp;
[~, MinHugoniotWindowIndex] = min(abs(HugoniotUp-4));
[~, MaxHugoniotWindowIndex] = min(abs(HugoniotUp-20));
[~, FoamHugoniotIndex] = min(abs(ReflectedHugoniotPressure - RayleighLine));
FoamApproximatePressure = RayleighLine(FoamHugoniotIndex);
FoamApproximateUp = HugoniotUp(FoamHugoniotIndex);
FoamApproximateRho = FoamDensity / (1 - FoamApproximateUp / AverageFoamShockFrontVelocity);
QuartzApproximatePressure = HugoniotPressure(QuartzApproximateUsIndex);
QuartzApproximateRho = 2.6 / (1 - QuartzApproximateUp / AverageQuartzShockFrontVelocity);
elseif Measured==1
%Plot Rayleigh line for Foam. Find intercept, and thus foam parameters
RayleighLine = FoamDensity * MeasuredFoamShockVelocity * HugoniotUp;
[~, MinHugoniotWindowIndex] = min(abs(HugoniotUp-4));
[~, MaxHugoniotWindowIndex] = min(abs(HugoniotUp-20));
[~, FoamHugoniotIndex] = min(abs(ReflectedHugoniotPressure - RayleighLine));
FoamApproximatePressure = RayleighLine(FoamHugoniotIndex);
FoamApproximateUp = HugoniotUp(FoamHugoniotIndex);
FoamApproximateRho = FoamDensity / (1 - FoamApproximateUp / MeasuredFoamShockVelocity);
QuartzApproximatePressure = HugoniotPressure(QuartzApproximateUsIndex);
QuartzApproximateRho = 2.6 / (1 - QuartzApproximateUp / AverageQuartzShockFrontVelocity);
else
    error('Use "Measured" to select a velocity to use')
end



Frequency1 =  3E8/(450E-9); % Frequency of laser light
CriticalDensity1 = (Frequency1/9)^2; % Electron density per m-3
CriticalDensity1 = CriticalDensity1/10^6 % Electron density per cm3
CritDensityThreshold = ElectronDensity>CriticalDensity1;
CritDensityIndex = arrayfun(@(x)find(CritDensityThreshold(:,x),1,'last'),1:size(CritDensityThreshold,2), 'UniformOutput',false);
tf = cellfun('isempty',CritDensityIndex); % true for empty cells
CritDensityIndex(tf) = {1};
CritDensityIndex= [CritDensityIndex{:}];
CritDensityIndexReshaped=sub2ind(size(Radius),CritDensityIndex,1:length(Time));
CritDensityRadius=Radius(CritDensityIndexReshaped);

Frequency2 =  3E8/(530E-9); % Frequency of laser light
CriticalDensity2 = (Frequency2/9)^2; % Electron density per m-3
CriticalDensity2 = CriticalDensity2/10^6 % Electron density per cm3
CritDensityThreshold2 = ElectronDensity>CriticalDensity2;
CritDensityIndex2 = arrayfun(@(x)find(CritDensityThreshold2(:,x),1,'first'),1:size(CritDensityThreshold2,2), 'UniformOutput',false);
tf = cellfun('isempty',CritDensityIndex2); % true for empty cells
CritDensityIndex2(tf) = {1};
CritDensityIndex2= [CritDensityIndex2{:}];
CritDensityIndexReshaped2=sub2ind(size(Radius),CritDensityIndex2,1:length(Time));
CritDensityRadius2=Radius(CritDensityIndexReshaped2);


if figures==1

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
try
xlim([0 Time(round(maxShockinFoam*1.1))*10^9]);
ylim([0 1.1*FoamMaxRadius.*10000])
catch
    xlim([0 Time(round(maxShockinFoam))*10^9]);
ylim([0 FoamMaxRadius.*10000])
end
hold on
plot3(Time(1:length(ShockFrontRadius))*10^9,ShockFrontRadius.*10000, max(Density(:, 1:length(ShockFrontRadius))), 'b');
plot3(Time(1:length(ShockRadius))*10^9,ShockRadius.*10000, max(Density(:, 1:length(ShockRadius))), 'w');
hold off
xline(Time(minShockinFoam)*10^9);
xline(Time(minShockinQuartz)*10^9);
xline(Time(maxShockinFoam)*10^9);

figure;
DensityPlot = surf(Time.*10^9,Radius(1:end-1,:).*10000, IonTempKelvin);
try
xlim([0 Time(round(maxShockinFoam*1.1))*10^9]);
ylim([0 1.1*FoamMaxRadius.*10000])
catch
    xlim([0 Time(round(maxShockinFoam))*10^9]);
ylim([0 FoamMaxRadius.*10000])
end
colorbar;
colormap(jet);
title('Ion Temp Plot');
xlabel('Time (ns)');
ylabel('Radius (cm)');
zlabel('Mass Density (g/cm^3)');
shading interp
set(gca,'colorscale','log')
view(2)

figure;
DensityPlot = surf(Time.*10^9,Radius(1:end-1,:).*10000, IonTempKelvin);
try
xlim([0 Time(round(maxShockinFoam*1.1))*10^9]);
ylim([0 1.1*FoamMaxRadius.*10000])
catch
    xlim([0 Time(round(maxShockinFoam))*10^9]);
ylim([0 FoamMaxRadius.*10000])
end
colorbar;
colormap(jet);
title('Ion Temp Plot');
xlabel('Time (ns)');
ylabel('Radius (cm)');
zlabel('Mass Density (g/cm^3)');
shading interp
set(gca,'colorscale','log')
view(2)
caxis([270 1500])


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
plot(Time, TotalLaserPower)
yyaxis right
plot(Time, AdjustedEnergy)
xline(Time(maxShockinFoam))
yline(AdjustedEnergy(maxShockinFoam))


figure
plot(Time*10^9, IonTemp(end,:))
xlabel('Time (ns)')
ylabel('Rear zone temperature (eV)')

figure
plot(Time*10^9, IonTempKelvin(end,:))
xlabel('Time (ns)')
ylabel('Rear zone temperature (K)')


figure
title('Breakout Temperature')
yyaxis left
plot(IonTemp(minQuartzIndex:end, maxShockinFoam))
ylabel('Temperatures at time of foam break out')
ylim([0 7])
yyaxis right
plot(IonTemp(minQuartzIndex:end, maxShockinQuartz))
ylabel('Temperatures at time of Al break out')
ylim([0 7])
xlabel('Zone Index')

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
   
figure
subplot(3,2,1)
plot(Time(ShockFrontValidQuartz).*10^9, ShockFrontPressure(ShockFrontValidQuartz).*10^-3)
yline(AverageQuartzShockFrontPressure, 'r');
yline(QuartzApproximatePressure, '--');
ylabel('Pressure (GPa)')
xlabel('Time (ns)')
title('Quartz')

subplot(3,2,3)
plot(Time(ShockFrontValidQuartz).*10^9, ShockFrontParticleVelocity(ShockFrontValidQuartz))
yline(AverageQuartzShockFrontParticleVelocity, 'r');
yline(QuartzApproximateUp, '--');
ylabel('Velocity (km/s)')
xlabel('Time (ns)')

subplot(3,2,5)
plot(Time(ShockFrontValidQuartz).*10^9, ShockFrontDensity(ShockFrontValidQuartz))
yline(AverageQuartzShockFrontDensity, 'r');
yline(QuartzApproximateRho, '--');
ylabel('Density (g/cm^3)')
xlabel('Time (ns)')

subplot(3,2,2)
plot(Time(ShockFrontValidFoam).*10^9, ShockFrontPressure(ShockFrontValidFoam).*10^-3)
yline(AverageFoamShockFrontPressure, 'r');
yline(FoamApproximatePressure, '--');
ylabel('Pressure (GPa)')
xlabel('Time (ns)')
title('Foam')

subplot(3,2,4)
plot(Time(ShockFrontValidFoam).*10^9, ShockFrontParticleVelocity(ShockFrontValidFoam))
yline(AverageFoamShockFrontParticleVelocity, 'r');
yline(FoamApproximateUp, '--');
ylabel('Velocity (km/s)')
xlabel('Time (ns)')

subplot(3,2,6)
plot(Time(ShockFrontValidFoam).*10^9, ShockFrontDensity(ShockFrontValidFoam))
yline(AverageFoamShockFrontDensity, 'r');
yline(FoamApproximateRho, '--');
ylabel('Density (g/cm^3)')
xlabel('Time (ns)')





figure
subplot(3,2,1)
plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)).*10000, Pressure(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)).*10^-3)
hold on
for i = 2:NumofIntervals
    plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)).*10000, Pressure(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)).*10^-3)
end
yline(AverageQuartzShockFrontPressure, 'r');
yline(QuartzApproximatePressure, '--');
hold off
ylabel('Pressure (GPa)')
xlabel('Distance (\mu m)')
title('Quartz')

subplot(3,2,3)
plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)).*10000, Velocity(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)))
hold on
for i = 2:NumofIntervals
    plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)).*10000, Velocity(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)))
end
yline(AverageQuartzShockFrontParticleVelocity, 'r');
yline(QuartzApproximateUp, '--');
hold off
ylabel('Velocity (km/s)')
xlabel('Distance (\mu m)')

subplot(3,2,5)
plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)).*10000, Density(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)))
hold on
for i = 2:NumofIntervals
    plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)).*10000, Density(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)))
end
yline(AverageQuartzShockFrontDensity, 'r');
yline(QuartzApproximateRho, '--');
hold off
ylabel('Density (g/cm^3)')
xlabel('Distance (\mu m)')


subplot(3,2,2)
plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)).*10000, Pressure(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)).*10^-3)
hold on
for i = 2:NumofIntervals
    plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)).*10000, Pressure(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)).*10^-3)
end
yline(AverageFoamShockFrontPressure, 'r');
yline(FoamApproximatePressure, '--');
hold off
ylabel('Pressure (GPa)')
xlabel('Distance (\mu m)')
title('Foam')

subplot(3,2,4)
plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)).*10000, Velocity(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)))
hold on
for i = 2:NumofIntervals
    plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)).*10000, Velocity(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)))
end
yline(AverageFoamShockFrontParticleVelocity, 'r');
yline(FoamApproximateUp, '--');
hold off
ylabel('Velocity (km/s)')
xlabel('Distance (\mu m)')

subplot(3,2,6)
plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)).*10000, Density(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)))
hold on
for i = 2:NumofIntervals
    plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)).*10000, Density(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)))
end
yline(AverageFoamShockFrontDensity, 'r');
yline(FoamApproximateRho, '--');
hold off
ylabel('Density (g/cm^3)')
xlabel('Distance (\mu m)')


figure
subplot(3,2,1)
plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)).*10000,  ElectronDensity(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)))
xline(CritDensityRadius(QuartzTimeIndex(1)).*10000);
hold on
for i = 2:NumofIntervals
plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)).*10000,  ElectronDensity(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)))
xline(CritDensityRadius(QuartzTimeIndex(i)).*10000);
end
yline(CriticalDensity1, '--');
hold off
ylabel('Electron Density (cm^{-3})')
xlabel('Distance (\mu m)')
title('Quartz')

subplot(3,2,2)
plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)).*10000,  ElectronDensity(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)))
hold on
xline(CritDensityRadius(FoamTimeIndex(1)).*10000);
for i = 2:NumofIntervals
plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)).*10000,  ElectronDensity(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)))
xline(CritDensityRadius(FoamTimeIndex(i)).*10000);
end
yline(CriticalDensity1, '--');
hold off
ylabel('Electron Density (cm^{-3})')
xlabel('Distance (\mu m)')
title('Foam')

subplot(3,2,3)
plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)).*10000,  IonTemp(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(1)))
xline(CritDensityRadius(QuartzTimeIndex(1)).*10000);
hold on
for i = 2:NumofIntervals
plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)).*10000,  IonTemp(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(i)))
xline(CritDensityRadius(QuartzTimeIndex(i)).*10000);
end
hold off
ylabel('Ion Temp (eV)')
xlabel('Distance (\mu m)')

subplot(3,2,4)
plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)).*10000,  IonTemp(minFoamIndex:maxFoamIndex, FoamTimeIndex(1)))
hold on
xline(CritDensityRadius(FoamTimeIndex(1)).*10000);
for i = 2:NumofIntervals
plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)).*10000,  IonTemp(minFoamIndex:maxFoamIndex, FoamTimeIndex(i)))
xline(CritDensityRadius(FoamTimeIndex(i)).*10000);
end
hold off
ylabel('Ion Temp (eV)')
xlabel('Distance (\mu m)')

subplot(3,2,5)
plot(Radius(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(4)).*10000,  IonTemp(minQuartzIndex:maxQuartzIndex, QuartzTimeIndex(4)))
xline(CritDensityRadius(QuartzTimeIndex(4)).*10000);
yline(IonTemp(CritDensityIndex(QuartzTimeIndex(4)), QuartzTimeIndex(4)));
ylabel('Ion Temp (eV)')
xlabel('Distance (\mu m)')
xlim([CritDensityRadius(QuartzTimeIndex(4)).*10000*0.99 CritDensityRadius(QuartzTimeIndex(4)).*10000*1.01])

subplot(3,2,6)
plot(Radius(minFoamIndex:maxFoamIndex, FoamTimeIndex(4)).*10000,  IonTemp(minFoamIndex:maxFoamIndex, FoamTimeIndex(4)))
xline(CritDensityRadius(FoamTimeIndex(4)).*10000);
yline(IonTemp(CritDensityIndex(FoamTimeIndex(4)), FoamTimeIndex(4)));
ylabel('Ion Temp (eV)')
xlabel('Distance (\mu m)')
xlim([CritDensityRadius(FoamTimeIndex(4)).*10000*0.99 CritDensityRadius(FoamTimeIndex(4)).*10000*1.01])



figure
plot(HugoniotUp, HugoniotPressure);
xlabel('Particle Velocity Up (km/s)')
ylabel('Pressure (GPa)')
hold on
plot(HugoniotUp, ReflectedHugoniotPressure);
plot(HugoniotUp, RayleighLine);
scatter(QuartzApproximateUp, QuartzApproximatePressure, 'k')
scatter(FoamApproximateUp, FoamApproximatePressure, 'r')
hold off
xline(HugoniotUp(MaxHugoniotWindowIndex), '--k');
xline(HugoniotUp(MinHugoniotWindowIndex), '--k');
xline(ReflectedHugoniotUp(MinHugoniotWindowIndex), '--r');
xline(ReflectedHugoniotUp(MaxHugoniotWindowIndex), '--r');
xlim([HugoniotUp(MinHugoniotWindowIndex)-2 ReflectedHugoniotUp(MinHugoniotWindowIndex)+2])
legend('Principle', 'Reflected', 'Rayleigh Line', 'Quartz shocked', 'Foam shocked', 'Quartz validity range', 'Quartz validity range', 'Foam validity range', 'Foam validity range')


figure
subplot(2,1,1)
scatter(Time(minShockinQuartz:maxShockinQuartz).*10^9, ShockFrontRadius(minShockinQuartz:maxShockinQuartz).*10000, 20, 'filled', 'r');
hold on
scatter(Time(ShockFrontValidQuartz).*10^9, ShockFrontRadius(ShockFrontValidQuartz).*10000, 20, 'filled', 'b');
plot(Time(ShockFrontValidQuartz).*10^9, ((Fit1Coeffs(1)*Time(ShockFrontValidQuartz))+ Fit1Coeffs(2)).*10000, 'g', 'LineWidth', 2);
plot(Time(minShockinQuartz:maxShockinQuartz).*10^9, (((MeasuredQuartzShockVelocity/10^-5)*Time(minShockinQuartz:maxShockinQuartz)+ MeasuredQuartzTrajectoryIntercept).*10000),'--k', 'LineWidth', 2);
hold off
xlabel('Time (ns)')
ylabel('Radius (\mu m)')
title('Quartz Shock Trajectory')
legend('Shock Tracking', 'Used for calculation', 'Fit', 'Measured from breakout', 'Location', 'northwest')
txt = ['Fit: ' num2str(AverageQuartzShockFrontVelocity) ' km/s, Measured: ' num2str(MeasuredQuartzShockVelocity) ' km/s'];
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
% text(2.6,50,txt)
text(SE(1), SE(2), txt, 'VerticalAlignment','bottom', 'HorizontalAlignment','right')

subplot(2,1,2)
scatter(Time(minShockinFoam:maxShockinFoam).*10^9, ShockFrontRadius(minShockinFoam:maxShockinFoam).*10000, 20, 'filled', 'r');
hold on
scatter(Time(ShockFrontValidFoam).*10^9, ShockFrontRadius(ShockFrontValidFoam).*10000, 20, 'filled', 'b');
plot(Time(ShockFrontValidFoam).*10^9, ((Fit2Coeffs(1)*Time(ShockFrontValidFoam))+ Fit2Coeffs(2)).*10000, 'g', 'LineWidth', 2);
% xline(Time(minShockinFoamforCalculation)*10^9, 'r');
plot(Time(minShockinFoam:maxShockinFoam).*10^9, (((MeasuredFoamShockVelocity/10^-5)*Time(minShockinFoam:maxShockinFoam)+ MeasuredFoamTrajectoryIntercept).*10000),'--k', 'LineWidth', 2);
hold off
xlabel('Time (ns)')
ylabel('Radius (\mu m)')
title('Foam Shock Trajectory')
legend('Shock Tracking', 'Used for calculation', 'Fit', 'Measured from breakout', 'Location', 'northwest')
txt = ['Fit: ' num2str(AverageFoamShockFrontVelocity) ' km/s, Measured: ' num2str(MeasuredFoamShockVelocity) ' km/s'];
SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
% text(2.6,50,txt)
text(SE(1), SE(2), txt, 'VerticalAlignment','bottom', 'HorizontalAlignment','right')
end


[FoamApproximateUp_RealIM, FoamApproximatePressure_RealIM, ~, ~, ~, ~, ~, QuartzApproximateUp_RealIM, QuartzApproximatePressure_RealIM] = ExactImpedanceMatch(MeasuredQuartzShockVelocity, MeasuredFoamShockVelocity)


% figure
% p=panel;
% p.pack(2,2)
% p(1,1).select()
% plot(Time(minTimeIndex:maxTimeIndex).*10^9, ShockFrontPressure(minTimeIndex:maxTimeIndex).*10^-3,  'Color', 'k', 'LineWidth',2)
% % yline(AverageQuartzShockFrontPressure, 'r');
% yline(QuartzApproximatePressure_RealIM, '--', 'Color', cm(11,:), 'LineWidth',2);
% ylabel('Pressure (GPa)')
% xlabel('Time (ns)')
% % ylim([0 400])
% % title('Quartz')
% box on
% p(2,1).select()
% plot(Time(minTimeIndex2:maxTimeIndex2).*10^9, ShockFrontPressure(minTimeIndex2:maxTimeIndex2).*10^-3,  'Color', 'k', 'LineWidth',2)
% % yline(AverageFoamShockFrontPressure, 'r');
% yline(FoamApproximatePressure_RealIM, '--', 'Color', cm(11,:), 'LineWidth',2);
% ylabel('Pressure (GPa)')
% % ylim([0 50])
% xlabel('Time (ns)')
% % title('Foam')
% box on
% p(1,2).select()
% plot(Time(minTimeIndex:maxTimeIndex).*10^9, ShockFrontParticleVelocity(minTimeIndex:maxTimeIndex),  'Color', 'k', 'LineWidth',2)
% % yline(AverageQuartzShockFrontParticleVelocity, 'r');
% yline(QuartzApproximateUp_RealIM, '--', 'Color', cm(11,:), 'LineWidth',2);
% ylabel('Particle vel. (km/s)')
% % ylim([0 10])
% xlabel('Time (ns)')
% box on
% p(2,2).select()
% plot(Time(minTimeIndex2:maxTimeIndex2).*10^9, ShockFrontParticleVelocity(minTimeIndex2:maxTimeIndex2),  'Color', 'k', 'LineWidth',2)
% % yline(AverageFoamShockFrontParticleVelocity, 'r');
% yline(FoamApproximateUp_RealIM, '--', 'Color', cm(11,:), 'LineWidth',2);
% ylabel('Particle vel. (km/s)')
% % ylim([0 15])
% xlabel('Time (ns)')
% box on
% p.de.margin = 6;
% disp(sprintf('p.margin is [ %i %i %i %i ]', p.margin));
% p.margin = [8 8 2 6];
% p(2).margintop = 14
% p(2,1).marginright = 10
% p(1,1).marginright = 10
% fig=gcf;
% fig.Units               = 'centimeters';
% fig.Position(3)         = 8.6
% fig.Position(4)         = 7;
% % and let's set the global font properties, also. we can do
% % this at any point, it doesn't have to be here.
% p.fontsize = 8;
% p.fontname = 'times';
% dim1 = [0.43,0.33,0.2,0.2];
% str1 = ['Foam'];
% annotation('textbox',dim1,'String',str1,'FitBoxToText','on',  'FontSize',     10, 'LineStyle', 'none', 'color', 'k',  'FontName',     'Times', 'HorizontalAlignment', 'center','FontWeight', 'bold');
% dim2 = [0.43,0.82,0.2,0.2];
% str2 = ['Quartz'];
% annotation('textbox',dim2,'String',str2,'FitBoxToText','on',  'FontSize',     10, 'LineStyle', 'none', 'color', 'k',  'FontName',     'Times', 'HorizontalAlignment', 'center','FontWeight', 'bold');



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


