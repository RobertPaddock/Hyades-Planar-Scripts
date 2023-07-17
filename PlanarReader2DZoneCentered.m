[FileName,PathName]=uigetfile('*.cdf','Select the INPUT DATA FILE(s)','MultiSelect','on');

Finished = 1;

for FileIndex=1:length(FileName)
    file=[char(PathName) char(FileName(FileIndex))]

Time(FileIndex) = ncreadatt(file,'/','Time');

%Import the key variables, and get in correct format. Do this in each
%sequence for each time file.
Radius(:,:,FileIndex)  = ncread(file,'R');
Axial(:,:,FileIndex)  = ncread(file,'Z');
Density(:,:,FileIndex)   = ncread(file,'Rho'); %g/cm^3
ElecTemp(:,:,FileIndex)  = ncread(file,'Te')*1000; %Converted to eV from keV
IonTemp(:,:,FileIndex)  = ncread(file,'Ti')*1000; %Converted to eV from keV
Pressure(:,:,FileIndex) = ncread(file,'Pres').*0.0000001; %Converted to J/cm^3 from dyn/cm^2
RadialVelocity(:,:,FileIndex) = ncread(file, 'Ur')./100000; %Converted to km/s from cm/s
AxialVelocity(:,:,FileIndex) = ncread(file, 'Uz')./100000; %Converted to km/s from cm/s
% IonThermalEnergy  = ncread(file,'Eion').*10^-7; %Ion Thermal Energy, converted from erg to J.
% ElectronThermalEnergy  = ncread(file,'Eelc').*10^-7; %Electron Thermal Energy, converted from erg to J.
% %KineticEnergy  = ncread(file,'Ekint').*10^-7;
DepositedLaserPowerZones(:,:,FileIndex) = ncread(file,'Deplas').*10^-7; %Deposited Laser Energy per zone, converted from erg to J.
DepositedLaserPower(FileIndex) = sum(sum(DepositedLaserPowerZones(:,:,FileIndex)));
Volume(:,:,FileIndex) =  ncread(file,'Vol'); %cm^3
Mass(:,:,FileIndex) = ncread(file,'Xmass');
LaserIntensityZones(:,:,FileIndex) = ncread(file,'Xlsint').*10^-7; %Laser Intensity per zone, converted from erg to J, and converted to account for instability.

IonThermalEnergy(:,:,FileIndex)  = ncread(file,'Eion').*10^-7; %Ion Thermal Energy, converted from erg to J.
ElectronThermalEnergy(:,:,FileIndex)  = ncread(file,'Eelc').*10^-7; %Electron Thermal Energy, converted from erg to J.
end

%Because these parameters are zonecentered, the first zone is given as a 0
%in each direction. So, remove the first row and column from these arrays.
Density(1, :, :) = []; Density(:, 1, :) = [];
ElecTemp(1, :, :) = []; ElecTemp(:, 1, :) = [];
IonTemp(1, :, :) = []; IonTemp(:, 1, :) = [];
Pressure(1, :, :) = []; Pressure(:, 1, :) = [];
Volume(1, :, :) = []; Volume(:, 1, :) = [];
% Mass(1, :, :) = []; Mass(:, 1, :) = [];
DepositedLaserPowerZones(1, :, :) = []; DepositedLaserPowerZones(:, 1, :) = [];
LaserIntensityZones(1, :, :) = []; LaserIntensityZones(:, 1, :) = [];

IonThermalEnergy(1, :, :) = []; IonThermalEnergy(:, 1, :) = [];
ElectronThermalEnergy(1, :, :) = []; ElectronThermalEnergy(:, 1, :) = [];


%Calculate deposited laser power
DepositedEnergy = cumtrapz(Time, DepositedLaserPower);

if Finished==1
    
    %Identify bulk regions, and type of quartz/quartz density
BulkRegions = Radius(1:end-1,1,1)<0.0123;
if max(abs(Density(:,1,1) - 2.2)<0.01)
    QuartzDensity = 2.2;
    disp('Amorphous quartz, 2.2 density')
elseif max(abs(Density(:,1,1) - 2.6)<0.01)
    QuartzDensity = 2.6;
    disp('Silica, 2.6 density')
    else
            QuartzDensity = 2.6;
    disp('Quartz type not identified')
end
    
%Identify the different regions based on density.
Quartzzones = (abs(Density(:,1,1) - QuartzDensity)<0.01);
Quartzzones = Quartzzones .* BulkRegions;
CHzones = abs(Density(:,1,1) - 0.253)<0.01;
CHzones = CHzones .* BulkRegions;
[~,minQuartzIndex] = find(Quartzzones.', 1, 'first');
[~,maxQuartzIndex] = find(Quartzzones.', 1, 'last');
[~,minFoamIndex] = find(CHzones.', 1, 'first');
[~,maxFoamIndex] = find(CHzones.', 1, 'last');

QuartzMaxRadius = Axial(maxQuartzIndex,1, 1);
QuartzMinRadius = Axial(minQuartzIndex,1, 1);
FoamMinRadius = Axial(minFoamIndex, 1,1);
FoamMaxRadius = Axial(maxFoamIndex, 1,1);
    
%Define shock front as last zone in increasing radius in the sim where
%there is a substantial difference between density now and original density
%- since density doesn't really change until shock reaches it. Doesn't work
%so well for gold/ start of al where there is preheat.
%Shock = (Density - Density(:,:,1);) > 0.1;
Shock = (Density ./ Density(:,:,1)) > 1.5;
BulkRegions = Axial(1:end-1,1:end-1,1)<FoamMaxRadius;
%Shock = Shock.*BulkRegions;

%For each time step, search to find the last axial position of the shock. Save
%index of this position. Take this axial position (which is in a 2D array
%corresponding to radial position and time, and convert into the correct
%index in radius/axial array for that position/time. Use this to find
%coordinates of the shock. Then, reshape this list of the coordinate so it corresponds
%to other coordinate and time.
for TimeIndex = 1:length(Time)
ShockIndex(:,TimeIndex)= arrayfun(@(x)find(Shock(:,x,TimeIndex),1,'last'),1:size(Shock,2), 'UniformOutput',false);
end
tf = cellfun('isempty',ShockIndex); % true for empty cells
ShockIndex(tf) = {1};
ShockIndex= [ShockIndex{:}];
ShockIndexReshaped=sub2ind(size(Radius),ShockIndex, repmat(1:size(Density,2), 1,size(Density,3)), repelem(1:length(Time),size(Density,2)));
ShockAxial=Axial(ShockIndexReshaped);
ShockAxial = reshape(ShockAxial, size(Density,2), size(Density,3));
ShockRadius=Radius(ShockIndexReshaped);
ShockRadius = reshape(ShockRadius, size(Density,2), size(Density,3));

%Once position of front has been identified, calculate the rate at which
%this moves through the target. Divide by 10^5 to change from cm/s to km/s.
ShockVelocity = (gradient(ShockAxial)./gradient(Time))/100000;
%SmoothShockVelocity = smoothdata(ShockVelocity, 'movmean', 10);


%Identify when the shock is in Al or Foam by comparing shock radius to
%above radii. Final condition ensures the shock has not yet reached the
%rear surface (as I had problems with the shock radius at times past this
%in some files).
ShockInQuartz=(QuartzMaxRadius>ShockAxial).*(ShockAxial>QuartzMinRadius);
ShockInFoam=(FoamMaxRadius>ShockAxial).*(ShockAxial>FoamMinRadius);

%Identify the indexes of when the shock enters and leaves the different
%regions.
for i = 1:size(ShockInQuartz,1)
    try
[~,minShockinQuartz(i)] = find(ShockInQuartz(i, :), 1, 'first');
[~,maxShockinQuartz(i)] = find(ShockInQuartz(i, :), 1, 'last');
[~,minShockinFoam(i)] = find(ShockInFoam(i, :), 1, 'first');
[~,maxShockinFoam(i)] = find(ShockInFoam(i, :), 1, 'last');
    catch 
minShockinQuartz(i) = 0;
maxShockinQuartz(i) = 0;
minShockinFoam(i) = 0;
maxShockinFoam(i) = 0;
    end
end



%Reshape Shock Index into the normal array format
ShockIndex2D = reshape(ShockIndex, size(Density,2), size(Density,3));

%Calculate pressure behind the shock front along each radial line within
%the laser spot, and at each time point.
for i = 1:size(ShockInQuartz,1)
for j=1:minShockinQuartz-1
    PressureTracked(i,j) = 0;
    PressureAtFrontTracked(i,j) = 0;
        PressureMaxTracked(i,j) = 0;

end
for j=minShockinQuartz:maxShockinQuartz
    PressureTracked(i,j) = mean(Pressure(minQuartzIndex:ShockIndex2D(i, j), i, j));
    PressureAtFrontTracked(i,j) = Pressure(ShockIndex2D(i,j), i,  j);
    try
    PressureMaxTracked(i,j) = max(Pressure(minQuartzIndex:ShockIndex2D(i, j), i, j));
    catch
       PressureMaxTracked(i,j) = nan; 
    end
    
end
for j=maxShockinQuartz+1:maxShockinFoam
    PressureTracked(i,j) = mean(Pressure(minFoamIndex:ShockIndex2D(i,j), i, j));
    PressureAtFrontTracked(i,j) = Pressure(ShockIndex2D(i,j),i,  j);
       try
           PressureMaxTracked(i,j) = max(Pressure(minFoamIndex:ShockIndex2D(i, j), i, j));
 catch
       PressureMaxTracked(i,j) = nan; 
    end
end
end

%For each radial line within the laser spot, smooth the shock velocity and
%then plot. Also take a single average, and find the time averaged velocity
%in the two regions.
[~, InLaserIndex] = max(Radius(1,:,1)>0.02);
for k = 1:InLaserIndex
    SmoothShockVelocity(k,:) = smoothdata(ShockVelocity(k,:), 'movmean', 30);
    SmoothPressure(k, :) = smoothdata(PressureTracked(k,:), 'movmean', 10);
end
AverageShockVelocity = mean(ShockVelocity(1:InLaserIndex, :));
TimeAveragedShockVelocityinFoam = mean(AverageShockVelocity(minShockinFoam:maxShockinFoam));
TimeAveragedShockVelocityinQuartz = mean(AverageShockVelocity(minShockinQuartz:maxShockinQuartz));
AveragePressure = mean(PressureTracked(1:InLaserIndex, :), 'omitnan');
TimeAveragedPressureinFoam = mean(AveragePressure(minShockinFoam:maxShockinFoam), 'omitnan');
TimeAveragedPressureinQuartz = mean(AveragePressure(minShockinQuartz:maxShockinQuartz), 'omitnan');

%Calculate time in each layer, and thus the average velocity within the
%layer
TimeInQuartz = Time(maxShockinQuartz(1:InLaserIndex)) - Time(minShockinQuartz(1:InLaserIndex))
TimeInFoam = Time(maxShockinFoam(1:InLaserIndex)) - Time(minShockinFoam(1:InLaserIndex))
AverageQuartzVelocity = (QuartzMaxRadius - QuartzMinRadius)./TimeInQuartz
AverageFoamVelocity = (FoamMaxRadius - FoamMinRadius)./TimeInFoam


%How do coordinates work? Coordinate 1 is the axial coordinate (length
%along the z axis - Z, or the k mesh point). Coordinate 2 is the radial
%coordinate (radius from axis, R, or l mesh point). Coordinate 3 is the
%time. Opening the variable, this means that the picture is 90 degrees
%rotated - the columns are axial position along z, the rows are radial
%position along R.

end

figure
title('Laser Power and input energy')
yyaxis left
plot(Time*10^9, DepositedLaserPower/10^12)
ylabel('Deposited Laser Power (TW)')
yyaxis right
plot(Time*10^9, DepositedEnergy)
ylabel('Deposited Energy (J)')
xlabel('Time (ns)')

figure;
DensityPlot = surf(Axial(1:end-1,1:end-1,4) ,Radius(1:end-1,1:end-1,4), Density(:,:,4));
colorbar;
colormap(jet);
xlabel('Z (cm)');
ylabel('R (cm)');
zlabel('Mass Density (g/cm^3)');
shading interp
set(gca,'colorscale','log')
view(2)
hold on
plot3(ShockAxial(:,4),ShockRadius(:,4), squeeze(max(Density(:,:,4))), '--w', 'LineWidth',1);
hold off
title('Density at a given time')




figure;
DensityPlot = surf(Axial(1:end-1,1:end-1,4) ,Radius(1:end-1,1:end-1,4), IonTemp(:,:,4));
colorbar;
colormap(jet);
xlabel('Z (cm)');
ylabel('R (cm)');
zlabel('Mass Density (g/cm^3)');
shading interp
view(2)
title('Ion Temp at a given time')


%I want how density changes as a function of time for a given R coodinate
%(i.e R = 0). This means I want to keep the R coordinate (a row, coordinate
%2) constant, and plot the density against time and Axial (since I'm
%interested in how the mesh moves along the z axis).
figure;
DensityPlot = surf(Time.*10^9 ,squeeze(Axial(1:end-1,1,:)).*10000, squeeze(Density(:,1,:)));
colorbar;
colormap(jet);
xlabel('Time (ns)');
ylabel('Axial position (\mum)');
zlabel('Mass Density (g/cm^3)');
shading interp
set(gca,'colorscale','log')
view(2)
hold on
plot3(Time.*10^9,(ShockAxial(1,:)).*10000, squeeze(max(Density(:,1,:))), 'w');
hold off
title('Density at a constant Radius as a fn of time')

% %Mass change moving along in z
% figure
% AxialMassDifference = 100*diff(Mass(:,1,1))./Mass(2:end,1,1);
% plot(AxialMassDifference)
% ylabel('Percentage Axial Mass Difference')
% yyaxis right;
% plot(Mass(:,1,1), ':');
% ylabel('Zone Mass')
% title('Axial Mass Difference')
% xlabel('Zone')
% 
% %Mass change moving along in r
% figure
% RadialMassDifference = 100*diff(Mass(1,:,1))./Mass(1,2:end,1);
% plot(RadialMassDifference)
% ylabel('Percentage Radial Mass Difference')
% yyaxis right;
% plot(Mass(1,:,1), ':');
% ylabel('Zone Mass')
% title('Radial Mass Difference')
% xlabel('Zone')

%Plot densities as sanity check
figure
AxialDensityDifference = 100*diff(Density(:,1,1))./Density(2:end,1,1);
plot(AxialDensityDifference)
ylabel('Percentage Axial Density Difference')
yyaxis right;
plot(Density(:,1,1), ':');
ylabel('Zone Density')
title('Axial Zone Density Plot')
xlabel('Zone')

figure
RadialDensityDifference = 100*diff(Density(1,:,1))./Density(1,2:end,1);
plot(RadialDensityDifference)
ylabel('Percentage Radial Mass Difference')
yyaxis right;
plot(Density(1,:,1), ':');
ylabel('Zone Density')
title('Radial Zone Density Plot')
xlabel('Zone')

%Plot zone thicknesses 
figure
AxialDifference = diff(Axial(:,1,1));
plot(AxialDifference)
ylabel('Axial Thickness')
title('Axial Zone Thickness')

figure
RadialDifference = diff(Radius(1,:,1));
plot(RadialDifference)
ylabel('Radial Thickness')
title('Radial Zone Thickness')

Area = Volume./AxialDifference;
Power = LaserIntensityZones.*Area;

%Plot shock velocity for each radial line
figure
plot(Time*10^9, ShockVelocity)
xlim([Time(minShockinQuartz(1)) Time(maxShockinFoam(1))]*10^9);
title('Shock velocity')
ylabel('Shock velocity (km/s)')
xlabel('Time (ns)')

figure;
DensityPlot = surf(Time(1:size(PressureTracked,2)).*10^9,ShockRadius(:, 1:size(PressureTracked,2)).*10000, PressureTracked.*10^-3);
colorbar;
colormap(jet);
title('Shock Pressure vs Time and Radius')
xlabel('Time (ns)');
ylabel('Radius (\mu m)')
zlabel('Pressure (GPa)');
shading interp
view(2)

figure;
DensityPlot = surf(Time(1:size(PressureTracked,2)).*10^9,ShockRadius(:, 1:size(PressureTracked,2)).*10000, PressureTracked.*10^-3);
colorbar;
colormap(jet);
title('Shock Pressure vs Time and Radius')
xlabel('Time (ns)');
ylabel('Radius (\mu m)')
zlabel('Pressure (GPa)');
shading flat
view(2)

figure;
DensityPlot = surf(Time(1:size(PressureMaxTracked,2)).*10^9,ShockRadius(:, 1:size(PressureMaxTracked,2)).*10000, PressureMaxTracked.*10^-3);
colorbar;
colormap(jet);
title('Shock Pressure vs Time and Radius')
xlabel('Time (ns)');
ylabel('Radius (\mu m)')
zlabel('Pressure (GPa)');
shading flat
view(2)
set(gca,'colorscale','log')


figure
plot(Time*10^9, SmoothShockVelocity)
hold on
yline(TimeAveragedShockVelocityinFoam);
yline(TimeAveragedShockVelocityinQuartz);
hold off
xlim([Time(minShockinQuartz(1)) Time(maxShockinFoam(1))]*10^9);
title('Shock velocity')
ylabel('Shock velocity (km/s)')
xlabel('Time (ns)')

figure
plot(Time*10^9, SmoothShockVelocity(1:InLaserIndex,:))
hold on
yline(TimeAveragedShockVelocityinFoam);
yline(TimeAveragedShockVelocityinQuartz);
hold off
xlim([Time(minShockinQuartz(1)) Time(maxShockinFoam(1))]*10^9);
title('Shock velocity')
ylabel('Shock velocity (km/s)')
xlabel('Time (ns)')

%Plot Pressures (10^-3 converts to GPa)
figure
plot(Time(1:size(SmoothPressure,2))*10^9, SmoothPressure.*10^-3)
hold on
yline(TimeAveragedPressureinFoam.*10^-3);
yline(TimeAveragedPressureinQuartz.*10^-3);
hold off
xlim([Time(minShockinQuartz(1)) Time(maxShockinFoam(1))]*10^9);
title('Pressure')
ylabel('Pressure (GPa)')
xlabel('Time (ns)')

%Plot Pressures (10^-3 converts to GPa)
figure
plot(Time(1:size(PressureMaxTracked,2))*10^9, PressureMaxTracked.*10^-3)
hold on
yline(TimeAveragedPressureinFoam.*10^-3);
yline(TimeAveragedPressureinQuartz.*10^-3);
hold off
xlim([Time(minShockinQuartz(1)) Time(maxShockinFoam(1))]*10^9);
title('Pressure')
ylabel('Pressure (GPa)')
xlabel('Time (ns)')

figure
plot(Time*10^9, AverageShockVelocity)
hold on
yline(TimeAveragedShockVelocityinFoam);
yline(TimeAveragedShockVelocityinQuartz);
hold off
xlim([Time(minShockinQuartz(1)) Time(maxShockinFoam(1))]*10^9);
title('Shock velocity')
ylabel('Shock velocity (km/s)')
xlabel('Time (ns)')

figure
plot(AverageQuartzVelocity)
hold on
plot(AverageFoamVelocity)
hold off
ylabel('Velocity')
legend('Measured Quartz', 'Measured Foam')

figure
scatter(Radius(1,1:InLaserIndex,1), Time(maxShockinQuartz(1:InLaserIndex)))
ylabel('Transition Time (ns)')
legend('Radius')
title('Quartz exit time')

figure
scatter(Radius(1,1:InLaserIndex,1), Time(maxShockinFoam(1:InLaserIndex)))
ylabel('Transition Time (ns)')
legend('Radius')
title('Foam exit time')

save([PathName, 'Data'])
