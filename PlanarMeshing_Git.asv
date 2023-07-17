  %Meshing script. This script meshes a 4 layer CH/Au/Quartz/Foam target,
  %with the dimensions listed below. The dimensions can be varied, but the
  %thicknesses listed below must also be varied. These will also need to be
  %varied if the densities are changed. The number of layers can be varied,
  %but this will require additional lines of code to inser the new layer in
  %the appropraite place, and mesh it as with the others. The meshing is
  %done accordign to the theroy described in my thesis.

close all
clear all

QuartzType = 1; %Quartz type, 1= alpha, 2 = amorphous

%Specify material dimensions
CHLength = 0.0040; 
GoldLength = 0.0001;
QuartzLength = 0.004; 
FoamLength = 0.0040; 
    
%Insert densities of the materials.
if QuartzType == 1
    DensityQuartz = 2.6; %Alpha quartz / silica density
else
DensityQuartz = 2.2; %Amorphous quartz density
end
DensityFoam=1.0000000e-01;
DensityCH = 1.040000e+00;
DensityGold=19.3;

%Specifies zone thicknesses in each material at the interfaces. This is key
%to the success of the meshing. The zone size at either edge of the target
%(CHThickness1 and FoamThickness2) should be fixed. Zones must be close in
%mass difference for accuracy, but big to reduce number of zones (and sim
%time). Here, we set the zone thickness for each material at each interface
%to accomplish this. The difference across a material should be sufficiently
%large to reduce number of zones, but if it is too big the mass difference
%may be too big, or the problem won't mesh (i.e. not possible to solve).
%The exact numbers are found through trial and error!
CHThickness1 = 0.000001;
CHThickness2 = 0.00008; 
GoldThickness1 = CHThickness2*DensityCH/DensityGold;
GoldThickness2 = GoldThickness1*1.5
QuartzThickness1 = GoldThickness2*DensityGold/DensityQuartz;
QuartzThickness2 = QuartzThickness1/15; 
FoamThickness1 = QuartzThickness2*DensityQuartz/DensityFoam; 
FoamThickness2 = 0.000001;

%We take the length of the material region, plus the desired thickness at
%each edge, and then use the theory from my thesis to attempt to calculate
%the number of layers and thickness ratio between subsequent layers that
%achieves this.
[RatioCH, LayersCH] = ThicknessSolver(0, CHLength, CHThickness1, CHThickness2);
[RatioGold, LayersGold] = ThicknessSolver(0, GoldLength, GoldThickness1, GoldThickness2);
[RatioQuartz, LayersQuartz] = ThicknessSolver(0, QuartzLength, QuartzThickness1, QuartzThickness2);
[RatioFoam, LayersFoam] = ThicknessSolver(0, FoamLength, FoamThickness1, FoamThickness2);

%Round number of layers so that it is a whole number
LayersCH = round(LayersCH);
LayersGold = round(LayersGold);
LayersQuartz = round(LayersQuartz);
LayersFoam = round(LayersFoam);

%Function takes these outputs, and determines the radii and density of each
%zone. Then, next function uses this to calculate the mass difference
%between each zone. The BoundarySolver function will need to be adapted if
%the number of regions is changed.
[Radii, Density] = BoundarySolver(CHLength, GoldLength, QuartzLength, FoamLength, RatioCH, RatioGold, RatioQuartz, RatioFoam, LayersCH,  LayersGold, LayersQuartz, LayersFoam, DensityCH, DensityGold, DensityQuartz, DensityFoam)
[ZoneMass, MassDiff] = MassDifference(Radii, Density);

        
      
    
%Calculate number of layers, mass difference
NumLayers = LayersCH+LayersGold+LayersQuartz+LayersFoam;
MaxDiff = max(abs(MassDiff(1:end)));
        
%Plot the zoning
figure
plot(MassDiff)
ylabel('Percentage Mass Difference')
yyaxis right;
plot(ZoneMass, ':');
ylabel('Zone Mass')
title('Zoning Plot')
xlabel('Zone')

        
        
      
    












%Function to recreate the Hyades mesh function - given the mesh boundaries
%j and radii to stretch between r and the ratio, it will create the mesh
%coordinates.
function Radii = RatioIncrement(j1, j2, r1, r2, Ratio)
Index = (1+j1-j1):(j2-j1);
Increment = Ratio.^(Index);
Radius = cumsum(Increment);
Scaling = (r2-r1)/Radius(end);
Radii = [r1+(Radius*Scaling)];
end

%Given the radii and density, will calculate the zone masses and return
%the percentage mass difference
function [ZoneMass, MassDiff] = MassDifference(Radii, Density)
ZoneMass = diff(Radii).*Density;
MassDiff = 100*(diff(ZoneMass)./ZoneMass(2:end));
end

%Calculates the mesh
function [Radii, Density] = BoundarySolver(CHLength, GoldLength, AlLength, FoamLength, RatioCH, RatioGold, RatioAl, RatioFoam, LayersCH,  LayersGold, LayersAl, LayersFoam, DensityCH, DensityGold, DensityAl, DensityFoam)

BoundaryCH = LayersCH+1;
BoundaryGold = BoundaryCH+LayersGold;
BoundaryAl = BoundaryGold+LayersAl;
BoundaryFoam = BoundaryAl + LayersFoam;

CHzonepositions = RatioIncrement(1, BoundaryCH, 0, CHLength, RatioCH);
Goldzonepositions = RatioIncrement(BoundaryCH, BoundaryGold, CHLength, CHLength+GoldLength, RatioGold);
Alzonepositions = RatioIncrement(BoundaryGold, BoundaryAl, CHLength+GoldLength, CHLength+GoldLength+AlLength, RatioAl);
Foamzonepositions = RatioIncrement(BoundaryAl, BoundaryFoam, CHLength+GoldLength+AlLength, CHLength+GoldLength+AlLength+FoamLength, RatioFoam);

Density = [ DensityCH*ones(1, LayersCH), DensityGold*ones(1, (LayersGold)), DensityAl*ones(1, (LayersAl)),  DensityFoam*ones(1, (LayersFoam))];

Radii = [0, CHzonepositions, Goldzonepositions, Alzonepositions, Foamzonepositions];


end




%Function to solve for optimal ratio and number of layers. Given the radii
%and desired thickness at either end of the region, this function will
%solve the two simultaneous equations to determine the optimal number of
%layers and ratio to use.
function [Ratio, Layers] = ThicknessSolver(LowerRadius, UpperRadius, LowerThickness, UpperThickness)
%Turn off numerical solver warning message
warning('off','symbolic:solve:FallbackToNumerical');

%Solve analytic formula to find ratio and layers for ice region .
syms a n
eqn1 = log(1 - (1-a)*(UpperRadius-LowerRadius)/LowerThickness)/log(a) == n;
eqn2 = a == (UpperThickness/LowerThickness)^(1/(n-1));
sol = solve([eqn1, eqn2], [a, n], 'Real', true);
aSol = sol.a;
nSol = sol.n;
Ratio=double(aSol);
Layers=double(nSol);

end

