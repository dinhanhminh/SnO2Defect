clear global;

%% GLOBAL VARIABLES
global k temperature dielectric_const E_MK ...
    dopingConcentration totalDopingDefects ...
    mu_O mu_Metal mu_H defect vbm cbm Eperf ...
    total_dos energy efermi cellsize dosCellsize;



%% CONSTANT
% CONSTANT
[ energy,total_dos,efermi,pdos ] = import_doscar('DOSCAR');
% vbm = energy(find(energy < efermi, 1, 'last' ));
vbm = 6.2704;
test = [energy total_dos];
% cbm = min(test(test(:,1) > efermi & test(:,2) > 0));
cbm = vbm + 3.6;
% defectFile = 'defects.txt';
%defectFile = 'defects.txt';

systemElement = ["O", "Sn", "H"];


%%
%%%%% MODULE SELECTION
% 1 = Kroger-Vink undoped
% 1.5 = Kroger-Vink H-doped
% 2 = Kroger-Vink with doping
% 3 = Scanned H vs F Doping level
% 4 = Formation vs mu_e at pO2 = 1atm
% 5 = Scan [C] of H at 1 atm across different element
% 6 = Scan [C] of H at 1 atm across different mu_E

module = 2; 

switch module
    case 1
        defectFile = "no-Hydrogen-bigCell.txt";
        outputFile = "K-V-no-H.txt";
        element = systemElement(1:2);
    case 1.5
        defectFile = "defects-bigCell.txt";
        outputFile = "K-V-H.txt";
        element = systemElement;
%    case 2
    case 4
        defectFile = "no-Hydrogen.txt";
        element = systemElement(1:2);
    otherwise
        defectElement = "F";
        element = [systemElement, defectElement];
        defectFile = 'defects-'+defectElement+'.txt';
        outputFile = 'K-V-'+defectElement+'.txt' ;
end

%%
%%%% INPUT %%%


latexSymbolCharge = ["'","\times","\bullet"];

allDefectElement = ["In", "Mn", "Nb", "Pd", "Rh", "Sb", "Ta", "Ti", "Zr"];



%%% test
% defectFile = "defect.txt";        
defect = importdata(defectFile);
defect = defect.data;
dopingConcentration = 10^-4 ;
totalDopingDefects = find(defect(:,5)>0, 1, 'last' )-find(defect(:,5)>0, 1 )+1;
if isempty(totalDopingDefects)
    totalDopingDefects = 0;
end

k = 8.6173303*10^-5; %eV/K
temperature = 1200; % K
dielectric_const = 24;
%E_MK =  2.337359/4; % eV
E_MK = 1.907001;

%% Energy constant:
% MxOy
xM = 1;
yO = 2;

E_Metalbulk = -3.8162; % eV. Bulk metal
E_MetalOxide = -19.6874; % eV. Energy per unit cell
E_water = -14.2176; % eV. One water molecule
E_O2_over = 1.22369; % eV. Over binding of O2 due to DFT
E_O2 = -9.8579+E_O2_over; % eV. Energy of O2 molecule
E_MetalFluoride = -23.854444; % eV


cellsize = 32 ; % formula units
dosCellsize = cellsize/2 ; % Cell Size in DOS
Eperf = E_MetalOxide*cellsize;

%% Standard chemical potential
%mu_O_p0_table = [[100,-0.08];[200,-0.17];[300,-0.27];[400,-0.38];...
%    [500,-0.50];[600,-0.61];[700,-0.73];[800,-0.85];[900,-0.98];...
%    [1000,-1.10]];
mu_H2O_p0_table = [[298.15, -0.48086], [300, -0.48448], [400, -0.68568],...
    [500, -0.89585], [600, -1.1134], [700, -1.3372], [800, -1.5666],...
    [900, -1.801], [1000, -2.04], [1100, -2.2833], [1200, -2.5306],...
    [1300, -2.7817], [1400, -3.0364], [1500, -3.2945], [1600, -3.5558],...
    [1700, -3.8203], [1800, -4.0877], [1900, -4.358], [2000, -4.6311]];
mu_O_p0_table = [[298.150000, -0.271962]; [300.000000, -0.273931];...
    [350.000000, -0.327734]; [400.000000, -0.382644];...
    [450.000000, -0.438528]; [500.000000, -0.495295];...
    [600.000000, -0.611189]; [700.000000, -0.729865];...
    [800.000000, -0.850985]; [900.000000, -0.974297];...
    [1000.000000, -1.09959]; [1100.000000, -1.22669];...
    [1200.000000, -1.35546]; [1300.000000, -1.48577];...
    [1400.000000, -1.61752]; [1500.000000, -1.7506];...
    [1600.000000, -1.88495]; [1700.000000, -2.02049];...
    [1800.000000, -2.15717]; [1900.000000, -2.29491];...
    [2000.000000, -2.43368]];
[M,I]=min(abs(mu_O_p0_table(:,1) - temperature));
mu_O_p0 = mu_O_p0_table(I,2);
[M,I]=min(abs(mu_H2O_p0_table(:,1) - temperature));
mu_H2O_p0 = mu_H2O_p0_table(I,2);
mu_H = 1/2*(E_water + mu_H2O_p0 - 1/2*E_O2 - mu_O_p0);



mu_O_min = (E_MetalOxide - xM*E_Metalbulk) / yO;
mu_O_max = (E_O2) / 2;


mu_e_range = [0 3];



%% test range %%
% mu_O = -5.5401;
% pO2 = exp((mu_O-E_O2/2-mu_O_p0)/k/temperature);
% mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
% mu_H = (E_water - mu_O)/2;
% 
% mu = [mu_O mu_Metal];
%  
% systemTotalCharge(1)
%systemTotalCharge(2)
%systemTotalCharge(3)

%%

%%% Latex legend for GLE plot %%%
latexLineLegend = ["$P_{O_2}$","$\mu_e$","e","h"];
for i = 1:length(defect)
    temp = defect(i,1:5) ;
    if temp(1) == 0 
        suffix = "^{"+latexSymbolCharge(2)+"}";
    else
        suffix = "^{"+...
            strjoin(repmat(latexSymbolCharge(temp(1)/abs(temp(1))+2),1,abs(temp(1))),'')+"}";
    end
    if sum(temp(2:5)) == -1
        latexName = "$V_{"+element(find(temp(2:5)))+"}"+suffix+ "$"; %#ok<FNDSB>
    elseif sum(temp(2:5)) == 1
        latexName =  "$" + element(find(temp(2:5)))+"_{i}" + suffix + "$"; %#ok<FNDSB>
    else
        latexName = "$" + element(find(temp(2:5) == 1)) + "_{" ...
            + element(find(temp(2:5) == -1 )) + "}" + suffix+ "$"; %#ok<FNDSB>
    end
    latexLineLegend = [latexLineLegend, latexName];             %#ok<AGROW>
end
fID = fopen('header.txt','w');
fprintf(fID,'%s\n',strjoin(latexLineLegend));
fclose(fID);
%latexLineLegendUndoped = latexLineLegend(1:20);

%% DEFECTS CALCULATION %%

%%%%%% KROGER-VINK No hydrogen

if module == 1
    fprintf("KROGER-VINK UNDOPED\n")
    result = []; %#ok<UNRCH>
    for i = E_O2/2+mu_O_p0:-0.05:-9.4
        mu_O = i;
        pO2 = exp(2*(mu_O-E_O2/2-mu_O_p0)/k/temperature);
        mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
        %mu_H = (E_water - mu_O)/2;
        %mu_e_result = double(vpasolve(systemTotalCharge(mu_e) == 0, mu_e, mu_e_range));
        mu_e_result = fzero(@(x)systemTotalCharge(x), 2);
        
        if(isempty(mu_e_result))
            temp = [pO2, 0];
            for j = 1:length(defect)
                temp = [temp,0, 0];
            end
        else
            
            temp = [pO2,mu_e_result, nc(mu_e_result), pv(mu_e_result)];
            mu = [mu_O mu_Metal];
            for j = 1:length(defect)
                charge = defect(j,1);
                nMu = defect(j,2:3);
                noSite = defect(j,end);
                E = defect(j,end-1);
                E_form = formationEnergy(E,charge,mu_e_result, nMu,mu);
                c = concentration(E_form, noSite);
                %temp = [temp,E_form, c];
                temp = [temp, c];
            end
        end
        result = [result;temp];
    end
    output = double(result);
    %output = [latexLineLegend;output];
    %save(outputFile,'output','-ascii')
    dlmwrite(char(outputFile),char(strjoin(latexLineLegend)),'delimiter','');
    dlmwrite(char(outputFile),output,'-append','delimiter',' ');
    %save temp.txt result -ASCII
    fprintf("Output at "+outputFile+"\n")
end

%%
%%%%%% KROGER-VINK H doped

if module == 1.5
    fprintf("KROGER-VINK H-DOPED\n")
    result = []; %#ok<UNRCH>
    for i = E_O2/2+mu_O_p0:-0.05:-9.4
        mu_O = i;
        pO2 = exp(2*(mu_O-E_O2/2-mu_O_p0)/k/temperature)
        mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
        %mu_H = (E_water - mu_O)/2;
        %mu_H = 1000;
        %mu_e_result = double(vpasolve(systemTotalCharge(mu_e) == 0, mu_e, mu_e_range));
        mu_e_result = fzero(@(x)systemTotalCharge(x), 2);
        
        if(isempty(mu_e_result))
            temp = [pO2, 0];
            for j = 1:length(defect)
                temp = [temp,0, 0];
            end
        else
            
            temp = [pO2,mu_e_result, nc(mu_e_result), pv(mu_e_result)];
            mu = [mu_O mu_Metal mu_H 0];
            for j = 1:length(defect)
                charge = defect(j,1);
                nMu = defect(j,2:5);
                noSite = defect(j,end);
                E = defect(j,end-1);
                E_form = formationEnergy(E,charge,mu_e_result, nMu,mu);
                c = concentration(E_form, noSite);
                %temp = [temp,E_form, c];
                temp = [temp, c];
            end
        end
        result = [result;temp];
    end
    output = double(result);
    %output = [latexLineLegend;output];
    %save(outputFile,'output','-ascii')
    dlmwrite(char(outputFile),char(strjoin(latexLineLegend)),'delimiter','');
    dlmwrite(char(outputFile),output,'-append','delimiter',' ');
    %save temp.txt result -ASCII
    fprintf("Output at "+outputFile+"\n")
end
%%
%%%% KROGER-VINK w Doping
syms mu_e mu_doping
format shortG

if module == 2
    fprintf("KROGER-VINK W DOPING\n")
    result = [];
    output = [];
    for i = E_O2/2+mu_O_p0:-0.05:-8.4
        mu_O = i;
        pO2 = exp(2*(mu_O-E_O2/2-mu_O_p0)/k/temperature)
        mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
        %mu_H = (E_water - mu_O)/2;
        mu_e_result = vpasolve([systemTotalChargeWithDoping(mu_e, mu_doping) == 0, ...
            totalDopingConcentration(mu_e, mu_doping) == dopingConcentration],[mu_e, mu_doping], mu_e_range);
        mu_doping_result = mu_e_result.mu_doping;
        mu_e_result = mu_e_result.mu_e;
        
        %mu_e_result = fzero(@(x)systemTotalChargeWithDoping(x), 2);
        if(isempty(mu_e_result))
            temp = [pO2,mu_O,0];
            temp2 = [pO2,0];
            for j = 1:length(defect)
                temp = [temp,0, 0]; %#ok<AGROW>
                temp2 = [temp2, 0]; %#ok<AGROW>
            end
        else
            temp = [pO2,mu_O,mu_e_result, nc(mu_e_result), pv(mu_e_result)];
            temp2 = [pO2,mu_e_result, nc(mu_e_result), pv(mu_e_result)];
            mu = [mu_O mu_Metal mu_H mu_doping_result];
            for j = 1:length(defect)
                charge = defect(j,1);
                nMu = defect(j,2:5);
                noSite = defect(j,end);
                E = defect(j,end-1);
                E_form = formationEnergy(E,charge,mu_e_result, nMu,mu);
                c = concentration(E_form, noSite);
                temp2 = [temp2, c]; %#ok<AGROW>
                temp = [temp,E_form, c]; %#ok<AGROW>
                %if (j == 19) | (j == 27) 
                %    temp = [temp,c];
                %end
            end
        end
        result = [result;temp]; %#ok<AGROW>
        output = [output;temp2]; %#ok<AGROW>
        
    end
    result = double(result);
    output = double(output);
    %output = [latexLineLegend;output];
    %save(outputFile,'output','-ascii')
    dlmwrite(char(outputFile),char(strjoin(latexLineLegend)),'delimiter','');
    dlmwrite(char(outputFile),output,'-append','delimiter',' ');
    save temp.txt result -ASCII;
    fprintf("Output at "+outputFile+"\n")
end

% for i = mu_O_min:0.1:mu_O_max
%     mu_O = i
%     mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
%     mu_H = (E_water - mu_O)/2;
%     vpasolve([systemTotalCharge(mu_e, mu_doping) == 0, ...
%     totalDopingConcentration(mu_e, mu_doping)],[mu_e, mu_doping])
% end
% result =[];

%%
%%% Scanned H vs F Doping level
if module == 3
    fprintf("Scanned H vs F Doping level\n")
    Fresult = [];  %#ok<UNRCH>
    for dopingConcentration = logspace(-10,0,11)
        
        result = [];
        
        mu_O = E_O2/2+mu_O_p0;
        pO2 = exp(2*(mu_O-E_O2/2-mu_O_p0)/k/temperature);
        mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
        %mu_H = (E_water - mu_O)/2;
        mu_e_result = vpasolve([systemTotalChargeWithDoping(mu_e, mu_doping) == 0, ...
            totalDopingConcentration(mu_e, mu_doping) == dopingConcentration],[mu_e, mu_doping], mu_e_range);
        mu_doping_result = mu_e_result.mu_doping;
        mu_e_result = mu_e_result.mu_e;
        
        %mu_e_result = fzero(@(x)systemTotalChargeWithDoping(x), 2);
        if(isempty(mu_e_result))
            temp = [pO2,mu_O,0];
            for j = 1:length(defect)
                temp = [temp,0, 0]; 
            end
        else
            temp = [pO2,mu_O,mu_e_result, nc(mu_e_result), pv(mu_e_result)];
            mu = [mu_O mu_Metal mu_H mu_doping_result];
            for j = 1:length(defect)
                charge = defect(j,1);
                nMu = defect(j,2:5);
                noSite = defect(j,end);
                E = defect(j,end-1);
                E_form = formationEnergy(E,charge,mu_e_result, nMu,mu);
                c = concentration(E_form, noSite);
                temp = [temp,E_form, c];
            end
            
            result = [result;temp];
            
            result = double(result);
            %save temp.txt result -ASCII
            [dopingConcentration/2,sum(result(end,43:2:61))]
    Fresult = [Fresult;[dopingConcentration/2,sum(result(end,43:2:61))]];
    
    Fresult
        end
    end
    [dopingConcentration/2,sum(result(end,43:2:61))]
    Fresult = [Fresult;[dopingConcentration/2,sum(result(end,43:2:61))]];
    
    Fresult
end

% for i = mu_O_min:0.1:mu_O_max
%     mu_O = i
%     mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
%     mu_H = (E_water - mu_O)/2;
%     vpasolve([systemTotalCharge(mu_e, mu_doping) == 0, ...
%     totalDopingConcentration(mu_e, mu_doping)],[mu_e, mu_doping])
% end
% result =[];

%%
%%%%%% Formation vs mu_e at pO2 = 1atm
if module == 4
    fprintf("Formation vs mu_e at pO2 = 1 atm\n")
    mu_O = E_O2/2+mu_O_p0; %#ok<UNRCH>
    mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
    %mu_H = (E_water - mu_O)/2;
    EformVsMuE = [];
    for i = 0:0.1:4.0
        mu_e_temp = i;
        mu = [mu_O mu_Metal mu_H 0];
        temp = mu_e_temp;
        for j = 1:length(defect)-totalDopingDefects
            charge = defect(j,1);
            nMu = defect(j,2:5);
            noSite = defect(j,end);
            E = defect(j,end-1);
            E_form = formationEnergy(E,charge,mu_e_temp, nMu,mu);
            temp = [temp, E_form];
        end
        EformVsMuE = [EformVsMuE; temp];
    end
    save Evsmu_E.txt EformVsMuE -ASCII
    fprintf("Output at Evsmu_E.txt\n")
end
% csvwrite('Eform-muE.txt',EformVsMuE);
%totalDefectsCharge(0.3,-5)

% vpasolve(totalDefectsCharge(1,mu_doping) == 10,mu_doping)

%%
%%%%%%% Scan [C] of H at 1 atm across different element
if module == 5
    fprintf("Scan [C] of H at 1 atm across different element\n");
    result = [];
    temperatureRange = [600,900,1200];
    for i=1:length(allDefectElement)
        defectElement = allDefectElement(i);
        fprintf(defectElement)
        element = [systemElement, defectElement];
        defectFile = 'defects-'+defectElement+'.txt';
        outputFile = 'K-V-'+defectElement+'.txt' ;
        defect = importdata(defectFile);
        defect = defect.data;
        dopingConcentration = 10^-4 ;
        totalDopingDefects = find(defect(:,5)>0, 1, 'last' )-find(defect(:,5)>0, 1 )+1;
        temp = [];
        for z = 1:length(temperatureRange)
            temperature = temperatureRange(z);
            [~,I]=min(abs(mu_O_p0_table(:,1) - temperature));
            mu_O_p0 = mu_O_p0_table(I,2);
            [~,I]=min(abs(mu_H2O_p0_table(:,1) - temperature));
            mu_H2O_p0 = mu_H2O_p0_table(I,2);
            mu_H = 1/2*(E_water + mu_H2O_p0 - 1/2*E_O2 - mu_O_p0);
            mu_O = E_O2/2+mu_O_p0;
            pO2 = 1;
            mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
            %mu_H = (E_water - mu_O)/2;
            mu_e_result = vpasolve([systemTotalChargeWithDoping(mu_e, mu_doping) == 0, ...
                totalDopingConcentration(mu_e, mu_doping) == dopingConcentration],[mu_e, mu_doping], mu_e_range);
            mu_doping_result = mu_e_result.mu_doping;
            mu_e_result = mu_e_result.mu_e;
            cSum = 0;
            %mu_e_result = fzero(@(x)systemTotalChargeWithDoping(x), 2);
            if(isempty(mu_e_result))
                mu_e_result = 0;
                cSum = 0;
            else
                mu = [mu_O mu_Metal mu_H mu_doping_result];
                for j = 17:28
                    charge = defect(j,1);
                    nMu = defect(j,2:5);
                    noSite = defect(j,end);
                    E = defect(j,end-1);
                    E_form = formationEnergy(E,charge,mu_e_result, nMu,mu);
                    c = concentration(E_form, noSite);
                    cSum = cSum + c;
                    %if (j == 19) | (j == 27)
                    %    temp = [temp,c];
                    %end
                end
                
            end
            temp2 = repmat("?",1,length(temperatureRange));
            temp2(z) = double(cSum);
            temp = [temp;strjoin([double(mu_e_result), defectElement, temp2])];
            dlmwrite('Hvsmu_E-GLEnames.txt',char(strjoin(["amove xg(",double(mu_e_result),") yg(",double(cSum*10),")"])),'-append','delimiter','');
            dlmwrite('Hvsmu_E-GLEnames.txt',char(strjoin(["text ",defectElement])),'-append','delimiter','');
        end
        result = [result;temp];
    end
    dlmwrite('Hvsmu_E.txt',char("$\mu_e$ elements "+strjoin(string(temperatureRange)+"K")),'delimiter','');
    dlmwrite('Hvsmu_E.txt',char(result),'-append','delimiter','');
    fprintf("Output at Hvsmu_E.txt and Hvsmu_E-GLEnames.txt\n")
end
%%
%%%%%%% Scan [C] of H at 1 atm across different mu_E
if module == 6
    fprintf("Scan [H] at 1 atm across different mu_E\n");
    result = [];
    temperatureRange = [600,900,1200];
    defectFile = 'w-o-doping.txt';
    defect = importdata(defectFile);
    defect = defect.data;
    for i=0:0.01:2.6
        temp = [];
        for z = 1:length(temperatureRange)
            temperature = temperatureRange(z);
            [~,I]=min(abs(mu_O_p0_table(:,1) - temperature));
            mu_O_p0 = mu_O_p0_table(I,2);
            [~,I]=min(abs(mu_H2O_p0_table(:,1) - temperature));
            mu_H2O_p0 = mu_H2O_p0_table(I,2);
            mu_H = 1/2*(E_water + mu_H2O_p0 - 1/2*E_O2 - mu_O_p0);
            mu_O = E_O2/2+mu_O_p0;
            pO2 = 1;
            mu_Metal = 1/xM*(E_MetalOxide - yO*mu_O);
            %mu_H = (E_water - mu_O)/2;
            %mu_e_result = vpasolve([systemTotalChargeWithDoping(mu_e, mu_doping) == 0, ...
            %    totalDopingConcentration(mu_e, mu_doping) == dopingConcentration],[mu_e, mu_doping], mu_e_range);
            %mu_doping_result = mu_e_result.mu_doping;
            %mu_e_result = mu_e_result.mu_e;
            mu_e_result = i;
            cSum = 0;
            %mu_e_result = fzero(@(x)systemTotalChargeWithDoping(x), 2);
            if(isempty(mu_e_result))
                mu_e_result = 0;
                cSum = 0;
            else
                mu = [mu_O mu_Metal mu_H];
                for j = 17:28
                    charge = defect(j,1);
                    nMu = defect(j,2:4);
                    noSite = defect(j,end);
                    E = defect(j,end-1);
                    E_form = formationEnergy(E,charge,mu_e_result, nMu,mu);
                    c = concentration(E_form, noSite);
                    cSum = cSum + c;
                    %if (j == 19) | (j == 27)
                    %    temp = [temp,c];
                    %end
                end
                
            end
            temp2 = repmat("?",1,length(temperatureRange));
            temp2(z) = double(cSum);
            temp = [temp;strjoin([double(mu_e_result), "?", temp2])];
        end
        result = [result;temp];
    end
    dlmwrite('Hvsmu_E-full.txt',char("$\mu_e$ elements "+strjoin(string(temperatureRange)+"K")),'delimiter','');
    dlmwrite('Hvsmu_E-full.txt',char(result),'-append','delimiter','');
    fprintf("Output at Hvsmu_E-full.txt\n")
end



%%%%% Supporting functions %%%%%%%

function y = totalDopingConcentration(mu_e, mu_doping)
global mu_O mu_Metal mu_H defect totalDopingDefects;
mu = [mu_O mu_Metal mu_H mu_doping];
y = 0;
for i = length(defect)-totalDopingDefects+1:length(defect)
    charge = defect(i,1);
    nMu = defect(i,2:5);
    noSite = defect(i,end);
    E = defect(i,end-1);
    E_form = formationEnergy(E,charge,mu_e, nMu,mu);
    c = concentration(E_form, noSite);
    y = y + c;
end
end

function y = systemTotalChargeWithDoping(mu_e, mu_doping)
y = pv(mu_e)-nc(mu_e)+totalDefectsChargeWithDoping(mu_e, mu_doping);
end

function y = systemTotalCharge(mu_e) %#ok<DEFNU>
y = pv(mu_e)-nc(mu_e)+totalDefectsCharge(mu_e);
end

function y = totalDefectsCharge(mu_e)
global mu_O mu_Metal mu_H defect;
mu = [mu_O mu_Metal mu_H 0];
y = 0;
for i = 1:length(defect)
    charge = defect(i,1);
    nMu = defect(i,2:5);
    noSite = defect(i,end);
    E = defect(i,end-1);
    E_form = formationEnergy(E,charge,mu_e, nMu,mu);
    c = charge*concentration(E_form, noSite);
    y = y + c;
end
end


function y = totalDefectsChargeWithDoping(mu_e, mu_doping)
global mu_O mu_Metal mu_H defect;
mu = [mu_O mu_Metal mu_H mu_doping];
y = 0;
for i = 1:length(defect)
    charge = defect(i,1);
    nMu = defect(i,2:5);
    noSite = defect(i,end);
    E = defect(i,end-1);
    E_form = formationEnergy(E,charge,mu_e, nMu,mu);
    c = charge*concentration(E_form, noSite);
    y = y + c;
end
end

function y = concentration(Eform, noSite)
global k temperature
y = noSite*exp(-Eform/(k*temperature));
end

function y = formationEnergy(Etot, charge, mu_e, nMu, mu)
global dielectric_const E_MK Eperf vbm;
y = Etot - Eperf ...
    + charge^2*E_MK/dielectric_const + charge*(vbm+mu_e)  ...
    - dot(nMu,mu);
end

function y = fermi_dirac(E, mu_e)
global temperature k
y = 1./(1+exp((E-mu_e)./(k*temperature)));
end

function y = nc(mu_e)
global total_dos energy cbm vbm dosCellsize
step = energy(2)-energy(1);
y = sum(total_dos(energy >= cbm).* step .* ...
    fermi_dirac(energy(energy >= cbm),mu_e+vbm))/dosCellsize;
end

function y = pv(mu_e)
global total_dos energy vbm dosCellsize
step = energy(2)-energy(1);
y = sum(total_dos(energy <= vbm).* step .* ...
    fermi_dirac(-1*energy(energy <= vbm),-mu_e-vbm))/dosCellsize;
end