clc;
close all;
clear;
%------------This is the main code for diffuser design of ramjet (considering self-start efficient)----------%
%Defining user input properties:

%First the input free stream properties (at Station 1):
P1= input("Enter the Input Free Stream Pressure: ");
T1=input("Enter the Input Free Stream Temperature: ");
M1=input("Enter a list of Mach values: ");

%Defining the input Shock strength for the ramjet diffuser design:
p2b_over_p1=input("Enter the Input Normal Shock Strength Pressure Ratio: ");

%Defining Burner Entry Mach value:
M3=input("Enter the Input Burner Entry Mach value: ");

%Defining input Required Thrust for the ramjet design:
F=input("Enter the Input Performance Thrust value for Ramjet: ");

%Defining the Gas Dynamic properties:
gamma=1.4;
R=287;

%|=========================== STATION 1 : ENTRY OF DIFFUSER ============================================================|

%Calculating the input free stream stagnation properties and other gas props:
P01=P1*(1+((gamma-1)/2).*M1.^2).^(gamma/(gamma-1)); %Isentropic relation

T01=T1*(1+((gamma-1)/2).*M1.^2); %Isentropic relation

V1=M1.*(gamma*R*T1)^(0.5);  

rho1= P1./(R*T1); % Ideal gas relation

%Calculating the Normal shock Mach upstream value Ms from the pressure ratio given as shock strength:
Ms_normal_upstream = inverse_normal_shock(p2b_over_p1);

assert (Ms_normal_upstream>1); %Checking if Upstream normal shock Mach is supersonic or not: should be supersonic, if otherwise ->error

%|================================== END OF STATION 1 =================================================================|

%If shock is after the throat

if (Ms_normal_upstream>=1.3)
    %|===================================== STATION 2 STARTS:THROAT REGION CALCS ===========================================|

    %|==============================MASS FLOW RATE AND AREA CALCULATIONS ================================================|

    m_dot=F./V1; %mass flow rate

    A1=m_dot./(rho1*V1); % Inlet area of crossection

    A_throat=A1./area_mach_func(M1); %Throat area of crossection

    %|======================================= STATE PROPERTIES AT THROAT ===========================================|

    % Calculate the stagnation pressure and temperature at the throat
    P02 = P01;
    T02 = T01;

    %Static properties at the throat:
    P2=P02/isentropic_press_ratio(M2);
    T2=T02/isentropic_temp_ratio(M2);

    %|============================= STATION 2B: NORMAL SHOCK ========================================================|

    %Post Normal Shock Analysis: Static and Stag Properties calculation and verification assertions for errors:

    M2b=((Ms_normal_upstream^2*(gamma-1)+2)/(2*gamma*Ms_normal_upstream^2-(gamma-1)))^(0.5); % Mach value post normal shock (at station 2)

    %Static props post shock

    P2b=P2.*p2b_over_p1;  %Static pressure after shock

    T2b_over_T2=norm_shock_temp_static(Ms_normal_upstream); %Static temp ratio calc post shock

    T2b=T2.*T2b_over_T2; %Static temperature calculations post shock

    %Stagantion props post shock

    T02b=T02; %Stagnation Temp constant across shock -> Adiabatic shock assumption

    P02b=P02.*norm_shock_press_stag(Ms_normal_upstream); %Stagnation Pressure calculations across shock -> Adiabatic shock assumption

    %Gas properties post normal shock:

    rho2b=P2b./(R*T2b);

    V2b=M2b.*(gamma*R*T2b)^(0.5);

    %|=============================== STATION 2b:NORMAL SHOCK ENDS ============================================================|
    %|=============================== STATION 3: ISENTROPIC SUBSONIC DIFFUSION STARTS ===================================|
    %Calculating Static and Stag values here for throat entry values (at station 3), also gas properties as well.

    assert ((M3>0) &&  all(M2b>0) && (M3<1) && all(M2b<1)); %Checking for the Mach values at burner and post shock to be subsonic and valid positive values.

    P03=P02b; % Isentropic process

    T03= T02b; %Isentropic process

    P3 = P03/isentropic_press_ratio(M3);

    T3 = T03/isentropic_temp_ratio(M3);

    %Gas properties:
    rho3=P3/(R*T3); % density at burner entry

    V3=M3.*(gamma*R*T3)^(0.5); %Velocity inlet at burner entry

    %|================================ STATION 3 END : BURNER STARTS ====================================================|
    A3=m_dot./(rho3*V3); % Burner Entry Area

end    

% If shock is at throat

if (Ms_normal_upstream<1.3)
   %Taking the post shock conditions here itself
   %|============================= STATION 2B: NORMAL SHOCK ========================================================|

    %Post Normal Shock Analysis: Static and Stag Properties calculation and verification assertions for errors:

    M2b=((Ms_normal_upstream^2*(gamma-1)+2)/(2*gamma*Ms_normal_upstream^2-(gamma-1)))^(0.5); % Mach value post normal shock (at station 2)

    %Static props post shock

    P2b=P1.*p2b_over_p1;  %Static pressure after shock

    T2b_over_T1=norm_shock_temp_static(Ms_normal_upstream); %Static temp ratio calc post shock

    T2b=T1.*T2b_over_T1; %Static temperature calculations post shock

    %Stagantion props post shock

    T02b=T01; %Stagnation Temp constant across shock -> Adiabatic shock assumption

    P02b=P01.*norm_shock_press_stag(Ms_normal_upstream); %Stagnation Pressure calculations across shock -> Adiabatic shock assumption

    %Gas properties post normal shock:

    rho2b=P2b./(R*T2b);

    V2b=M2b.*(gamma*R*T2b)^(0.5);

    %|=============================== STATION 2b:NORMAL SHOCK ENDS ============================================================|
    %|=============================== STATION 3: ISENTROPIC SUBSONIC DIFFUSION STARTS ===================================|
    %Calculating Static and Stag values here for throat entry values (at station 3), also gas properties as well.

    assert ((M3>0) &&  all(M2b>0) && (M3<1) && all(M2b<1)); %Checking for the Mach values at burner and post shock to be subsonic and valid positive values.

    P03=P02b; % Isentropic process

    T03= T02b; %Isentropic process

    P3 = P03/isentropic_press_ratio(M3);

    T3 = T03/isentropic_temp_ratio(M3);

    %Gas properties:
    rho3=P3/(R*T3); % density at burner entry

    V3=M3.*(gamma*R*T3)^(0.5); %Velocity inlet at burner entry

    %|================================ STATION 3 END : BURNER STARTS ====================================================|
    A3=m_dot./(rho3*V3); % Burner Entry Area

end

% Stagnation Pressure loss metric --> Important for diffuser efficiency optimization:

n_diff=P03./P01;

%|================================= PLOTTING SECTION ================================================================|

figure ('Name','Stagnation Pressure Loss V/S Mach Value');
plot(M1,n_diff,"-",LineWidth=2,Color="b");
xlabel ('Mach No.');
ylabel ('Stagnation Pressure Loss Ratio');

%|==================== OUTPUT RESULTS FOR DIFFUSER DESIGN ===========================================================|
fprintf('\n');
fprintf('The stagnation pressure loss ratio is: %.4f\n',n_diff);
fprintf('\n');
fprintf('The diffuser throat area is: %.2f\n',A_throat);
fprintf('\n');
fprintf('The inlet area for the diffuser is: %.2f\n',A1);
fprintf('\n');
fprintf('The burner entry area for the diffuser is: %.2f\n',A3);
fprintf('\n');
fprintf('The normal shock area of cross section for the diffuser is: %.2f\n',A2);
fprintf('\n');
fprintf('The Burner entry pressure is: %.2f\n',P3);
fprintf('\n');
fprintf('The Burner entry temperature is: %.2f\n',T3);
fprintf('\n');
fprintf('The Burner entry Mach is: %.2f\n',M3);
fprintf('\n');


   



