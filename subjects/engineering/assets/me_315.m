%% Parameters
n1 = 3000;      % Shaft 1 speed (RPM)
H1 = 40000;     % Power at Input (W)
H2 = 16000;     % Combined Power of Shaft 3 Outputs (W)
H3 = 24000;     % Combined Power of Shaft 2 Outputs (W)
phi = 20;       % Pressure Angle (degrees)
psi = 0;        % Helical Angle (degrees)

phi = phi*pi/180;   % Converting to radians
psi = psi*pi/180;   % Converting to radians

m = .004;       % Module (m)
N1A = 17;       % # of Teeth for Gear 1A
N1B = 17;       % # of Teeth for Gear 1B
N2 = 73;        % # of Teeth for Gear 2
N3 = 64;        % # ofTeeth for Gear 3

L = .2;          % Length of shaft (m)
A = (2/3)*L;     % Distance of Gears 1B and 2 down the shaft
B = (1/3)*L;     % Distance of Gears 1A and 3 down the shaft

% Shaft Diameters
D1 = .02;
D2 = .02;
D3 = .03;

Life = 10400;   % Life of system (hours)

%Gear Hardness
HB_1A = 300;
HB_1B = 300;
HB_2 = 250;
HB_3 = 250;

%Shaft Material Parameters
%AISI 1020 Steel
Sy = 295e6;
Sut = 395e6;
Kf = 3;
Kfs = 2.8;

%% Secondary Parameters
%Parameters which depend on other parameters and should be adjusted as
%needed

% Geometry factors for bending
    % Varies with N1A, N1B, N2, N3
Yj_1A = .3;
Yj_1B = .3;
Yj_2 = .43;
Yj_3 = .42;

% Load distribution factors
    % Varies with bw, d1A, d1B
Km_1A = 1.19;
Km_1B = 1.19;
Km_2 = Km_1B;
Km_3 = Km_1A;

%% Output Speeds

i_AB = N1B/N2;  %Gear ratios
i_CB = N1A/N3;

n2 = i_AB*n1;
n3 = i_CB*n1;
fprintf('\nn2 = %6.2f RPM\nn3 = %6.2f RPM\n\n',n2,n3);

%% Gear Forces

d1A = m*N1A;        %Pitch diameters
d1B = m*N1B;
d2 = m*N2;
d3 = m*N3;

w1 = n1*2*pi/60;    %Angular shaft speeds
w2 = n2*2*pi/60;
w3 = n3*2*pi/60;

T1 = H1/w1;         %Torque at input
T2 = H2/w2;         %Torque at Gear 2
T3 = H3/w3;         %Torque at Gear 3

%Forces on Gear 2
Wt2 = 2*T2/d2;
W2 = Wt2/(cos(phi)*cos(psi));
Wr2 = W2*sin(phi);
Wa2 = W2*cos(phi)*sin(psi);
fprintf('Wt2 = %6.4f kN\nWr2 = %6.4f kN\nWa2 = %6.4f kN\n\n',Wt2/1000,Wr2/1000,Wa2/1000);

%Forces on Gear 3
Wt3 = 2*T3/d3;
W3 = Wt3/(cos(phi)*cos(psi));
Wr3 = W3*sin(phi);
Wa3 = W3*cos(phi)*sin(psi);
fprintf('Wt3 = %6.4f kN\nWr3 = %6.4f kN\nWa3 = %6.4f kN\n\n',Wt3/1000,Wr3/1000,Wa3/1000);

%Forces on Gear 1A
fprintf('Wt1A = %6.4f kN\nWr1A = %6.4f kN\nWa1A = %6.4f kN\n\n',Wt3/1000,Wr3/1000,Wa3/1000);

%Forces on Gear 1B
fprintf('Wt1B = %6.4f kN\nWr1B = %6.4f kN\nWa1B = %6.4f kN\n\n',Wt2/1000,Wr2/1000,Wa2/1000);

T1A = T1;               %Torque at Gear 1A
T1B = T1 - Wt3*d1A/2;   %Torque at Gear 1B

%% Bearing Reactions

%SHAFT 1
%Shaft 1, y-dir
N11y = Wr2*(L-A)/L - Wr3*(L-B)/L;   %Bearing closer to motor
N12y = Wr2 - Wr3 - N11y;    %Bearing farther from motor

%Shaft 1, z-dir
N11z = -Wt2*(L-A)/L + Wt3*(L-B)/L;  %Bearing closer to motor
N12z = -Wt2 + Wt3 - N11z;   %Bearing farther from motor

%SHAFT 2
%Shaft 2, y-dir
N22y = -Wr2*A/L;            %Bearing farther from motor
N21y = -Wr2 - N22y;         %Bearing closer to motor

%Shaft 2, z-dir
N22z = Wt2*A/L;             %Bearing farther from motor
N21z = Wt2 - N22z;          %Bearing closer to motor

%SHAFT 3
%Shaft 3, y-dir
N32y = Wr3*B/L;             %Bearing farther from motor
N31y = Wr3 - N32y;          %Bearing closer to motor

%Shaft 3, z-dir
N32z = -Wt3*B/L;            %Bearing farther from motor
N31z = -Wt3 - N32z;         %Bearing closer to motor

%% Gear Geometry

bw = 12*m;   %Gear width
pc = pi*m;   %Circular pitch
t = pc/2;    %Tooth thickness
a = m;       %Addendum
b = 1.25*m;  %Dedendum
l = 2.25*m;  %Whole depth

%% Materials
% Find these by hand using graphs

%% Allowable Gear Stress

KT = 1;   % No special temperature requirements
KR = 1;   % 99 percent reliability

%Life of individual gears (in revolutions)
Life_1 = Life*60*n1;    %Gear 1A/1B life (revs)
Life_2 = Life*60*n2;    %Gear 2 life (revs)
Life_3 = Life*60*n3;    %Gear 3 life (revs)

if Life_1 < 3*10^6 || Life_2 < 3*10^6 || Life_3 < 3*10^6
    fprintf('Life cycle out of bounds.\nSome allowable stress values may not be accurate')
end

%Bending Stress Cycle Factors
YN_1 = 1.6831*(Life_1 ^ -.0323);
YN_2 = 1.6831*(Life_2 ^ -.0323);
YN_3 = 1.6831*(Life_3 ^ -.0323);

%Contact Stress Cycle Factors
ZN_1 = 2.466*(Life_1 ^ -.056);
ZN_2 = 2.466*(Life_2 ^ -.056);
ZN_3 = 2.466*(Life_3 ^ -.056);

%Hardness Ratio Factors
CH_2 = 1 + ((8.89*10^-3)*(HB_1B/HB_2) - 8.29*10^-3)*(N2/N1B - 1);
CH_3 = 1 + ((8.89*10^-3)*(HB_1A/HB_3) - 8.29*10^-3)*(N3/N1A - 1);

%Allowable Bending Stress (before modifications)
Sb_1A_all = (.102*HB_1A + 16.4)*6.895e+6;
Sb_1B_all = (.102*HB_1B + 16.4)*6.895e+6;
Sb_2_all = (.0773*HB_2 + 12.8)*6.895e+6;
Sb_3_all = (.0773*HB_3 + 12.8)*6.895e+6;

%Allowable Contact Stress (before modifications)
Sc_1A_all = (.349*HB_1A + 34.3)*6.895e+6;
Sc_1B_all = (.349*HB_1B + 34.3)*6.895e+6;
Sc_2_all = (.322*HB_2 + 29.1)*6.895e+6;
Sc_3_all = (.322*HB_3 + 29.1)*6.895e+6;

%Allowable Bending Stress (modified)
Sb_1A_all = Sb_1A_all*YN_1/(KT*KR);
Sb_1B_all = Sb_1B_all*YN_1/(KT*KR);
Sb_2_all = Sb_2_all*YN_2/(KT*KR);
Sb_3_all = Sb_3_all*YN_3/(KT*KR);

%Allowable Contact Stress (modified)
Sc_1A_all = Sc_1A_all*ZN_1/(KT*KR);
Sc_1B_all = Sc_1B_all*ZN_1/(KT*KR);
Sc_2_all = Sc_2_all*ZN_2*CH_2/(KT*KR);
Sc_3_all = Sc_3_all*ZN_3*CH_3/(KT*KR);

%% Predicted Gear Stress
%Factors omitted from analysis if == 1 by project definition

KE = 473713;  %Elastic factor (Pa^.5)

%Geometric factors for contact stress)
I_1A = pi*cos(phi)*sin(phi)/(1+d1A/d3);
I_1B = pi*cos(phi)*sin(phi)/(1+d1B/d2);
I_2 = I_1B;
I_3 = I_1A;

Ks = 1;  %Size factor

%Dynamic factor
Qv = 6;
Bv = ((12 - Qv) ^ (2/3))/4;
Av = 50 + 56*(1-Bv);

V_1A = pi*n1*d1A;
V_1B = pi*n1*d1B;
V_2 = V_1B;
V3 = V_1A;

Kv_1A = ((Av + sqrt(V_1A))/Av)^Bv;
Kv_1B = ((Av + sqrt(V_1B))/Av)^Bv;
Kv_2 = Kv_1B;
Kv_3 = Kv_1A;

%Predicted Bending Stress
Sb_1A = Wt3*Ks*Km_1A*Kv_1A/(m*bw*Yj_1A);
Sb_1B = Wt2*Ks*Km_1B*Kv_1B/(m*bw*Yj_1B);
Sb_2 = Wt2*Ks*Km_2*Kv_2/(m*bw*Yj_2);
Sb_3 = Wt3*Ks*Km_3*Kv_3/(m*bw*Yj_3);

%Predicted Contact Stress
Sc_1A = KE*sqrt(Wt3*Ks*Km_1A*Kv_1A/(bw*d1A*I_1A));
Sc_1B = KE*sqrt(Wt2*Ks*Km_1B*Kv_1B/(bw*d1B*I_1B));
Sc_2 = Sc_1B;
Sc_3 = Sc_1A;

%% Factors of Safety

%Bending
nb_1A = Sb_1A_all/Sb_1A;
nb_1B = Sb_1B_all/Sb_1B;
nb_2 = Sb_2_all/Sb_2;
nb_3 = Sb_3_all/Sb_3;

%Contact
nc_1A = Sc_1A_all/Sc_1A;
nc_1B = Sc_1B_all/Sc_1B;
nc_2 = Sc_2_all/Sc_2;
nc_3 = Sc_3_all/Sc_3;

%% Shaft Moments

%Shaft 1 Moment about z-axis
M1Az = N11y*B;
M1Bz = N12y*(L-A);
%Shaft 1 Moment about y-axis
M1Ay = N11z*B;
M1By = N12z*(L-A);
%Shaft 1 Total Moment
M1A = sqrt((M1Az^2)+(M1Ay^2));
M1B = sqrt((M1Bz^2)+(M1By^2));
if M1A > M1B
    M1_max = M1A;   %Max moment is @ gear 1A
else
    M1_max = M1B;   %Max moment is @ gear 1B
end

M2z = N21y*A;   %Shaft 2 Moment about z-axis
M2y = N21z*A;   %Shaft 2 Moment about y-axis
M2_max = sqrt((M2z^2)+(M2y^2)); %Max moment (@ gear 2)

M3z = N31y*B;   %Shaft 3 Moment about z-axis
M3y = N31z*B;   %Shaft 3 Moment about y-axis
M3_max = sqrt((M3z^2)+(M3y^2)); %Max moment (@ gear 3)

%% Shaft Bending Stresses

%Sbs_1 = ...
%Sbs_2 = ...
Sbs_3 = 32*M3_max/(pi*D3^3);    %Shaft 3 (@ gear 3)


%% Shaft Torque Shears

%Tts_1 = ...
%Tts_1_a = ...
%Tts_2_m = ...

%Tts_2 = ...
%Tts_2_a = ...
%Tts_2_m = ...

Tts_3 = 16*T3/(pi*D3^3);    %Shaft 3 (@ gear 3, motor side)
Tts_3_a = Tts_3/2;              %Alternating
Tts_3_m = Tts_3/2;              %Mean


%% Shaft Von Mises

%Svm_1_a = ...
%Svm_1_m = ...

%Svm_2_a = ...
%Svm_2_m = ...

Svm_3_a = sqrt(((Kf*Sbs_3)^2)+3*(Kfs*Tts_3_a)^2);
Svm_3_m = sqrt(3)*Tts_3_m;

%% Goodman and Yield Lines

Se = .5*Sut;

kf = .95;   %Ground surface
kr = .82;   %99% survival

ks_1 = 1.189*D1^-.112;
ks_2 = 1.189*D2^-.112;
ks_3 = 1.189*D3^-.112;

%ns_1 = ...
%ny_1 = ...

%ns_2 = ...
%ns_2 = ...

ns_3 = (Svm_3_a/(kf*ks_3*kr*Se) + Svm_3_m/(Sut))^-1;        %Goodman
ny_3 = (sqrt((Sbs_3^2) + 3*(Tts_3_a + Tts_3_m)^2)/Sy)^-1;   %Yield

%% Bearing Selection

%L10_1 = ...
%C_11 = ...
%C_12 = ...

%L10_2 = ...
%C_21 = ...
%C_22 = ...

L10_3 = Life_3/10^6;
C_31 = sqrt((N31y^2) + (N31z^2))*L10_3^(1/3);
C_32 = sqrt((N32y^2) + (N32z^2))*L10_3^(1/3);

fprintf('Bearing 3.1: C = %6.4f kN\nBearing 3.2: C = %6.4f kN\n\n',C_31/1000,C_32/1000);
