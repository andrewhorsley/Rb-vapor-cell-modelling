function output=RbDlines_func(Dline,Bdc)
% RbDlines_func(Dline,Bdc)
%
% Andrew Horsley 25/1/2016
%
% Calculates Rb D line energy levels, optical transition strengths and
% optical transition frequencies
% Can accept Bdc as a vector
% all SI units
%
% To debug, use RbDlines_test.m

% The script breaks if Bdc=0 and the mF level degeneracy is unbroken
% (ElecSus has the same problem)
% we get around this by setting Bdc=0 to Bdc=1uT=0.01G
Bdc(Bdc==0)=1e-6;

% Threshold for transition strengths. Strengths below are generally due 
% to numerical errors, and are set to zero
TransStrengthThreshold=1e-8;

hplanck=6.62606896e-34; % Planck's Constant
muB=9.27400915e-24; % Bohr Magneton

% angular momenta
I87 = 3/2; % 87Rb nuclear
I85= 5/2; % 85Rb nuclear
J_5S1_2 = 1/2; % 5S1/2 electron
J_5P1_2 = 1/2; % 5P1/2 electron
J_5P3_2 = 3/2; % 5P3/2 electron

% Ahfs87_5S1_2 = hplanck*3.41734130545215e9; % 87Rb 5S1/2 hyperfine constant; Steck
% Ahfs87_5P1_2 = hplanck*407.24e6; % 87Rb 5P1/2 hyperfine constant; Steck
% Ahfs87_5P3_2 = hplanck*84.7185e6; % 87Rb 5P3/2 hyperfine constant; Steck
% Bhfs87_5P3_2 = hplanck*12.4965e6; % 87Rb 5P3/2 hyperfine constant; Steck
% 
% Ahfs85_5S1_2 = hplanck*1.0119108130e9; % 85Rb 5S1/2 hyperfine constant; Steck
% Ahfs85_5P1_2 = hplanck*120.527e6; % 85Rb 5P1/2 hyperfine constant; Steck
% Ahfs85_5P3_2 = hplanck*25.0020e6; % 85Rb 5P3/2 hyperfine constant; Steck
% Bhfs85_5P3_2 = hplanck*25.790e6; % 85Rb 5P3/2 hyperfine constant; Steck

Ahfs87_5S1_2 = hplanck*3.41734130545215e9; % 87Rb 5S1/2 hyperfine constant; ElecSus
Ahfs87_5P1_2 = hplanck*406.147e6; % 87Rb 5P1/2 hyperfine constant; ElecSus
Ahfs87_5P3_2 = hplanck*84.7185e6; % 87Rb 5P3/2 hyperfine constant; ElecSus
Bhfs87_5P3_2 = hplanck*12.4965e6; % 87Rb 5P3/2 hyperfine constant; ElecSus

Ahfs85_5S1_2 = hplanck*1.0119108130e9; % 85Rb 5S1/2 hyperfine constant; ElecSus
Ahfs85_5P1_2 = hplanck*120.640e6; % 85Rb 5P1/2 hyperfine constant; ElecSus
Ahfs85_5P3_2 = hplanck*25.038e6; % 85Rb 5P3/2 hyperfine constant; ElecSus
Bhfs85_5P3_2 = hplanck*26.011e6; % 85Rb 5P3/2 hyperfine constant; ElecSus

% g-factors
gJ_5S1_2 = 2.002331; % Steck
gJ_5P1_2 = 0.666;  % Steck
gJ_5P3_2 = 1.3362;  % Steck
gI_85 = - 0.00029364000; % Steck
gI_87 = - 0.00099514; % Steck

% isotope shifts
Rb85_D1shift = 21.624e6; % Elecsus
Rb85_D2shift = 21.734e6; % Elecsus
Rb87_D1shift = -56.077e6; % Elecsus
Rb87_D2shift = -56.361e6; % Elecsus

% optical transition matrices:
% First term comes from J->L reduced density matrix conversion (Steck):
% <J||er||J'>=<L||er||L'>*(-1)^(J'+L+S+1)*sqrt((2J'+1)(2L+1))*Wigner6J({L L' 1};{J' J S})
% Gives a factor of sqrt(1/3) for D1 line, sqrt(2/3) for D2 line. 
% Matrix term is
% <J',mJ'|q|J,mJ>/<J||er||J'>=(-1)^(J'+mJ-1)*sqrt((2J+1)*Wigner3J({J 1 J'};{mJ q -mJ'})
% See Steck, or bk.5 p.26 for details

% 87Rb D1
J87D1p = sqrt(1/3)*[zeros(4,4),-sqrt(2/3)*diag([1,1,1,1]);zeros(4,8)];
J87D1m = sqrt(1/3)*[zeros(4,8);sqrt(2/3)*diag([1,1,1,1]),zeros(4,4)];
J87D1z = sqrt(1/3)*[sqrt(1/3)*eye(4),zeros(4);zeros(4),-sqrt(1/3)*eye(4)];
% 87Rb D2: checked at Bdc=0: OD factor of 2 too large
J87D2p = sqrt(2/3)*[diag([1,1,1,1]),zeros(4,4);...
    zeros(4,4),sqrt(1/3)*diag([1,1,1,1]); zeros(8,8)];
J87D2m = sqrt(2/3)*[zeros(8,8);sqrt(1/3)*diag([1,1,1,1]),zeros(4,4);...
    zeros(4,4),diag([1,1,1,1])];
J87D2z = sqrt(2/3)*[zeros(4,8);...
    sqrt(2/3)*eye(4),zeros(4);zeros(4),sqrt(2/3)*eye(4);...
    zeros(4,8)];
% 85Rb D1: checked at Bdc=0: good +/-0.5%, larger deviations at larger Bdc
J85D1p = sqrt(1/3)*[zeros(6,6),-sqrt(2/3)*diag([1,1,1,1,1,1]);zeros(6,12)];
J85D1m = sqrt(1/3)*[zeros(6,12);sqrt(2/3)*diag([1,1,1,1,1,1]),zeros(6,6)];
J85D1z = sqrt(1/3)*[sqrt(1/3)*eye(6),zeros(6);zeros(6),-sqrt(1/3)*eye(6)];
% 85Rb D2: checked at Bdc=0: OD factor of 2 too large
J85D2p = sqrt(2/3)*[diag([1,1,1,1,1,1]),zeros(6,6);...
    zeros(6,6),sqrt(1/3)*diag([1,1,1,1,1,1]); zeros(12,12)];
J85D2m = sqrt(2/3)*[zeros(12,12);sqrt(1/3)*diag([1,1,1,1,1,1]),zeros(6,6);...
    zeros(6,6),diag([1,1,1,1,1,1])];
J85D2z = sqrt(2/3)*[zeros(6,12);...
    sqrt(2/3)*eye(6),zeros(6);zeros(6),sqrt(2/3)*eye(6);...
    zeros(6,12)];

% prepare arrays for transition strengths and frequencies:
TransFreq87D1_p=zeros(12,length(Bdc))*NaN;
TransFreq87D1_m=zeros(12,length(Bdc))*NaN;
TransFreq87D1_z=zeros(14,length(Bdc))*NaN;
% TransFreq87D2_p=zeros(22,length(Bdc))*NaN;
% TransFreq87D2_m=zeros(22,length(Bdc))*NaN;
% TransFreq87D2_z=zeros(24,length(Bdc))*NaN;

TransStrength87D1_p=zeros(12,length(Bdc))*NaN;
TransStrength87D1_m=zeros(12,length(Bdc))*NaN;
TransStrength87D1_z=zeros(14,length(Bdc))*NaN;
% TransStrength87D2_p=zeros(22,length(Bdc))*NaN;
% TransStrength87D2_m=zeros(22,length(Bdc))*NaN;
% TransStrength87D2_z=zeros(24,length(Bdc))*NaN;

TransFreq87D1_p_array=zeros(8,8,length(Bdc))*NaN;
TransFreq87D1_m_array=zeros(8,8,length(Bdc))*NaN;
TransFreq87D1_z_array=zeros(8,8,length(Bdc))*NaN;
% TransFreq87D2_p_array=zeros(16,8,length(Bdc))*NaN;
% TransFreq87D2_m_array=zeros(16,8,length(Bdc))*NaN;
% TransFreq87D2_z_array=zeros(16,8,length(Bdc))*NaN;
TransStrength87D1_p_array=zeros(8,8,length(Bdc))*NaN;
TransStrength87D1_m_array=zeros(8,8,length(Bdc))*NaN;
TransStrength87D1_z_array=zeros(8,8,length(Bdc))*NaN;
% TransStrength87D2_p_array=zeros(16,8,length(Bdc))*NaN;
% TransStrength87D2_m_array=zeros(16,8,length(Bdc))*NaN;
% TransStrength87D2_z_array=zeros(16,8,length(Bdc))*NaN;

TransFreq85D1_p=zeros(20,length(Bdc))*NaN;
TransFreq85D1_m=zeros(20,length(Bdc))*NaN;
TransFreq85D1_z=zeros(22,length(Bdc))*NaN;
% TransFreq85D2_p=zeros(26,length(Bdc))*NaN;
% TransFreq85D2_m=zeros(27,length(Bdc))*NaN;
% TransFreq85D2_z=zeros(28,length(Bdc))*NaN;

TransStrength85D1_p=zeros(20,length(Bdc))*NaN;
TransStrength85D1_m=zeros(20,length(Bdc))*NaN;
TransStrength85D1_z=zeros(22,length(Bdc))*NaN;
% TransStrength85D2_p=zeros(26,length(Bdc))*NaN;
% TransStrength85D2_m=zeros(27,length(Bdc))*NaN;
% TransStrength85D2_z=zeros(28,length(Bdc))*NaN;

TransFreq85D1_p_array=zeros(12,12,length(Bdc))*NaN;
TransFreq85D1_m_array=zeros(12,12,length(Bdc))*NaN;
TransFreq85D1_z_array=zeros(12,12,length(Bdc))*NaN;
% TransFreq85D2_p_array=zeros(16,12,length(Bdc))*NaN;
% TransFreq85D2_m_array=zeros(16,12,length(Bdc))*NaN;
% TransFreq85D2_z_array=zeros(16,12,length(Bdc))*NaN;
TransStrength85D1_p_array=zeros(12,12,length(Bdc))*NaN;
TransStrength85D1_m_array=zeros(12,12,length(Bdc))*NaN;
TransStrength85D1_z_array=zeros(12,12,length(Bdc))*NaN;
% TransStrength85D2_p_array=zeros(16,12,length(Bdc))*NaN;
% TransStrength85D2_m_array=zeros(16,12,length(Bdc))*NaN;
% TransStrength85D2_z_array=zeros(16,12,length(Bdc))*NaN;

clear TransFreq87D2_p TransFreq87D2_m TransFreq87D2_z
clear TransStrength87D2_p TransStrength87D2_m TransStrength87D2_z
clear TransFreq87D2_p_array TransFreq87D2_m_array TransFreq87D2_z_array
clear TransStrength87D2_p_array TransStrength87D2_m_array TransStrength87D2_z_array
clear TransFreq85D2_p TransFreq85D2_m TransFreq85D2_z
clear TransStrength85D2_p TransStrength85D2_m TransStrength85D2_z
clear TransFreq85D2_p_array TransFreq85D2_m_array TransFreq85D2_z_array
clear TransStrength85D2_p_array TransStrength85D2_m_array TransStrength85D2_z_array


for b=1:length(Bdc)
%% Rb 5S1/2

%%%%%% For 85Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice of basis: psi(1)=|mJ=1/2,mI=5/2>, psi(2)=|mJ=1/2,mI=3/2>,
% psi(3)=|mJ=1/2,mI=1/2>,... psi(8)=|mJ=-1/2,mI=-5/2>
Jz = 1/2*[eye(6),zeros(6);zeros(6),-eye(6)];
Jp = [zeros(6),eye(6);zeros(6),zeros(6)];
Jm = [zeros(6),zeros(6);eye(6),zeros(6)];
Jx = ( Jp + Jm ) / 2;
Jy = ( Jp - Jm ) / (2*1i);

Iz = diag([5/2 3/2,1/2,-1/2,-3/2,-5/2,5/2,3/2,1/2,-1/2,-3/2,-5/2]);
Ip = [diag([sqrt(5), sqrt(8),3,sqrt(8), sqrt(5)],1),zeros(6);zeros(6),diag([sqrt(5), sqrt(8),3,sqrt(8), sqrt(5)],1)];
Im = [diag([sqrt(5), sqrt(8),3,sqrt(8), sqrt(5)],-1),zeros(6);zeros(6),diag([sqrt(5), sqrt(8),3,sqrt(8), sqrt(5)],-1)];
Ix = ( Ip + Im ) / 2;
Iy = ( Ip - Im ) / (2*1i);

Hz = muB*( gJ_5S1_2*Jz + gI_85*Iz )* Bdc(b); % Zeeman hamiltonian in this  basis
Hhfs = Ahfs85_5S1_2*( Ix*Jx + Iy*Jy + Iz*Jz ); % Hyperfine hamiltonian
H = Hhfs + Hz; % full hamiltonian
[ev,ew] = eig(H); % diagonalize hamiltonian
E85_5S1_2(:,b)=diag(ew);
% Matrices of the vector descriptions of each state at each field.
% A85_5S1_2(1,:) = lowest energy state (approaches mI=5/2, mJ=-1/2 at high field)
% A85_5S1_2(12,:) = highest energy state (approaches mI=5/2, mJ=1/2)
A85_5S1_2(:,:,b) = ev(:,:);

%%%%%% For 87Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice of basis: psi(1)=|mJ=1/2,mI=3/2>, psi(2)=|mJ=1/2,mI=1/2>,
% psi(3)=|mJ=1/2,mI=-1/2>,... psi(8)=|mJ=-1/2,mI=-3/2>
Jz = 1/2*[eye(4),zeros(4);zeros(4),-eye(4)];
Jp = [zeros(4),eye(4);zeros(4),zeros(4)];
Jm = [zeros(4),zeros(4);eye(4),zeros(4)];
Jx = ( Jp + Jm ) / 2;
Jy = ( Jp - Jm ) / (2*1i);

Iz = diag([3/2,1/2,-1/2,-3/2,3/2,1/2,-1/2,-3/2]);
Ip = [diag([sqrt(3),2,sqrt(3)],1),zeros(4);zeros(4),diag([sqrt(3),2,sqrt(3)],1)];
Im = [diag([sqrt(3),2,sqrt(3)],-1),zeros(4);zeros(4),diag([sqrt(3),2,sqrt(3)],-1)];
Ix = ( Ip + Im ) / 2;
Iy = ( Ip - Im ) / (2*1i);

Hz = muB*( gJ_5S1_2*Jz + gI_87*Iz )* Bdc(b); % Zeeman hamiltonian in this  basis
Hhfs = Ahfs87_5S1_2*( Ix*Jx + Iy*Jy + Iz*Jz ); % Hyperfine hamiltonian
H = Hhfs + Hz; % full hamiltonian
[ev,ew] = eig(H); % diagonalize hamiltonian
E87_5S1_2(:,b)=diag(ew);
% Matrices of the vector descriptions of each state at each field.
% A87_5S1_2(1,:) = lowest energy state (approaches mI=3/2, mJ=-1/2 at high field)
% A87_5S1_2(8,:) = highest energy state (approaches mI=3/2, mJ=1/2)
A87_5S1_2(:,:,b) = ev(:,:);

if strcmp(Dline,'D1')
%% Rb 5P1/2

%%%%%% For 85Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice of basis: psi(1)=|mJ=1/2,mI=5/2>, psi(2)=|mJ=1/2,mI=3/2>,
% psi(3)=|mJ=1/2,mI=1/2>,... psi(8)=|mJ=-1/2,mI=-5/2>
Jz = 1/2*[eye(6),zeros(6);zeros(6),-eye(6)];
Jp = [zeros(6),eye(6);zeros(6),zeros(6)];
Jm = [zeros(6),zeros(6);eye(6),zeros(6)];
Jx = ( Jp + Jm ) / 2;
Jy = ( Jp - Jm ) / (2*1i);

Iz = diag([5/2 3/2,1/2,-1/2,-3/2,-5/2,5/2,3/2,1/2,-1/2,-3/2,-5/2]);
Ip = [diag([sqrt(5), sqrt(8),3,sqrt(8), sqrt(5)],1),zeros(6);zeros(6),diag([sqrt(5), sqrt(8),3,sqrt(8), sqrt(5)],1)];
Im = [diag([sqrt(5), sqrt(8),3,sqrt(8), sqrt(5)],-1),zeros(6);zeros(6),diag([sqrt(5), sqrt(8),3,sqrt(8), sqrt(5)],-1)];
Ix = ( Ip + Im ) / 2;
Iy = ( Ip - Im ) / (2*1i);

Hz = muB*( gJ_5P1_2*Jz + gI_85*Iz )* Bdc(b); % Zeeman hamiltonian in this  basis
Hhfs = Ahfs85_5P1_2*( Ix*Jx + Iy*Jy + Iz*Jz ); % Hyperfine hamiltonian
H = Hhfs + Hz; % full hamiltonian
[ev,ew] = eig(H); % diagonalize hamiltonian
E85_5P1_2(:,b)=diag(ew);
% Matrices of the vector descriptions of each state at each field.
% A85_5S1_2(1,:) = lowest energy state (approaches mI=5/2, mJ=-1/2 at high field)
% A85_5S1_2(12,:) = highest energy state (approaches mI=5/2, mJ=1/2)
A85_5P1_2(:,:,b) = ev(:,:);

%%%%%% For 87Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice of basis: psi(1)=|mJ=1/2,mI=3/2>, psi(2)=|mJ=1/2,mI=1/2>,
% psi(3)=|mJ=1/2,mI=-1/2>,... psi(8)=|mJ=-1/2,mI=-3/2>

Jz = 1/2*[eye(4),zeros(4);zeros(4),-eye(4)];
Jp = [zeros(4),eye(4);zeros(4),zeros(4)];
Jm = [zeros(4),zeros(4);eye(4),zeros(4)];
Jx = ( Jp + Jm ) / 2;
Jy = ( Jp - Jm ) / (2*1i);

Iz = diag([3/2,1/2,-1/2,-3/2,3/2,1/2,-1/2,-3/2]);
Ip = [diag([sqrt(3),2,sqrt(3)],1),zeros(4);zeros(4),diag([sqrt(3),2,sqrt(3)],1)];
Im = [diag([sqrt(3),2,sqrt(3)],-1),zeros(4);zeros(4),diag([sqrt(3),2,sqrt(3)],-1)];
Ix = ( Ip + Im ) / 2;
Iy = ( Ip - Im ) / (2*1i);

Hz = muB*( gJ_5P1_2*Jz + gI_87*Iz )* Bdc(b); % Zeeman hamiltonian in this  basis
Hhfs = Ahfs87_5P1_2*( Ix*Jx + Iy*Jy + Iz*Jz ); % Hyperfine hamiltonian
H = Hhfs + Hz; % full hamiltonian
[ev,ew] = eig(H); % diagonalize hamiltonian
E87_5P1_2(:,b)=diag(ew);
% Matrices of the vector descriptions of each state at each field.
% A87_5P1_2(1,:,:) = lowest energy state (approaches mI=3/2, mJ=-1/2 at high field)
% A87_5P1_2(8,:,:) = highest energy state (approaches mI=3/2, mJ=1/2)
A87_5P1_2(:,:,b) = ev(:,:);


%% D1 transition strengths and frequencies:
% transition strengths calculated using <e|J+|g>, 
% with J+ modified using appropriate 3j and 6j matrices

%%%%%% For 85Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D1 Transition strengths. TransFreq85D1_X_array gives 
% the frequencies in the mI,mJ basis, which is useful for checking
% the validity of the results. For TransFreq85D1_X, the strengths are
% placed into a single vector, for export to Esus_func

% sigma+ transitions:
TransStrength85D1_p_array(:,:,b)=A85_5P1_2(:,:,b)'*J85D1p*A85_5S1_2(:,:,b);
TransStrength85D1_p_array(abs(TransStrength85D1_p_array)<TransStrengthThreshold)=NaN;
temp=TransStrength85D1_p_array(:,:,b);
select_p=~isnan(temp);
TransStrength85D1_p(:,b)=temp(select_p);
% sigma- transitions:
TransStrength85D1_m_array(:,:,b)=A85_5P1_2(:,:,b)'*J85D1m*A85_5S1_2(:,:,b);
TransStrength85D1_m_array(abs(TransStrength85D1_m_array)<TransStrengthThreshold)=NaN;
temp=TransStrength85D1_m_array(:,:,b);
select_m=~isnan(temp);
TransStrength85D1_m(:,b)=temp(select_m);
% pi transitions:
TransStrength85D1_z_array(:,:,b)=A85_5P1_2(:,:,b)'*J85D1z*A85_5S1_2(:,:,b);
TransStrength85D1_z_array(abs(TransStrength85D1_z_array)<TransStrengthThreshold)=NaN;
temp=TransStrength85D1_z_array(:,:,b);
select_z=~isnan(temp);
TransStrength85D1_z(:,b)=temp(select_z);

% D1 Transition frequencies:
% For a given Bdc, we scan over the transition strength matrix, and
% calculate the energies of g->e transitions with nonzero strength
% TransFreq85D1_X_array gives the frequencies in the mI,mJ basis, which is
% useful for checking the validity of the results
for ii=1:12
    for f=1:12
        if abs(TransStrength85D1_p_array(f,ii,b))>TransStrengthThreshold
            TransFreq85D1_p_array(f,ii,b) = E85_5P1_2(f,b)-E85_5S1_2(ii,b);
        end
        if abs(TransStrength85D1_m_array(f,ii,b))>TransStrengthThreshold
            TransFreq85D1_m_array(f,ii,b) = E85_5P1_2(f,b)-E85_5S1_2(ii,b);
        end
        if abs(TransStrength85D1_z_array(f,ii,b))>TransStrengthThreshold
%             sprintf('%g->%g',ii,f)
            TransFreq85D1_z_array(f,ii,b) = E85_5P1_2(f,b)-E85_5S1_2(ii,b);       
        end
    end
end
% Extracting frequencies from nonzero transitions and
% placing into a single vector for each Bdc value
temp_p=TransFreq85D1_p_array(:,:,b); 
select = isfinite(1./temp_p);
TransFreq85D1_p(:,b)=temp_p(select); 

temp_m=TransFreq85D1_m_array(:,:,b);
select = isfinite(1./temp_m);
TransFreq85D1_m(:,b)=temp_m(select);

temp_z=TransFreq85D1_z_array(:,:,b);
select = isfinite(1./temp_z);
TransFreq85D1_z(:,b)=temp_z(select); 

%%%%%% For 87Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D1 Transition strengths.
% sigma+ transitions:
TransStrength87D1_p_array(:,:,b)=A87_5P1_2(:,:,b)'*J87D1p*A87_5S1_2(:,:,b);
TransStrength87D1_p_array(abs(TransStrength87D1_p_array)<TransStrengthThreshold)=NaN;
temp=TransStrength87D1_p_array(:,:,b);
select_p=~isnan(temp);
TransStrength87D1_p(:,b)=temp(select_p);
% sigma- transitions:
TransStrength87D1_m_array(:,:,b)=A87_5P1_2(:,:,b)'*J87D1m*A87_5S1_2(:,:,b);
TransStrength87D1_m_array(abs(TransStrength87D1_m_array)<TransStrengthThreshold)=NaN;
temp=TransStrength87D1_m_array(:,:,b);
select_m=~isnan(temp);
TransStrength87D1_m(:,b)=temp(select_m);


% pi transitions:
TransStrength87D1_z_array(:,:,b)=A87_5P1_2(:,:,b)'*J87D1z*A87_5S1_2(:,:,b);
% TransStrength87D1_z_array(:,:,b)
TransStrength87D1_z_array(abs(TransStrength87D1_z_array)<TransStrengthThreshold)=NaN;
temp=TransStrength87D1_z_array(:,:,b);
select_z=~isnan(temp);
TransStrength87D1_z(:,b)=temp(select_z);

% D1 Transition frequencies:
for ii=1:8
    for f=1:8
        if abs(TransStrength87D1_p_array(f,ii,b))>TransStrengthThreshold
            TransFreq87D1_p_array(f,ii,b) = E87_5P1_2(f,b)-E87_5S1_2(ii,b);
        end
        if abs(TransStrength87D1_m_array(f,ii,b))>TransStrengthThreshold
            TransFreq87D1_m_array(f,ii,b) = E87_5P1_2(f,b)-E87_5S1_2(ii,b); 
        end
        if abs(TransStrength87D1_z_array(f,ii,b))>TransStrengthThreshold
%             sprintf('%g->%g',ii,f)        
            TransFreq87D1_z_array(f,ii,b) = E87_5P1_2(f,b)-E87_5S1_2(ii,b);
        end
    end
end
% Extracting nonzero transition strengths and frequencies,
% placing into a single vector for each Bdc value
temp_p=TransFreq87D1_p_array(:,:,b);
select = isfinite(1./temp_p); TransFreq87D1_p(:,b)=temp_p(select);
temp_m=TransFreq87D1_m_array(:,:,b);
select = isfinite(1./temp_m); TransFreq87D1_m(:,b)=temp_m(select);
temp_z=TransFreq87D1_z_array(:,:,b);       
select = isfinite(1./temp_z); TransFreq87D1_z(:,b)=temp_z(select); 


elseif strcmp(Dline,'D2')
%% 87 Rb 5P3/2 - D2 line
%%%%%% For 85Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice of basis: psi(1)=|mJ=3/2,mI=5/2>, psi(2)=|mJ=3/2,mI=3/2>,
% psi(3)=|mJ=3/2,mI=1/2>,... psi(24)=|mJ=-3/2,mI=-5/2>

Jz = kron(diag([3/2,1/2,-1/2,-3/2]),eye(6));
Jp = kron(diag([sqrt(3),2,sqrt(3)],1),eye(6));
Jm = kron(diag([sqrt(3),2,sqrt(3)],-1),eye(6));
% Jm=Jp';
Jx = ( Jp + Jm ) / 2;
Jy = ( Jp - Jm ) / (2*1i);

Iz = kron(eye(4),diag([5/2,3/2,1/2,-1/2,-3/2,-5/2]));
Ip = kron(eye(4),diag([sqrt(5),sqrt(8),3,sqrt(8),sqrt(5)],1));
Im = kron(eye(4),diag([sqrt(5),sqrt(8),3,sqrt(8),sqrt(5)],-1));
Ix = ( Ip + Im ) / 2;
Iy = ( Ip - Im ) / (2*1i);
IdJ=Ix*Jx+Iy*Jy+Iz*Jz; % I dot J

Hz = muB*( gJ_5P3_2*Jz + gI_85*Iz )* Bdc(b); % Zeeman hamiltonian in this  basis
Hhfs = Ahfs85_5P3_2*( Ix*Jx + Iy*Jy + Iz*Jz ) +...
    Bhfs85_5P3_2*(3*(Ix*Jx+Iy*Jy+Iz*Jz)^2+3/2*(Ix*Jx+Iy*Jy+Iz*Jz)-eye(24)*I85*(I85+1)*J_5P3_2*(J_5P3_2+1))/...
    (2*I85*(2*I85-1)*J_5P3_2*(2*J_5P3_2-1)); % Hyperfine hamiltonian
H = Hhfs + Hz; % full hamiltonian
[ev,ew] = eig(H); % diagonalize hamiltonian
E85_5P3_2(:,b)=diag(ew);
% Matrices of the vector descriptions of each state at each field.
% A87_5S1_2(1,:,:) = lowest energy state (approaches mI=3/2, mJ=-1/2 at high field)
% A87_5S1_2(16,:,:) = highest energy state (approaches mI=3/2, mJ=1/2)
% for k=1:16
A85_5P3_2(:,:,b) = ev(:,:);
% end


%%%%%% For 87Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choice of basis: psi(1)=|mJ=3/2,mI=3/2>, psi(2)=|mJ=3/2,mI=1/2>,
% psi(3)=|mJ=3/2,mI=-1/2>,... psi(16)=|mJ=-3/2,mI=-3/2>

Jz = kron(diag([3/2,1/2,-1/2,-3/2]),eye(4));
Jp = kron(diag([sqrt(3),2,sqrt(3)],1),eye(4));
Jm = kron(diag([sqrt(3),2,sqrt(3)],-1),eye(4));
Jx = ( Jp + Jm ) / 2;
Jy = ( Jp - Jm ) / (2*1i);

Iz = kron(eye(4),diag([3/2,1/2,-1/2,-3/2]));
Ip = kron(eye(4),diag([sqrt(3),2,sqrt(3)],1));
Im = kron(eye(4),diag([sqrt(3),2,sqrt(3)],-1));
Ix = ( Ip + Im ) / 2;
Iy = ( Ip - Im ) / (2*1i);
IdJ=Ix*Jx+Iy*Jy+Iz*Jz; % I dot J

Hz = muB*( gJ_5P3_2*Jz + gI_87*Iz )* Bdc(b); % Zeeman hamiltonian in this  basis
Hhfs = Ahfs87_5P3_2*IdJ +...
    Bhfs87_5P3_2*...
    (3*IdJ^2+3/2*IdJ-eye(16)*I87*(I87+1)*J_5P3_2*(J_5P3_2+1))/...
    (2*I87*(2*I87-1)*J_5P3_2*(2*J_5P3_2-1)); % Hyperfine hamiltonian
H = Hhfs + Hz; % full hamiltonian
[ev,ew] = eig(H); % diagonalize hamiltonian
E87_5P3_2(:,b)=diag(ew);
% Matrices of the vector descriptions of each state at each field.
% A87_5S1_2(1,:,:) = lowest energy state (approaches mI=3/2, mJ=-1/2 at high field)
% A87_5S1_2(16,:,:) = highest energy state (approaches mI=3/2, mJ=1/2)
for k=1:16
A87_5P3_2(:,k,b) = ev(:,k);
end


%% D2 transition strengths and frequencies:
% transition strengths calculated using <e|J+|g>, 
% with J+ modified using appropriate 3j and 6j matrices

%%%%%% For 85Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D2 Transition strengths. TransFreq85D2_X_array gives 
% the frequencies in the mI,mJ basis, which is useful for checking
% the validity of the results. For TransFreq85D2_X, the strengths are
% placed into a single vector, for export to Esus_func

% sigma+ transitions:
TransStrength85D2_p_array(:,:,b)=A85_5P3_2(:,:,b)'*J85D2p*A85_5S1_2(:,:,b);
TransStrength85D2_p_array(abs(TransStrength85D2_p_array)<TransStrengthThreshold)=NaN;
temp=TransStrength85D2_p_array(:,:,b);
select_p=~isnan(temp);
TransStrength85D2_p(:,b)=temp(select_p);
% sigma- transitions:
TransStrength85D2_m_array(:,:,b)=A85_5P3_2(:,:,b)'*J85D2m*A85_5S1_2(:,:,b);
TransStrength85D2_m_array(abs(TransStrength85D2_m_array)<TransStrengthThreshold)=NaN;
temp=TransStrength85D2_m_array(:,:,b);
select_m=~isnan(temp);
TransStrength85D2_m(:,b)=temp(select_m);
% pi transitions:
TransStrength85D2_z_array(:,:,b)=A85_5P3_2(:,:,b)'*J85D2z*A85_5S1_2(:,:,b);
TransStrength85D2_z_array(abs(TransStrength85D2_z_array)<TransStrengthThreshold)=NaN;
temp=TransStrength85D2_z_array(:,:,b);
select_z=~isnan(temp);
TransStrength85D2_z(:,b)=temp(select_z);

% D2 Transition frequencies:
% For a given Bdc, we scan over the transition strength matrix, and
% calculate the energies of g->e transitions with nonzero strength
% TransFreq85D2_X_array gives the frequencies in the mI,mJ basis, which is
% useful for checking the validity of the results
for ii=1:size(TransStrength85D2_p_array,2)
    for f=1:size(TransStrength85D2_p_array,1)
        if abs(TransStrength85D2_p_array(f,ii,b))>TransStrengthThreshold
            TransFreq85D2_p_array(f,ii,b) = E85_5P3_2(f,b)-E85_5S1_2(ii,b);
        end
        if abs(TransStrength85D2_m_array(f,ii,b))>TransStrengthThreshold
            TransFreq85D2_m_array(f,ii,b) = E85_5P3_2(f,b)-E85_5S1_2(ii,b);
        end
        if abs(TransStrength85D2_z_array(f,ii,b))>TransStrengthThreshold
%             sprintf('%g->%g',ii,f)
            TransFreq85D2_z_array(f,ii,b) = E85_5P3_2(f,b)-E85_5S1_2(ii,b);       
        end
    end
end
% Extracting frequencies from nonzero transitions and
% placing into a single vector for each Bdc value
temp_p=TransFreq85D2_p_array(:,:,b); 
select = isfinite(1./temp_p);
TransFreq85D2_p(:,b)=temp_p(select); 

temp_m=TransFreq85D2_m_array(:,:,b);
select = isfinite(1./temp_m);
TransFreq85D2_m(:,b)=temp_m(select);

temp_z=TransFreq85D2_z_array(:,:,b);
select = isfinite(1./temp_z);
TransFreq85D2_z(:,b)=temp_z(select); 

%%%%%% For 87Rb %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D2 Transition strengths.
% sigma+ transitions:
% TransStrength87D2_p_array(:,:,b)
TransStrength87D2_p_array(:,:,b)=A87_5P3_2(:,:,b)'*J87D2p*A87_5S1_2(:,:,b);
TransStrength87D2_p_array(abs(TransStrength87D2_p_array)<TransStrengthThreshold)=NaN;
temp=TransStrength87D2_p_array(:,:,b);
select_p=~isnan(temp);
TransStrength87D2_p(:,b)=temp(select_p);
% sigma- transitions:
TransStrength87D2_m_array(:,:,b)=A87_5P3_2(:,:,b)'*J87D2m*A87_5S1_2(:,:,b);
TransStrength87D2_m_array(abs(TransStrength87D2_m_array)<TransStrengthThreshold)=NaN;
temp=TransStrength87D2_m_array(:,:,b);
select_m=~isnan(temp);
TransStrength87D2_m(:,b)=temp(select_m);
% pi transitions:
TransStrength87D2_z_array(:,:,b)=A87_5P3_2(:,:,b)'*J87D2z*A87_5S1_2(:,:,b);
TransStrength87D2_z_array(abs(TransStrength87D2_z_array)<TransStrengthThreshold)=NaN;
temp=TransStrength87D2_z_array(:,:,b);
select_z=~isnan(temp);
TransStrength87D2_z(:,b)=temp(select_z);

% D2 Transition frequencies:
for ii=1:8
    for f=1:16
        if abs(TransStrength87D2_p_array(f,ii,b))>TransStrengthThreshold
            TransFreq87D2_p_array(f,ii,b) = E87_5P3_2(f,b)-E87_5S1_2(ii,b);
        end
        if abs(TransStrength87D2_m_array(f,ii,b))>TransStrengthThreshold
            TransFreq87D2_m_array(f,ii,b) = E87_5P3_2(f,b)-E87_5S1_2(ii,b); 
        end
        if abs(TransStrength87D2_z_array(f,ii,b))>TransStrengthThreshold
%             sprintf('%g->%g',ii,f)        
            TransFreq87D2_z_array(f,ii,b) = E87_5P3_2(f,b)-E87_5S1_2(ii,b);
        end
    end
end
% Extracting nonzero transition strengths and frequencies,
% placing into a single vector for each Bdc value
temp_p=TransFreq87D2_p_array(:,:,b);
% select = ~isnan(temp_p)
select = isfinite(1./temp_p);
TransFreq87D2_p(:,b)=temp_p(select);
temp_m=TransFreq87D2_m_array(:,:,b); select = isfinite(1./temp_m); 
TransFreq87D2_m(:,b)=temp_m(select);
temp_z=TransFreq87D2_z_array(:,:,b); select = isfinite(1./temp_z);
TransFreq87D2_z(:,b)=temp_z(select); 

end
end

%%%%%%%%% OUTPUT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5S1/2 eigenvectors and energies:
output.E85_5S1_2=E85_5S1_2; output.A85_5S1_2=A85_5S1_2;
output.E87_5S1_2=E87_5S1_2; output.A87_5S1_2=A87_5S1_2;


if strcmp(Dline,'D1')
    % 5P1/2 eigenvectors and energies:
    output.E85_5P1_2=E85_5P1_2; output.A85_5P1_2=A85_5P1_2;
    output.E87_5P1_2=E87_5P1_2; output.A87_5P1_2=A87_5P1_2;
    
    % Rb87 transition strengths and frequecies
    output.TransFreq87_p=TransFreq87D1_p/hplanck-Rb87_D1shift;
    output.TransFreq87_m=TransFreq87D1_m/hplanck-Rb87_D1shift;
    output.TransFreq87_z=TransFreq87D1_z/hplanck-Rb87_D1shift;
    output.TransStrength87_p=TransStrength87D1_p;
    output.TransStrength87_m=TransStrength87D1_m;
    output.TransStrength87_z=TransStrength87D1_z;
    % Arrays of Rb87 transition strengths and frequecies, in the mI,mJ basis
    output.TransFreq87_p_array=TransFreq87D1_p_array/hplanck-Rb87_D1shift;
    output.TransFreq87_m_array=TransFreq87D1_m_array/hplanck-Rb87_D1shift;
    output.TransFreq87_z_array=TransFreq87D1_z_array/hplanck-Rb87_D1shift;
    output.TransStrength87_p_array=TransStrength87D1_p_array;
    output.TransStrength87_m_array=TransStrength87D1_m_array;
    output.TransStrength87_z_array=TransStrength87D1_z_array;
    
    % Rb85 transition strengths and frequecies
    output.TransFreq85_p=TransFreq85D1_p/hplanck-Rb85_D1shift;
    output.TransFreq85_m=TransFreq85D1_m/hplanck-Rb85_D1shift;
    output.TransFreq85_z=TransFreq85D1_z/hplanck-Rb85_D1shift;
    output.TransStrength85_p=TransStrength85D1_p;
    output.TransStrength85_m=TransStrength85D1_m;
    output.TransStrength85_z=TransStrength85D1_z;
    % Arrays of Rb85 transition strengths and frequecies, in the mI,mJ basis
    output.TransFreq85_p_array=TransFreq85D1_p_array/hplanck-Rb85_D1shift;
    output.TransFreq85_m_array=TransFreq85D1_m_array/hplanck-Rb85_D1shift;
    output.TransFreq85_z_array=TransFreq85D1_z_array/hplanck-Rb85_D1shift;
    output.TransStrength85_p_array=TransStrength85D1_p_array;
    output.TransStrength85_m_array=TransStrength85D1_m_array;
    output.TransStrength85_z_array=TransStrength85D1_z_array;
 
    
elseif strcmp(Dline,'D2')
    % 5P3/2 eigenvectors and energies:
    output.E85_5P3_2=E85_5P3_2; output.A85_5P3_2=A85_5P3_2;
    output.E87_5P3_2=E87_5P3_2; output.A87_5P3_2=A87_5P3_2;
   
    % Rb87 transition strengths and frequecies
    output.TransFreq87_p=TransFreq87D2_p/hplanck-Rb87_D2shift;
    output.TransFreq87_m=TransFreq87D2_m/hplanck-Rb87_D2shift;
    output.TransFreq87_z=TransFreq87D2_z/hplanck-Rb87_D2shift;
    output.TransStrength87_p=TransStrength87D2_p;
    output.TransStrength87_m=TransStrength87D2_m;
    output.TransStrength87_z=TransStrength87D2_z;
    % Arrays of Rb87 transition strengths and frequecies, in the mI,mJ basis
    output.TransFreq87_p_array=TransFreq87D2_p_array/hplanck-Rb87_D2shift;
    output.TransFreq87_m_array=TransFreq87D2_m_array/hplanck-Rb87_D2shift;
    output.TransFreq87_z_array=TransFreq87D2_z_array/hplanck-Rb87_D2shift;
    output.TransStrength87_p_array=TransStrength87D2_p_array;
    output.TransStrength87_m_array=TransStrength87D2_m_array;
    output.TransStrength87_z_array=TransStrength87D2_z_array;
    
    % Rb85 transition strengths and frequecies
    output.TransFreq85_p=TransFreq85D2_p/hplanck-Rb85_D2shift;
    output.TransFreq85_m=TransFreq85D2_m/hplanck-Rb85_D2shift;
    output.TransFreq85_z=TransFreq85D2_z/hplanck-Rb85_D2shift;
    output.TransStrength85_p=TransStrength85D2_p;
    output.TransStrength85_m=TransStrength85D2_m;
    output.TransStrength85_z=TransStrength85D2_z;
    % Arrays of Rb85 transition strengths and frequecies, in the mI,mJ basis
    output.TransFreq85_p_array=TransFreq85D2_p_array/hplanck-Rb85_D2shift;
    output.TransFreq85_m_array=TransFreq85D2_m_array/hplanck-Rb85_D2shift;
    output.TransFreq85_z_array=TransFreq85D2_z_array/hplanck-Rb85_D2shift;
    output.TransStrength85_p_array=TransStrength85D2_p_array;
    output.TransStrength85_m_array=TransStrength85D2_m_array;
    output.TransStrength85_z_array=TransStrength85D2_z_array;
end
