% Andrew Horsley 25/1/2016
% test script for debugging RbDlines_func

clear
constants; units

Dline='D2';
% static magnetic field Bdc along z is the quantization axis
Bdc = linspace(0,8000,500)*Gauss;
% Bdc=0;
output=RbDlines_func(Dline,Bdc);

E85_5S1_2=output.E85_5S1_2;
A85_5S1_2=output.A85_5S1_2;
E87_5S1_2=output.E87_5S1_2;
A87_5S1_2=output.A87_5S1_2;
if strcmp(Dline,'D1')
E85_5P=output.E85_5P1_2;
A85_5P=output.A85_5P1_2;
E87_5P=output.E87_5P1_2;
A87_5P=output.A87_5P1_2;
else
E85_5P=output.E85_5P3_2;
A85_5P=output.A85_5P3_2;
E87_5P=output.E87_5P3_2;
A87_5P=output.A87_5P3_2;
end

% Rb87 transition strengths and frequencies
TransFreq87_p=output.TransFreq87_p;
TransFreq87_m=output.TransFreq87_m;
TransFreq87_z=output.TransFreq87_z;
TransStrength87_p=output.TransStrength87_p;
TransStrength87_m=output.TransStrength87_m;
TransStrength87_z=output.TransStrength87_z;
% Arrays of Rb87 transition strengths and frequecies, in the mI,mJ basis
TransFreq87_p_array=output.TransFreq87_p_array;
TransFreq87_m_array=output.TransFreq87_m_array;
TransFreq87_z_array=output.TransFreq87_z_array;
TransStrength87_p_array=output.TransStrength87_p_array;
TransStrength87_m_array=output.TransStrength87_m_array;
TransStrength87_z_array=output.TransStrength87_z_array;

% Rb85 transition strengths and frequencies
TransFreq85_p=output.TransFreq85_p;
TransFreq85_m=output.TransFreq85_m;
TransFreq85_z=output.TransFreq85_z;
TransStrength85_p=output.TransStrength85_p;
TransStrength85_m=output.TransStrength85_m;
TransStrength85_z=output.TransStrength85_z;
% Arrays of Rb85 transition strengths and frequecies, in the mI,mJ basis
TransFreq85_p_array=output.TransFreq85_p_array;
TransFreq85_m_array=output.TransFreq85_m_array;
TransFreq85_z_array=output.TransFreq85_z_array;
TransStrength85_p_array=output.TransStrength85_p_array;
TransStrength85_m_array=output.TransStrength85_m_array;
TransStrength85_z_array=output.TransStrength85_z_array;
%%
fontsize=15;

figure(10)
hPlot = plot(Bdc/Gauss,E87_5P./hplanck/1e6);
xlabel('B_0 [G]');
ylabel('E/h [MHz]');
title('Level energies','FontWeight','Normal','FontSize',fontsize);
% xlim([0,500]);
% ylim([-20,20]);

figure(11)
plot(Bdc/Gauss,squeeze(A87_5P(7,:,:)));
% plot(Bdc(2:end)/Gauss,A1(:,2:end));
xlabel('B_0 [G]');
ylabel('Matrix element');
title('Eigenvectors','FontWeight','Normal','FontSize',fontsize);
% xlim([0,500]);
%%

figure(12)
data=TransFreq87_z;
plot(Bdc/Gauss,data(1,:)/GHz);
hold on
for a=2:size(data,1)
plot(Bdc/Gauss,data(a,:)/GHz);
end
hold off
xlabel('B_0 [G]');
ylabel('Transition Detuning (GHz)');
title('87Rb Transition Frequencies','FontWeight','Normal','FontSize',fontsize);
% xlim([0,500]);

figure(13)
data=TransStrength87_z.^2;
plot(Bdc/Gauss,data(1,:));
hold on
for a=2:size(data,1)
plot(Bdc/Gauss,data(a,:));
end
hold off    
xlabel('B_0 [G]');
ylabel('Transition Strength');
title('87Rb Transition Strengths','FontWeight','Normal','FontSize',fontsize);
% xlim([0,500]);

%%
figure(14)
data=TransFreq85_z;
plot(Bdc/Gauss,data(1,:)/GHz);
hold on
for a=2:size(data,1)
plot(Bdc/Gauss,data(a,:)/GHz);
end
hold off
xlabel('B_0 [G]');
ylabel('Transition Detuning (GHz)');
title('85Rb Transition Frequencies','FontWeight','Normal','FontSize',fontsize);
% xlim([0,500]);

figure(15)
data=TransStrength85_z;
plot(Bdc/Gauss,abs(data(1,:)));
hold on
for a=2:size(data,1)
plot(Bdc/Gauss,abs(data(a,:)));
end
hold off    
xlabel('B_0 [G]');
ylabel('Transition Strength');
title('85Rb Transition Strengths','FontWeight','Normal','FontSize',fontsize);
% xlim([0,500]);

%%
% 
% figure(12)
% plot(Bdc/Gauss,squeeze(TransFreq87_z_array(8,1,:))/GHz);
% hold on
% for b=1:8
%     for a=1:8
%     plot(Bdc/Gauss,squeeze(TransFreq87_z_array(b,a,:))/GHz);
%     end
% end
% hold off
% xlabel('B_0 [G]');
% ylabel('Transition Detuning (GHz)');
% title('Transition Frequencies','FontWeight','Normal','FontSize',fontsize);
% % xlim([0,500]);
% 
% figure(13)
% plot(Bdc/Gauss,squeeze(TransStrength87_z_array(3,5,:)));
% hold on
% for b=1:8
%     for a=1:8
%     plot(Bdc/Gauss,squeeze(TransStrength87_z_array(b,a,:)));
%     end
% end
% hold off    
% xlabel('B_0 [G]');
% ylabel('Transition Strength');
% title('Transition Strengths','FontWeight','Normal','FontSize',fontsize);
% % xlim([0,500]);
% 
% %%
% figure(14)
% plot(Bdc/Gauss,squeeze(TransFreq85_z_array(8,1,:))/GHz);
% hold on
% for b=1:12
%     for a=1:12
%     plot(Bdc/Gauss,squeeze(TransFreq85_z_array(b,a,:))/GHz);
%     end
% end
% hold off
% xlabel('B_0 [G]');
% ylabel('Transition Detuning (GHz)');
% title('Transition Frequencies','FontWeight','Normal','FontSize',fontsize);
% % xlim([0,500]);
% 
% figure(15)
% plot(Bdc/Gauss,squeeze(TransStrength85_z_array(3,5,:)));
% hold on
% for b=1:12
%     for a=1:12
%     plot(Bdc/Gauss,squeeze(TransStrength85_z_array(b,a,:)));
%     end
% end
% hold off    
% xlabel('B_0 [G]');
% ylabel('Transition Strength');
% title('Transition Strengths','FontWeight','Normal','FontSize',fontsize);
% % xlim([0,500]);
