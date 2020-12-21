clear; close all;beep off;
clc;

%% load datas
testbruit=0;
testbiais=0;
plot_BOP=0;

n_var=1
biais_tab=0;
sigma_tab=0;
if testbruit==1
sigma_tab=0:1000:5000;
n_var=length(sigma_tab);
biais_tab=zeros(1,n_var);

end

if testbiais==1
biais_tab=[100 500 ];
n_var=length(biais_tab);
sigma_tab=zeros(1,n_var);
nb_satellites_biaised=0;
end

nb_tot_sat_biaised=1;
for satellite_baise=1:nb_tot_sat_biaised
nb_satellites_biaised=satellite_baise;
MSE_tab=zeros(1,n_var);
for i = 1:n_var
data=load('donnees_GPS_TP.mat');
trajectoire=load('trajectoire_TP.mat');

% 8 Satellites, 954 instants
XYZsat=data.XYZsat;  %contient les positions satellites (8 mesures GPS à chaque instant = Y(t)
[nb_sat,nb_coord,nb_instants]=size(XYZsat);
PRN=data.PRN;          %PRN contient les pseudos-distances [CANAUX DE POURSUITE X INSTANTS]
trajectoire_vehicule=trajectoire.Xloc; % CE QU'ON DOIT RETROUVER
% figure,plot(trajectoire_vehicule(1,:),trajectoire_vehicule(2,:))

%% Paramètres pour ajouter des perturbations
biais_added=biais_tab(i);
sigma2_bruit=sigma_tab(i); %graph plot 100 et 1000 et plot RMSE avec (1,10,10,100,1000)

%On ajoute du bruit si y'en a
if sigma2_bruit~=0
    for n_sat=1:nb_sat
        for t=1:nb_instants
            if PRN(n_sat,t) ~= 0
                noise=sqrt(sigma2_bruit)*randn;
                PRN(n_sat,t) = PRN(n_sat,t) + noise;
            end
        end
    end
end

%On ajoute un biais si y'en a
if biais_added~=0
    for n_sat=1:nb_satellites_biaised
        for t=1:nb_instants
            if PRN(n_sat,t) ~= 0
                PRN(n_sat,t) = PRN(n_sat,t) + biais_added;
            end
        end
    end    
end

%% 2. Po repère 2 => Po repère 1
%latitude, longitude
lat_lambda=(44+48/60)*(pi/180);
long_phi= (-35/60)*(pi/180);
hauteur=0;

Po2_rad=[lat_lambda;long_phi;hauteur]; % Po repère 2 NED en deg

%%coordonnées cartésienne Po dans repere 1 (en m)
[Po1_cart] = llh2xyz(Po2_rad);

%Matrice de passage du repère 1 Terre au repère 2 Talence
[M] = matrice_passage(lat_lambda,long_phi);

%On complete la premiere estimation avec P0
X_estime = zeros(4,nb_instants+1);
X_estime(1) = Po1_cart(1);
X_estime(2) = Po1_cart(2);
X_estime(3) = Po1_cart(3);
biais_o_clock=0;
X_estime(4) = biais_o_clock;

%Vecteur pour stocker les mesures d'erreurs
DOP_tab = zeros(1,nb_instants);
PDOP_tab = zeros(1,nb_instants);
VDOP_tab = zeros(1,nb_instants);
HDOP_tab = zeros(1,nb_instants);

%% On retrouve la trajectoire à partir des observations des GPS

for t = 1:nb_instants
    
    observations_disponibles = PRN(:,t);
    i_sat = find(observations_disponibles);
    
    H=zeros(length(i_sat),nb_coord);
    r=zeros(length(i_sat),1);
    
    H_rep_2 = zeros(length(i_sat),nb_coord);

    
    for n_sat = 1:length(i_sat)
        
        r(n_sat) = sqrt((X_estime(1,t)-XYZsat(i_sat(n_sat),1,t)).^2 + (X_estime(2,t)-XYZsat(i_sat(n_sat),2,t)).^2 + (X_estime(3,t)-XYZsat(i_sat(n_sat),3,t)).^2 )+biais_o_clock;
        
        H(n_sat,1) = (X_estime(1,t)-XYZsat(i_sat(n_sat),1,t))/r(n_sat);
        H(n_sat,2) = (X_estime(2,t)-XYZsat(i_sat(n_sat),2,t))/r(n_sat);
        H(n_sat,3) = (X_estime(3,t)-XYZsat(i_sat(n_sat),3,t))/r(n_sat);
        H(n_sat,4) = 1;
        
        %H repère 2 pour les calculs mesures erreur
        Talence_loc = M'*(([XYZsat(i_sat(n_sat),1,t) XYZsat(i_sat(n_sat),2,t) XYZsat(i_sat(n_sat),3,t)])'-Po1_cart');
        
        X_estime_dans_Talence = M'*(([X_estime(1,t) X_estime(2,t) X_estime(3,t)])'-Po1_cart');
        
        r_rep_2 = sqrt((X_estime_dans_Talence(1)-Talence_loc(1)).^2 + (X_estime_dans_Talence(2)-Talence_loc(2)).^2 + (X_estime_dans_Talence(3)-Talence_loc(3)).^2 );
        
        H_rep_2(n_sat,1) = (X_estime_dans_Talence(1)-Talence_loc(1))/r_rep_2;
        H_rep_2(n_sat,2) = (X_estime_dans_Talence(2)-Talence_loc(2))/r_rep_2;
        H_rep_2(n_sat,3) = (X_estime_dans_Talence(3)-Talence_loc(3))/r_rep_2;
        H_rep_2(n_sat,4) = 1;
        
        
        
    end
    
    H_espagnol= inv(H'*H);
    H_espagnol_Talence = inv(H_rep_2'*H_rep_2);
    
    % BOP à la DaBaby
    DOP_tab(t) = sqrt(trace(H_espagnol));
    PDOP_tab(t) = sqrt(H_espagnol(1,1)+H_espagnol(2,2)+H_espagnol(3,3));
    VDOP_tab(t) = sqrt(H_espagnol_Talence(3,3));
    HDOP_tab(t) = sqrt(H_espagnol_Talence(1,1)+H_espagnol_Talence(2,2));
    
    biais_o_clock = X_estime(4,t); 
    Z = PRN(i_sat',t) - r+  + H*X_estime(:,t);    
    X_estime(:,t+1) = pinv(H)*Z;
    
end

%% Trajectoire qu'on obtient dans le repère 1 Terre ECEF
figure,plot3(X_estime(1,:),X_estime(2,:),X_estime(3,:));
% hold on;plot3(X_estime(1,587:684),X_estime(2,587:684),X_estime(3,587:684),'-k');

title('Trajectoire repere 1 Terre ECEF')

%% Matrice de passage M: On repasse dans le repère 2 Talence NED
X_estime_repere2 = zeros(3,nb_instants+1);
for j = 1:nb_instants+1
    X_estime_repere2(:,j) = M\(X_estime(1:3,j)-Po1_cart');
end

%% On compare la vrai trajectoire et celle qu'on a obtenue dans le repère 2 Talence NED
figure,
plot(X_estime_repere2(1,:),X_estime_repere2(2,:));
% title("Trajectoire repere 2 Talence NED pour sigma^2_b_r_u_i_t = "+sigma2_bruit)
title("Trajectoire NED "+nb_satellites_biaised+" satellites biased avec biais = "+biais_added)
hold on;plot(trajectoire_vehicule(1,:),trajectoire_vehicule(2,:));
% hold on;plot(trajectoire_vehicule(1,587:684),trajectoire_vehicule(2,587:684),'-k');
legend('trajectoire estimée','Trajectoire réelle');

%% Erreur quadratique moyenne
X_estime_repere2(:,1)=[];
MSE=mean(mean((abs(X_estime_repere2-trajectoire_vehicule)).^2)); %faire varier bruit de 1 10 100
MSE_tab(i)=MSE;

end

%% BIAIS et MSE et tout
% figure,
% hold on;(plot(biais_tab,MSE_tab));
% figure,(plot(sigma_tab,MSE_tab));
% hold on,scatter(sigma_tab,MSE_tab,'r+');
% title('MSE en fonction de la variance du bruit');
% title('MSE en fonction du biais ajouté et du nombre de satellites biaisé');
% ylabel('MSE')
% xlabel('biais ajouté')
% legend("sat biaised="+nb_satellites_biaised+"sat biaised="+nb_satellites_biaised+"sat biaised="+nb_satellites_biaised);

end

% legend("satellites biaisés = 1","satellites biaisés = 2","satellites biaisés = 3","satellites biaisés = 4","satellites biaisés = 5","satellites biaisés = 6");


if plot_BOP==1
%% Question 6 : Facteur DOP

% DOP
figure(3)
plot(DOP_tab)
title("Variation du facteur DOP à chaque instant")
xlabel("Numero d'échantillon")
ylabel('Valeur du facteur DOP')
grid on

% HDOP
figure(4)
plot(HDOP_tab)
title("Variation du facteur HDOP à chaque instant")
xlabel("Numero d'échantillon")
ylabel('Valeur du facteur HDOP')
grid on


% PDOP
figure(5)
plot(PDOP_tab)
title("Variation du facteur PDOP à chaque instant")
xlabel("Numero d'échantillon")
ylabel('Valeur du facteur PDOP')
grid on

% VDOP
figure(6)
plot(VDOP_tab)
title("Variation du facteur VDOP à chaque instant")
xlabel("Numero d'échantillon")
ylabel('Valeur du facteur VDOP')
grid on

end
