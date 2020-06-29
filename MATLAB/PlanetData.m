% data for all celestial bodies in the Kerbol system (relative to celestial
% longitude and orbit of primary body

%% star data %%
%star = struct('name','Kerbol','mu',1.1723328*10^18,'eqr',261600000,'color',[0.9290,0.6940,0.1250]);
star = Body('Kerbol',261600000,1.1723328*10^18,[],Orbit(),[0.9290,0.6940,0.1250]);

%% planet data %%

name = {'Moho','Eve','Kerbin','Duna','Dres','Jool','Eeloo'};
a       = [5263138304, 9832684544, 13599840256, 20726155264, ... 
           40839348203, 68773560320, 90118820000];
e       = [0.2, 0.01, 0, 0.051, 0.145, 0.05, 0.26];
inc     = [7, 2.1, 0, 0.06, 5, 1.304, 6.15] * pi/180;
w       = [15, 0, 0, 0, 90, 0, 260] * pi/180;
Omega   = [70, 15, 0, 135.5, 280, 52, 50] * pi/180;
Mo      = [3.14, 3.14, 3.14, 3.14, 3.14, 0.1, 3.14];
eqr     = [250000, 700000, 600000, 320000, 138000, 6000000, 210000];
mu      = [1.6860938*10^11, 8.1717302*10^12, 3.5315984*10^12,...
           3.0136321*10^11, 2.1484489*10^10, 2.8252800*10^14,...
           7.4410815*10^10];
SOI     = [9646663.0, 85109365, 84159271, 47921949, 32832840,...
           2.4559852*10^9, 1.1908294*10^8];
color   = [[0.4,  0.4,  0.4 ]; 
           [0.6,  0.2,  0.4 ];
           [0,    0.5,  1   ];
           [0.9,  0.36, 0.27];
           [0.5,  0.5,  0.5 ];
           [0.2,  0.8,  0.2 ];
           [0.8,  0.8,  0.9 ]]';
              
planet(1:length(name)) = Body();     
for ii = 1:length(name)
    planet(ii) = Body(name{ii},eqr(ii),mu(ii),SOI(ii),Orbit(a(ii),e(ii),inc(ii),Omega(ii),w(ii),Mo(ii),star),color(:,ii));
end       

%% moon data %%

name = {'Gilly','Mun','Minmus','Ike','Laythe','Vall','Tylo','Bop','Pol'};
a       = [31500000, 12000000, 47000000, 3200000, 27184000, 43152000, ...
           68500000, 128500000, 179890000];
e       = [0.55, 0, 0, 0.03, 0, 0, 0, 0.235, 0.171];
inc     = [12, 0, 6, 0.2, 0, 0, 0.025, 15, 4.25] * pi/180;
w       = [10, 0, 38, 0, 0, 0, 0, 25, 15] * pi/180;
Omega   = [80, 0, 78, 0, 0, 0, 0, 10, 2] * pi/180;
Mo      = [0.9, 1.7, 0.9, 1.7, 3.14, 0.9, 3.14, 0.9, 0.9];
eqr     = [13000, 200000, 60000, 130000, 500000, 300000, 600000, ...
           65000, 44000];
mu      = [8289449.8, 6.5138398*10^10, 1.7658000*10^9, 1.8568369*10^10, ...
           1.9620000*10^12, 2.0748150*10^11, 2.8252800*10^12, ...
           2.4868349*10^9, 7.2170208*10^8];
SOI     = [126123.27, 2429559.1, 2247428.4, 1049598.9, 3723645.8, ...
           2406401.4, 10856518, 1221060.9, 1042138.9];
primary = [planet(2),planet(3),planet(3),planet(4),planet(6),planet(6),...
           planet(6),planet(6),planet(6)];
color   = [[0.6,  0.55,  0.5 ]; 
           [0.45, 0.45,  0.45]; 
           [0.75,   1,   0.9 ]; 
           [0.4,  0.4,   0.4 ];
           [0.3,  0.3,   0.55];
           [0.6,  0.6,   0.6 ];
           [0.9,  0.9,   0.9 ];
           [0.55, 0.5,   0.45];
           [0.9,  0.75,  0.55]]';

moon(1:length(name)) = Body();     
for ii = 1:length(name)
    moon(ii) = Body(name{ii},eqr(ii),mu(ii),SOI(ii),Orbit(a(ii),e(ii),inc(ii),Omega(ii),w(ii),Mo(ii),primary(ii)),color(:,ii));
end  
%% clear temp variables
clear name a e inc w Omega Mo eqr mu SOI color primary ii

%% test plotting planets and moons
% figure
% for ii=1:length(planet)
%     planet(ii).plotTrajectory(0,1000000,'o','none')
% end
% xlim([-1.2*planet(7).orbit.a,1.2*planet(7).orbit.a])
% ylim([-1.2*planet(7).orbit.a,1.2*planet(7).orbit.a])
% zlim([-1.2*planet(7).orbit.a,1.2*planet(7).orbit.a])
% figure
% for jj=1:length(moon)
%     if isequal(moon(jj).primary,planet(6))
%         moon(jj).plotTrajectory(0,100000,'o','none')
%     end
% end
% xlim([-1.2*moon(9).orbit.a,1.2*moon(9).orbit.a])
% ylim([-1.2*moon(9).orbit.a,1.2*moon(9).orbit.a])
% zlim([-1.2*moon(9).orbit.a,1.2*moon(9).orbit.a])
