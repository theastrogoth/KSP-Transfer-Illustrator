clear
clc
close all
PlanetData

%%%%% inputs
startParkOrbit = 100000;
endParkOrbit = 100000;
startPlanet = planet(3);
endPlanet = planet(4);
minStartTime = 0;
type = 'ballistic';
%%%%%

startOrbit = Orbit(startParkOrbit+startPlanet.eqr,0,0,0,0,0,startPlanet);
endOrbit = Orbit(endParkOrbit+endPlanet.eqr,0,0,0,0,0,endPlanet);
maxStartTime = minStartTime + 2*startPlanet.period;
maxFlightTime = sqrt((startPlanet.period^2 + endPlanet.period^2)/4);
minFlightTime = maxFlightTime/3;
%% generate transfer table
trTable = TransferTable(startOrbit,endOrbit,minStartTime,maxStartTime,minFlightTime,maxFlightTime,type);
%% prompt for transfer selection and draw trajectories
figure
chosenTransfer = trTable.getBestTransfer;
while true
    % porkchop plot
    ax1 = subplot(1,3,1);
    cla(ax1,'reset');
    trTable.porkchop(chosenTransfer)
    % ejection trajectory
    ax2 = subplot(1,3,2);
    cla(ax2,'reset');
    chosenTransfer.plotEjection('w')
    % transfer trajectory
    ax3 = subplot(1,3,3);
    cla(ax3,'reset');
    chosenTransfer.plotTransferOrbit('w')
    % link rotation in 3D plots
    hlink = linkprop([ax2,ax3],{'CameraPosition','CameraUpVector'});
    rotate3d on
    pause()
    % prompt for next transfer
    [x,y] = getpts(ax1);
    tStart = x*6*3600; % days (6h) to seconds
    tFlight = y*6*3600;
    if strcmp(trTable.type,'ballistic')
        chosenTransfer = Transfer(trTable.startOrbit,trTable.endOrbit,tStart,tFlight,false);
    elseif strcmp(trTable.type,'plane change')
        chosenTransfer = Transfer(trTable.startOrbit,trTable.endOrbit,tStart,tFlight,true);
    elseif strcmp(trTable.type,'optimal')
        bTransfer = Transfer(trTable.startOrbit,trTable.endOrbit,tStart,tFlight,false);
        bDV = norm(bTransfer.ejectionDV)+norm(bTransfer.insertionDV)+norm(bTransfer.planeChangeDV);
        pTransfer = Transfer(trTable.startOrbit,trTable.endOrbit,tStart,tFlight,true);
        pDV = norm(pTransfer.ejectionDV)+norm(pTransfer.insertionDV)+norm(pTransfer.planeChangeDV);
        if bDV<pDV
            chosenTransfer = bTransfer;
        else
            chosenTransfer = pTransfer;
        end
    else
        error('Invalid transfer type');
    end
end
