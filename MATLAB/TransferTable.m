classdef TransferTable < handle
   properties (Access = public)
      startOrbit        % starting parking orbit
      endOrbit          % ending parking orbit (assume circular, not necessarily equatorial)
      minStartTime      % earliest departure time (s)
      maxStartTime      % latest departure time (s)
      minFlightTime     % fastest flight interval (s)
      maxFlightTime     % longest flight interval (s)
      transfers         % orbit between ejection and insertion
      type              % string denoting ballistic, plane change, or optimal transfers
      deltaV            % total delta v required for burns (m/s)
     
   end
   methods (Access = public, Static)
       function trtb = TransferTable(startOrbit,endOrbit,minStartTime,maxStartTime,minFlightTime,maxFlightTime,type)
           if nargin > 0
              trtb.startOrbit = startOrbit;
              trtb.endOrbit = endOrbit;
              trtb.minStartTime = minStartTime;
              trtb.maxStartTime = maxStartTime;
              trtb.minFlightTime = minFlightTime;
              trtb.maxFlightTime = maxFlightTime;
              if nargin>6
                trtb.type = type;
              else
                trtb.type = 'ballistic';
              end
              trtb.fillTable
           end
       end
   end
   methods (Access = public)
       function fillTable(trtb)
           tStart = linspace(trtb.minStartTime,trtb.maxStartTime,201);
           tFlight = linspace(trtb.minFlightTime,trtb.maxFlightTime,201);
           table(1:length(tStart),1:length(tFlight)) = Transfer();
           dvTable = zeros(length(tStart),length(tFlight));
           for ii = 1:length(tStart)
               for jj = 1:length(tFlight)
                   if strcmp(trtb.type,'ballistic')
                       table(ii,jj) = Transfer(trtb.startOrbit,trtb.endOrbit,tStart(ii),tFlight(jj),false);
                       dvTable(ii,jj) = norm(table(ii,jj).ejectionDV)+norm(table(ii,jj).insertionDV)+norm(table(ii,jj).planeChangeDV);
                   elseif strcmp(trtb.type,'plane change')
                       table(ii,jj) = Transfer(trtb.startOrbit,trtb.endOrbit,tStart(ii),tFlight(jj),true);
                       dvTable(ii,jj) = norm(table(ii,jj).ejectionDV)+norm(table(ii,jj).insertionDV)+norm(table(ii,jj).planeChangeDV);                       
                   elseif strcmp(trtb.type,'optimal')
                       bTransfer = Transfer(trtb.startOrbit,trtb.endOrbit,tStart(ii),tFlight(jj),false);
                       bDV = norm(bTransfer.ejectionDV)+norm(bTransfer.insertionDV)+norm(bTransfer.planeChangeDV);
                       pTransfer = Transfer(trtb.startOrbit,trtb.endOrbit,tStart(ii),tFlight(jj),true);
                       pDV = norm(pTransfer.ejectionDV)+norm(pTransfer.insertionDV)+norm(pTransfer.planeChangeDV);
                       if bDV<pDV
                           table(ii,jj) = bTransfer;
                           dvTable(ii,jj) = bDV; 
                       else
                           table(ii,jj) = pTransfer;
                           dvTable(ii,jj) = pDV;                            
                       end
                   else
                       error('Unrecognized transfer type');
                   end
               end
           end
           trtb.transfers = table;
           trtb.deltaV = dvTable;
       end
       function bestTransfer = getBestTransfer(trtb)
           [~,tStartIndex] = min(min(trtb.deltaV'));
           [~,tFlightIndex] = min(min(trtb.deltaV));
           bestTransfer = trtb.transfers(tStartIndex,tFlightIndex);
       end
       function porkchop(trtb,chosenTransfer)
           dvMin = min(min(trtb.deltaV));
           dvMax = max(max(trtb.deltaV));
           if dvMax > 20*dvMin
              dvMax = 20*dvMin; 
           end
           tStart = linspace(trtb.minStartTime,trtb.maxStartTime,size(trtb.deltaV,1));
           tFlight = linspace(trtb.minFlightTime,trtb.maxFlightTime,size(trtb.deltaV,2));
           tStartChosen = chosenTransfer.startTime;
           tFlightChosen = chosenTransfer.flightTime;
           % convert all times to days (6 hours)
           tStartDays = tStart/21600;
           tStartChosenDays = tStartChosen/21600;
           tFlightDays = tFlight/21600;
           tFlightChosenDays = tFlightChosen/21600;
           % set contour levels and draw plot
           levels = 10.^(linspace(log10(dvMin*1.01),log10(dvMax),25)); %logarithmic
           contour(tStartDays,tFlightDays,trtb.deltaV',levels);
           % indicate chosen point
           hold on
           plot(tStartChosenDays*ones(length(tFlightDays),1),tFlightDays,'-k',...
               tStartDays,tFlightChosenDays*ones(length(tStartDays)),'-k')
           plot(tStartChosenDays,tFlightChosenDays,'xk','LineWidth',2,'MarkerSize',10)
           % scale colormap
           cMap = jet(256);
           scalingIntensity = 4;
           xx = 1:length(cMap);
           xx = scalingIntensity * xx/max(xx);
           xx = exp(abs(xx)); %exponentiate
           xx = xx - min(xx); %set bottom back to 1 and max to 511
           xx = xx*511/max(xx)+1;
           newMap = interp1(xx, cMap, 1:512);
           colormap(newMap)
           % add labels
           caxis([dvMin levels(end-3)])
           grid on
           xlabel('Departure Date (day #)')
           ylabel('Time of Flight (days)')
           c = colorbar('Ticks',ceil(linspace(dvMin,levels(end-3)-1,8)));
           c.Label.String = '\Deltav (m/s)';
           title(sprintf('Flight from %s to %s',trtb.startOrbit.primary.name,trtb.endOrbit.primary.name));
           axis square
           hold off
       end
   end
end
