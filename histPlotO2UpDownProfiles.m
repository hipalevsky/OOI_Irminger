function histPlotO2UpDownProfiles(struct_in,fignum)

indup = find(struct_in.profile_direction > 0);
inddown = find(struct_in.profile_direction < 0);

O2grad_up = -5*diff(struct_in.O2conc(:,indup));
O2grad_down = 5*diff(struct_in.O2conc(:,inddown));
O2conc_up = struct_in.O2conc(:,indup);
O2conc_down = struct_in.O2conc(:,inddown);
        tol = 10;
    upcut = find(abs(O2grad_up) < tol); 
    downcut = find(abs(O2grad_down) < tol);
        minval = 240; maxval = 290;
    upcutconc = find(O2conc_up > minval & O2conc_up < maxval); 
    downcutconc = find(O2conc_down > minval & O2conc_down < maxval); 
    
figure(fignum); clf
    subplot(211)
histogram(O2conc_up(upcutconc)); hold on;
histogram(O2conc_down(downcutconc)); hold on;
title('O_2 concentration, \muM (raw)'); legend('Up','Down')
    subplot(212)
histogram(O2grad_up(upcut)); hold on;
histogram(O2grad_down(downcut)); hold on;
title('O_2 gradient (dO_2/dz), \muM/m'); legend('Up','Down')

end