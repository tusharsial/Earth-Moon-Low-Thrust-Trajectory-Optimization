function [rM, thetaM, vrM, vthetaM] = ECRF_to_MCRF(rE, thetaE, vrE, vthetaE, D, omg)
    % Step 1: Convert Earth-centered polar to Cartesian
    xE = rE.*cos(thetaE);
    yE = rE.*sin(thetaE);

    vxM = vrE .*cos(thetaE) - vthetaE.*sin(thetaE);
    vyM = vrE.*sin(thetaE) + vthetaE.*cos(thetaE);

    % Step 2: Translate to Moon-centered Cartesian
    xM = xE + D;
    yM = yE;
    
    % Step 3: Convert Moon-centered Cartesian to polar
    rM = sqrt(xM.^2 + yM.^2);
    thetaM = atan2(yM, xM);

    vrM = vxM.*cos(thetaM) + vyM.*sin(thetaM);
    vthetaM = -vxM.*sin(thetaM) + vyM.*cos(thetaM);
   