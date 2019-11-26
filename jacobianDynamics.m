function dpdot_dp = jacobianDynamics(saturate, index)
    Nw=2;
    f=0.01;
    Iz=2667;
    a=1.35;
    b=1.45;
    By=0.27;
    Cy=1.2;
    Dy=0.7;
    Ey=-1.6;
    Shy=0;
    Svy=0;
    m=1400;
    g=9.806;

    syms x u y v psi_ r delta_f F_x
    assume(x, 'real')
    assume(u, 'real')
    assume(y, 'real')
    assume(v, 'real')
    assume(psi_, 'real')
    assume(r, 'real')
    assume(delta_f, 'real')
    assume(F_x, 'real')
    sym_vars = [x u y v psi_ r delta_f F_x];

    %slip angle functions in degrees
    a_f=(180/pi)*(delta_f-atan2(v+a*r,u));
    a_r=(180/pi)*(-atan2((v-b*r),u));

    %Nonlinear Tire Dynamics
    phi_yf=(1-Ey)*(a_f+Shy)+(Ey/By)*atan(By*(a_f+Shy));
    phi_yr=(1-Ey)*(a_r+Shy)+(Ey/By)*atan(By*(a_r+Shy));

    F_zf=b/(a+b)*m*g;
    F_yf=F_zf*Dy*sin(Cy*atan(By*phi_yf))+Svy;

    F_zr=a/(a+b)*m*g;
    F_yr=F_zr*Dy*sin(Cy*atan(By*phi_yr))+Svy;

    
    F_total=sqrt((Nw*F_x)^2+(F_yr^2));
    F_max=0.7*m*g;

    if saturate == 1
        F_x=F_max/F_total*F_x;

        F_yr=F_max/F_total*F_yr;
    end

    %vehicle dynamics
    dzdt_dvar= [u*cos(psi_)-v*sin(psi_);...
              (-f*m*g+Nw*F_x-F_yf*sin(delta_f))/m+v*r;...
              u*sin(psi_)+v*cos(psi_);...
              (F_yf*cos(delta_f)+F_yr)/m-u*r;...
              r;...
              (F_yf*a*cos(delta_f)-F_yr*b)/Iz];
          
    if index == 1 %Grabbing state jacobian
        dpdot_dp = [diff(dzdt_dvar, sym_vars(1)),diff(dzdt_dvar, sym_vars(2)),diff(dzdt_dvar, sym_vars(3)),diff(dzdt_dvar, sym_vars(4)),diff(dzdt_dvar, sym_vars(5)),diff(dzdt_dvar, sym_vars(6))];
    end
    if index == 2 %Grabbing input jacobian
        dpdot_dp = [diff(dzdt_dvar, sym_vars(7)),diff(dzdt_dvar, sym_vars(8))];
    end
    if index == 3 %Grabbing saturation condition
        dpdot_dp = F_total;
    end
end