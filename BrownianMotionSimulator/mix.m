function XYZ = mix(N,ratio,cc_param)
N_surf = round(N*ratio);
N_in = round(N*(1-ratio));
XYZ_surf = surface_pt(cc_param,N_surf);
XYZ_in = internal_pt(cc_param,N_in);
XYZ = [XYZ_surf;XYZ_in];

%{
plot3(XYZ_surf(:,1),XYZ_surf(:,2),XYZ_surf(:,3),'r.')
hold on
plot3(XYZ_in(:,1),XYZ_in(:,2),XYZ_in(:,3),'b.')
daspect([1 1 1])

%}

end

    
