function x_norm = xNorm(XYZ,cc_param,factor)

R = cc_param.R;
L = cc_param.L;

R_0 = factor*R;
L_0 = factor*L;

len = size(XYZ,1);
x_norm = zeros(len,1);

for ii = 1:len
    
    z = XYZ(ii,3);
    if abs(z)< (L_0/2 + R_0)
    x_norm(ii) = XYZ(ii,1)./(l(z));
    else 
        x_norm(ii) = XYZ(ii,1);
    end
    
end

    function r = l(z)
        if abs(z)<= L_0/2
            r = R_0;
            return
        end
        
        if z>L_0/2
            r = sqrt(R_0.^2-(z-L_0/2).^2);
            return
        end
        
        if z<-L_0/2 
            r = sqrt(R_0.^2-(z+L_0/2).^2);
            return
        end
        
    end


end