function rand_pt = surface_pt(cc_param,N)

R = cc_param.R;
L = cc_param.L;
shift = (L/2);

%{
syms d z
syms l(z)
l(z) = piecewise(-(shift+R)<=z<-L/2, sqrt(R.^2-(z+shift).^2),...
    -L/2<=z<L/2, R,...
    L/2<=z<=(shift+R),sqrt(R.^2-(z-shift).^2));
%}

rand_num = unifrnd(-(shift+R),(shift+R),N,1);
rand_angle = unifrnd(0,2*pi,N,1);
rand_pt = zeros(N,3);
%rand_pt = double([l(rand_num).*cos(rand_angle),l(rand_num).*sin(rand_angle),rand_num]);

for ii = 1:N
    r = l(rand_num(ii));
    rand_pt(ii,:) = [r*cos(rand_angle(ii)),r*sin(rand_angle(ii)),rand_num(ii)];
end

%{

plot3(rand_pt(:,1),rand_pt(:,2),rand_pt(:,3),'k.')
daspect([1 1 1])
%}

    function r = l(z)
        if abs(z)<= L/2
            r = R;
            return
        end
        
        if z>L/2
            r = sqrt(R.^2-(z-L/2).^2);
            return
        end
        
        if z<-L/2
            r = sqrt(R.^2-(z+L/2).^2);
            return
        end
        
    end



end