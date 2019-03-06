function a_z_arr = a_z(V,c,g,zu,za,zt,zq,zug,zfug,Lg,zag,zfag,u,alpha,theta,q,u_g,a_g,a_fg,zd,sigmaug_V,sigmaag,elev,w1,w3)
    N = length(u); %length of measured data
    a_z_arr = zeros(N,1);
    for i = 1:N
        a_z_arr(i) = (V/c)*q(i) -(zu*u(i)+za*alpha(i)+zt*theta(i)+zq*q(i)+(zug-zfug*V/Lg*(c/V))*u_g(i)+zag*a_g(i)+zfag*(c/V)*a_fg(i)) - (zd*elev(i)+ zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg)*w1(i)+ zfag*(c/V)*sigmaag*sqrt(3*V/Lg)*w3(i)) ;
        a_z_arr(i) = V*a_z_arr(i)/g;
    end
end
