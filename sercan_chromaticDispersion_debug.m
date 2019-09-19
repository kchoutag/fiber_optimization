function [ng,nph,dispersion,dispersionerror,gd,gderror] =dispersion(alpha,horind,index_profile,dlamda,ro_core,ro_dip_inner,ro_dip_outer,n_center,n_dip,n_outer,m,lambda) 

% dlamda=0.01e-9; % m 
% lamda1=lamda-dlamda/2; 
% lamda2=lamda+dlamda/2; 
spl=299792458; % speed of light 
omega_center=2*pi*spl/lambda; 
s1=1; 
s2=1; 

% lamda_low=lamda_center-dlamda*s1; 
% lamda_high=lamda_center+dlamda*s2; 
% lamda=linspace(lamda_low,lamda_high,s1+s2+1); 

domega=2*pi*spl/lambda^2*dlamda; 
omega_low=omega_center-domega*s1; 
omega_high=omega_center+domega*s2; 
omega=linspace(omega_low,omega_high,s1+s2+1); 
nphase=zeros(1,s1+s2+1); 
beta=zeros(1,s1+s2+1); 
plotvar=0; 
% plotlabel='o'; 
refine=1; 
% [~,~,beta1]=normefield(beta,index_profile,lamda1,m,ro_core,ro_dip_inner,ro_dip_outer,n_center,n_dip,n_outer,refine,degree,plotvar,plotlabel,choosesign_center); 
% [~,~,beta2]=normefield(beta,index_profile,lamda2,m,ro_core,ro_dip_inner,ro_dip_outer,n_center,n_dip,n_outer,refine,degree,plotvar,plotlabel,choosesign_center); 
for k1=1:s1+s2+1 
    %     [~,~,beta(k1)]=normefield(beta0,index_profile,lamda(k1),m,ro_core,ro_dip_inner,ro_dip_outer,... 
    %         sellmeier(lamda(k1),n_center),sellmeier(lamda(k1),n_dip),sellmeier(lamda(k1),n_outer),... 
    %         refine,degree,plotvar,plotlabel,choosesign0); 
    %     nphase(k1)=phaseindex(beta(k1),lamda(k1)); 
    betatemp=solve_radial_mode_exact(index_profile,alpha,2*pi*spl/omega(k1),m,ro_core,ro_dip_inner,ro_dip_outer,... 
        sellmeier_ge_silica(2*pi*spl/omega(k1),n_center,0),sellmeier_fl_silica(2*pi*spl/omega(k1),n_dip,0),sellmeier_ge_silica(2*pi*spl/omega(k1),n_outer,0),... 
        refine,plotvar); 
    betatemp=sort(betatemp,'descend'); 
    betatemp 
    %     if(~isempty(betatemp)) 
    %         [~,horind]=min(abs(betatemp-beta0)); 
   if(horind<length(betatemp)) 
    beta(k1)=betatemp(horind); 
   end 
     
     
     
    %     else 
    %         beta(k1)=beta0; 
    %     end 
    nphase(k1)=phaseindex(beta(k1),2*pi*spl/omega(k1)); 
end 
% dbetaoverdlamda=(beta2-beta1)/(lamda2-lamda1); 
% dnoverdlamda=beta/2/pi+lamda/2/pi*dbetaoverdlamda; 
% ng=n-lamda*dnoverdlamda; 


% figure,plot(lamda,nphase,'o'); 
% figure,plot(lamda,beta,'o'); 
% p=polyfit(lamda,beta,3); 
% d=-1/2/pi/spl*lamda.^2.*(6*p(1)*lamda+2*p(2))-1/pi/spl*lamda.*(3*p(1)*lamda.^2+2*p(2)*lamda+p(3)); 
% figure,plot(lamda,d*1e6,'o'); 

% p1=4.89117496467341e-009; 
% p2=-79370.6452251412; 
% beta=beta-p1*omega-p2; 
% figure,plot(2*pi*spl*omega.^-1,beta); 
% figure,plot(omega,beta); 
[p,~,mu]=polyfit(omega,beta,1); 
% p=polyfit(omega,beta,1) 

dbetaoverdomega=p(1)/mu(2); 
dnoverdomega=-beta(s1+1)*spl/omega_center^2+spl/omega_center*dbetaoverdomega; 
nph=nphase(s1+1); 
ng=nph+omega_center*dnoverdomega; 

% d=-1/2/pi/spl*omega.^2.*(12*p(1)*omega.^2+6*p(2)*omega+2*p(3)); 
% figure,plot(2*pi*spl*omega.^-1,d*1e6,'o'); 

[pdispersion,~,mudispersion]=polyfit(omega,beta,1); 

figure; 
gd=pdispersion(1)/mudispersion(2); 
disp(['Group delay=',num2str(gd*1e12),' ps/m']); 
betaminusline=beta-gd*omega; 
plot(omega,betaminusline,'o'); 
title('\beta-(constant line) vs \omega'); 
xlabel('\omega'); 
ylabel('\beta-(constant line)'); 
hold on; 

[pdispersion2,~,mudispersion2]=polyfit(omega,betaminusline,2); 
omegacont=linspace(omega_low,omega_high,300); 
plot(omegacont,pdispersion2(1)*((omegacont-mudispersion2(1))/mudispersion2(2)).^2+pdispersion2(2)*(omegacont-mudispersion2(1))/mudispersion2(2)+pdispersion2(3),'r-'); 
dispersion=-omega_center/lambda*2*pdispersion2(1)/(mudispersion2(2))^2; 
gderror=abs(dispersion/omega_center*lambda)*(omega_high-omega_low); 
disp(['Group delay error=',num2str(gderror*1e12),' ps']); 
disp(['Dispersion=',num2str(dispersion*1e6),' ps/{nm km}']); 

[pdispersion3,~,mudispersion3]=polyfit(omega,betaminusline,3); 
d3betaoverdomega3=6*pdispersion3(1)/(mudispersion3(2))^3; 
dispersionerror=abs(omega_center/lambda*d3betaoverdomega3*(omega_high-omega_low)); 
disp(['Dispersion error=',num2str(dispersionerror*1e6),' ps/{nm km}']); 
disp(['Third derivative of beta=',num2str(d3betaoverdomega3)]); 
end