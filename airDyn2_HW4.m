clc; close all; 
g = 1.4;                                 
beta = rad2deg(0:(pi/180):(pi/2));            
m = 0;    
                                         
for M1 = 1:0.2:5                        
    m = m+1;
    Nr = ((M1^2)*((sind(beta)).^2))-1;    
    Dr = ((g+(cosd(2*beta)))*M1^2)+2;
    theta = atand(2*cotd(beta).*Nr./Dr);
    a(m) = max(theta);                  
    b(m) = beta(theta==a(m));     
    plot(theta,beta,'-b')
    hold on
end
plot(a,b,'-r','Linewidth',1.5)
xlabel('\theta')
ylabel('\beta')
axis([0 58 0 90])