function main

% x = 0:.01:1;
% 
% plot(x,f(x))
% 
% I = integrate_trap(@f,0,1,100)


x=0.001:.01:1.011;
hold on

%area(2*(.1-x)./x,2*(1-x)./x)

plot(x,2*(.1-x)./x,x,2*(1-x)./x,x,2*(.1-x)./x)
plot(x,.1./x)

axis([0 1 0 1]);

end


function [z] = f(y)

a=0.1;
w=0.2;
beta = 4.5708;

z = 2*y.*(2*(1-a)+y).^2 - log(2*a./(2+y)).*(4*a*w*(2+y) + 2*w*(w+y).^2).*(1-y).^(beta-1);
z = z./((4-y).^2.*(2+y)*beta);



end

function [I] = integrate_trap(f,a,b,N)

    x = linspace(a,b,N);
    I = (b-a)*sum(f(x(1:end-1)) + f(x(2:end)))/2/N;
    
end

