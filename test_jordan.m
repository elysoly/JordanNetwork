function [err,z,E]=test_jordan(error_type,a,a_r,b,bias_out,T,z_0,y_0)
f_h=@(m) (1/(1+exp(-m)));
f_prime_h=@(m) (f_h(m))*(1-f_h(m));
f_out=@(m) (1/(1+exp(-m)));
f_prime_out=@(m) (f_out(m))*(1-f_out(m));
% X=(X-min(X))/(max(X)-min(X));

hidden_neuron=size(b,1);
z(1)=z_0;
y(1:hidden_neuron,1)=y_0;
for t=1:length(T)-1
        for i=1:hidden_neuron
            U(i)=a(i)*T(t) + a_r(i)*z(t); 
            y(i,t+1)=f_h(U(i));
        end
        v(t+1)=sum(b(:).*y(:,t+1)) + bias_out(i);
        z(t+1)=f_out(v(t+1));

        E(t+1)=(z(t+1)-T(t+1))^2/2;

end

z=z.';
switch error_type
    case 'MAE'
        err=(sum(abs(z-T)))/length(T);
    case 'RMAE'
        err=((sum(abs(z-T)))/length(T))/mean(T);
    case 'MSE'
        err=sum((z-T).^2)/length(T);
    case 'PI'
        nom=sum((z-T).^2);
        denom=0;
        for jj=2:length(T)
            denom=denom+((T(jj)-T(jj-1))^2 );
        end
        err=1- (nom/denom);

end
z=z.';

