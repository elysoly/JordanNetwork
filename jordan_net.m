function  [selected_a,selected_ar,selected_b,selected_bout,selected_z,selected_y]=jordan_net(X,val_data,ts_data,epsilon1,epsilon2,hidden_neuron,handles)
% clear;
% close all;
% clc;
% 
% [data,timestamps] = xlsread('daily-minimum-temperatures-in-me.xlsx', 'Daily minimum temperatures in M', 'A:B');
% l=length(data);
% data=(data-min(data))/(max(data)-min(data));
% 
% X=data(1:round(0.7*l),:);
% val_data=data(length(X)+1:length(X)+round(0.15*l),:);
% l_prime=length(X)+length(val_data)+1;
% ts_data=data(l_prime:l);



% epsilon1=0.5;
% epsilon2=0.5;
% hidden_neuron=3;

target_mean=mean(X);
threshold=str2num(handles.lr_threshold.String);
pause_time=str2num(handles.pause_time.String);
min_error=inf;

if(handles.mse.Value==1)
    measure_type='MSE';
end
if (handles.mae.Value==1)
        measure_type='MAE';
end
if (handles.rmae.Value==1)
    measure_type='RMAE';
end
if (handles.pi.Value==1)
    measure_type='PI';
    min_error=-inf;
end

        
a=ones(hidden_neuron,1).*rand(hidden_neuron,1); %a
a_r=ones(hidden_neuron,1).*rand(hidden_neuron,1);%a_r
b=zeros(hidden_neuron,2).*rand(hidden_neuron,2);%b
bias_hidden=ones(hidden_neuron,1);
bias_out=ones(hidden_neuron,1);

y=zeros(hidden_neuron,1).*rand(hidden_neuron,1);%y
z(1)=X(1);
selected_a=a(:,size(a,2));
selected_ar=a_r(:,size(a_r,2));
selected_b=b(:,size(b,2));
selected_bout=bias_out;
selected_bhidden=bias_hidden;
selected_z=z(1);
selected_y=y;
for j=1:hidden_neuron
    b(j,1)=0;
end



f_h=@(m) (1/(1+exp(-m)));
f_prime_h=@(m) (f_h(m))*(1-f_h(m));
f_out=@(m) (1/(1+exp(-m)));
f_prime_out=@(m) (f_out(m))*(1-f_out(m));
for i=1:hidden_neuron
    phi(i,1)=0;
end
epoch=1;
while(epoch<=threshold)
    for t=1:length(X)-1
        for i=1:hidden_neuron
            U(i)=a(i,t)*X(t) + a_r(i,t)*z(t); 
            y(i,t+1)=f_h(U(i));
        end
        v(t+1)=sum(b(:,t+1).*y(:,t+1)) + bias_out(i,t);
        z(t+1)=f_out(v(t+1));

        E(t+1)=(z(t+1)-X(t+1))^2/2;

        %calculate delta weights 
        p(t+1)=z(t+1)*(z(t+1)-X(t+1)) * (1-z(t+1));
        if(isnan(z(t+1)))
            fprintf('yes');
        end
        for i=1:hidden_neuron
            Q(i,t+1)=z(t+1)*(1-z(t+1))*b(i,t+1)*y(i,t+1)*(1-y(i,t+1));
            phi(i,t+1)=Q(i,t+1)*(z(t)+ a_r(i,t)*phi(i,t));
            delta_b(i,t+1) = -epsilon1*p(t+1)*y(i,t+1);
            delta_a(i,t)=-epsilon2*p(t+1)* b(i,t+1)*f_prime_h(U(i))*X(t);
             if(isnan(delta_a(i,t)))
                fprintf('yes');
            end
            delta_ar(i,t)=-epsilon2*(z(t+1)-X(t+1))*phi(i,t+1);
        end

        %%% do update
        b(:,t+2)=b(:,t+1)+delta_b(:,t+1);
        a(:,t+1)=a(:,t)+delta_a(:,t);
         if(isnan(a(:,t+1)))
            fprintf('yes');
        end
        a_r(:,t+1)=a_r(:,t)+delta_ar(:,t);
%         bias_hidden(:,t+1)=-epsilon2*p(t+1)* b(i,t+1)*f_prime_h(U(i));
        bias_out(:,t+1)=-epsilon1*p(t+1)+bias_out(:,t);
    end
%     tr_error(epoch)=sum(E);
    z=z.';
    switch measure_type
        case 'MAE'
            tr_error(epoch)=(sum(abs(z-X)))/length(X);
        case 'RMAE'
            tr_error(epoch)=((sum(abs(z-X)))/length(X))/target_mean;
        case 'MSE'
            tr_error(epoch)=sum((z-X).^2)/length(X);
        case 'PI'
            nom=sum((z-X).^2);
            denom=0;
            for jj=2:length(X)
                denom=denom+((X(jj)-X(jj-1))^2 );
            end
            tr_error(epoch)=1- (nom/denom);

    end
        
        
    z=z.';
    ll=length(a);
    [err,z_val,E_val]=test_jordan(measure_type,a(:,ll),a_r(:,ll),b(:,ll+1),bias_out(:,ll),val_data,z(ll),y(:,ll));
    val_error(epoch)=err;
    
    if( ( err<min_error && strcmp(measure_type,'PI')==0) || ( err>min_error && strcmp(measure_type,'PI')==1) )
        min_error=err;
        selected_a=a(:,size(a,2));
        selected_ar=a_r(:,size(a_r,2));
        selected_b=b(:,size(b,2));
        selected_bout=bias_out(:,size(bias_out,2));
        selected_bhidden=bias_hidden(:,size(bias_hidden,2));
        selected_z=z(ll);
        selected_y=y(:,ll);
    end
    %plot
    axes(handles.axes2);
    plot(X,'-r');
    hold on 
    plot(z,'-b');
    indx=(length(X)+1: length(X)+length(val_data));
    plot(indx,val_data,'--g')
    plot(indx,z_val,'--c');
    legend('tr data','tr prediction','val data','val prediction');
    title(strcat(' epoch ',num2str(epoch)));
    ylabel('output and target');
    hold off;
    
    axes(handles.axes3);
    plot(E,'-r');
    hold on
    plot(indx,E_val,'--g');
    xlabel('Time');
    ylabel('Error');
    legend('tr data error','val data error')
    hold off;
    

    axes(handles.axes6);
    plot(delta_b(1,:),'--r');
    hold on 
    axes(handles.axes5);
    plot(delta_a(1,:),'-r');
    hold on
    plot(delta_ar(1,:),'-y');
    legend_a(1,:)='delta a1';
    legend_a(2,:)='delta r1';
    legend_b(1,:)='delta b1';
    if(size(b,1)>1)
        axes(handles.axes6);
        plot(delta_b(2,:),'--b');
        axes(handles.axes5);
        plot(delta_a(2,:),'-b');
        plot(delta_ar(2,:),'-c');

        legend_a(3,:)='delta a2';
        legend_a(4,:)='delta r2';
        legend_b(2,:)='delta b2';
        if(size(b,1)>2)
            axes(handles.axes6);
            plot(delta_b(3,:),'--g');
            axes(handles.axes5);
            plot(delta_a(3,:),'-g');
            legend_a(5,:)='delta a3';
            legend_a(6,:)='delta r3';
            plot(delta_ar(3,:),'-m');
            legend_b(3,:)='delta b3';
        end
    end
    axes(handles.axes5);
    legend(legend_a);
    hold off;
    axes(handles.axes6);
    legend(legend_b);
    hold off;
    
    %%%jadidi
    epoch_a(:,epoch)=a(:,size(a,2));
    epoch_ar(:,epoch)=a_r(:,size(a_r,2));
    epoch_b(:,epoch)=b(:,size(b,2));
    epoch_bout(:,epoch)=bias_out(:,size(bias_out,2));
    epoch_bhidden(:,epoch)=bias_hidden(:,size(bias_hidden,2));
    
    
    pause(pause_time);
    
    
    
    epoch=epoch+1;
    tmp=size(b,2);
    b(:,2)=b(:,tmp);
    a(:,1)=a(:,tmp-1);
    a_r(:,1)=a_r(:,tmp-1);
    y(:,1)=y(:,tmp-1);
    phi(1)=phi(tmp-1);
    bias_out(:,1)=bias_out(:,tmp-1);
end
   
% axes(handles.axes2);
% plot(X,'-r');
% hold on 
% plot(z,'-b');
% legend('data','prediction');
% ylabel('output and target');
% title('Training data');

axes(handles.axes4);
% [~,ll]=min(val_error);
[ts_err,z_test,e_test]=test_jordan(measure_type,selected_a,selected_ar,selected_b,selected_bout,ts_data,selected_z,selected_y);
plot(ts_data,'-r');
hold on 
plot(z_test,'-b');
legend('data','prediction');
ylabel('output and target');
title('Test data');
label=handles.test_board.String;
handles.test_board.String=strcat(label, num2str(ts_err));
axes(handles.axes7);
 plot(e_test,'-r');
xlabel('Time');
ylabel('Error test');



axes(handles.axes1);
plot(tr_error,'-b');
hold on 
plot(val_error,'-r');
legend('train','validation');
xlabel('epoch');
ylabel(measure_type);



axes(handles.axes8);
Y=repmat([1:1:hidden_neuron],epoch-1,1);
X=repmat([1:1:epoch-1],size(Y,2),1).';
Z=epoch_a.';
if(hidden_neuron~=1)
    surf(X,Y,Z,'FaceColor','g');
else
    scatter3(X,Y,Z,'FaceColor','g')
end
hold on
Z=epoch_b.';
if(hidden_neuron~=1)
    surf(X,Y,Z,'FaceColor','b');
else
    scatter3(X,Y,Z,'FaceColor','b')
end

hold on
Z=epoch_ar.';
if(hidden_neuron~=1)
    surf(X,Y,Z,'FaceColor','r');
else
    scatter3(X,Y,Z,'FaceColor','r');
end
xlabel('x-epoch');
ylabel('y-neuron');
zlabel('z-weight');

legend('a','b','r');
hold off;

fprintf('end');
