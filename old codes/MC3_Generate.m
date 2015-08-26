function [out] = MC3_Generate(mu1, mu2, Sigma11, Sigma12, p)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
minbeta=0.2;
beta_set=Select_temp(minbeta,mu1,mu2,Sigma11,Sigma12,p);
d_beta=length(beta_set);
MC3=cell(1,d_beta);
d=length(mu1);

x0=zeros(1,d);

for k=1:d_beta
    MC3{1,k}(1,:)=x0;
end
nsamples=4000;
i=2;
propvar=power(2.4,2)/d;
Sigma1=propvar*eye(d);
mu1=zeros(1,d);

counter=0;

while i<=nsamples
    %smpl=zeros(nsamples,d);
    for k=1:d_beta
        epsilon= mvnrnd(mu1, Sigma1,1);
        y= MC3{1,k}((i-1),:)+epsilon;
        hand=@(x) pdfmix(x,mu1,mu2,Sigma11,Sigma12,p);
        beta_hand= @(x) power(hand(x),k);
        accept_rate= min(1,(beta_hand(y)/beta_hand(MC3{1,k}((i-1),:))));
        w=rand(1);
        if w<accept_rate
            MC3{1,k}(i,:)=y;
            counter=counter+1;
        else
            MC3{1,k}(i,:)=MC3{1,k}((i-1),:);
        end
    end
    
    for k=2:d_beta
        beta_hand_1= @(x) power(hand(x),k);
        beta_hand_2= @(x) power(hand(x),k-1);
        u= MC3{1,k}(i,:);
        v= MC3{1,k-1}(i,:);
        swap_rate= min(1,(beta_hand_1(v)*beta_hand_2(u)/beta_hand_1(u)*beta_hand_2(v)));
        w=rand(1);
        if w<swap_rate
            temp=MC3{1,k}(i,:);
            MC3{1,k}(i,:)=MC3{1,k-1}(i,:);
            MC3{1,k-1}(i,:)=temp;
        end
    end
    
    i=i+1;
end

out=MC3{1,1};

    
end

