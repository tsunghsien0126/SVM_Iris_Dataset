clear;
close all;
clc;

color = lines(6);

s=readmatrix("iris.txt");
Tr=s([51:75,101:125],3:5); %training data

for i=1:50
    if(Tr(i,3)==2)
        Tr(i,3)=1;
    else
        Tr(i,3)=-1;
    end
end

Te=s([76:100,126:150],3:5); %test data
for i=1:50
    if(Te(i,3)==2)
        Te(i,3)=1;
    else
        Te(i,3)=-1;
    end
end

%% Sigma = 5
C = 10;
sigma = 5;
f = ones(1,50).*(-1);
lb = zeros(50,1);
ub = ones(50,1).*C;
Aeq = Tr(:,3)';
beq = 0;
A = [];
b = [];

H = zeros(50,50);
for i = 1:50
    for j = 1:50
        H(j,i) = Tr(i,3)*Tr(j,3)*(exp(((norm(Tr(j,1:2)-Tr(i,1:2))/sigma)^2)/(-2)));
    end
end

alphaOri = quadprog(H,f,A,b,Aeq,beq,lb,ub);
alpha = alphaOri;
alpha(alpha<   sqrt(eps) ) = 0;
alpha(alpha>(C-sqrt(eps))) = C;
alpha = round(alpha,6);

AlphaSum(1,1) = sum(alpha);

for i = 1:50
    if (alpha(i)<C) && (alpha(i)>0)
        ak = [i,alpha(i)];
        break
    end
end

Kb = 0;
for i = 1:50
    Kb = Kb+(alpha(i)*Tr(i,3)*exp(((norm(Tr(i,1:2)-Tr(ak(1),1:2))/sigma)^2)/(-2)));
end

bias = (1/Tr(ak(1),3))-Kb;
Bias(1,1) = round(bias,4);

 % testing 
Dt = zeros(50,1);
for i = 1:50
    for j = 1:50
        Dt(i) = Dt(i)+(alpha(j)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Te(i,1:2))/sigma)^2)/(-2)));
    end
end
Dt = Dt+bias;
for i = 1:50
    if Dt(i)>=0
        Dt(i) = 1;
    else
        Dt(i) = -1;
    end
end

CR1 = 0;
for i = 1:50
    if Dt(i) == Te(i,3)
        CR1 = CR1+1;
    end
end
CR1 = CR1/50;
Alpha(:,1) = alpha;

%% Sigma = 1
C = 10;
sigma = 1;
f = ones(1,50).*(-1);
lb = zeros(50,1);
ub = ones(50,1).*C;
Aeq = Tr(:,3)';
beq = 0;
A = [];
b = [];

H = zeros(50,50);
for i = 1:50
    for j = 1:50
        H(j,i) = Tr(i,3)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Tr(i,1:2))/sigma)^2)/(-2));
    end
end

alphaOri = quadprog(H,f,A,b,Aeq,beq,lb,ub);
alpha = alphaOri;
alpha(alpha<   sqrt(eps) ) = 0;
alpha(alpha>(C-sqrt(eps))) = C;
alpha = round(alpha,6);

AlphaSum(1,2) = sum(alpha);

for i = 1:50
    if (alpha(i)<C) && (alpha(i)>0)
        ak = [i,alpha(i)];
        break
    end
end

Kb = 0;
for i = 1:50
    Kb = Kb+(alpha(i)*Tr(i,3)*exp(((norm(Tr(i,1:2)-Tr(ak(1),1:2))/sigma)^2)/(-2)));
end

bias = (1/Tr(ak(1),3))-Kb;
Bias(1,2) = round(bias,4);

 % testing 
Dt = zeros(50,1);
for i = 1:50
    for j = 1:50
        Dt(i) = Dt(i)+(alpha(j)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Te(i,1:2))/sigma)^2)/(-2)));
    end
end
Dt = Dt+bias;
for i = 1:50
    if Dt(i)>=0
        Dt(i) = 1;
    else
        Dt(i) = -1;
    end
end

CR2 = 0;
for i = 1:50
    if Dt(i) == Te(i,3)
        CR2 = CR2+1;
    end
end
CR2 = CR2/50;
Alpha(:,2) = alpha;

%% Sigma = 0.5
C = 10;
sigma = 0.5;
f = ones(1,50).*(-1);
lb = zeros(50,1);
ub = ones(50,1).*C;
Aeq = Tr(:,3)';
beq = 0;
A = [];
b = [];

H = zeros(50,50);
for i = 1:50
    for j = 1:50
        H(j,i) = Tr(i,3)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Tr(i,1:2))/sigma)^2)/(-2));
    end
end

alphaOri = quadprog(H,f,A,b,Aeq,beq,lb,ub);
alpha = alphaOri;
alpha(alpha<   sqrt(eps) ) = 0;
alpha(alpha>(C-sqrt(eps))) = C;
alpha = round(alpha,6);

AlphaSum(1,3) = sum(alpha);

for i = 1:50
    if (alpha(i)<C) && (alpha(i)>0)
        ak = [i,alpha(i)];
        break
    end
end

Kb = 0;
for i = 1:50
    Kb = Kb+(alpha(i)*Tr(i,3)*exp(((norm(Tr(i,1:2)-Tr(ak(1),1:2))/sigma)^2)/(-2)));
end

bias = (1/Tr(ak(1),3))-Kb;
Bias(1,3) = round(bias,4);

 % testing 
Dt = zeros(50,1);
for i = 1:50
    for j = 1:50
        Dt(i) = Dt(i)+(alpha(j)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Te(i,1:2))/sigma)^2)/(-2)));
    end
end
Dt = Dt+bias;
for i = 1:50
    if Dt(i)>=0
        Dt(i) = 1;
    else
        Dt(i) = -1;
    end
end

CR3 = 0;
for i = 1:50
    if Dt(i) == Te(i,3)
        CR3 = CR3+1;
    end
end
CR3 = CR3/50;
Alpha(:,3) = alpha;

%% Sigma = 0.1
C = 10;
sigma = 0.1;
f = ones(1,50).*(-1);
lb = zeros(50,1);
ub = ones(50,1).*C;
Aeq = Tr(:,3)';
beq = 0;
A = [];
b = [];

H = zeros(50,50);
for i = 1:50
    for j = 1:50
        H(j,i) = Tr(i,3)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Tr(i,1:2))/sigma)^2)/(-2));
    end
end

alphaOri = quadprog(H,f,A,b,Aeq,beq,lb,ub);
alpha = alphaOri;
alpha(alpha<   sqrt(eps) ) = 0;
alpha(alpha>(C-sqrt(eps))) = C;
alpha = round(alpha,6);

AlphaSum(1,4) = sum(alpha);

for i = 1:50
    if (alpha(i)<C) && (alpha(i)>0)
        ak = [i,alpha(i)];
        break
    end
end

Kb = 0;
for i = 1:50
    Kb = Kb+(alpha(i)*Tr(i,3)*exp(((norm(Tr(i,1:2)-Tr(ak(1),1:2))/sigma)^2)/(-2)));
end

bias = (1/Tr(ak(1),3))-Kb;
Bias(1,4) = round(bias,4);

 % testing 
Dt = zeros(50,1);
for i = 1:50
    for j = 1:50
        Dt(i) = Dt(i)+(alpha(j)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Te(i,1:2))/sigma)^2)/(-2)));
    end
end
Dt = Dt+bias;
for i = 1:50
    if Dt(i)>=0
        Dt(i) = 1;
    else
        Dt(i) = -1;
    end
end

CR4 = 0;
for i = 1:50
    if Dt(i) == Te(i,3)
        CR4 = CR4+1;
    end
end
CR4 = CR4/50;
Alpha(:,4) = alpha;

%% Sigma = 0.05
C = 10;
sigma = 0.05;
f = ones(1,50).*(-1);
lb = zeros(50,1);
ub = ones(50,1).*C;
Aeq = Tr(:,3)';
beq = 0;
A = [];
b = [];

H = zeros(50,50);
for i = 1:50
    for j = 1:50
        H(j,i) = Tr(i,3)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Tr(i,1:2))/sigma)^2)/(-2));
    end
end

alphaOri = quadprog(H,f,A,b,Aeq,beq,lb,ub);
alpha = alphaOri;
alpha(alpha<   sqrt(eps) ) = 0;
alpha(alpha>(C-sqrt(eps))) = C;
alpha = round(alpha,6);

AlphaSum(1,5) = sum(alpha);

for i = 1:50
    if (alpha(i)<C) && (alpha(i)>0)
        ak = [i,alpha(i)];
        break
    end
end

Kb = 0;
for i = 1:50
    Kb = Kb+(alpha(i)*Tr(i,3)*exp(((norm(Tr(i,1:2)-Tr(ak(1),1:2))/sigma)^2)/(-2)));
end

bias = (1/Tr(ak(1),3))-Kb;
Bias(1,5) = round(bias,4);

 % testing 
Dt = zeros(50,1);
for i = 1:50
    for j = 1:50
        Dt(i) = Dt(i)+(alpha(j)*Tr(j,3)*exp(((norm(Tr(j,1:2)-Te(i,1:2))/sigma)^2)/(-2)));
    end
end
Dt = Dt+bias;
for i = 1:50
    if Dt(i)>=0
        Dt(i) = 1;
    else
        Dt(i) = -1;
    end
end

CR5 = 0;
for i = 1:50
    if Dt(i) == Te(i,3)
        CR5 = CR5+1;
    end
end
CR5 = CR5/50;
Alpha(:,5) = alpha;

%% Table
CR = [CR1 CR2 CR3 CR4 CR5];
t1 = [round(Alpha(1:50,:),4) ; AlphaSum ; Bias;CR ];
T1 = table(t1(:,1),t1(:,2),t1(:,3),t1(:,4),t1(:,5), ...
    'VariableNames',{'sigma=5','sigma=1','sigma=0.5','sigma=0.1','sigma=0.05'}, ...
    'RowNames',{'alpha_1','alpha_2','alpha_3','alpha_4','alpha_5','alpha_6','alpha_7','alpha_8','alpha_9','alpha_10','alpha_11','alpha_12','alpha_13','alpha_14','alpha_15','alpha_16','alpha_17','alpha_18','alpha_19','alpha_20','alpha_21','alpha_22','alpha_23','alpha_24','alpha_25','alpha_26','alpha_27','alpha_28','alpha_29','alpha_30','alpha_31','alpha_32','alpha_33','alpha_34','alpha_35','alpha_36','alpha_37','alpha_38','alpha_39','alpha_40','alpha_41','alpha_42','alpha_43','alpha_44','alpha_45','alpha_46','alpha_47','alpha_48','alpha_49','alpha_50','Total','bias','CR'});
format short;
disp(T1);

