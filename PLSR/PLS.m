clc;
clear;

%装载数据
load('FitnessData.txt');

%计算均值和标准差
mu = mean(FitnessData);
sig = std(FitnessData);

%求相关系数和数据标准化
rr = corrcoef(FitnessData);
data = zscore(FitnessData);

%因变量个数和自变量个数
n = 3;
m = 3;

X0 = FitnessData(:,1:n);
Y0 = FitnessData(:,n+1:end);

E0 = data(:,1:n);
F0 = data(:,n+1:end);

num = size(E0, 1);

chg = eye(n);

for i = 1 : n
    matrix = (E0')*(F0)*(F0')*E0;
    
    [vec, val] = eig(matrix);
    val = diag(val);
    [val, ind] = sort(val, 'descend');
    w(:,i) = vec(:,ind(1));
    w_star(:,i) = chg * w(:, i);
    t(:,i) = E0*w(:,i);
   alpha = E0'*t(:,i)/(t(:,i)'*t(:,i));
   chg = chg*(eye(n) - w(:,i)*alpha');
   E = E0 - t(:,i)*alpha';
   E0 = E;
   
   %计算ss(i)的值
   beta = [t(:,1:i),ones(num,1)]\F0;
   beta(end,:) = [];
   cancha = F0 - t(:,1:i)*beta;
   ss(i) = sum(sum(cancha.^2));
   
   %计算press(i)
   for j = 1 : num
       t1 = t(:,1:i);
       F1 = F0;
       she_t = t1(j,:);
       she_f = F1(j,:);
       t1(j,:) = [];
       F1(j,:) = [];
       beta1 = [t1, ones(num-1, 1)]\F1;
       beta1(end,:) = [];
       cancha = she_f - she_t*beta1;
       press_i(j) = sum(cancha.^2);
   end
   press(i) = sum(press_i);
   
   if i > 1
       Q_h2(i) = 1 - press(i)/ss(i-1);
   else
       Q_h2(1) = 1;
   end
   
   if Q_h2(i) < 0.0975
       fprintf('提出的成分个数r=%d',i);
       r = i;
       break;
   end
end

beta_z = [t(:,1:r),ones(num,1)]\F0;
beta_z(end,:) = [];
xishu = w_star(:,1:r)*beta_z;
mu_x = mu(1:n);
mu_y = mu(n+1:end);
sig_x = sig(1:n);
sig_y = sig(n+1:end);

for i = 1 : m
    ch0(i) = mu_y(i) - mu_x ./ sig_x*sig_y(i)*xishu(:,i);
end

for i = 1 : m
    xish(:,i) = xishu(:,i)./sig_x'*sig_y(i);
end


sol = [ch0;xish];

ch0 = repmat(ch0, num, 1);
yhat = ch0 + X0*xish;
y1max = max(yhat);
y2max = max(Y0);
ymax =  max([y1max;y2max]);
cancha = yhat - Y0;
subplot(2,2,1);
plot(0:ymax(1), 0:ymax(1), yhat(:,1), Y0(:,1),'*');
subplot(2,2,2);
plot(0:ymax(2), 0:ymax(2), yhat(:,2), Y0(:,2), 'O');
subplot(2,2,3);
plot(0:ymax(3), 0:ymax(3), yhat(:,3), Y0(:,3), 'H');