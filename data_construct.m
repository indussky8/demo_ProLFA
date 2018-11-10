clc, clear all
load image_0001_sift.mat
Data = features.data(1:20,:);
num(1) = size(features.data,1);
load image_0002_sift.mat
Data = [Data;features.data(1:20,:)];
num(2) = size(features.data,1);
load image_0003_sift.mat
Data = [Data;features.data(1:20,:)];
num(3) = size(features.data,1);
load image_0004_sift.mat
Data = [Data;features.data(1:20,:)];
num(4) = size(features.data,1);

num = [20,20,20,20];

n = length(num);
X = Data;
[N,d] = size(X);
G = zeros(n,N);
for i =1:1:n
    if i == 1
        G(i,1:num(i)) = 1;
    else
        G(i, sum(num(1:i-1))+1 : sum(num(1:i))) = 1;
    end
end

Y = [1,0;1,0;0,1;0,1];


load image_0001_sift.mat
Data = features.data;
num(1) = size(features.data,1);
load image_0002_sift.mat
Data = [Data;features.data];
num(2) = size(features.data,1);
load image_0003_sift.mat
Data = [Data;features.data];
num(3) = size(features.data,1);
load image_0004_sift.mat
Data = [Data;features.data];
num(4) = size(features.data,1);

n = length(num);
X = Data;
[N,d] = size(X);
G = zeros(n,N);
for i =1:1:n
    if i == 1
        G(i,1:num(i)) = 1;
    else
        G(i, sum(num(1:i-1))+1 : sum(num(1:i))) = 1;
    end
end

Y = [1,0;1,0;0,1;0,1];
