a = input('Tuliskan Nama Anda =  ','s');
b = input('Tuliskan Nim Anda = ');
d = fix(b/1000000000);
d1 = mod(b,100);
n = d*100+d1;
for x = 1:n
c1 = fix(rand(5,5)*100);
c2 = fix(rand(5,5)*100);
c3 = fix(rand(5,5)*100);
end
R =c1
G = c2
B = c3
save(a,'b','R','G','B')

%Nomor 2
[kolom, baris] = size(R);
hasil = zeros(kolom,baris);

for x = 1 : kolom
    for y = 1 : baris
        hasil(x,y) = R(x,y) * 0.4 + G(x,y) * 0.32 + B(x,y) * 0.28;
    end
end

%% Histogram array
myhist = zeros(256,1);
for k = 0:255
    myhist(k+1) = numel(find(I==k)); % number of elements where I has gray level equal to 'k'     
end
figure, stem(myhist,'r');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Number of Elements')
title('Histogram of Image')
grid on
% End of Histogram array
% Calculate cdf
cdf = zeros(256,1);
cdf(1) = myhist(1);
for k = 2:256
    cdf(k) = cdf(k-1)+myhist(k);
end
figure, stem(cdf,'r');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Cumulative Distribution')
title('CDF of Image')
grid on
% End of Calculate cdf

%% Find Equalized histogram array
J = I;
cumprob = cdf/(rows.*cols);
equalizedhist = floor((cumprob).*255);
 
for i = 1:rows
    for j = 1:cols
        for m = 0:255
            if (I(i,j)==m)
                J(i,j) = equalizedhist(m+1);
            end
        end
    end
end
 
figure, imshow(J);
title('Equalized Image')

%% Equalized Histogram array
myeqhist = zeros(256,1);
for k = 0:255
    myeqhist(k+1) = numel(find(J==k));
end
figure, stem(myeqhist,'r');
set(gca,'XLim',[0 255])
xlabel('Gray Level')
ylabel('Number of Elements')
title('Histogram of Equalized Image')
grid on
% End of Equalized Histogram array

%NOMOR 5
h=[1 1 1;1 4 1;1 1 1];
[kolom, baris]=size(hasil);
z=zeros(kolom,baris);
[kolom_h , baris_h]=size(h);

for x = 1 : kolom
    for y=1:baris
        for k1=1:kolom_h
            for k2=1:baris_h
                aa=x-2+k1;
                bb=y-2+k2;

                if aa == 00 || bb==0 || aa==kolom+1 || bb==kolom+1
                   z(x,y)=z(x,y)+(h(k1,k2)*0);
                else
                   z(x,y)=z(x,y)+h(k1,k2)*hasil(aa,bb);
                end
            end
        end
    end
end