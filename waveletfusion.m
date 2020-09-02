close all;clear;clc
%% Read Images
%% the size of images must be equal
a=imread("CT_Image2.jpg");
a=rgb2gray(a);
a=imresize(a,[276,276]);
b=imread("MRI_Image.jpg");
b=rgb2gray(b);
b=imresize(b,[276,276]);
%%   Wavelet Transform 
[a1,b1,c1,d1]=dwt2(a,'db2');
%a1 is approximate coefficient matrix
%b1,c1,d1 are detail horizontal,vertical and diagonal coefficient matrix respectively
[a2,b2,c2,d2]=dwt2(b,'db2');
%a2 is approximate coefficient matrix
%b2,c2,d2 are detail horizontal,vertical and diagonal coefficient matrix respectively
[k1,k2]=size(a1);
%% Fusion Rules
a3 = zeros([k1 k2]);%preallocating with zeros
%% Average Rule
for i=1:k1
    for j=1:k2
        a3(i,j)=(a1(i,j)+a2(i,j))/2;
   end
end
%% Max Rule
b3 = zeros([k1 k2]);%preallocating with zeros
c3 = zeros([k1 k2]);%preallocating with zeros
d3 = zeros([k1 k2]);%preallocating with zeros
for i=1:k1
    for j=1:k2
        b3(i,j)=max(b1(i,j),b2(i,j));
        c3(i,j)=max(c1(i,j),c2(i,j));
        d3(i,j)=max(d1(i,j),d2(i,j));
    end
end
%% Inverse Wavelet Transform 
c=idwt2(a3,b3,c3,d3,'db2');
figure,imshow(a)
title('CT Scan of Brain')
figure,imshow(b)
title('MRI Scan of Brain')
figure,imshow(c,[])
title('Fused Image')
%% Showing the tumour part
imageRGB=imread('final.jpg');
brain_bw=im2bw(rgb2gray(imageRGB));
figure,imshow(brain_bw)
title('Fused Image in gray scale')
result=bwareafilt(brain_bw,1);
figure,imshow(result) 
title('Tumor Segmented')
%% Performance Criteria
CR1=corr2(a,c);
CR2=corr2(b,c);
S1=snrr(double(a),double(c));
S2=snrr(double(b),double(c));
fprintf('Correlation between first image and fused image =%f \n\n',CR1);
fprintf('Correlation between second image and fused image =%f \n\n',CR2);
fprintf('Signal to Noise Ratio between first image and fused image =%4.2f db\n\n',S1);
fprintf('Signal to Noise Ratio between second image and fused image =%4.2f db \n\n',S2);
function r = snrr(in, est)
    error = in - est;
    r = 10 * log10((255^2)/ mean(error(:).^2));
end
