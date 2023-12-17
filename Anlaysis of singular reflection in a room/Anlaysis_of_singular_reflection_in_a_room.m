clear all; close all; clc;

%parametri
bbbb=144;
gggg=2020;
hz=2-mod(bbbb,10)/10;
hm=1.5-1.4*mod(bbbb,8)/7;
xm=3+mod(bbbb,5)/4;
Uz=1+mod(bbbb,3);
Um=3-mod(bbbb,3);
M=1+mod(bbbb,3);

c=340;
fs=48000;
%hz=0.1; hm=0.1
%M=3;
%Um=4;
%Uz=2;
switch M
    case 1
        Alfa=[0.01 0.01 0.01 0.02 0.02 0.02 0.05 0.05 0.05 0.05]; 
    case 2
        Alfa=[0.18 0.18 0.12 0.10 0.09 0.08 0.07 0.07 0.07 0.07]; 
    case 3
        Alfa=[0.46 0.46 0.93 1.00 1.00 1.00 1.00 1.00 1 1]; 
end

f=[0 125 250 500 1000 2000 4000 8000 16000 24000]/24000;

D=sqrt((hz-hm)^2+(xm)^2);
R=sqrt((hz+hm)^2+(xm)^2); %praougli trougao iz teoreme likova

x1=(hz*xm)/(hz+hm);
x2=xm-x1;

%Prikaz položaja predajnika i prijemnika, refleksione ravni, direktne i reflektovane putanje
figure, plot([0 xm],[hz hm]), hold on, plot([0 x1 xm],[hz 0 hm], '--'), hold on,
plot(0,hz,'*',LineWidth=5), hold on,plot(xm,hm,'x', LineWidth=5),hold on,plot([-2 xm+2],[0 0], 'k',LineWidth=2)
xlim([-2 xm+2]); ylim([-0.5 max(hm,hz)+1]); legend('direktna komponenta','reflektovana komponenta', 'predajnik', 'prijemnik', 'refleksiona ravan');


A=sqrt(1-Alfa(1:length(Alfa)));
h=fir2(20,f,A);

td=ceil(D/c *fs);
tr=ceil(R/c *fs);


teta_md=atan((hz-hm)/xm);
teta_zd=teta_md;

teta_zr=atan(hz/x1);
teta_mr=teta_zr;

%direktna komponenta
switch Uz 
    case 1
        Gzd=cos(teta_zd); 
        Gzr=cos(teta_zr);
    case 2
        Gzd=1/2*(1+cos(teta_zd));
        Gzr=1/2*(1+cos(teta_zr));
    case 3
        Gzd=1/4*(1+3*cos(teta_zd)); 
        Gzr=1/4*(1+3*cos(teta_zr));
    case 4
        Gzd=1; 
        Gzr=1;

end
switch Um 
    case 1
        Gmd=cos(teta_md); 
        Gmr=cos(teta_mr);
    case 2
        Gmd=1/2*(1+cos(teta_md));
        Gmr=1/2*(1+cos(teta_mr));
    case 3
        Gmd=1/4*(1+3*cos(teta_md)); 
        Gmr=1/4*(1+3*cos(teta_mr));
    case 4
        Gmd=1; 
        Gmr=1;
end

d=zeros(1,fs);
r=zeros(1,fs);
pd=(sqrt(413/4*pi)/D)*Gzd*Gmd;
pr=(sqrt(413/4*pi)/R)*Gzd;
d(td)=pd;
r(tr)=pr;
r_a=filter(h,1,r)*Gmd;
%reflektovana komponenta

y=r_a+d;
y=y./max(abs(y));
Y=fft(y);
%Prikaz impulsnog odziva u prijemnoj tački.
figure, stem((0:length(Y)-1),y); title('impulsni odziv');
%Normalizovane prikaze frekvencijskog odziva i odziva po oktavama na mestu prijema u dB
figure, semilogx((0:length(Y)-1),20*log10(Y./(max(abs(Y))))); title('frekvencijski odziv u dB'); subtitle('logaritamska skala');
N=20*log10(Y./(max(abs(Y))));
figure, plot((0:length(Y)-1),N);title('frekvencijski odziv u dB '); subtitle('komb filtar');
%sound(y,fs);

br_oktava=11;
for br=1:br_oktava
    fc(br)=15.625*2^(br-1);
end
fg=fc/sqrt(2); fg(br_oktava+1)=fc(br_oktava)*sqrt(2);
for br=1:br_oktava
    O(br)=sum(abs(Y(fg(br):fg(br+1))).^2);
end
O=20*log10(O./(max(abs(O))));
for br=1:br_oktava
    y_rect(br*2-1)=O(br);
    y_rect(br*2)=O(br);
end
x_rect(1)=fg(1);x_rect(br_oktava*2)=fg(br_oktava+1);
for br=2:br_oktava
    x_rect(br*2-2)=fg(br);
    x_rect(br*2-1)=fg(br);
end
figure, semilogx(x_rect,y_rect); title('frekvencijski odziv po oktavama');
set(gca, 'Xtick',[16 31.5 63 125 250 500 1000 2000 4000 8000 16000 24000]);
set(gca, 'XtickLabel',{'16' '31.5' '63' '125' '250' '500' '1000' '2000' '4000' '8000' '16000' '24000'});
%figure, semilogx(fc,O); title('frekvencijski  odziv po oktavama');