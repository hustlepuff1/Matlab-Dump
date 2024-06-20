clc; clear all; close all;

figure(1);
subplot(2,1,1);
plot_fourier(3);
subplot(2,1,2);
plot_fourier(4);

figure(2);
subplot(2,1,1);
plot_fourier(3.5);
subplot(2,1,2);
plot_fourier(4.5);

figure(3);
subplot(2,1,1);
plot_fourier(10);
subplot(2,1,2);
plot_fourier(10.5);

figure(4)
period_rT = [3, 4, 3.5, 4.5, 10, 10.5];
for i=1:6
    subplot(6,1,i);
    plot_fourier(period_rT(i));
end

function plot_fourier(r)
    cn=[];

    for n=1:10*r
        temp1=0; temp2=0;

        for k = 1:2:2*r
            tmp_odd = exp(-1i*(k/r)*n*pi);
            temp1=temp1+tmp_odd;
        end

        for k = 2:2:2*r-1
            tmp_even = -exp(-1i*(k/r)*n*pi);
            temp2=temp2+tmp_even;
        end

        temp = -1/2 + temp1 + temp2 - 1/2*exp(-1i*2*n*pi);
        cn = [cn; 1i*temp/(pi*n)];

    end

    stem([0:1/r:n/r],[0; abs(cn)], 'o', 'filled')

    title(['Fourier coefficients of a square wave (r = ', num2str(r), ')'])
    xlabel('Frequency (Hz)')
    ylabel('Modulus(\mid\itc_n\rm\mid)')
end


