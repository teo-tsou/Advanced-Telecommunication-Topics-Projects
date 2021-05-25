%Theodoros Tsourdinis 2303   1978.352485 sec 1364.857234  %
tic;
M = [1 2 3 4];
SNRdb = [0 2 4 6 8 10 12 14 16 18 20];
os = 5;
a = 1;
gtx = zeros(1,os);
grx = gtx;
N = 1000;

while 1
    isOk = 1;
    for m = 1:length(M)
       if mod(N,M(m)) ~= 0
           isOk = 0;
           break
       end
    end
       
   if isOk == 1
       break
   else
       N = N + 1;
   end
end

SNRlin = zeros(length(SNRdb),1);
noise_power = zeros(length(SNRdb),1);
output_bits = zeros(N,1);
errors = zeros(length(M),length(SNRdb));
%BER = zeros(length(M),length(SNRdb));
bits = rand(1,N) < 0.5;

%Main - Function%
for m = 1:length(M)
    xn = zeros(round(N/M(m)),1);
    sni = zeros(round(N/M(m)), os);
    snq = zeros(round(N/M(m)), os);
    zni = zeros(round(N/M(m)) ,(2*os)-1);
    znq = zeros(round(N/M(m)),(2*os)-1);
    max_zni = zeros(round(N/M(m)),1);
    max_znq = zeros(round(N/M(m)),1);
    
    j = 1;
    for k = 1:M(m):N
        xni = cosd(bi2de(bits(k:k-1+M(m)), 'left-msb')*180/(2^(M(m)-1)));
        xnq = sind(bi2de(bits(k:k-1+M(m)), 'left-msb')*180/(2^(M(m)-1)))*1i;
        xni = sqrt(log2(2^M(m))) * xni;
        xnq = sqrt(log2(2^M(m))) * xnq;
        xn(j) = xni + xnq;
        j = j + 1;
    end
    
    
%     Tx Filtering

      for j = 1:os
        gtx(j) = sqrt(1/os) ;
      end 

      for i=1:length(xn)
          sni(i,:) = conv(real(xn(i)),gtx);
          snq(i,:) = conv(imag(xn(i)),gtx);
      end
      
      sni=sni';
      sni = sni(:);
      
      snq = snq';
      snq = snq(:);
      
    
    for i = 1:length(SNRdb)
        SNRlin(i) = 10.^(SNRdb(i)./10);
        noise_power = 1./SNRlin(i);
        real_noise = sqrt(noise_power) .* randn(1,length(sni));
        imag_noise = sqrt(noise_power) .* randn(1,length(snq));
        yni = sni + real_noise(:);
        ynq = snq + imag_noise(:);
        
        yni = yni(:);
        ynq = ynq(:);
        
        %Matched Filtering
        for j = 1:os
        grx(j) = sqrt(1/os) ;
        end 
        
          j=1;
          for k=1:os:length(yni)
              zni(j,:) = conv(yni(k:k+os-1),grx);
              znq(j,:) = conv(ynq(k:k+os-1),grx);
              j = j + 1;
          end
          
          new_zni = zni';
          new_zni = new_zni(:);
         
          new_znq = znq';
          new_znq = new_znq(:);
          
          j=1;
          for k=1:(2*os)-1:length(new_zni)
            max_zni_temp = new_zni(k:k+(2*os)-2);
            max_zni(j) = max_zni_temp(round(((2*os)-1)/2));
            max_znq_temp = new_znq(k:k+(2*os)-2);
            max_znq(j) = max_znq_temp(round(((2*os)-1)/2));
            j = j + 1;
          end
          
          zn = max_zni + 1i*max_znq;
         
        %Decision Bit Maker%
          p = 1;
          dist = zeros(length(xn),1);
          output_symbol = zeros(length(zn),1);
          for j = 1:length(zn)
            for k = 1:length(xn)   
                dist(k) = norm(xn(k)-zn(j));
            end
            [minv,idx] = min(dist);
            output_symbol(j) = xn(idx);
            phi= radtodeg(angle(output_symbol(j)));
            
            if (phi < 0)
                decoded = dec2bin((360+phi)*(2^(M(m)-1))/180)-'0';
            else
                decoded = dec2bin(phi*(2^(M(m)-1))/180)-'0';
            end
            output_bits(p:p+M(m)-1) = [zeros(1, M(m)-length(decoded)) decoded];
            dist = zeros(length(xn),1);
            p = p+M(m);
          end    
                
        for j = 1:N
            if(output_bits(j) ~= bits(j))
                errors(m,i) = errors(m,i) + 1;
            end 
        end    
    end
          figure(m);
          scatter(max_zni , max_znq,'*r')
          title(num2str(2^M(m),'Constellation Diagram of %-d-PSK'));
end

BER = errors ./ N;

figure(m+1);
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; 
semilogy(SNRdb,qfunc(sqrt(SNRlin)),'y-','linewidth',2.0); 
hold on;
semilogy(SNRdb,BER(1,:),'b-','color', C{1}, 'linewidth' ,2.0);
legend('Theoretical BPSK', 'Simulation BPSK');
hold off;

figure(6);
for i=1:length(M)
    semilogy(SNRdb,BER(i,:),'b-','color', C{i}, 'linewidth' ,2.0);
    hold on;
end

 for l = 1:length(M)
      legendCell{l} = num2str(2^M(l),'%-d-PSK');
 end
legend(legendCell);
hold off;
toc;