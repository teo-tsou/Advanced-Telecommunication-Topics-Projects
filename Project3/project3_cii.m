%   Name: Theodoros Tsourdinis 
%   AEM: 2303
%   Project 3c: MRC m-PSK - 16-QAM BER simulation 
%   Only tested on MATLAB (version: 2016b)
%   Estimated execution-elapsed time:  N = 10^4 ----> 21 seconds (7 seconds after 1st execution)
%                                      N = 10^5 ----> 1-2 minutes 
%                                      N = 10^6 ----> 11 minutes
%   N: Number of Bits
tic;
T = 1000;
N = 10^5;
M = [1 4];
SNRdb = [0 2 4 6 8 10 12 14 16 18 20 22 24 26 28];
os = 1;
a = 1;
gtx = zeros(1,os);
grx = gtx;

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
bits = rand(1,N) < 0.5;
correct_packets = zeros(length(M),length(SNRdb));
total_packets = 0 ;

for j = 1:length(bits)
    if mod(j,T) == 0
            total_packets = total_packets + 1;
    end
end

%Looping through modulations%
for m = 1:length(M)
    xn = zeros(round(N/M(m)),1);
    sni = zeros(round(N/M(m)), os);
    snq = zeros(round(N/M(m)), os);
    zni = zeros(round(N/M(m)) ,(2*os)-1);
    znq = zeros(round(N/M(m)),(2*os)-1);
    max_zni = zeros(round(N/M(m)),1);
    max_znq = zeros(round(N/M(m)),1);
    m_symbols = zeros(M(m),1);
    y1 = zeros(round(N/M(m))*os,1);
    y2 = zeros(round(N/M(m))*os,1);
    y = zeros(round(N/M(m))*os,1);
    
    if M(m) ~= 4
    %Symbol m-psk encoder
    j = 1;
    for k = 1:M(m):N
        xni = cosd(bi2de(bits(k:k-1+M(m)), 'left-msb')*180/(2^(M(m)-1)));
        xnq = sind(bi2de(bits(k:k-1+M(m)), 'left-msb')*180/(2^(M(m)-1)))*1i;
        xni = sqrt(log2(2^M(m))) * xni;
        xnq = sqrt(log2(2^M(m))) * xnq;
        xn(j) = xni + xnq;
        j = j + 1;
    end
    m_symbols = unique(xn);
    end
    
    if M(m) == 4
       % 16 - QAM Implementation
        
        % Radius of inner and outer circle
        r1 = 1;
        r2 = 2;
        % Define lookup - mapping table with Gray mapping
        lookup(1) = r1 * exp(1i* 0);
        lookup(2) = r1 * exp(1i* pi/4);
        lookup(3) = r1 * exp(1i* 3*pi/4);
        lookup(4) = r1 * exp(1i* pi/2);
        lookup(5) = r1 * exp(1i* 7*pi/4);
        lookup(6) = r1 * exp(1i* 3*pi/2);
        lookup(7) = r1 * exp(1i* pi);
        lookup(8) = r1 * exp(1i* 5*pi/4);
        lookup(9:16) = lookup(1:8) ./ r1 .* r2;
   
        % Bits to Symbols
        for k = 1:4:length(bits)

            bit_map = bits(k:k+3);

            bit_index = 2^3 * bit_map(1) + 2^2 * bit_map(2) + 2^1 * bit_map(3) + 2^0 * bit_map(4);

        % Mapping
            xn((k - 1)/4 + 1) = lookup(bit_index + 1);
        end
    end
    
%     Tx Filtering

      for j = 1:os
        gtx(j) = sqrt(1/os) ;
      end 

      for i=1:length(xn)
          sni(i,:) = conv(real(xn(i)),gtx);
          snq(i,:) = conv(imag(xn(i)),gtx);
      end
      
      new_sni=sni';
      new_sni = new_sni(:);
      
      new_snq = snq';
      new_snq = new_snq(:);
      
      sn = new_sni + 1i*new_snq;
      
      %Flat Fading Channel - Two Different Channels%
    
        h1i = sqrt(1/sqrt(2))*randn(round(length(new_sni)/T)+1,1);
        h1q = sqrt(1/sqrt(2))*randn(round(length(new_sni)/T)+1,1);
        h1 = h1i + 1i*h1q ;

        h2i = sqrt(1/sqrt(2))*randn(round(length(new_sni)/T)+1,1);
        h2q = sqrt(1/sqrt(2))*randn(round(length(new_sni)/T)+1,1);
        h2 = h2i + 1i*h2q ;
      
    %Looping through SNR values
    for i = 1:length(SNRdb)
        SNRlin(i) = 10.^(SNRdb(i)./10);
        noise_power = 1./SNRlin(i);
        real_noise = sqrt(noise_power) .* randn(length(new_sni),1);
        imag_noise = sqrt(noise_power) .* randn(length(new_snq),1);
        w = real_noise + 1i*imag_noise;
        
      
% %         Flat Fading Channel
% 
%         hi = sqrt(1/sqrt(2))*randn(1,length(new_sni)/Rb);
%         hq = sqrt(1/sqrt(2))*randn(1,length(new_snq)/Rb);
%         h = hi + 1i*hq ; 
%         
% 
% %       Multiplying sn with different values of h every T*Rb and Adding noise%
% 
%         h_index=1;
%         for j=1:length(sn)
%             y(j) = h(h_index).*sn(j) + w(j);    
%             if mod(j,T*Rb)==0
%                 h_index = h_index + 1;
%             end
%         end
%         
%         
% %       Equalization
% 
%         h_index=1;
%         for j=1:length(y)
%             y(j) = (conj(h(h_index)) * y(j))/(abs(h(h_index)^2));
%             if mod(j,T*Rb)==0 
%                 h_index = h_index + 1;
%             end
%         end
        
       
        
%       Multiplying sn with different values of h every T*Rb and Adding noise%

        h_index=1;
        for j=1:length(sn)
            y1(j) = h1(h_index) *sn(j) + w(j);
            y2(j) = h2(h_index) *sn(j) + w(j);
            if mod(j,T)==0
                h_index = h_index + 1;
            end
        end
        
%       MRC       

        h_index=1;
        for j=1:length(y1)
            y_vector = [y1(j);y2(j)];
            h_vector = [h1(h_index);h2(h_index)];
            y(j) = (h_vector'/(norm(h_vector)^2)) * y_vector ;
            if mod(j,T)==0 
                h_index = h_index + 1;
            end
        end

        
        yni = real(y);
        ynq = imag(y);
        
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
         if M(m) ~= 4
              p = 1;
              dist = zeros(length(m_symbols),1);
              output_symbol = zeros(length(zn),1);
              for j = 1:length(zn)
                for k = 1:length(m_symbols)   
                    dist(k) = norm(m_symbols(k)-zn(j));
                end
                [minv,idx] = min(dist);
                output_symbol(j) = m_symbols(idx);
                phi= radtodeg(angle(output_symbol(j)));

                if (phi < 0)
                    decoded = dec2bin((360+phi)*(2^(M(m)-1))/180)-'0';
                else
                    decoded = dec2bin(phi*(2^(M(m)-1))/180)-'0';
                end
                output_bits(p:p+M(m)-1) = [zeros(1, M(m)-length(decoded)) decoded];
                p = p+M(m);
              end
         end
         
          
         if M(m) == 4
                % 16-QAM Decision
            output_bits = zeros(1, N / 4);
            for j = 1:length(zn)
              [minv, idx] = min(zn(j) - lookup);
              bit_index = idx - 1;
              bit_str = dec2bin(bit_index, 4);
              output_bits((j-1)*4 + 1) = str2double(bit_str(1));
              output_bits((j-1)*4 + 2) = str2double(bit_str(2));
              output_bits((j-1)*4 + 3) = str2double(bit_str(3));
              output_bits((j-1)*4 + 4) = str2double(bit_str(4));
            end
          end
          
                %Calculating errors - Packets
        correct_packet = 1;
        for j = 1:N
            if(output_bits(j) ~= bits(j))
                errors(m,i) = errors(m,i) + 1;
                correct_packet = 0;
            end    
            if mod(j,T) == 0
                if correct_packet == 1
                    correct_packets(m,i) = correct_packets(m,i) + 1;
                end
                correct=1;
            end
            
        end 
        
    end
        %Constellation Diagrams
          figure(m);
          scatter(max_zni , max_znq,'*r')
          title(num2str(2^M(m),'Constellation Diagram of %-d-PSK'));
end

%Calculating BER
BER = errors ./ N;

%BER Diagrams
C = {'r','g','y',[.5 .6 .7],[.8 .2 .6]}; 


figure(m+1);
for i=1:length(M)
    semilogy(SNRdb,BER(i,:),'b-','color', C{i}, 'linewidth' ,2.0);
    title('BER Diagrams over SNR values after the MRC');
    hold on;
end

for l = 1:length(M)
    if M(l) ~= 4
        legendCell{l} = num2str(2^M(l),'%-d-PSK');
    end
    if M(l) == 4
        legendCell{l} = num2str(2^M(l),'%-d-QAM');
    end
 end
legend(legendCell);
hold off;

goodput = correct_packets ./ total_packets;
 figure(m+2);
for i=1:length(M)
    plot(SNRdb,goodput(i,:),C{i});
    title('Goodput - MRC');
    hold on;
end

 for l = 1:length(M)
    if M(l) ~= 4
        legendCell{l} = num2str(2^M(l),'%-d-PSK');
    end
    if M(l) == 4
        legendCell{l} = num2str(2^M(l),'%-d-QAM');
    end
 end
 
legend(legendCell);
hold off;

toc;