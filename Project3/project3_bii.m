%   Name: Theodoros Tsourdinis 
%   AEM: 2303
%   Project 2a: m-PSK - 16-QAM Flat Fading with Equalization system simulation 
%   Only tested on MATLAB (version: 2016b)
%   Estimated execution-elapsed time:  N = 10^4 ----> 21 seconds (7 seconds after 1st execution)
%                                      N = 10^5 ----> 2-3 minutes 
%                                      N = 10^6 ----> 11 minutes
%   N: Number of Bits
tic;
T = 1;
Rb = 10^3;
N = 10^5;
M = [4];
SNRdb = [0 2 4 6 8 10 12 14 16 18 20];
os = 5;
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
        % Define mapping table applying Gray mapping
        mappingTable(1) = r1 * exp(1i* 0);
        mappingTable(2) = r1 * exp(1i* pi/4);
        mappingTable(3) = r1 * exp(1i* 3*pi/4);
        mappingTable(4) = r1 * exp(1i* pi/2);
        mappingTable(5) = r1 * exp(1i* 7*pi/4);
        mappingTable(6) = r1 * exp(1i* 3*pi/2);
        mappingTable(7) = r1 * exp(1i* pi);
        mappingTable(8) = r1 * exp(1i* 5*pi/4);
        mappingTable(9:16) = mappingTable(1:8) ./ r1 .* r2;
   
        % Map bits to symbols
        for k = 1:4:length(bits)

            symbolBits = bits(k:k+3);

            symbolIndex = 2^3 * symbolBits(1) + 2^2 * symbolBits(2) + 2^1 * symbolBits(3) + 2^0 * symbolBits(4);

        % Mapping
            xn((k - 1)/4 + 1) = mappingTable(symbolIndex + 1);
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
      
      % Flat Fading Channel

        hi = sqrt(1/sqrt(2))*randn(round(length(new_sni)/Rb)+1,1);
        hq = sqrt(1/sqrt(2))*randn(round(length(new_sni)/Rb)+1,1);
        h = hi + 1i*hq ; 
      
    %Looping through SNR values
    for i = 1:length(SNRdb)
        SNRlin(i) = 10.^(SNRdb(i)./10);
        noise_power = 1./SNRlin(i);
        real_noise = sqrt(noise_power) .* randn(length(new_sni),1);
        imag_noise = sqrt(noise_power) .* randn(length(new_snq),1);
        w = real_noise + 1i*imag_noise;
        
     
        

%       Multiplying sn with different values of h every T*Rb and Adding noise%

        h_index=1;
        for j=1:length(sn)
            y(j) = h(h_index) *sn(j) + w(j);    
            if mod(j,T*Rb)==0
                h_index = h_index + 1;
            end
        end
        
        
        
%       Equalization

        h_index=1;
        for j=1:length(y)
            y(j) = (conj(h(h_index))/abs(h(h_index)^2)) * y(j);
            if mod(j,T*Rb)==0 
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
        if M(m) ~= 4
        %Decision Bit Maker%
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
             % Decision and demapping
            output_bits = zeros(1, N / 4);
            for j = 1:length(zn)
              [minv, idx] = min(zn(j) - mappingTable);
              symbolIndex = idx - 1;
              bitString = dec2bin(symbolIndex, 4);
              output_bits((j-1)*4 + 1) = str2double(bitString(1));
              output_bits((j-1)*4 + 2) = str2double(bitString(2));
              output_bits((j-1)*4 + 3) = str2double(bitString(3));
              output_bits((j-1)*4 + 4) = str2double(bitString(4));
            end
        end
          
        %Calculating errors
        for j = 1:N
            if(output_bits(j) ~= bits(j))
                errors(m,i) = errors(m,i) + 1;
            end 
        end    
    end
        %Constellation Diagrams
          figure(m);
          scatter(max_zni , max_znq,'*r')
          if M(m) ~= 4
            title(num2str(2^M(m),'Constellation Diagram of %-d-PSK'));
          end
          if M(m) == 4
            title(num2str(2^M(m),'Constellation Diagram of %-d-QAM'));
          end
end

%Calculating BER
BER = errors ./ N;

%BER Diagrams
C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]}; 


figure(m+1);
for i=1:length(M)
    semilogy(SNRdb,BER(i,:),'b-','color', C{i}, 'linewidth' ,2.0);
    title(num2str('BER Diagrams over SNR values After Equalization (1 Antenna)'));
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