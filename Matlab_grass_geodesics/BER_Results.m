%% CONSTELLATION PARAMETERS
L = 16; % constellation size
M = 2; % transmit antennas
N = 2; % receive antennas
T = 4; % time slots > max(M,N);

C = Cg16; %Adjust this line after running one of the Code_T_M scripts
bit_labels = Bg16;%This is a matrix of L rows and log2(L) columns being K the number of codewords
%already defined by running an example (e.g. Codes_T4_M2)

%% SIMULATION PARAMETERS
SNR = 0:2:20;
ns = length(SNR);
% NumSim = 1e5.*ones(1,ns);
NumSim = [1e3 1e3 1e3 1e4 1e4 1e4 1e4 1e4 1e5 1e5 1e5].*10; %This can be adjusted manually

NoiseVar = (M/T)*10.^(-SNR/10); % noise variance


%% SIMULATION    
BER = zeros(1,length(SNR)); % symbol error rate
SER = zeros(1,length(SNR)); % symbol error rate

for cc = 1:ns % SNR loop
    disp(['SNR =  ' int2str(SNR(cc))])
    NoiseVarIter = NoiseVar(cc);
    BERc = 0; % BER counter
    SERc = 0; % symbol-error-rate counter

    for ss = 1:NumSim(cc) % Here parfor can be used for paralell computing
        
        if ~mod(ss,floor(NumSim(cc)/10)), fprintf('.'); end
        
        % TX Signal
        TX_codeword = randi(L);
        X = C(:,:,TX_codeword);
        
        % Rayleigh MIMO channel (fixed during the coherence block)
        H = randn(M,N)/sqrt(2) + 1i*randn(M,N)/sqrt(2); 
        
        % AWGN 
        Noise = sqrt(NoiseVarIter/2)*randn(T,N) + 1i*sqrt(NoiseVarIter/2)*randn(T,N);
        
        % RX Signal
        Y = X*H + Noise;
        
        % ML Detector
        RX_codeword = MLGrass(C,Y);

        % Error counter
        SERc = SERc + (RX_codeword ~= TX_codeword);

        %Bits tx and rx
        tx_bits = bit_labels(TX_codeword, :);
        rx_bits = bit_labels(RX_codeword, :);

        % Bit error counter
        BERc = BERc + sum(tx_bits ~= rx_bits);

    end
    
    SER(cc) = SERc / NumSim(cc); % SER computation
    BER(cc) = BERc / (NumSim(cc)*log2(L)); % BER computation
    

    fprintf('\n')

end

%% PLOT RESULTS
fs = 11;
lw = 1.5;
ms = 8;

% BER curve
figure();clf;
% hold on
semilogy(SNR,SER,'r--v','MarkerSize',ms,'LineWidth',lw);hold on;
semilogy(SNR,BER,'r--o','MarkerSize',ms,'LineWidth',lw);hold on;
xlabel('SNR (dB)');
ylabel('SER/BER');
legend('SER','BER');
title(['T = ' num2str(T) ', ' 'M = ' num2str(M) ', ' 'N = ' num2str(N) ', ' 'L = ' int2str(L)])
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
grid on

