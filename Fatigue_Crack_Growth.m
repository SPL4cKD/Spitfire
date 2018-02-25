clear all;
clc;
%% Fatigue Crack Growth - Center Crack
choice = menu('Choose a case ','Infinite Plate','Finite Plate');

%Calculations for either infinite or finite plate fatigue crack growth
if choice==1
    
    prompt = {'Enter a Max Normal Stress', 'Enter a c Constant','Enter a m Constant',...
        'Enter a Fracture Toughness','Enter a Initial Crack Length'};
    answer = inputdlg(prompt);
    
    syms a_c a;
    %Infinite Plate
    NS = str2num(answer{1}); %normal stress
    c = str2num(answer{2}); %material constant from da/dN vs. K plot
    m = str2num(answer{3}); %material constant from da/dN vs. K plot
    K_IC = str2num(answer{4}); %critical stress intensity factor
    a_0 = str2num(answer{5}); %initial crack length
    
    f_aw = 1; %geometric correction factor
    K = f_aw*NS*sqrt(pi*a_c); %stress intensity factor
    g = solve(K == K_IC,a_c); %solving for critical crack length
    Critical_Crack = double(g) %converting to double precision value
    I = 1/(c*(f_aw*NS*sqrt(pi*a))^m); %Paris's law
    N = int(I,a,a_0,Critical_Crack); %integral = Paris' law for number of cycles
    Cycles = double(N) %converting to double precision value
    
elseif choice==2
    
    prompt = {'Enter Max Normal Stress', 'Enter c Constant','Enter m Constant',...
        'Enter Fracture Toughness','Enter Initial Crack Length','Enter Plate Width'};
    title = 'Givens';
    answer = inputdlg(prompt);
    
    syms a_c a;
    %Finite Plate
    NS = str2num(answer{1}); %normal stress
    c = str2num(answer{2}); %material constant from da/dN vs. K plot
    m = str2num(answer{3}); %material constant from da/dN vs. K plot
    K_IC = str2num(answer{4}); %critical stress intensity factor
    a_0 = str2num(answer{5}); %initial crack length
    
    f_aw = (1-(a_c/(2*b))+(0.326*(a_c/b)^2))/sqrt(1-(a_c/b)); %geometric correction factor
    K = f_aw*NS*sqrt(pi*a_c);  %stress intensity factor
    g = solve(K == K_IC,a_c); %solving for critical crack length
    Critical_Crack = double(g) %converting to double precision value
    I = 1/(c*(f_aw*NS*sqrt(pi*a))^m); %Paris's law
    N = int(I,a,a_0,g); %integral = Paris' law for number of cycles
    Cycles = double(N) %converting to double precision value
end