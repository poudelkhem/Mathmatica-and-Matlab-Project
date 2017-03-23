%%****************************************************
% Middle Tennesse State University
%  Date: 02/25/2017, @Khem Narayan Poudel M01385540
%  Project 3:Planetary Motion part Part1
% Input:
% N=1000000
% Methods 1. Euler 2. Cromer and 3. RK 
% Interval=5000
%Output:  Diplay all Plannet Data in tabular fom and display graphical
%views to
%Output looks like:
%######################### Output file ######################################
% Units are:Distance:10^6 km, Masses:kg, Speed:km/s and period:days
% 
% UsefulPlanetaryData = 
% 
%                   Mass       PerihelionDistance    PerihelionSpeed    AphelionDistance    AphelionSpeed    SemimajorAxis    OrbitalPeriod    OrbitEccentricity    MeanorbitVelocity
%                __________    __________________    _______________    ________________    _____________    _____________    _____________    _________________    _________________
% 
%     Mercury     3.302e+23    4.6028e+07            58.942             6.9817e+07          38.859           5.7922e+07        87.97             0.20536            47.354           
%     Venus      4.8685e+24    1.0749e+08            35.254               1.09e+08          34.765           1.0825e+08        224.7           0.0071115            35.008           
%     Earth      5.9722e+24     1.471e+08            30.292             1.5234e+08          29.252           1.4972e+08       365.26            0.017565            29.764           
%     Mars       6.4185e+23    2.0667e+08            26.496             2.4939e+08          21.957           2.2803e+08       686.98             0.09368            24.065           
%     Jupiter    1.8981e+27    7.4063e+08            13.715             8.1871e+08          12.407           7.7967e+08       4332.6            0.050036            13.034           
%     Saturn     5.6832e+26    1.3539e+09            10.173             1.5149e+09          9.0916           1.4344e+09        10759             0.05611            9.6061           
%     Uranus      8.681e+25     2.751e+09            7.0959             3.0047e+09          6.4966           2.8778e+09        30685             0.04411            6.7874           
%     Neptune    1.0241e+26    4.4534e+09            5.4884             4.5541e+09          5.3671           4.5038e+09        60189            0.011131            5.4272           
%###########################################################################################################################################################################################
%%*****************************************************

function [  ] = project3_part1mat( N, method, interval )
    
     % initialization of variables
    m = zeros(1,8);
    r = zeros(3,N,8);
    v = zeros(3,N,8);
    G = 1e-9*6.67e-11;
    M = 1.989e30;      
     % planets years in seconds     
   T = 24*60*60*[87.97, 224.701, 365.256, 686.980, 4332.59, 10759.22,30685.4, 60189];
           
    % Store all the  planet data from txt file
    [ m(1), r(:,1,1), v(:,1,1) ] = MercuryDatafile;
    [ m(2), r(:,1,2), v(:,1,2) ] = VenusDatafile;
    [ m(3), r(:,1,3), v(:,1,3) ] = EarthDatafile;
    [ m(4), r(:,1,4), v(:,1,4) ] = MarsDatafile;
    [ m(5), r(:,1,5), v(:,1,5) ] = JupiterDatafile;
    [ m(6), r(:,1,6), v(:,1,6) ] = SaturnDatafile;
    [ m(7), r(:,1,7), v(:,1,7) ] = UranusDatafile;
    [ m(8), r(:,1,8), v(:,1,8) ] = NeptuneDatafile;
    
    str = [ 'Percentage complete: ', num2str(0) ];
    disp(str)
      % Use different methods 1. Euler 2. Cromer and 3. RK     
      
      for i = 1:8
        
        if( method == 1 )
            [ r(:,:,i), v(:,:,i) ] =  Euler( G, M, r(:,:,i), v(:,:,i), N, T(i));
        elseif( method == 2 )
            [ r(:,:,i), v(:,:,i) ] = Cromer( G, M, r(:,:,i), v(:,:,i), N, T(i));
        else
            [ r(:,:,i), v(:,:,i) ] =    RK( G, M, r(:,:,i), v(:,:,i), N, T(i));
        end
        
        str = [ 'Percentage complete: ', num2str(100*i/8) ];
        disp(str)
        
      end
    
    figure
    hold on
    for i = 1:4
        plot(r(1,:,i),r(2,:,i))
    end
    hold off
    title('Terrestiral Planets ','color','r','fontsize', 16);

legend('Mercury', 'Venus','Earth', 'Mars');
figure
    hold on
    for i = 5:8
        plot(r(1,:,i),r(2,:,i))
    end
    hold off
  title('Jovian Planets ','color','r','fontsize', 16);

legend('Jupiter', 'Saturn','Uranus', 'Neptune');  
    figure
    hold on
    for i = 1:8
        plot(r(1,:,i),r(2,:,i))
    end
    hold off
  title('All Planets ','color','r','fontsize', 16);

legend('Mercury', 'Venus','Earth', 'Mars','Jupiter', 'Saturn','Uranus', 'Neptune');    
    planetdatatable( r, v, m, T, N )
  
 end

%Reading the mass radius andd velocity data from Mercury.txt 
function [ mass, r, v ] = MercuryDatafile( )    
    filename = 'Mercury.txt';
    file = fopen(filename, 'r');    
    % Skip  first 5 lines
    textscan(file,'%s',5,'Delimiter','\n');
    % Read mass data
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline),{' ', '=', '(', ')', '^'});
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(3)),'^');
    value3 = strsplit(char(headerline(5)),'+-');    
     % find  mass
    mass = str2num(char(value3(1)))*str2num(char(value1(1)))^str2num(char(value2(1)));
    
    %     disp(mass)    
    % skip remaining 47 lines and search for radius x,y and z
    textscan(file,'%s',47,'Delimiter','\n');
    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    x = str2num(char(value1));
    y = str2num(char(value2));
    z = str2num(char(value3));    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    vx = str2num(char(value1));
    vy = str2num(char(value2));
    vz = str2num(char(value3));    
    r = [x;y;z];
    v = [vx;vy;vz];    
end

%Reading the mass radius andd velocity data from Venus.txt 
function [ mass, r, v ] = VenusDatafile( )
    
    filename = 'Venus.txt';
    file = fopen(filename, 'r');    
    % Skip  first 5 lines
    textscan(file,'%s',5,'Delimiter','\n');
    % Read mass data
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline),{' ', '=', '(', ')', '^'});
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(3)),'^');
    value3 = strsplit(char(headerline(5)),'+-');
    
     % find  mass
    mass = str2num(char(value3(1)))*str2num(char(value1(1)))^str2num(char(value2(1)));
    
    %     disp(mass)    
    % skip remaining 47 lines and search for radius x,y and z
    textscan(file,'%s',47,'Delimiter','\n');
    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    x = str2num(char(value1));
    y = str2num(char(value2));
    z = str2num(char(value3));    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    vx = str2num(char(value1));
    vy = str2num(char(value2));
    vz = str2num(char(value3));    
    r = [x;y;z];
    v = [vx;vy;vz];    
    
end
%Reading the mass radius andd velocity data from Earth.txt 
function [ mass, r, v ] = EarthDatafile( )
    
    filename = 'Earth.txt';
    file = fopen(filename, 'r');    
    % Skip  first 4 lines
    textscan(file,'%s',4,'Delimiter','\n');
    
   %Read mass data
    firstline = fscanf(file, '%s %s %s %s %f%s%f', [1, inf]);   
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline),' ');
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(5)),'+-');    
    % Find mass
    mass = str2num(char(value2(1)))*str2num(char(value1(1)))^str2num(char(value1(2)));
    
  %     disp(mass)    
    % skip remaining 47 lines and search for radius x,y and z
    textscan(file,'%s',52,'Delimiter','\n');
    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    x = str2num(char(value1));
    y = str2num(char(value2));
    z = str2num(char(value3));    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    vx = str2num(char(value1));
    vy = str2num(char(value2));
    vz = str2num(char(value3));    
    r = [x;y;z];
    v = [vx;vy;vz];    
end
%Reading the mass radius andd velocity data from Mars.txt 

function [ mass, r, v ] = MarsDatafile( )
    
    filename = 'Mars.txt';
     file = fopen(filename, 'r');    
    % Skip  first 5 lines
    textscan(file,'%s',5,'Delimiter','\n');
    % Read mass data
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline),{' ', '=', '(', ')', '^'});
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(3)),'^');
    value3 = strsplit(char(headerline(5)),'+-');
    
     % find  mass
    mass = str2num(char(value3(1)))*str2num(char(value1(1)))^str2num(char(value2(1)));
    
    %     disp(mass)    
    % skip remaining 47 lines and search for radius x,y and z
    textscan(file,'%s',47,'Delimiter','\n');
    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    x = str2num(char(value1));
    y = str2num(char(value2));
    z = str2num(char(value3));    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    vx = str2num(char(value1));
    vy = str2num(char(value2));
    vz = str2num(char(value3));    
    r = [x;y;z];
    v = [vx;vy;vz];        
end
%Reading the mass radius andd velocity data from Jupiter.txt 
function [ mass, r, v ] = JupiterDatafile( )
    
    filename = 'Jupiter.txt';
     file = fopen(filename, 'r');    
    % Skip  first 4 lines
    textscan(file,'%s',4,'Delimiter','\n');
    % Read mass data
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline),{' ', '=', '(', ')', '^'});
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(3)),'^');
    value3 = strsplit(char(headerline(5)),'+-');
    
     % find  mass
    mass = str2num(char(value3(1)))*str2num(char(value1(1)))^str2num(char(value2(1)));
    
    %     disp(mass)    
    % skip remaining 51 lines and search for radius x,y and z
    textscan(file,'%s',51,'Delimiter','\n');
    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    x = str2num(char(value1));
    y = str2num(char(value2));
    z = str2num(char(value3));    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    vx = str2num(char(value1));
    vy = str2num(char(value2));
    vz = str2num(char(value3));    
    r = [x;y;z];
    v = [vx;vy;vz];    
end
%Reading the mass radius andd velocity data from Saturn.txt 

function [ mass, r, v ] = SaturnDatafile( )
    
    filename = 'Saturn.txt';
    file = fopen(filename, 'r');    
    % Skip  first 4 lines
    textscan(file,'%s',4,'Delimiter','\n');
    % Read mass data
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline),{' ', '=', '(', ')', '^'});
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(3)),'^');
    value3 = strsplit(char(headerline(5)),'+-');
    
     % find  mass
    mass = str2num(char(value3(1)))*str2num(char(value1(1)))^str2num(char(value2(1)));
    
    %     disp(mass)    
    % skip remaining 51 lines and search for radius x,y and z
    textscan(file,'%s',51,'Delimiter','\n');
    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    x = str2num(char(value1));
    y = str2num(char(value2));
    z = str2num(char(value3));    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    vx = str2num(char(value1));
    vy = str2num(char(value2));
    vz = str2num(char(value3));    
    r = [x;y;z];
    v = [vx;vy;vz];      
end
%Reading the mass radius andd velocity data from Uranus.txt 

function [ mass, r, v ] = UranusDatafile( )
    
    filename = 'Uranus.txt';
    file = fopen(filename, 'r');    
    % Skip  first 4 lines
    textscan(file,'%s',4,'Delimiter','\n');
    % Read mass data
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline),{' ', '=', '(', ')', '^'});
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(3)),'^');
    value3 = strsplit(char(headerline(5)),'+-');
    
     % find  mass
    mass = str2num(char(value3(1)))*str2num(char(value1(1)))^str2num(char(value2(1)));
    
    %     disp(mass)    
    % skip remaining 51 lines and search for radius x,y and z
    textscan(file,'%s',51,'Delimiter','\n');
    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    x = str2num(char(value1));
    y = str2num(char(value2));
    z = str2num(char(value3));    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    vx = str2num(char(value1));
    vy = str2num(char(value2));
    vz = str2num(char(value3));    
    r = [x;y;z];
    v = [vx;vy;vz];  
    
end
%Reading the mass radius andd velocity data from Neptune.txt 

function [ mass, r, v ] = NeptuneDatafile( )
    
    filename = 'Neptune.txt';
    file = fopen(filename, 'r');    
    % Skip  first 4 lines
    textscan(file,'%s',4,'Delimiter','\n');
    % Read mass data
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline),{' ', '=', '(', ')', '^'});
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(3)),'^');
    value3 = strsplit(char(headerline(5)),'+-');
    
     % find  mass
    mass = str2num(char(value3(1)))*str2num(char(value1(1)))^str2num(char(value2(1)));
    
    %     disp(mass)    
    % skip remaining 51 lines and search for radius x,y and z
    textscan(file,'%s',51,'Delimiter','\n');
    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    x = str2num(char(value1));
    y = str2num(char(value2));
    z = str2num(char(value3));    
    head = textscan(file,'%s',1,'Delimiter','\n');
    headerline = head{1};
    headerline = strsplit(char(headerline), {' ', '='});    
    value1 = strsplit(char(headerline(2)),'^');
    value2 = strsplit(char(headerline(4)),'+-');
    value3 = strsplit(char(headerline(6)),'+-');    
    vx = str2num(char(value1));
    vy = str2num(char(value2));
    vz = str2num(char(value3));    
    r = [x;y;z];
    v = [vx;vy;vz];  
    
end



%      1. Euler Method
function [ r, v ] = Euler( G, M, r, v, N, T)
    
    dt = 1.0*T/N;
    a = zeros(3,1);
    
    for i = 1:N-1
        r(:,i+1) = r(:,i) + dt*v(:,i);
        a(:) = -G*M/norm(r(:,i+1))^3*r(:,i+1);
        v(:,i+1) = v(:,i) + dt*a(:);
    end
    
end
% 2. Cromer Method
 function [ r, v ] = Cromer( G, M, r, v, N, T)
    
    dt = 1.0*T/N;
    a = zeros(3,1);
    
    for i = 1:N-1
        a(:) = -G*M/norm(r(:,i))^3*r(:,i);
        v(:,i+1) = v(:,i) + dt*a(:);
        r(:,i+1) = r(:,i) + dt*v(:,i+1);
    end
    
end
% 3. RK Method
function [ r, v ] = RK( G, M, r, v, N, T)
    
    dt = 1.0*T/N;
    kr2 = zeros(3,1);
    kv1 = zeros(3,1);
    kv2 = zeros(3,1);
    
    for i = 1:N-1
        kv1(:) = -G*M/norm(r(:,i))^3*r(:,i);
        
        kr2(:) = v(:,i) + 0.5*dt*kv1(:);
        kv2(:) = -G*M/norm(r(:,i) + 0.5*dt*v(:,i))^3*(r(:,i) + 0.5*dt*v(:,i));
        
        r(:,i+1) = r(:,i) + dt*kr2(:);
        v(:,i+1) = v(:,i) + dt*kv2(:);
    end
    
end

function [  ] = planetdatatable( r, v, m, T, N )    
    distance   = zeros(N,8);
    normvel = zeros(N,8);
    meanvel  = zeros(8,1);
    Apheldis    = zeros(8,1);
    Aphelvel    = zeros(8,1);
    peridist  = zeros(8,1);
    perivel  = zeros(8,1);
    majoraxis  = zeros(8,1);
    minaxis  = zeros(8,1);
    eccen    = zeros(8,1);
   center = zeros(3,1);    
    for i = 1:8
        for j = 1:N
            distance(j,i)   = norm(r(:,j,i));
            normvel(j,i) = norm(v(:,j,i));
        end
    end
    str = [ 'Percentage complete: ', num2str(0) ];
    disp(str)    
    for i = 1:8        
        meanvel(i) = mean(normvel(:,i));    
        [ Apheldis(i), indexAph ] = max(distance(:,i));
        Aphelvel(i) = norm(v(:,indexAph,i));
        [ peridist(i), indexperi ] = min(distance(:,i));
        perivel(i) = norm(v(:,indexperi,i));        
        APH_val = r(:,indexperi,i)-r(:,indexAph,i);
        majoraxis(i) = 0.5*norm(APH_val);
        center = 0.5*(r(:,indexperi,i)+r(:,indexAph,i));        
        n = 1;
        minimum_val(i) = n;
        minimum_init = 0;
        tempdata = 0;
        val_norm = 1;
        while( n <= N && minimum_val(i) >= minimum_init)            
            val_norm = norm(cross(APH_val,r(:,n,i)));
            minimum_init = minimum_val(i);
            if( val_norm > tempdata )
                minimum_val(i) = n;
                tempdata = val_norm;
            end
            n = n + 1;
        end
        minaxis(i) = norm(r(:,minimum_val(i),i)-center);
        eccen(i) = abs(sqrt(1-(minaxis(i)/majoraxis(i))^2));        
        str = [ 'Percentage complete: ', num2str(100*i/8) ];
        disp(str)
    end
    %Diplay all Plannet Data in tabular fom
    col = {'Mass', 'PerihelionDistance', 'PerihelionSpeed', 'AphelionDistance', 'AphelionSpeed', 'SemimajorAxis','OrbitalPeriod',  'OrbitEccentricity', 'MeanorbitVelocity'};
    row = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'};
    disp('Units are:Distance:10^6 km, Masses:kg, Speed:km/s and period:days') 
   
    
    UsefulPlanetaryData = table(m',peridist,perivel,Apheldis,Aphelvel,majoraxis,T'/60/60/24,eccen,meanvel,'RowNames',row,'VariableNames',col)
    
end






