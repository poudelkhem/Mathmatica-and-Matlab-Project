
%%****************************************************
% Middle Tennesse State University
%  Date: 02/25/2017, @Khem Narayan Poudel M01385540
%  Project 3:Planetary Motion part Part2:Graphical Rotation
% Input:
% N=1000000
% methods 1. Euler 2. Cromer and 3. RK 
% Interval=5000
%Output:  Graphical Rotation of all plannet
%%*****************************************************
function [  ] = project3_part2( N, method, interval )
    
     % initialization of variables
    m = zeros(1,8);
    r = zeros(3,N,8);
    v = zeros(3,N,8);
    G = 1e-9*6.67e-11;
    M = 1.989e30;      
     % planet years in seconds
     % 60189 is Nuptune's orbital period in seconds
       T = 60189*24*60*60;
           
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
    if( method == 1 )
        [ r, v ] =  Euler( G, M, m, r, v, N, T);
    elseif( method == 2 )
        [ r, v ] = Cromer( G, M, m, r, v, N, T);
    else
        [ r, v ] =    RK( G, M, m, r, v, N, T);
    end
           % Call the plotting function
    graphicalmotionplanet( r, m, N, interval )    
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
% Graphical Display Function
function [  ] = graphicalmotionplanet( r, m, N, interval )
n = 10;
X = zeros(n,n,8);
Y = zeros(n,n,8);
Z = zeros(n,n,8);
%calculate the Radii
R = log10(m).^6/2;
colors = [  [0 1 0],
            [1 0.4 0.6],
            [0 0 1],
            [1 0 0],
            [0 0 0],
            [0.4 0.2 0],
            [0.3 0 0.3],
            [0 0 1]];
figure(1)
for i = 1:interval:N    
    for j = 1:8
        [X(:,:,j),Y(:,:,j),Z(:,:,j)] = sphere(n-1);
        X(:,:,j) = r(1,i,j) + R(j)*X(:,:,j);
        Y(:,:,j) = r(2,i,j) + R(j)*Y(:,:,j);
        Z(:,:,j) = r(3,i,j) + R(j)*Z(:,:,j);
    end    
    for j = 1:8
        surf(X(:,:,j),Y(:,:,j),Z(:,:,j),'FaceColor',colors(j,:),'EdgeColor','none');
        axis([ -3e9 3e9 -3e9 3e9 -3e9 3e9])
        lighting gouraud
        camlight left        
        hold on
    end
    hold off       
    pause(0.00000001)
    
end

end


%      1. Euler Method
function [ r, v ] = Euler( G, M, m, r, v, N, T)    
    dt = T/N;
    a = zeros(3,1);    
    for i = 1:N-1
        for j = 1:8
            r(:,i+1,j) = r(:,i,j) + dt*v(:,i,j);
            %Using equation 75 of class material
            a(:) = -G*M/norm(r(:,i+1,j))^3*r(:,i+1,j);            
            for k = 1:8
                if k ~= j
                    a(:) = a(:) - (r(:,i+1,j)-r(:,i+1,k))*G*m(k)/norm(r(:,i+1,j)-r(:,i+1,k))^3;
                end
            end            
            v(:,i+1,j) = v(:,i,j) + dt*a(:);
        end        
        if mod(i,N/100) == 0
            str = [ 'Percentage complete: ', num2str(100*i/N) ];
            disp(str)
        end
    end
    str = [ 'Percentage complete: ', num2str(100) ];
    disp(str)
    
end
% 2. Cromer Method
function [ r, v ] = Cromer( G, M, m, r, v, N, T)    
    dt = T/N;
    a = zeros(3,1);    
    for i = 1:N-1
        for j = 1:8
            a(:) = -G*M/norm(r(:,i,j))^3*r(:,i,j);            
            for k = 1:8
                if k ~= j
                    a(:) = a(:) - (r(:,i,j)-r(:,i,k))*G*m(k)/norm(r(:,i,j)-r(:,i,k))^3;
                end
            end
            
            v(:,i+1,j) = v(:,i,j) + dt*a(:);
            r(:,i+1,j) = r(:,i,j) + dt*v(:,i+1,j);
        end
        
        if mod(i,N/100) == 0
            str = [ 'Percentage complete: ', num2str(100*i/N) ];
            disp(str)
        end
    end
    str = [ 'Percentage complete: ', num2str(100) ];
    disp(str)
    
end
% 3. RK Method
function [ r, v ] = RK( G, M, m, r, v, N, T)    
    dt = T/N;
    kr1 = zeros(3,8);
    kv1 = zeros(3,8);
    kr2 = zeros(3,8);
    kv2 = zeros(3,8);    
    for i = 1:N-1
        for j = 1:8
            kr1(:,j) = v(:,i,j);
            kv1(:,j) = -G*M/norm(r(:,i,j))^3*r(:,i,j);
            for k = 1:8
                if k ~= j
                    kv1(:,j) = kv1(:,j) - (r(:,i,j)-r(:,i,k))*G*m(k)/norm(r(:,i,k)-r(:,i,j))^3;
                end
            end
        end
           
        for j = 1:8
            kr2(:,j) = kr1(:,j) + 0.5*dt*kv1(:,j);
            kv2(:,j) = -G*M/norm(r(:,i,j) + 0.5*dt*kr1(:,j))^3*(r(:,i,j) + 0.5*dt*kr1(:,j));
            for k = 1:8
                if k ~= j
                    kv2(:,j) = kv2(:,j) - ((r(:,i,j) + 0.5*dt*kr1(:,j))-(r(:,i,k) + 0.5*dt*kr1(:,k)))*G*m(k)/norm((r(:,i,k) + 0.5*dt*kr1(:,k))-(r(:,i,j) + 0.5*dt*kr1(:,j)))^3;
                end
            end
            
            r(:,i+1,j) = r(:,i,j) + dt*kr2(:,j);
            v(:,i+1,j) = v(:,i,j) + dt*kv2(:,j);
        end
        
        if mod(i,N/100) == 0
            str = [ 'Percentage complete: ', num2str(100*i/N) ];
            disp(str)
        end
    end
    str = [ 'Percentage complete: ', num2str(100) ];
    disp(str)
    
end





