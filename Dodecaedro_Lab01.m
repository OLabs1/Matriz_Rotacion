%%%%%%%Graficando el Dodecaedro(12 caras)%%%%%%%%%
%Sólido platonico porque todas su caras son regulares iguales entre sí
%Todas sus coordenadas quedan definidos aplicando la propoción Aurea y toma base las coordenadas enteras de un cubo.

clc; clear all

l=10 %Factor de escala, recordar longitud del cubo interno
P=(1+sqrt(5))/2  %Proporción Aurea
q=1/P;

%DIBUJANDO LOS VERTICES

A=[P;0;q]*l;
B=[-P;0;q]*l;
C=[-P;0;-q]*l;
D=[P;0;-q]*l;

E=[q;P;0]*l;
F=[q;-P;0]*l;
G=[-q;-P;0]*l;
H=[-q;P;0]*l;

I=[0;q;P]*l;
J=[0;q;-P]*l;
K=[0;-q;-P]*l;
L=[0;-q;P]*l;

M=[l;l;l];
N=[l;-l;l];
O=[-l;-l;l];
PV=[-l;l;l];
Q=[-l;l;-l];
R=[l;l;-l];
S=[l;-l;-l];
T=[-l;-l;-l];
na=NaN(3,1);

%DIBUJANDO LA CARA DEL DODECAEDRO
%CARA1-AMILN
CA1=horzcat(A,M,I,L,N,A);
%CARA2- IPHEM
CA2=horzcat(M,E,H,PV,I);
%CARA3- LIPBO
CA3=horzcat(PV,B,O,L);
%CARA4- LNFGO
CA4=horzcat(N,F,G,O);
%CARA5- NFSDA
CA5=horzcat(F,S,D,A);
%CARA6- ADREM
CA6=horzcat(D,R,E);
%CARA7- ERJQH
CA7=horzcat(R,J,Q,H);
%CARA8- DRJKS
CA8=horzcat(S,K,J);
%CARA9- KJQTC
CA9=horzcat(Q,C,T,K);
%CARA10- PHCGO
CA10=horzcat(B,C);
%CARA11- 
CA11=horzcat(T,G);
%CARA12 Por defecto despues de haber llenado el resto.

%Armando la matriz del dodecaedro
M0C=horzcat(CA1,CA2,CA3,CA4,na,CA5,na,CA6,na,CA7,na,CA8,CA9,na,CA10,na,CA11)
XC0=M0C(1,:);
YC0=M0C(2,:);
ZC0=M0C(3,:);


%% 
%Realizando rotaciones
figure
subplot(2,2,1);plot3(XC0,YC0,ZC0);axis equal;grid on;title('Dodecaedro dibujado');
set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
xlabel('eje x');ylabel('eje y');zlabel('eje z')
axis on

%Arreglo para graficar el vector de origen Z
scale=1;%Para establecer la escala del vector habilitado.
az = [0 0 -2*l];% Punto inicio para graficar el vector del origen Z
bz = [0 0 2*l];% Punto fin para graficar el vector del origen Z
inicioz=az;% "az" es el origen de mi vector para dibujar del punto a hacia bz y se crea un nuevo punto "cz"
cz = bz-az;% cz es la nueva coordinada de b con origen en a, este arreglo es para crear un nuevo origen.
finz=cz;

%Arreglo para graficar el vector de origen y
scale=1;%Para establecer la escala del vector habilitado.
ay = [0 -2*l  0];% Punto inicio para graficar el vector del origen y
by = [0 2*l 0];% Punto fin para graficar el vector del origen y
inicioy=ay;% "ay" es el origen de mi vector para dibujar del punto a hacia bz y se crea un nuevo punto "cy"
cy = by-ay;% cy es la nueva coordinada de b con origen en a, este arreglo es para crear un nuevo origen.
finy=cy;

%Arreglo para graficar el vector de origen X
scale=1;%Para establecer la escala del vector habilitado.
ax = [-2*l 0 0];% Punto inicio para graficar el vector del origen x
bx = [2*l 0 0];% Punto fin para graficar el vector del origen x
iniciox=ax;% "ax" es el origen de mi vector para dibujar del punto a hacia bz y se crea un nuevo punto "cx"
cx = bx-ax;% cx es la nueva coordinada de b con origen en a, este arreglo es para crear un nuevo origen.
finx=cx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Realizando rotaciones sobre el eje Z
for Rota=0:30:360;
    M1C=rotz(Rota)*M0C;  %Rotando 15 grados en cada evaluación
    X1C=M1C(1,:)
    Y1C=M1C(2,:)
    Z1C=M1C(3,:)
    subplot(2,2,2);plot3(X1C,Y1C,Z1C);title(sprintf('Rotación sobre el eje Z de %d grados',Rota));axis equal;grid on;
    hold on
    quiver3(inicioz(:,1), inicioz(:,2), inicioz(:,3), finz(:,1), finz(:,2), finz(:,3),scale);% Dibujando el vector
    subplot(2,2,2);plot3(X1C(4),Y1C(4),Z1C(4),'*r'); %Graficando un punto de referencia
    subplot(2,2,2);plot3(X1C(3),Y1C(3),Z1C(3),'*g'); %Graficando un punto de referencia
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')
    axis on
    pause(1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Realizando rotaciones sobre el eje y
for Rota=0:30:360;
    M2C=roty(Rota)*M0C;  %Rotando 15 grados en cada evaluación
    X2C=M2C(1,:)
    Y2C=M2C(2,:)
    Z2C=M2C(3,:)
    subplot(2,2,3);plot3(X2C,Y2C,Z2C);title(sprintf('Rotación sobre el eje Y de %d grados',Rota));axis equal;grid on;
    hold on
    quiver3(inicioy(:,1), inicioy(:,2), inicioy(:,3), finy(:,1), finy(:,2), finy(:,3),scale);% Dibujando el vector
    subplot(2,2,3);plot3(X2C(4),Y2C(4),Z2C(4),'*r'); %Graficando un punto de referencia
    subplot(2,2,3);plot3(X2C(3),Y2C(3),Z2C(3),'*g'); %Graficando un punto de referencia
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')
    axis on
    pause(1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Realizando rotaciones sobre el eje X
for Rota=0:30:360;
    M3C=rotx(Rota)*M0C;  %Rotando 15 grados en cada evaluación
    X3C=M3C(1,:)
    Y3C=M3C(2,:)
    Z3C=M3C(3,:)
    subplot(2,2,4);plot3(X3C,Y3C,Z3C);title(sprintf('Rotación sobre el eje X de %d grados',Rota));axis equal;grid on;
    hold on
    quiver3(iniciox(:,1), iniciox(:,2), iniciox(:,3), finx(:,1), finx(:,2), finx(:,3),scale);% Dibujando el vector
    subplot(2,2,4);plot3(X3C(4),Y3C(4),Z3C(4),'*r'); %Graficando un punto de referencia
    subplot(2,2,4);plot3(X3C(3),Y3C(3),Z3C(3),'*g'); %Graficando un punto de referencia
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')
    axis on
    pause(1)
end
%% 
%%%% GIRANDO EL OCTAEDRO SOBRE UN EJE ARBITRARIO UN DETERMINADO ANGULO PHI

%Estableciendo las variables, rx ry rz como coordenadas del vector y phi angulo de rotación radianes.
% syms rx ry rz phi
%Creando las coordenadas de x y z del vector.
x=15;y=15;z=15;
%Arreglo para graficar el vector a partir de dos puntos.
scale=1;%Para establecer la escala del vector habilitado.
a = [0 0 0];% Punto 1 desde origen, para partir de otro punto elegir otras coordenadas
b = [x y z];% Los valores de la coordenada de "R"

inicio=a;% "a" es el origen de mi vector para dibujar del punto a hacia b y se crea un nuevo punto "c"
c = b-a;% c es la nueva coordinada de b con origen en a, este arreglo es para crear un nuevo origen.
fin=c;

%Hallando los valores del vector unitario [rx ry rz] a partir de los puntos x y z
%Magnitud del vector desde el origen
Mr=norm(b)
rx=(x/Mr);ry=(y/Mr);rz=(z/Mr); %Obteniendo los valores del vector unitario

%Estableciendo los valores de rotación en grados.
%Recordar que la matriz del octaedro es MOC
for Rotb=0:30:360;
    %Convirtiendo a radiales
    phi=Rotb*pi/180;
    %Creando la matriz de rotación "MR"
    V=1-cos(phi);
    format short
    MR=[rx^2*V+cos(phi) rx*ry*V-rz*sin(phi) rx*rz*V+ry*sin(phi);...
        rx*ry*V+rz*sin(phi) ry^2*V+cos(phi) ry*rz*V-rx*sin(phi);...
        rx*rz*V-ry*sin(phi) ry*rz*V+rx*sin(phi) rz^2*V+cos(phi)]
    
    MR1=MR*M0C  %Rotación de phi grados

    MRX=MR1(1,:);
    MRY=MR1(2,:);
    MRZ=MR1(3,:);

    %Graficando las rotaciones
    figure(2)
    subplot(2,2,1);plot3(XC0,YC0,ZC0);axis equal;grid on;title('Dodecaedro sin rotar');
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')

    subplot(2,2,2);plot3(XC0,YC0,ZC0);axis equal;grid on;title('Dodecaedro y vector arbitrario');
    hold on
    plot3(a(1),a(2),a(3),'*b')  %Dibujando el punto de inicio
    plot3(b(1),b(2),b(3),'*g')  %Dibujando el punto de fin
    quiver3(inicio(:,1), inicio(:,2), inicio(:,3), fin(:,1), fin(:,2), fin(:,3),scale);% Dibujando el vector
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')

    subplot(2,2,3);plot3(MRX,MRY,MRZ);axis equal;grid on;title(sprintf('Octaedro girado sobre un vector arbitrario y %d grados',Rotb));
    hold on
    plot3(a(1),a(2),a(3),'*b')  %Dibujando el punto de inicio
    plot3(b(1),b(2),b(3),'*b')  %Dibujando el punto de fin
    plot3(MRX(4),MRY(4),MRZ(4),'*r'); %Graficando un punto de referencia
    plot3(MRX(3),MRY(3),MRZ(3),'*g'); %Graficando un punto de referencia
    plot3(MRX(2),MRY(2),MRZ(2),'*k'); %Graficando un punto de referencia
    quiver3(inicio(:,1), inicio(:,2), inicio(:,3), fin(:,1), fin(:,2), fin(:,3),scale);% Dibujando el vector
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')
    pause(1)
end
%% 
%Creando las coordenadas de x1 y1 z1 del vector.
x1=-4;y1=7;z1=-8;
%Arreglo para graficar el vector a partir de dos puntos.
scale=1;%Para establecer la escala del vector habilitado.
a1 = [0 0 0];% Punto 1 desde origen, para partir de otro punto elegir otras coordenadas
b1 = [x1 y1 z1];% Los valores de la coordenada de "R"

inicio1=a1;% "a" es el origen de mi vector para dibujar del punto a hacia b y se crea un nuevo punto "c"
c1 = b1-a1;% c es la nueva coordinada de b con origen en a, este arreglo es para crear un nuevo origen.
fin1=c1;

%Hallando los valores del vector unitario [rx ry rz] a partir de los puntos x y z
%Magnitud del vector desde el origen
Mr1=norm(b)
rx1=(x1/Mr);ry1=(y1/Mr1);rz1=(z1/Mr1); %Obteniendo los valores del vector unitario

%Estableciendo los valores de rotación en grados.
%Recordar que la matriz del octaedro es MOC
for Rotb=0:30:360;
    %Convirtiendo a radiales
    phi=Rotb*pi/180;
    %Creando la matriz de rotación "MR"
    V=1-cos(phi);
    format short
    MR3=[rx^2*V+cos(phi) rx*ry*V-rz*sin(phi) rx*rz*V+ry*sin(phi);...
        rx*ry*V+rz*sin(phi) ry^2*V+cos(phi) ry*rz*V-rx*sin(phi);...
        rx*rz*V-ry*sin(phi) ry*rz*V+rx*sin(phi) rz^2*V+cos(phi)]
    
    MR4=MR3*M0C  %Rotación de phi grados

    MRX4=MR4(1,:);
    MRY4=MR4(2,:);
    MRZ4=MR4(3,:);

    %Graficando las rotaciones
    subplot(2,2,4);plot3(MRX4,MRY4,MRZ4);axis equal;grid on;title(sprintf('Octaedro girado sobre un vector arbitrario y %d grados',Rotb));
    hold on
    plot3(a1(1),a1(2),a1(3),'*b')  %Dibujando el punto de inicio
    plot3(b1(1),b1(2),b1(3),'*b')  %Dibujando el punto de fin
    plot3(MRX4(4),MRY4(4),MRZ4(4),'*r'); %Graficando un punto de referencia
    plot3(MRX4(3),MRY4(3),MRZ4(3),'*g'); %Graficando un punto de referencia
    quiver3(inicio1(:,1), inicio1(:,2), inicio1(:,3), fin1(:,1), fin1(:,2), fin1(:,3),scale);% Dibujando el vector
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')
    pause(1)
end
