%%%%%%%Graficando el Octaedro%%%%%%%%%
%Sólido platonico porque todas su caras son regulares iguales entre sí
%Todas sus coordenadas quedan definidos con coordenadas enteras

clc; clear all

L=10
XC0=[0 0 2*L 2*L 0 L 2*L NaN L 0 NaN L 2*L NaN 0 L 0 NaN 2*L L 2*L]
YC0=[0 2*L 2*L 0 0 L 2*L NaN L 2*L NaN L 0 NaN 2*L L 0 NaN 0 L 2*L]
ZC0=[L L L L L 0 L NaN 0 L NaN 0 L NaN L 2*L L NaN L 2*L L]
M0C=[XC0;YC0;ZC0]      %Organizando la matriz del octaedro
%% 
%Realizando rotaciones
figure
subplot(2,2,1);plot3(XC0,YC0,ZC0);axis equal;grid on;title('Octaedro dibujado');
set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
xlabel('eje x');ylabel('eje y');zlabel('eje z')


for Rota=0:30:360;
    M1C=rotz(Rota)*M0C;  %Rotando 15 grados en cada evaluación
    X1C=M1C(1,:)
    Y1C=M1C(2,:)
    Z1C=M1C(3,:)
    subplot(2,2,2);plot3(X1C,Y1C,Z1C);title(sprintf('Rotación sobre el eje Z de %d grados',Rota));axis equal;grid on;
    hold on
    subplot(2,2,2);plot3(X1C(20),Y1C(20),Z1C(20),'*r');  %Graficando un punto de referencia
    hold off
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
    subplot(2,2,3);plot3(X2C(20),Y2C(20),Z2C(20),'*r'); %Graficando un punto de referencia
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')
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
    subplot(2,2,4);plot3(X3C(20),Y3C(20),Z3C(20),'*r'); %Graficando un punto de referencia
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')
    pause(1)
end
%% 
%%%% GIRANDO EL OCTAEDRO SOBRE UN EJE ARBITRARIO UN DETERMINADO ANGULO PHI

%Estableciendo las variables, rx ry rz como coordenadas del vector y phi angulo de rotación radianes.
% syms rx ry rz phi
%Creando las coordenadas de x y z del vector.
x=10;y=9;z=15;

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
    subplot(2,2,1);plot3(XC0,YC0,ZC0);axis equal;grid on;title('Octaedro sin rotar');
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')

    subplot(2,2,2);plot3(XC0,YC0,ZC0);axis equal;grid on;title('Octaedro y vector arbitrario');
    hold on
    plot3(a(1),a(2),a(3),'*b')  %Dibujando el punto de inicio
    plot3(b(1),b(2),b(3),'*g')  %Dibujando el punto de fin
    quiver3(inicio(:,1), inicio(:,2), inicio(:,3), fin(:,1), fin(:,2), fin(:,3),scale);% Dibujando el vector
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')

    subplot(2,2,[3 4]);plot3(MRX,MRY,MRZ);axis equal;grid on;title(sprintf('Octaedro girado sobre un vector arbitrario y %d grados',Rotb));
    hold on
    plot3(a(1),a(2),a(3),'*b')  %Dibujando el punto de inicio
    plot3(b(1),b(2),b(3),'*g')  %Dibujando el punto de fin
    plot3(MRX(20),MRY(20),MRZ(20),'*r'); %Graficando un punto de referencia
    quiver3(inicio(:,1), inicio(:,2), inicio(:,3), fin(:,1), fin(:,2), fin(:,3),scale);% Dibujando el vector
    hold off
    set(gca, 'XDir','reverse') %Invirtiendo la dirección del eje x
    set(gca, 'YDir','reverse') %Invirtiendo la dirección del eje y
    xlabel('eje x');ylabel('eje y');zlabel('eje z')
    pause(1.5)
end
