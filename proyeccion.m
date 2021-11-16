% Esta funcion calcula la matriz H, y la imagen de proyeccion
% Ver.1.5   03-09-2004
%*********************************************************

function [t_proy] = proyeccion(input_points, dm, I)


t_proy=0;
tam=size(I);
b=[];
A=[];
Matriz_H=[];

base_points = [ 0    0
                dm   0
                dm   dm
                0    dm];
            
            
uv=input_points;
xy=base_points;

if size(uv,1)>=4
        
    for c_fil=1:size(uv,1)
        u= uv(c_fil, 1);v= uv(c_fil, 2); 
        x= xy(c_fil, 1);y= xy(c_fil, 2);
        A([c_fil*2-1 c_fil*2],:)=[x y 1.0000 0 0 0 -u*x -u*y;
                                0 0 0 x y 1.0000 -v*x -v*y];
        b([c_fil*2-1 c_fil*2],1)=[u;v];
	end
    
    if size(uv,1)==4
        % En este caso hay solucion directa 
        % Suponemos que H33=1
        
		Vector_H= A \ b;
		Vector_H([9])=1; 
		Matriz_H=[Vector_H([1:3]).'; Vector_H([4:6]).'; Vector_H([7:9]).'];   
		Inv_Matriz_H= inv(Matriz_H);
       
    else
        % Hay una solucion a traves de los minimos cuadrados
        % Suponemos que H33=1
        Vector_H= (A'*A)\(A'*b);
        Vector_H([9])=1; 
        Matriz_H=[Vector_H([1:3]).'; Vector_H([4:6]).'; Vector_H([7:9]).']   
		Inv_Matriz_H= inv(Matriz_H);
    end
    
    %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    % Proceso de generacion de la imagen de salida a partir de la matriz H
    
	ima = []; %---- Borra la imagen de salida
    
    % Calculo las dimensiones finales de la imagen de salida
    %dimen = dimension(tam, Inv_Matriz_H);
    dimen= [max(base_points(:)), max(base_points(:))];
    delta_XY=[];
	u=0;v=0;	
	for x=1:(dimen(1)-1)
        for y=1:(dimen(2)-1) 
            Vector=[x;y;1];
            new_posicion=(Matriz_H*Vector); %Inv_Matriz_H
            if new_posicion(3)~=0
                
                % Pasando a coordenadas homogeneas
                v = new_posicion(1)/new_posicion(3);
                u = new_posicion(2)/new_posicion(3);
                
                int_u = fix(u); int_v = fix(v);
                delta_X = u - int_u; 
                delta_Y = v - int_v;
                
            else
                u = 0; v = 0;
            end
            %:::::::::::::::::::::::::::::::::::::::::::::::::::
            % Calculo de la Interpolacion Bilineal
            
            if  int_u>0 & int_v>0 & int_u<tam(1) & int_v<tam(2)
                ima (y, x)= [double(I(int_u+1,int_v))- double(I(int_u,int_v))]*delta_X+...
                            [double(I(int_u,int_v))- double(I(int_u,int_v+1))]*delta_Y+...
                            [double(I(int_u+1,int_v+1))+double(I(int_u,int_v))-double(I(int_u+1,int_v))-double(I(int_u,int_v+1))]*delta_X*delta_Y+...
                            double(I(int_u,int_v));
            end           

        end
	end
   % Salida de la funcion, t_proy::funcion de proyeccion.
    t_proy = uint8(ima);
    
else 
   % Hay menos de 4 puntos 
   t_proy = 1;
   
end

    


