%--------------------------------------------------------------------------
%Descripci�n: Matriz rigidez viga 3D recta Coordenadas locales
%--------------------------------------------------------------------------
%Input:     - 
%--------------------------------------------------------------------------
%Output:    - 
%--------------------------------------------------------------------------
%Notas:     - Kirchoff
%--------------------------------------------------------------------------
% Written      May 05 2017 by V. Araya      victor.araya@alumnos.usm.cl
% edited       May 05 2017 by V. Araya      victor.araya@alumnos.usm.cl
% Edited       May 05 2018 by F. Mayorga    franco.mayorga@gmail.com
% - It is corrected the bz vector, It must only change the components 
% related to rotation respect to by vector.
%--------------------------------------------------------------------------

function k = k_frame3D_ki(EE,nu,A,I2,I3,J,kap,nod)
    %modulo de corte
        G = EE/(2*(1+nu));
    %Nodos de integraci�n (3)
        ri = [-sqrt(3/5),        0  ,   sqrt(3/5)];
        wi = [5/9       ,      8/9  ,         5/9];
    %Largo
        xi = nod(1,:);
        xj = nod(2,:);
        L  = sqrt((xj - xi)*(xj - xi)');
%--------------------------------------------------------------------------
% Integraci�n Gaussiana (3)
%--------------------------------------------------------------------------
    k = zeros(12,12);
    for i = 1:3;
        %Evaluar funci�n de forma
            [~,MdN_dr] = f_N(ri(i));
        %Evaluar funciones b
            by = f_by(MdN_dr,L)*(2/L);
            bz = f_bz(MdN_dr,L)*(2/L);
            byv = f_byv(L);
            byw = f_byw(L);
            ba = f_ba(MdN_dr)*(2/L);
            bt = f_bt(MdN_dr)*(2/L);
        %Suma a rigidez    
            k = k + wi(i)*(by'*EE*I2*by + ...
                           bz'*EE*I3*bz + ...
                           byv'*G*A*kap*byv + ...
                           byw'*G*A*kap*byw + ...
                           ba'*EE*A*ba + ...
                           bt'*G*J*bt)*(L/2);    
    end

% END Main Function
%--------------------------------------------------------------------------    
    
%--------------------------------------------------------------------------
% Funciones de forma
%--------------------------------------------------------------------------
    function [MN,MdN_dr] = f_N(r)
        %Se evaluan las funciones de forma
        %   r     : vector columna
        %   MN    : [r x 3]
        %   MNh_dr: [r x 3]

        % Shape function
            N1 = (1/2)*(1 - r);
            N2 = (1/2)*(1 + r);
            N3 = 1 -r.^2;

            MN = [N1,N2,N3];

        % Derivadas
            dN1_dr = -(1/2);
            dN2_dr = (1/2);
            dN3_dr = -2*r;

            MdN_dr = [dN1_dr,dN2_dr,dN3_dr]; 
    return
%--------------------------------------------------------------------------
% flexi�n
%--------------------------------------------------------------------------
    function by = f_by(MdN_dr,L)
        % flexion y
            by = [0 ,0,-3*MdN_dr(3)/(2*L) ,0,-MdN_dr(1) + (3/4)*MdN_dr(3),0,...
                  0 ,0,+3*MdN_dr(3)/(2*L),0,-MdN_dr(2) + (3/4)*MdN_dr(3),0];
    return

    function bz = f_bz(MdN_dr,L)
        % flexion z
            bz = [0 ,-3*MdN_dr(3)/(2*L),0,0,0,MdN_dr(1) - (3/4)*MdN_dr(3),...
                  0 ,3*MdN_dr(3)/(2*L) ,0,0,0,MdN_dr(2) - (3/4)*MdN_dr(3)];
    return
%--------------------------------------------------------------------------
% Corte
%--------------------------------------------------------------------------
    function byv = f_byv(L)
        % Corte v
            byv = [0 ,-1/L + (1/L),0,0,0,-1/2 + 1/2,...
                   0 ,(1/L -(1/L)),0,0,0,-1/2 + 1/2];
    return

    function byw = f_byw(L)
        % Corte w
            byw = [0,0,-1/L + 1/L,0,1/2 - 1/2,0,...
                   0,0,1/L - 1/L ,0,1/2 - 1/2,0];
    return
%--------------------------------------------------------------------------
% Axial
%--------------------------------------------------------------------------
    function ba = f_ba(MdN_dr)
        % Axial
            ba = [MdN_dr(1),0,0,0,0,0,   MdN_dr(2),0,0,0,0,0];
    return
%--------------------------------------------------------------------------
% Torque
%--------------------------------------------------------------------------
    function bt = f_bt(MdN_dr)
        % Axial
            bt = [0,0,0,MdN_dr(1),0,0,   0,0,0,MdN_dr(2),0,0];
    return