function Mtest

%Structural Mechanics
%Matrix analysis of frame
%Program created by Dr. Antonio Gil
%Program Developed by Dr. Oubay Hassan and Shen Ma  2024.6 - 2024.11


clc;close all;

datafilename = input('Enter the input file name:');

% =========================================================================
% Data input phase. Read information from input data file
% =========================================================================
[ndof,X,connectivity,Ndof,Nelem,totaldof,E,A,I,H,F,U,alldof,Elastic,Bnode,m,forceunits,lengthunits,UDLFinternal]=readframeputdata1011(datafilename);

%assemble global stiffness matrix
[K]=compute_stiffness_matrix(ndof,X,Ndof,connectivity,E,A,I,totaldof,Nelem);


% =========================================================================
%  Modify matrix phase. Modify matrix according boundary conditions
% =========================================================================

%index fix and free degrees of freedom

fixdof  = find(alldof == 1);
freedof = find(alldof == 0);

%modify the K matrix according to fixed support
K_mod = K;
F_mod = F;
%zero the coefficients of the  corresponding line and column
K_mod(fixdof, :) = 0;
K_mod(:, fixdof) = 0;
%place 1 in the diagonal location of the same equation.
for iq = 1:length(fixdof)
    K_mod(fixdof(iq), fixdof(iq)) = 1;
end

%modify the F matrix according to prescribed displacement


% Now U1 is all prescribed displacement
U1 = U;

dpredof = find(U == 0);
predof  = find(U ~= 0);

if predof ~= 0
    for iu = 1:length(predof)
        for iux = 1:length(freedof)
            % F_mod(dpredof(iux)) is degree of freedom no prescribed displacement
            F_mod(freedof(iux)) = F_mod(freedof(iux)) - U(predof(iu))*K(freedof(iux),predof(iu));
            F_mod(predof)  = U(predof);
        end
    end
else
end

%Modify the K matrix according to elastic supports
El = Elastic(fixdof);
% Place elastic support value in the diagonal location of the same equation.
for iel = 1:length(fixdof)
    K_mod(fixdof(iel), fixdof(iel)) = K_mod(fixdof(iel), fixdof(iel))+El(iel);
end


% =========================================================================
%  Computing unknown value phase. Computing unknown displacements, final
%  displacements, unknown reactions, internal forces
% =========================================================================

%computing unknown displacements
U = K_mod \ F_mod;
U(fixdof) = 0;
% unknown displacements + prescribed displacement
U2 = U+U1;
F_mod = K_mod * U2;
R = K(fixdof,:) * U2 - F(fixdof);
U = U2;
F = F_mod;
F(fixdof) = R;


%computing internal forces
[Finternal]=compute_internal_forces(Nelem,X,connectivity,E,A,I,U2,UDLFinternal);
%print results in a tabular format
% printresults(X,U,F,H,Bnode,Ndof,Finternal,connectivity,Nelem,forceunits,lengthunits);
%plot structure only
% plotFrameStructure(X, connectivity, H)
%plot structure displacements
% Ud = U(Ndof(:,1:2)); scale = 100000;
% plotDeformedStructure(X, connectivity, Ud, scale)
%plot internal forces
factor=1;offset=0;
plot_internal_forces(X,connectivity,Finternal,factor,offset,m)


% =========================================================================
%  Sub Function
% =========================================================================


    function [K]=compute_stiffness_matrix(ndof,X,Ndof,connectivity,E,A,I,totaldof,Nelem)

        K = zeros(totaldof,totaldof);
        %loop over elements
        for ielem=1:size(connectivity,1)
            %extract connectivity information
            i=connectivity(ielem,1);
            j=connectivity(ielem,2);
            idof = Nelem(:,1:3,ielem);
            jdof = Nelem(:,4:6,ielem);
            %compute geometric information
            Xij=X(j,:)-X(i,:);
            L=norm(Xij,2); c=Xij(1)/L; s=Xij(2)/L;
            %frame transformation matrix
            T=[c s 0; -s c 0; 0 0 1];
            %compute element local stiffness blocks
            kii=[E(ielem)*A(ielem)/L               0                            0;...
                0             12*E(ielem)*I(ielem)/L^3     6*E(ielem)*I(ielem)/L^2;...
                0              6*E(ielem)*I(ielem)/L^2     4*E(ielem)*I(ielem)/L];

            kij=[-E(ielem)*A(ielem)/L,               0,                            0;...
                0,             -12*E(ielem)*I(ielem)/L^3,     6*E(ielem)*I(ielem)/L^2;...
                0,             - 6*E(ielem)*I(ielem)/L^2,     2*E(ielem)*I(ielem)/L];

            kji=[-E(ielem)*A(ielem)/L,               0,                            0;...
                0,             -12*E(ielem)*I(ielem)/L^3,    -6*E(ielem)*I(ielem)/L^2;...
                0,               6*E(ielem)*I(ielem)/L^2,     2*E(ielem)*I(ielem)/L];

            kjj=[E(ielem)*A(ielem)/L,               0,                             0;...
                0,             12*E(ielem)*I(ielem)/L^3,     -6*E(ielem)*I(ielem)/L^2;...
                0,             -6*E(ielem)*I(ielem)/L^2,      4*E(ielem)*I(ielem)/L];

            %compute element global stiffness blocks
            Kii=T'*kii*T;Kij=T'*kij*T;Kji=T'*kji*T;Kjj=T'*kjj*T;
            %assemble global stiffness matrix
            K(idof,idof)=K(idof,idof)+Kii; K(idof,jdof)=K(idof,jdof)+Kij;
            K(jdof,idof)=K(jdof,idof)+Kji; K(jdof,jdof)=K(jdof,jdof)+Kjj;
        end
    end


    function [Finternal]=compute_internal_forces(Nelem,X,connectivity,E,A,I,U2,UDLFinternal)
        %initial allocation
        Finternal=zeros(size(connectivity,1),2,3);

        %loop over elements
        for ielem=1:size(connectivity,1)


            %extract connectivity information
            i=connectivity(ielem,1);
            j=connectivity(ielem,2);


            idof = Nelem(:,1:3,ielem);
            jdof = Nelem(:,4:6,ielem);

            %compute geometric information
            Xij=X(j,:)-X(i,:);
            L=norm(Xij,2); c=Xij(1)/L; s=Xij(2)/L;
            %frame transformation matrix
            T=[c s 0; -s c 0; 0 0 1];

            %local displacements
            ui=U2(idof,1);uj=U2(jdof,1);


            %compute element local stiffness blocks
            kii=[E(ielem)*A(ielem)/L               0                            0;...
                0             12*E(ielem)*I(ielem)/L^3     6*E(ielem)*I(ielem)/L^2;...
                0              6*E(ielem)*I(ielem)/L^2     4*E(ielem)*I(ielem)/L];

            kij=[-E(ielem)*A(ielem)/L,               0,                            0;...
                0,             -12*E(ielem)*I(ielem)/L^3,     6*E(ielem)*I(ielem)/L^2;...
                0,             - 6*E(ielem)*I(ielem)/L^2,     2*E(ielem)*I(ielem)/L];

            kji=[-E(ielem)*A(ielem)/L,               0,                            0;...
                0,             -12*E(ielem)*I(ielem)/L^3,    -6*E(ielem)*I(ielem)/L^2;...
                0,               6*E(ielem)*I(ielem)/L^2,     2*E(ielem)*I(ielem)/L];

            kjj=[E(ielem)*A(ielem)/L,               0,                             0;...
                0,             12*E(ielem)*I(ielem)/L^3,     -6*E(ielem)*I(ielem)/L^2;...
                0,             -6*E(ielem)*I(ielem)/L^2,      4*E(ielem)*I(ielem)/L];


            kii=T'*kii*T;kij=T'*kij*T;kji=T'*kji*T;kjj=T'*kjj*T;

            %compute local internal forces

            fi=kii*ui+kij*uj;fj=kji*ui+kjj*uj;
            fi=T*fi;fj=T*fj;
            %extract axial, shear and bending moment for Node 1
            Finternal(ielem,1,1)=fi(1,1);
            Finternal(ielem,1,2)=fi(2,1);
            Finternal(ielem,1,3)=fi(3,1);
            %extract axial, shear and bending moment for Node 2
            Finternal(ielem,2,1)=fj(1,1);
            Finternal(ielem,2,2)=fj(2,1);
            Finternal(ielem,2,3)=fj(3,1);
        end

        Finternal = Finternal + UDLFinternal;
    end


    function printresults(X,U,F,H,Bnode,Ndof,Finternal,connectivity,Nelem,forceunits,lengthunits)


        ouputfilename = input('Enter the output file name:');

        fid = fopen(ouputfilename,'wt');
        fprintf(fid,'                 Table results \n');
        fprintf(fid,'       ------------------------------   \n');
        fprintf(fid,'   \n');
        if H == 0
            fprintf(fid,['\n  Displacements results '  'mm' '\n\n']);

            fprintf(fid,'   \n');

            fprintf(fid,['  Node             U' 'mm' '       V' 'mm' '   Theta' ' \n']);

            for i=1:size(X,1)
                idof=Ndof(i,:);
                fprintf(fid,'   %2g       %10.5g  %10.5g  %10.5g\n', i, U(idof,:)*1e3);
            end
            fprintf(fid,'   \n');
            fprintf(fid,'   \n');




            fprintf(fid,['\n  Node Forces '  forceunits '\n\n']);
            fprintf(fid,['\n  Node             FX' forceunits '      FY' forceunits '    M' forceunits '*' lengthunits '\n']);
            for i=1:size(Bnode,1)
                idof=Ndof(Bnode(i),:);
                fprintf(fid,'   %2g       %10.5g  %10.5g  %10.5g\n', Bnode(i), F(idof,:));
            end
            fprintf(fid,'   \n');
            fprintf(fid,['\n  Internal Forces '  forceunits '\n\n']);

            fprintf(fid,'   \n');
        else
            fprintf(fid,['\n  Displacements results '  'mm' '\n\n']);

            fprintf(fid,'   \n');

            fprintf(fid, '  Element       U (mm)     V (mm)     Theta      U (mm)       V (mm)       Theta\n');
            for i = 1:size(connectivity, 1)
                idof = Nelem(1, 1:3, i);
                jdof = Nelem(1, 4:6, i);
                fprintf(fid, '   %2g       %10.5g  %10.5g  %10.5g  %10.5g  %10.5g  %10.5g\n', i, U(idof,:) * 1e3, U(jdof,:) * 1e3);
            end


            fprintf(fid,['\n  Node Forces '  forceunits '\n\n']);
            fprintf(fid,['\n  Node             FX' forceunits '      FY' forceunits '    M' forceunits '*' lengthunits '\n']);
            for i=1:size(Bnode,1)
                idof=Ndof(Bnode(i),:);
                fprintf(fid,'   %2g       %10.5g  %10.5g  %10.5g\n', Bnode(i), F(idof,:));
            end
            fprintf(fid,'   \n');
            fprintf(fid,['\n  Internal Forces '  forceunits '\n\n']);

            fprintf(fid,'   \n');

        end



        %Axial force results
        k=1;
        fprintf(fid,['  Axial Force'  forceunits '\n']);
        fprintf(fid,'  Element        Node1        Node2 \n');
        for ielem=1:size(Finternal,1)
            % Finternal(ielem,2,k) = -Finternal(ielem,2,k);
            fprintf(fid,'   %2g       %10.5g    %10.5g\n', ielem, Finternal(ielem,1,k), Finternal(ielem,2,k));
        end

        %Shear force results
        k=2;
        fprintf(fid,['  Shear Force'  forceunits '\n']);
        fprintf(fid,'  Element        Node1        Node2 \n');
        for ielem=1:size(Finternal,1)
            % Finternal(ielem,2,k) = -Finternal(ielem,2,k);
            fprintf(fid,'   %2g       %10.5g    %10.5g\n', ielem, Finternal(ielem,1,k), Finternal(ielem,2,k));
        end

        %Bending moment results
        k=3;
        fprintf(fid,['  Bending moment' forceunits '*' lengthunits '\n']);
        fprintf(fid,'  Element        Node1        Node2 \n');
        for ielem=1:size(Finternal,1)
            % Finternal(ielem,2,k) = -Finternal(ielem,2,k);
            fprintf(fid,'   %2g       %10.5g    %10.5g\n', ielem, Finternal(ielem,1,k), Finternal(ielem,2,k));
        end

    end

% Function Only plot Frame Structure
    function plotFrameStructure(X, connectivity, H)
        % X: Node coordinates matrix (Nx2)
        % connectivity: Element connectivity matrix (Mx2), each row represents two connected nodes
        % H: Hinge connection vector (Mx2), each element has two values indicating if each end is hinged (1) or fixed (0)

        figure;
        hold on;

        % Loop through each frame element and plot
        for i = 1:size(connectivity, 1)
            node1 = X(connectivity(i, 1), :);  % Get coordinates of the first node
            node2 = X(connectivity(i, 2), :);  % Get coordinates of the second node

            % Calculate the midpoint coordinates of the element
            midpoint = (node1 + node2) / 2;

            if H(i, 1) == 1 && H(i, 2) == 1
                % Both ends are hinged - plot the entire element with blue dashed line
                plot([node1(1), node2(1)], [node1(2), node2(2)], '--b', 'LineWidth', 2);
            elseif H(i, 1) == 1 && H(i, 2) == 0
                % One end (near node1) is hinged - plot the half near node1 with blue dashed line, the other half with blue solid line
                plot([node1(1), midpoint(1)], [node1(2), midpoint(2)], '--b', 'LineWidth', 2);
                plot([midpoint(1), node2(1)], [midpoint(2), node2(2)], '-b', 'LineWidth', 2);
            elseif H(i, 1) == 0 && H(i, 2) == 1
                % One end (near node2) is hinged - plot the half near node2 with blue dashed line, the other half with blue solid line
                plot([node1(1), midpoint(1)], [node1(2), midpoint(2)], '-b', 'LineWidth', 2);
                plot([midpoint(1), node2(1)], [midpoint(2), node2(2)], '--b', 'LineWidth', 2);
            else
                % Both ends are fixed - plot the entire element with blue solid line
                plot([node1(1), node2(1)], [node1(2), node2(2)], '-b', 'LineWidth', 2);
            end
        end

        % Plot node positions
        plot(X(:,1), X(:,2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

        % Set plot parameters
        axis equal;
        grid on;
        title('2D Frame Structure with Hinge and Fixed Connections');
        xlabel('X');
        ylabel('Y');
        hold off;
    end


% Function plot structure displacements
    function plotDeformedStructure(X, connectivity, displacements, scale)
        % X: Matrix of original node coordinates (Nx2)
        % connectivity: Element connectivity matrix (Mx2)
        % displacements: Matrix of node displacements (Nx2)
        % scale: Scale factor for displacements

        % Calculate the deformed node positions
        X_deformed = X + displacements * scale;

        figure;
        hold on;

        % Plot original structure
        for i = 1:size(connectivity, 1)
            node1 = X(connectivity(i, 1), :);
            node2 = X(connectivity(i, 2), :);
            plot([node1(1), node2(1)], [node1(2), node2(2)], '-k', 'LineWidth', 1);
        end

        % Plot deformed structure
        for i = 1:size(connectivity, 1)
            node1 = X_deformed(connectivity(i, 1), :);
            node2 = X_deformed(connectivity(i, 2), :);
            plot([node1(1), node2(1)], [node1(2), node2(2)], '-r', 'LineWidth', 2); % Deformed structure in red
        end

        % Plot original nodes
        plot(X(:,1), X(:,2), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k'); % Original nodes

        % If there are prescribed displacements
        % plot(X_deformed(:,1), X_deformed(:,2), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); % Deformed nodes

        % Set plot parameters
        axis equal;
        grid on;
        title(['Deformed Structure with Scale Factor = ', num2str(scale)]);
        xlabel('X');
        ylabel('Y');
        legend('Original Structure', 'Deformed Structure');
        hold off;
    end


% Function plot Internal
    function plot_internal_forces(X,connectivity,Int_forces,factor,tol,m)
        %
        coodinate1 = zeros (22,2);
        M          = zeros (22,1);
        figure;
        hold on;
        for i = 1:size(connectivity, 1)
            node1 = connectivity(i, 1);
            node2 = connectivity(i, 2);
            nodei = X(connectivity(i, 1), :);
            nodej = X(connectivity(i, 2), :);

            Xij=X(node2,:)-X(node1,:);
            n=Xij/norm(Xij);v=[-n(2) n(1)];
            XT = norm(Xij);
            for im = 1:21
    
                x = XT/21 * (im-1);
                xn = XT/21 * (im);
                % create Line
                coodinate1(im,1) = X(node1,1) + (X(node2,1) - X(node1,1))/21 * (im-1);
                coodinate1(im,2) = X(node1,2) + (X(node2,2) - X(node1,2))/21 * (im-1);
                coodinate1(im+1,1) = X(node1,1) + (X(node2,1) - X(node1,1))/21 * (im) ;
                coodinate1(im+1,2) = X(node1,2) + (X(node2,2) - X(node1,2))/21 * (im) ;
                x1 = coodinate1(im,1);
                y1 = coodinate1(im,2);
                x2 = coodinate1(im+1,1);
                y2 = coodinate1(im+1,2);
                % create Moment
                k = 3;
                M(im)   = Int_forces(i,1,k) - m(i)*x/2+m(i)*x^2/2;
                M(im+1) = Int_forces(i,1,k) - m(i)*(xn)/2+m(i)*(xn)^2/2;

                
                plot([x1, x2], [y1, y2], '-k', 'LineWidth', 1);
                hold on;
                plot([x1, x2], [M(im)*factor, M(im+1)*factor], '-k', 'LineWidth', 1);


            end



        end
        %








    end






end





