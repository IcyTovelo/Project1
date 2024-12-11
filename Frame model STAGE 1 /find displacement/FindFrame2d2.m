function FindFrame2d2

%Structural Mechanics
%Matrix analysis of frame
%Program created by Dr. Antonio Gil
%Program Developed by Dr. Oubay Hassan and Shen Ma  2024.6 - 2024.11


clc;close all;

datafilename = input('Enter the input file name:');

% =========================================================================
% Data input phase. Read information from input data file
% =========================================================================
[ndof,X,connectivity,Ndof,Nelem,totaldof,E,A,I,H,F,U,alldof,Elastic,Bnode,w,forceunits,lengthunits,UDLFinternal]=readframeputdata1011(datafilename);


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
printresults(X,U,F,H,Bnode,Ndof,Finternal,connectivity,Nelem,forceunits,lengthunits);

%plot structure only
% plotFrameStructure(X, connectivity, H)

%plot structure displacements

%plot internal forces
scale1 = 0.1;
plot_BM1_forces(X, connectivity, Finternal, w, scale1)

plot_deformation(X,connectivity,U,Nelem,w,E,I,A,Finternal)
%plot structure displacements
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
            Finternal(ielem,1,2)= fi(2,1);
            Finternal(ielem,1,3)=fi(3,1);
            %extract axial, shear and bending moment for Node 2
            Finternal(ielem,2,1)= fj(1,1);
            Finternal(ielem,2,2)= fj(2,1);
            Finternal(ielem,2,3)= fj(3,1);
        end

        Finternal = Finternal + UDLFinternal;

        Finternal(:,2,3) = - Finternal(:,2,3);

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

% Function to plot Bending Moment and Shear Force
    function plot_BM_forces(X, connectivity, Finternal, w)

        coo = zeros(size(X,1),2);
        M   = zeros(size(connectivity,1),11);
        V   = zeros(size(connectivity, 1),11);
        % Create unified figures for bending moment and shear force
        figure('Name', 'Bending Moment Distribution', 'NumberTitle', 'off'); % Create bending moment figure
        hold on; % Allow overlaying plots
        xlabel('Length (m)');
        ylabel('Bending Moment (kN·m)');
        title('Bending Moment Distribution for All Elements');
        grid on;

        figure('Name', 'Shear Force Distribution', 'NumberTitle', 'off'); % Create shear force figure
        hold on;
        xlabel('Length (m)');
        ylabel('Shear Force (kN)');
        title('Shear Force Distribution for All Elements');
        grid on;

        % Loop through each element
        for ie = 1:size(connectivity,1)

            i = connectivity(ie,1);
            j = connectivity(ie,2);
            Xa = X((connectivity(ie,1)),1);
            Xb = X((connectivity(ie,2)),1);
            Ya = X((connectivity(ie,1)),2);
            Yb = X((connectivity(ie,2)),2);

            VA = Finternal(ie,1,2);
            VB = Finternal(ie,2,2);
            MA = Finternal(ie,1,3);
            MB = Finternal(ie,2,3);

            Xij=X(j,:)-X(i,:);
            L=norm(Xij,2);

            % Calculate bending moment and shear force at 11 points
            for i = 1:11
                coo(i,1) = Xa + (Xb - Xa) / 10 * (i - 1);
                coo(i,2) = Ya + (Yb - Ya) / 10 * (i - 1);

                x = (i - 1) / (11-1)*L;

                % Shear force calculation (derivative of bending moment with respect to x)
                V(ie, i) = VA + w(ie) * x;

                % Bending moment calculation
                M(ie,i) = -MA + VA * x + (w(ie) * x^2) / 2;
            end

            M1 = -M;
            M1(ie,11) = MB;

            disp('Bending Moments:');
            disp(M1);
            disp('Shear Force:');
            disp(V);

            % Plot bending moment
            figure(1); % Switch to bending moment figure
            % Plot element coordinates
            plot(coo(:, 1), coo(:, 2), '-bo', 'DisplayName', ['Coordinates (Element ' num2str(ie) ')']);
            % Plot bending moment curve
            plot(coo(:, 1), M1(ie, :), '-x', 'DisplayName', ['Element ' num2str(ie)]);
            % Draw vertical lines from element axis to bending moment points
            for ii = 1:11
                line([coo(ii, 1), coo(ii, 1)], [coo(ii, 2), M1(ie, ii)], 'Color', 'k', 'LineStyle', '--');
            end

            % Plot shear force
            figure(2); % Switch to shear force figure
            % Plot element coordinates
            plot(coo(:, 1), coo(:, 2), '-bo', 'DisplayName', ['Coordinates (Element ' num2str(ie) ')']);
            % Plot shear force curve
            plot(coo(:, 1), V(ie, :), '-o', 'DisplayName', ['Element ' num2str(ie)]);
            % Draw vertical lines from element axis to shear force points
            for ii = 1:11
                line([coo(ii, 1), coo(ii, 1)], [coo(ii, 2), V(ie, ii)], 'Color', 'k', 'LineStyle', '--');
            end
        end

        % Add legends to the figures
        figure(1);
        hold off;

        figure(2);
        hold off;
    end

% Function to plot Bending Moment and Shear Force with scaling
    function plot_BM1_forces(X, connectivity, Finternal, w, scale1)

        % Create unified figures for bending moment and shear force
        figure(1);
        hold on;
        xlabel('Local X-axis (m)');
        ylabel('Local Y-axis (m)');
        title('Bending Moment Distribution with Original Structure');
        grid on;

        figure(2);
        hold on;
        xlabel('Local X-axis (m)');
        ylabel('Local Y-axis (m)');
        title('Shear Force Distribution with Original Structure');
        grid on;

        % Loop through each element
        for ie = 1:size(connectivity, 1)
            % Extract node indices
            i = connectivity(ie, 1);
            j = connectivity(ie, 2);

            % Extract node coordinates
            Xa = X(i, 1); Ya = X(i, 2);
            Xb = X(j, 1); Yb = X(j, 2);

            % Internal forces
            VA = Finternal(ie, 1, 2);
            VB = Finternal(ie, 2, 2);
            MA = Finternal(ie, 1, 3);
            MB = Finternal(ie, 2, 3);

            % Element length and direction
            L = norm([Xb - Xa, Yb - Ya]);
            dx = (Xb - Xa) / L;
            dy = (Yb - Ya) / L;

            % Local coordinates
            x_local = linspace(0, L, 11); % 11 points along the element

            % Calculate bending moment and shear force
            M_local = arrayfun(@(x) -MA + VA * x + (w(ie) * x^2) / 2, x_local);
            V_local = arrayfun(@(x) VA + w(ie) * x, x_local);

            % Calculate global coordinates
            x_global = Xa + dx * x_local;
            y_global = Ya + dy * x_local;

            % Perpendicular offset direction
            offset_x = -dy; % Perpendicular x component
            offset_y = dx;  % Perpendicular y component
            scale = scale1 * L; % Scaling factor for visualization

            % Plot original structure
            % In bending moment figure
            figure(1);
            plot([Xa, Xb], [Ya, Yb], 'k-', 'LineWidth', 1.5, 'DisplayName', 'Original Structure'); % Element axis

            % Plot bending moment offsets
            for k = 1:length(x_local)
                bx = x_global(k) + offset_x * -M_local(k) * scale;
                by = y_global(k) + offset_y * -M_local(k) * scale;
                line([x_global(k), bx], [y_global(k), by], 'Color', 'r', 'LineStyle', '--');
            end
            plot(x_global + offset_x * -M_local * scale, y_global + offset_y * -M_local * scale, 'r-', ...
                'DisplayName', ['Element ' num2str(ie)]); % Bending moment curve

            % In shear force figure
            figure(2);
            plot([Xa, Xb], [Ya, Yb], 'k-', 'LineWidth', 1.5, 'DisplayName', 'Original Structure'); % Element axis

            % Plot shear force offsets
            for k = 1:length(x_local)
                sx = x_global(k) + offset_x * V_local(k) * scale;
                sy = y_global(k) + offset_y * V_local(k) * scale;
                line([x_global(k), sx], [y_global(k), sy], 'Color', 'b', 'LineStyle', '--');
            end
            plot(x_global + offset_x * V_local * scale, y_global + offset_y * V_local * scale, 'b-', ...
                'DisplayName', ['Element ' num2str(ie)]); % Shear force curve
        end

        % Add legends to the figures
        figure(1);
        hold off;

        figure(2);
        hold off;
    end



    function plot_deformation(X,connectivity,U,Nelem,w,E,I,A,Finternal)

        for ielem=1:size(connectivity,1)
            %extract connectivity information
            i=connectivity(ielem,1);
            j=connectivity(ielem,2);
            idof = Nelem(:,1:3,ielem);
            jdof = Nelem(:,4:6,ielem);
            % global y and rotation
            y1   = U(Nelem(:,2,ielem));
            r1   = U(Nelem(:,3,ielem));
            y2   = U(Nelem(:,5,ielem));
            r2   = U(Nelem(:,6,ielem));


            %compute geometric information
            Xij=X(j,:)-X(i,:);
            L=norm(Xij,2); c=Xij(1)/L; s=Xij(2)/L;

             s1 = asind(s);
             c1 = acosd(c); 
             s2 = sin((90-s1)*pi/180);  
             % local y and rotation
             ya = y1 / s2;
             ra = r1;
             yb = y2 / s2;
             rb = r2;

            d = ya;
            c = ra;

         % a  = (((w(ielem) * L^4) / 24)  + ((b * L^2) / 2) + c * L + d - yb) * (6 / L^3);
         % a  = (((w(ielem) * L^3) / 6)  + b * L + c - rb) * (2 / L^2);
         %    %

% 构造矩阵 A 和常数项 B
MA  = [L^3/6, L^2/2; L^2/2, L];

Mb  = [yb - (w(ielem) * L^4 / (24 * E(ielem)*I(ielem)));rb - (w(ielem) * L^3 / (6 * E(ielem)*I(ielem)))];


% 解方程组
Rab = MA \ Mb;

% 输出结果
a = Rab(1);
b = Rab(2);


Yd = zeros(11,1);
Rd = zeros(11,1);
Nd = zeros(11,1);
Xxd = zeros(11,1);
for i = 1:11
  
    Xd = ((i-1)/10*L);

    Yd(i) = ((w(ielem) * Xd^4) / (24 * E(ielem)*I(ielem))) + (a * Xd^3 / 6) + ((b * Xd^2) / 2) +c * Xd + d;
    Rd(i) = ((w(ielem) * Xd^3) / (6 * E(ielem)*I(ielem))) + (a * Xd^2 / 2)  + b * Xd + c;
    

end

Yd(1)  = ya;
Rd(1)  = ra;
Yd(11) = yb;
Rd(11) = rb;


% 添加  Finternal
Nd(1) = Finternal(1,1,ielem);
for i = 1:11
  
    Xd = ((i-1)/10*L);
    if i == 1
       Xxd(i) = U(Nelem(:,1,ielem));
        Nd(i) = Finternal(1,1,ielem);
    else
        Xxd(i) = Nd(i - 1) / (E(ielem) * A(ielem) / Xd) + Xxd(i-1);
       Nd(i-1) = (E(ielem) * A(ielem) / Xd) * (Xxd(i) - Xxd(i-1));
    end
   
end

disp(Yd);
disp(Rd);
disp(Xxd)


        end
    end



end