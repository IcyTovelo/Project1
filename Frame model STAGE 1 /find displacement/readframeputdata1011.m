function  [ndof,X,connectivity,Ndof,Nelem,totaldof,E,A,I,H,F,U,alldof,Elastic,Bnode,w,forceunits,lengthunits,UDLFinternal]=readframeputdata1011(datafilename)

%Structural Mechanics
%Input data Function of Matrix analysis of frame
%Program created by Dr.Oubay Hassan and Shen Ma

% Open the file
fid = fopen(datafilename, 'r');

% Check if the file was successfully opened
if fid == -1
    error('Cannot open file: %s', datafilename);
end

% Read and display the title
AnalysisTitle = fscanf(fid, 'TITLE = %[^\n]', 1);


% Move the file pointer to the next line
fgetl(fid);

% Read the subtitle

AnalysisStuff = fscanf(fid, 'SUBTITLE = %[^\n]', 1);

% Read the input list
stufflist = fscanf(fid, '\n%d %d %d %d %d %d %d', [7,1]);


% =========================================================================
% Data input phase. Read information from input data file 'datafilename.txt'
% =========================================================================


ndof  = 3;
% Unit section
lengthunits = fscanf(fid, '\nLENGTHUNITS = %s', 1);
forceunits  = fscanf(fid, '\nAXIALUNITS = %s', 1);



%-------------------------------------------------------
% Table of node coordinates
%
% 1.X  2.HN  3.HCT  4.Ndof
%-------------------------------------------------------

npoin = fscanf(fid, '\nNODE_COORDINATES = %d', 1);

if npoin ~= stufflist(1,1)
    disp('error: Number of node is wrong')
end

coord = fscanf(fid, '\n%d %f %f', [3,npoin]);
coord = coord';
coordinate = coord(:,1);

coordRows = size(coord,1);
%X: node coordinate
X = coord(:,2:3);





%-------------------------------------------------------
% Table of connectivities
%
% 1.emp  2.A  3.I  4.connectivity
%-------------------------------------------------------


% Total number of elements in the mesh
nelem = fscanf(fid, '\nELEMENTS = %d', 1);

if nelem ~= stufflist(2,1)

    disp('error: Number of elements is wrong')
end

lnods = fscanf(fid, '\n%d %d %d %d %f %f %d %d', [8,nelem]);
lnods = lnods';

%element metrials type
emp = lnods(:,4);
%cross-sectional areas
A = lnods(:,5);
%cross-sectional areas
I = lnods(:,6);

H = lnods(:,7:8);
%store only nodal connectivities in lnods
connectivity = lnods(:,2:3);
% Initialize an array to store the DOF for each node (3 DOFs per node)
nodeDOF = zeros(size(X,1), 3);
Ndof    = zeros(size(X,1), 3);
% Initialize next available DOF and a counter for the number of DOFs assigned
nextDOF = 1;
Nelem = zeros(1,6,size(nelem,1));
ncnode = zeros(size(X,1),3);



% Loop through all elements
for i = 1:nelem
    node1 = lnods(i, 2);  % First node of the element
    node2 = lnods(i, 3);  % Second node of the element
    modifyNode1 = lnods(i, 7);  % Modification flag for node1 (second-to-last column)
    modifyNode2 = lnods(i, 8);

    % Check and assign DOFs for node1


    if nodeDOF(node1, 1) == 0  % Node1 has not been assigned DOFs
        % Assign the first DOF to node1
        nodeDOF(node1, 1) = nextDOF;  % Assign current nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the first DOF

        % Assign the second DOF to node1
        nodeDOF(node1, 2) = nextDOF;  % Assign updated nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the second DOF

        % Assign the third DOF to node1
        nodeDOF(node1, 3) = nextDOF;  % Assign updated nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the third DOF

        Ndof(node1, 1:3) = nodeDOF(node1, 1:3);
        Nelem(1,1:3,i) = nodeDOF(node1, :);



        if modifyNode1 == 1  % Modify node1's third DOF

            ncnode(node1,1:2) =  nodeDOF(node1,1:2);
            ncnode(node1,3)   = -nodeDOF(node1,3);

            Nelem(1,1:3,i)  = abs(ncnode(node1,:));
            Ndof(node1, :)  = ncnode(node1,:);
        else
            ncnode(node1,1:3) =  nodeDOF(node1,1:3);
            Nelem(1,1:3,i)  = nodeDOF(node1, :);
            Ndof(node1, :)  = ncnode(node1,:);
        end

    elseif nodeDOF(node1, 1) ~= 0

        if modifyNode1 == 0 && Ndof(node1, 3) <= 0 % Modify node1's third DOF

            nodeDOF(node1, 3) = nextDOF;  % Assign updated nextDOF
            nextDOF = nextDOF + 1;
            ncnode(node1,:) =  nodeDOF(node1,:);
            Ndof(node1, 3) = ncnode(node1, 3);
            Nelem(1,1:3,i)  = ncnode(node1,:);
        elseif modifyNode1 == 0 && Ndof(node1, 3) >= 0
            ncnode(node1,:) =  nodeDOF(node1,:);
            Nelem(1,1:3,i)  = Ndof(node1, :);
        elseif modifyNode1 == 1
            nodeDOF(node1, 3) = nextDOF;  % Assign updated nextDOF
            ncnode(node1,:) =  nodeDOF(node1,:);
            nextDOF = nextDOF + 1;
            Nelem(1,1:3,i)  = ncnode(node1,:);
        end


    end

    if nodeDOF(node2, 1) == 0  % Node1 has not been assigned DOFs
        % Assign the first DOF to node1
        nodeDOF(node2, 1) = nextDOF;  % Assign current nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the first DOF

        % Assign the second DOF to node1
        nodeDOF(node2, 2) = nextDOF;  % Assign updated nextDOF
        nextDOF = nextDOF + 1;        % Increment nextDOF after assigning the second DOF

        % Assign the third DOF to node1
        nodeDOF(node2, 3) = nextDOF;  % Assign updated nextDOF
        nextDOF = nextDOF + 1;
        Ndof(node2, 1:3) = nodeDOF(node2, 1:3);
        Nelem(1,4:6,i) = nodeDOF(node2, :);

        if modifyNode2 == 1  % Modify node1's third DOF

            ncnode(node2,1:2) =  nodeDOF(node2,1:2);
            ncnode(node2,3)   = -nodeDOF(node2,3);

            Nelem(1,4:6,i)  = abs(ncnode(node2,:));
            Ndof(node2, :)  = ncnode(node2,:);
        else
            ncnode(node2,1:3) =  nodeDOF(node2,1:3);
            Nelem(1,4:6,i)  = nodeDOF(node2, :);
            Ndof(node2, :)  = ncnode(node2,:);
        end


    elseif nodeDOF(node2, 1) ~= 0


        if modifyNode2 == 0 && Ndof(node2, 3) <= 0 % Modify node1's third DOF

            nodeDOF(node2, 3) = nextDOF;  % Assign updated nextDOF
            nextDOF = nextDOF + 1;
            ncnode(node2,:) =  nodeDOF(node2,:);
            Ndof(node2, 3) = ncnode(node2, 3);
            Nelem(1,4:6,i)  = ncnode(node2,:);
        elseif modifyNode2 == 0 && Ndof(node2, 3) >= 0
            ncnode(node2,:) =  nodeDOF(node2,:);

            Nelem(1,4:6,i)  = Ndof(node2, :);
        elseif modifyNode2 == 1
            nodeDOF(node2, 3) = nextDOF;  % Assign updated nextDOF
            ncnode(node2,:) =  nodeDOF(node2,:);
            nextDOF = nextDOF + 1;
            Nelem(1,4:6,i)  = ncnode(node2,:);
        end

        % else
        %       ncnode(node2,1:3) =  nodeDOF(node2,1:3);
        %
        %       Nelem(1,4:6,i)  = nodeDOF(node2, :);
        %       Ndof(node2, :)  = ncnode(node2,:);
    end

end

% disp(Nelem)
% disp(Ndof)



% disp(Nelem)
totaldof = (nextDOF-1);


alldof  = zeros(3*size(Ndof,1),1);
U       = zeros(3*size(Ndof,1),1);
Elastic = zeros(3*size(Ndof,1),1);
% disp(size(alldof))



%-------------------------------------------------------
% Table of material properties
%
% 1.E  2.alpha  3.rho
%-------------------------------------------------------


totalemp = fscanf(fid, '\nMATERIAL PROPERTIES = %d', 1);

if totalemp ~= stufflist(3,1)
    disp('error: Number of material properties is wrong')
end

%empdata是MATERIAL PROPERTIES的data矩阵
empdata = fscanf(fid, '\n%d %e %e %f',[4,totalemp]);

empdata = empdata';

%element material properties

%read all type of material properties
Elist     = empdata(:,2);
alphalist = empdata(:,3);
rholist   = empdata(:,4);

%initial allocation
E     = zeros(size(emp));
alpha = zeros(size(emp));
rho   = zeros(size(emp));
% loop over elements
for i = 1:length(emp)
    if emp(i) == 1
        E(i) = Elist(1);
        alpha(i) = alphalist(1);
        rho(i) = rholist(1);
    elseif emp(i) == 2
        E(i) = Elist(2);
        alpha(i) = alphalist(2);
        rho(i) = rholist(2);
    elseif emp(i) == 3
        E(i) = Elist(3);
        alpha(i) = alphalist(3);
        rho(i) = rholist(3);
    end
end








%-------------------------------------------------------
% Table of boundary conditions
%
% 1.fix boundary conditions 2.prescribed displacements
%-------------------------------------------------------

nnodefix = fscanf(fid,'\nNODES_WITH_PRESCRIBED_DISPLACEMENTS = %d',1);
if nnodefix ~= stufflist(4,1)
    disp('error: Number of BC is wrong')
end
ipresc   = fscanf(fid, '\n%d %d %d %d %f %f %f %f %f %f', [10,nnodefix]);
ipresc   = ipresc';
Bnode    = ipresc(:,1);


for ibc=1:size(ipresc,1)

    bnode = Ndof(Bnode(ibc),:);
    alldof(bnode)  = ipresc(ibc, 2:4);
    U(bnode)       = ipresc(ibc, 5:7);
    Elastic(bnode) = ipresc(ibc, 8:10);

end





%-------------------------------------------------------
% Table of load
%
% 1.P1  2.P2  3.UDLFinternal 4.F
%-------------------------------------------------------

%----  ----  ----  ----
% 1. Table of point load
%----  ----  ----  ----
nnodP = fscanf(fid,'\nNODES_WITH_POINT_LOAD = %d',1);

%initial allocation
F  = zeros(totaldof,1);
P1 = zeros(totaldof,1);
P2 = zeros(totaldof,1);

Pndof = zeros(nnodP,3);
Pnode = zeros(nnodP,1);

UDLFinternal = zeros(size(connectivity, 1), 2, 3);


if nnodP ~= stufflist(5,1)
    disp('error: Number of point load is wrong')
end
if nnodP == 0
    disp('No point load')
else
    Pdata1 = fscanf(fid,'\n%d %f %f %f', [4, nnodP]);
    Pdata1 = Pdata1';
    dnode = Pdata1(:,2);

    for iP = 1:nnodP


        Pnode(iP) = Pdata1(iP,1);

        P1(Ndof(Pnode(iP),1:2)) = Pdata1(iP,2:3);
        if Ndof(Pnode(iP),3) >= 0
            P1(Ndof(Pnode(iP),3)) = Pdata1(iP,4);
        end

    end
    F=P1;
end



%----  ----  ----  ----
% 2. Table of udl load
%----  ----  ----  ----

nelep = fscanf(fid,'\nELEMENT_WITH_UDL_LOAD = %d',1);
if nelep ~= stufflist(6,1)
    disp('error: Number of UDL load is wrong')
end
w = zeros(size(connectivity,1),1);
if nelep ~= 0
    Pdata2 = fscanf(fid,'\n%d %d %f', [3, nelep]);
    Pdata2 = Pdata2';
    %所有有udl的element1的集合
    elegather    = Pdata2(:,1);

    udldirection = Pdata2(:,2);

    udl          = Pdata2(:,3);


    % use for print BM figure



    for i = 1:size(elegather,1)

        w(elegather(i),1) = w(elegather(i),1) + udl(i);
    end

else

end

for ielem=1:nelep
    %extract connectivity information
    nodei=connectivity(elegather(ielem),1);
    idof=Nelem(1,1:3,elegather(ielem));
    nodej=connectivity(elegather(ielem),2);
    jdof=Nelem(1,4:6,elegather(ielem));

    %compute geometric information
    Xij=X(nodej,:)-X(nodei,:);
    L=norm(Xij,2); c=Xij(1)/L; s=Xij(2)/L;

    Ppoint   = udl(ielem)*L/2;
    Pbending = udl(ielem)*L^2/12;


    % udldirection = 0 means udl  apply in local element direction

    if udldirection(ielem) == 0 && s <= 0

        Px       =  Ppoint*s;
        Py       =  Ppoint*c;
        Pz       =  Pbending;
        Ptotal1   = [Px;Py;Pz];
        Ptotal2   = [Px;Py;-Pz];
        %assemble global stiffness matrix
        P2(idof,1)=P2(idof,1)+Ptotal1; P2(jdof,1)=P2(jdof,1)+Ptotal2;




        F=P1+P2;

    elseif udldirection(ielem) == 0 && udl <= 0 && s >= 0
        Px       =  -Ppoint*s;
        Py       =  Ppoint*c;
        Pz       =  Pbending;
        Ptotal1   = [Px;Py;Pz];
        Ptotal2   = [Px;Py;-Pz];
        %assemble global stiffness matrix
        P2(idof,1)=P2(idof,1)+Ptotal1; P2(jdof,1)=P2(jdof,1)+Ptotal2;



        F=P1+P2;

        % udldirection = 1 means udl only apply in global x direction
    elseif udldirection(ielem) == 1
        Px       = Ppoint;
        Py       = Ppoint;
        Pz       = Pbending;
        Ptotal1   = [Px;Py;-Pz];
        Ptotal2   = [Px;Py;Pz];
        %assemble global stiffness matrix
        P2(idof,1)=P2(idof,1)+Ptotal1; P2(jdof,1)=P2(jdof,1)+Ptotal2;

        F=P1+P2;

        % udldirection = 2 means udl apply in global y direction
    elseif udldirection(ielem) == 2
        Px       = Ppoint;
        Py       = Ppoint;
        Pz       = Pbending;
        Ptotal1   = [Px;Py;Pz];
        Ptotal2   = [Px;Py;-Pz];

        %assemble global stiffness matrix
        P2(idof,1)=P2(idof,1)+Ptotal1; P2(jdof,1)=P2(jdof,1)+Ptotal2;
        % P2(alldof == 1)=0;

        F=P1+P2;

    end



    %----  ----  ----  ----  ----
    % 3. add UDLFinternal effect
    %----  ----  ----  ----  ----

    % Loop over elements affected by UDL
    for ie = 1:size(elegather, 1)




        % Extract connectivity information
        nodei = connectivity(elegather(ie), 1);
        nodej = connectivity(elegather(ie), 2);

        % Compute geometric information
        Xij = X(nodej, :) - X(nodei, :);
        L = norm(Xij, 2);
        c = Xij(1) / L;
        s = Xij(2) / L;

        % Calculate UDL effects
        Ppoint = udl(ie) * L / 2;
        Pbending = udl(ie) * L^2 / 12;


        if udldirection(ie) == 0
            Px = Ppoint * s;
            Py = Ppoint;
            Pz = Pbending;
            % Fill UDLFinternal matrix
            UDLFinternal(elegather(ie), 1, 1) = 0; % Axial force at Node 1 (assuming no axial force due to UDL)
            UDLFinternal(elegather(ie), 1, 2) = -Py; % Shear force at Node 1
            UDLFinternal(elegather(ie), 1, 3) = -Pz; % Bending moment at Node 1

            UDLFinternal(elegather(ie), 2, 1) = 0; % Axial force at Node 2 (assuming no axial force due to UDL)
            UDLFinternal(elegather(ie), 2, 2) = -Py; % Shear force at Node 2
            UDLFinternal(elegather(ie), 2, 3) = Pz; % Bending moment at Node 2


            % udldirection = 1 means udl only apply in global x direction
        elseif  udldirection(ie) == 1
            Px = Ppoint;
            Pz = Pbending;
            % Fill UDLFinternal matrix
            UDLFinternal(elegather(ie), 1, 1) = -Px; % Axial force at Node 1 (assuming no axial force due to UDL)
            UDLFinternal(elegather(ie), 1, 2) = 0; % Shear force at Node 1
            UDLFinternal(elegather(ie), 1, 3) = -Pz; % Bending moment at Node 1

            UDLFinternal(elegather(ie), 2, 1) = -Px; % Axial force at Node 2 (assuming no axial force due to UDL)
            UDLFinternal(elegather(ie), 2, 2) = 0; % Shear force at Node 2
            UDLFinternal(elegather(ie), 2, 3) = Pz; % Bending moment at Node 2

            % udldirection = 2 means udl apply in global y direction

        elseif  udldirection(ie) == 2
            Py = Ppoint;
            Pz = Pbending;
            % Fill UDLFinternal matrix
            UDLFinternal(elegather(ie), 1, 1) = 0; % Axial force at Node 1 (assuming no axial force due to UDL)
            UDLFinternal(elegather(ie), 1, 2) = -Py; % Shear force at Node 1
            UDLFinternal(elegather(ie), 1, 3) = -Pz; % Bending moment at Node 1

            UDLFinternal(elegather(ie), 2, 1) = 0; % Axial force at Node 2 (assuming no axial force due to UDL)
            UDLFinternal(elegather(ie), 2, 2) = -Py; % Shear force at Node 2
            UDLFinternal(elegather(ie), 2, 3) = Pz; % Bending moment at Node 2

        end



    end

end

end
